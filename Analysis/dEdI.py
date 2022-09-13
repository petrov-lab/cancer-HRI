'''
	 Functions used to calculate dE/dI via breakpoint frequency and fractional overlap of CNs
	 by patient sample and an optionally specified gene track.
'''

#################
### Libraries ###
#################

import pandas as pd
import numpy as np
from MapGenes import create_gene_map
from bootstrap import sample
import warnings
warnings.simplefilter("ignore", UserWarning)

#################
### Functions ###
#################


def ReadAndAnnotateCNAs():
    ''' 
        Reads in raw COSMIC CN calls and returns a dataframe of CNs annotated by length, 
        and the CN change (amp/del). CNAs must have `Chromosome`,  `Start`, `End`, `CN` and 
        `Sample` in index.
    '''
    DataDir= "/labs/ccurtis2/tilk/scripts/hri-data/Raw/CNV/"
    raw_CNAs = pd.read_csv(DataDir + "CosmicCompleteCNA.tsv.gz", sep='\t')
    discarded_studies = {619}     # This COSMIC Study ID is cell line data, not cancers
    expanded = raw_CNAs['Chromosome:G_Start..G_Stop'].str.extract("(?P<Chromosome>[0-9]+):(?P<Start>\d+)..(?P<End>\d+)") # weird Chr:Start..Stop annotation to extract
    CNAs = expanded[['Start', 'End']].astype(int).assign(
        Chromosome = expanded['Chromosome'].map(lambda x: {'23':'X', '24':'Y'}.get(x, x)),   # Map Chr 23 -> X, 24 -> Y
        CN=raw_CNAs['MUT_TYPE'].map({'LOSS':-1, 'GAIN':+1, np.nan:0}).astype(int),
        Sample=raw_CNAs['ID_SAMPLE'],
		SAMPLE_NAME=raw_CNAs['SAMPLE_NAME'], # Used to get TCGA Barcode IDs
        Study=raw_CNAs['ID_STUDY']).query('End > Start and Study not in @discarded_studies')    # Remove CNAs with idiosyncratic lengths
    return(CNAs)



def GetGenomeTrack():
	'''
		Returns a dictionary of gene categories like passengers/drivers (which are the keys) and the
		copy number alterations fall within those groups of genes (which are the items). By default, 
		the dictionary contains one key ('All') of all mappabable genes.
	''' 
	AnnotationsDir=os.getcwd() + "/dEdI_Calculator/annotations/"
	chromosome_lengths = pd.read_csv(AnnotationsDir + 'biomart_chromosomes.csv', 
				index_col=['Chromosome'])['Length'].drop('MT') ### normal chromosome lengths
	all_genes = pd.read_pickle(AnnotationsDir + 'biomart_gene_annotations.pkl.gz')
	all_tracks = all_genes.assign(
	Chromosome  = all_genes.pop('Chromosome/scaffold name'),
	Start       = all_genes.pop('Transcript start (bp)'),
	End         = all_genes.pop('Transcript end (bp)')).query('Chromosome in @chromosome_lengths.index')
	track_sets = {'All': all_tracks}
	return(track_sets)


def MapAdditionalGeneSetsToGenomeTrack(track_sets, gene_set_file_name):
	'''
		Takes an input file path @gene_set_file_name and reads it in as a dataframe, @gene_sets, 
		which should contain a column of gene names ('GeneName') and 
		a column of unique string identifiers for a set of genes ('Type'). 
		Iterates through all available gene name identifiers (ENSG, ENST, Hugo Symbol, etc.) to return
		a list of positions associated with the gene set. Adds unique string identifier 'Type' to 
		the dictionary @track_sets.
	'''
	all_tracks=track_sets['All']
	gene_mappings=dict( 
	prev_symbol     = 'Gene name',
	alias_symbol    = 'Gene name',
	ensembl_gene_id = 'Gene stable ID')
	gene_sets = pd.read_csv(gene_set_file_name, sep=',')
	for gene_type in gene_sets['Type'].unique():    # Map each gene set separately
		print('Finding coordinate positions for gene set: ' + str(gene_type))
		unique_gene_names = gene_sets['GeneName'].unique()
		# Begin mapping `Gene names`, i.e. Hugo symbols
		mapped_exons = all_tracks['Gene name'].isin(unique_gene_names)
		gene_group = all_tracks.loc[mapped_exons]
		mapped = gene_group['Gene name'].unique()
		print('{:.2%} of {:} were mapped using Gene names.'.format(len(mapped)/len(unique_gene_names), gene_type))
		remainder = set(unique_gene_names) - set(mapped)
		if len(remainder) != 0: # There are still genes that haven't mapped
			for output_type, exon_type in gene_mappings.items():    # Iterate through all other possible nomenclatures
				transcripts = create_gene_map(np.array(list(remainder)), output_type, input_type='symbol', verbose=False)
				new_tracks = all_tracks.loc[all_tracks[exon_type].astype(str).isin(transcripts.astype(str))]
				gene_group = pd.concat([gene_group, new_tracks])
				mapped_output = frozenset(new_tracks[exon_type].unique().tolist())
				mapped = {transcript for transcript, mapping in transcripts.items() if mapping in mapped_output}
				print('{:.2%} of remainder were mapped via `{:}` -> `{:}`'.format(
					len(mapped)/len(remainder), output_type, exon_type))
				remainder -= set(mapped)
				print('{:.2%} could not be mapped.'.format(len(remainder)/len(unique_gene_names)))
		track_sets[ gene_type] = gene_group
	return(track_sets)



def GetFunctionalRegions():
	'''
		Returns a dataframe of regions that are considered "functional" in the human genome.
		'Functional' is the set of all genomic tracks that will fall into dE numerator.
		Everything that is not "functional" will fall into the dI (intergenic/intronic) regions of the genome, 
		used for the denominator. 
	'''
	all_tracks = GetGenomeTrack()['All']
	types = all_tracks['Gene type'].unique()
	meta_categories = dict(
	pseudogenes = {t for t in types if 'pseudogene' in t},
	ribosomal = {t for t in types if 'rRNA' in t or 'tRNA' in t} | {'ribozyme', 'vaultRNA', 'snoRNA', 'snRNA', 'scaRNA'},
	non_coding_RNAs = {t for t in types if 'ncRNA' in t} | {'non_coding'},
	regulatory = {'scRNA', 'miRNA', 'antisense', 'sense_intronic', 'sense_overlapping'} )
	functional_types = {'protein_coding'} | meta_categories['ribosomal'] | meta_categories['regulatory']  
	functional = all_tracks.loc[all_tracks['Gene type'].isin(functional_types)]
	return(functional)

def CalculateOverlaps(track, CNAs):
	'''
		Calculates dE/dI of every CNA with a genomic track 
		(both by fractional overlap and whether the breakpoint overlaps).
	'''
	functional = GetFunctionalRegions()
	mask = functional.groupby('Chromosome').apply(GetCentromeresAndTelomeres) 
	unmapped = mask.eval('start + cen_end - cen_start')
	chr_in_track = track['Chromosome'].iloc[0]
	CNAs['Chromosome'] = CNAs['Chromosome'].astype(str) ### maybe not necesary
	chromosomes_in_scnas = CNAs.groupby('Chromosome').groups.keys()
	if chr_in_track in chromosomes_in_scnas:
		cna_set = CNAs.groupby('Chromosome').get_group(chr_in_track)
		track_vector = np.zeros(max(mask.loc[chr_in_track, 'end'], cna_set['End'].max()) + 1, dtype=np.uint32)
		for start, stop in track[['Start', 'End']].drop_duplicates().values:
			track_vector[start:stop] = 1
		cumulative_track = track_vector.cumsum()
		normalization = cumulative_track[-1]/(len(track_vector) - unmapped[chr_in_track])
		return cna_set.assign(fractional_overlap = (cumulative_track[cna_set.End] - cumulative_track[cna_set.Start])/(cna_set.End - cna_set.Start),
							breakpoint_frequency = (track_vector[cna_set.Start] + track_vector[cna_set.End])/2).set_index(
				['Sample','SAMPLE_NAME', 'Study', 'Start', 'End', 'CN'], append=True).drop('Chromosome', axis=1)*(1/normalization)
	else:
		return(pd.DataFrame())

def GetCentromeresAndTelomeres(df):
	'''
		Returns a pandas Series of centromeric and telomeric regions of every chromosome.
		These regions will be masked (i.e. removed) from the analysis since we can't properly 
		characterize SCNAS in highly-repetitive regions of chromosomes.
		Telomeres = Start of first gene & End of last gene
		Centromere = Boundaries of largest gap in genes
		This methods works better than using the actual probe locations in an Affymetrix SNP6 array.
	'''
	V = np.r_[df['Start'].values, df['End'].values]
	V.sort()
	centromere = np.diff(V).argmax()
	return pd.Series({
		'start'     : df['Start'].min(), 
		'cen_start' : V[centromere],
		'cen_end'   : V[centromere+1],
		'end'       : df['End'].max()})

def GetMeanOverlapsByLength(S):
	'''
		Calculates mean dE/dI of samples by length of CNA.
	'''
	return S.groupby(['Track','Length_Category','xbin','NumberOfTumorsInBin']).mean()

def GetMeanOverlapsByType(S):
	'''
		Calculates mean dE/dI of samples by cancer type.
	'''
	return S.groupby(['Track','type','xbin','NumberOfTumorsInBin']).mean()

def GetMeanOverlapsByBroadType(S):
	'''
		Calculates mean dE/dI of samples by broad cancer groupings.
	'''
	return S.groupby(['Track','broad_type','xbin','NumberOfTumorsInBin']).mean()

def GetMeanOverlaps(S):
	'''
		Calculates mean dE/dI of samples by cancer type..
	'''
	return S.groupby(['Track','xbin','NumberOfTumorsInBin']).mean()

def AnnotateByLength(stats):
	'''
		Adds an additional column 'Length Category' that stratifies
		CNAs as either >100Kb in length or <=100Kb. 
	'''
	stats['Length'] = stats['End'] - stats['Start'] 
	stats['Length_Category'] = np.where(stats['Length'] < int(1e5), '<100Kb','>100Kb')
	return(stats)

def AddBinsByMutRate(df):
	'''
		Adds additional columns 'bin' and 'binStart', which dictate what mutational burden
		bin of CNAs that patient falls into and how many CNAs that bin contains.
	'''
	MutRate = df.groupby(['Sample']).size().reset_index().rename(columns={0:"MutRate"})
	MutRate['xbin'] = pd.cut(MutRate['MutRate'], bins=[0,3,10,30,100,300,1000,3000,10000,30000], 
                    labels=['1-3','3-10','10-30','30-100','100-300','300-1000','1000-3000','3000-10000','10000-30000'])
	MutRate = MutRate.merge(MutRate.groupby(['xbin']).size().reset_index(), left_on = 'xbin', right_on = 'xbin').rename(
		columns={0 : 'NumberOfTumorsInBin'})
	df = df.merge(MutRate, left_on = 'Sample', right_on = 'Sample')
	return(df)

def AddNumberOfMutsInBins(df):
	NumMuts = df.groupby(df.index.names).size().reset_index().rename(columns={0:'NumMutsInBin'})
	NewIndex = df.index.names + ['NumMutsInBin']
	if 'type' in NewIndex:
		df = df.reset_index().merge(NumMuts, left_on = ['xbin','type','NumberOfTumorsInBin','Track'],
			right_on = ['xbin','type','NumberOfTumorsInBin','Track'])
	if 'broad_type' in NewIndex:
		df = df.reset_index().merge(NumMuts, left_on = ['xbin','broad_type','NumberOfTumorsInBin','Track'],
			right_on = ['xbin','broad_type','NumberOfTumorsInBin','Track'])
	elif 'Length_Category' in NewIndex:
		df = df.reset_index().merge(NumMuts, left_on = ['xbin','Length_Category','NumberOfTumorsInBin','Track'],
			right_on = ['xbin','Length_Category','NumberOfTumorsInBin','Track'])
	else:
		df = df.reset_index().merge(NumMuts, left_on = ['xbin','NumberOfTumorsInBin','Track'],
			right_on = ['xbin','NumberOfTumorsInBin','Track'])
	return(df.set_index(NewIndex))

def AddCancerType(muts):
	'''Reads in a dataframe of mutations in TCGA and annotates the cancer type of each patient barcode.'''
	AnnotationsDir="/labs/ccurtis2/tilk/scripts/hri-data/Annotations/"
	SampleSourceCodes=pd.read_csv(AnnotationsDir + 'TCGA/TCGA_TSS_SourceSite.txt', sep='\t', dtype=str, keep_default_na=False)
	CancerTypeCodes=pd.read_csv(AnnotationsDir + 'TCGA/TCGA_CancerTypeNames.txt', sep='\t', header=None)
	MergedCodes = SampleSourceCodes.merge(CancerTypeCodes, left_on='Study Name', right_on=1, how='left')[['TSS Code',0]]
	MergedCodes['TSS Code'] = MergedCodes['TSS Code'].astype(str).str.zfill(2) ### these are not getting imported as 01, but as 1
	TCGABarcodesInCOSMIC = muts[muts['SAMPLE_NAME'].str.contains('TCGA', na=False)]['SAMPLE_NAME'].drop_duplicates()
	Barcodes=pd.DataFrame(TCGABarcodesInCOSMIC.str.split("-", 2).tolist(), columns = ['TCGA','sampleSource','restOfBarcode'])
	Barcodes['FullBarcode'] = TCGABarcodesInCOSMIC.tolist()
	Barcodes = Barcodes.merge(MergedCodes, right_on='TSS Code', left_on='sampleSource', how='left').rename(columns={0:'type'})
	return(muts.merge(Barcodes[['FullBarcode','type']], left_on='SAMPLE_NAME', right_on='FullBarcode', how='left'))

def GetCancerGroups():
    CancerGroups = pd.Series(dict(
        LAML='Circulatory', ACC='Endocrine', BLCA='Urinary', LGG='Nervous',
        BRCA='Reproductive', CESC='Reproductive', CHOL='Digestive', COAD='Digestive',
        ESCA='Digestive', GBM='Nervous', HNSC='Respiratory', KICH='Urinary', KIRC='Urinary',
        KIRP='Urinary', LIHC='Digestive', LUAD='Respiratory', LUSC='Respiratory', DLBC='Circulatory',
        MESO='Skeletal', OV='Reproductive', PAAD='Digestive', PCPG='Endocrine',PRAD='Reproductive',
        READ='Digestive', SARC='Skeletal', SKCM='Skin', STAD='Digestive', TGCT='Reproductive',
        THYM='Endocrine', THCA='Endocrine', UCS='Reproductive', UCEC='Reproductive', UVM='Skin',
        BTCA='Digestive', BOCA='Skeletal', CLLE='Circulatory', CMDI='Circulatory', EOPC='Reproductive',
        ESAD='Digestive', GACA='Digestive', LINC='Digestive', LIRI='Digestive', MALY='Circulatory',
        MELA='Skin', ORCA='Digestive', PAEN='Endocrine', PBCA='Nervous', RECA='Urinary'
    ))
    return(pd.DataFrame(CancerGroups).reset_index().rename(columns={'index':'type',0:'broad_type'}))


def GetdEdIByMutationalBin(CNAs=ReadAndAnnotateCNAs(), ByGroup = 'Length'):
	'''
		Takes in a list of CNAs and outputs a table of dE/dI values
		(using fractional overlap and breakpoint frequencies).
	'''
	CNAs[['Start','End']] = CNAs[['Start','End']].astype(float).astype(int)
	TrackSets = GetGenomeTrack()
	AllTracks = TrackSets['All'] # Track of where all genes in the genome are
	GeneSets="/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/All_CN_Drivers.csv"
	if len(GeneSets) != 0:
		TrackSets = MapAdditionalGeneSetsToGenomeTrack(TrackSets, GeneSets) # Add in where driver and passenger genes are
		TrackSets['Passengers'] = AllTracks[~(AllTracks.index.isin(TrackSets['Drivers'].index))]
	Iter = {(TrackName, Chr):TrackSubset for TrackName, Track in TrackSets.items() for Chr, TrackSubset in Track.groupby('Chromosome')}
	EmptyList=[]
	for Types in Iter.keys():
		if len(Iter.get(Types)) != 0: 
			Output = CalculateOverlaps(Iter.get(Types), CNAs)
			Output['Track'] = list(Types)[0]
			Output['Chromosome'] = list(Types)[1]
			EmptyList.append(Output)
	Final = pd.concat(EmptyList).fillna(0).reset_index()
	Final = AnnotateByLength(Final)
	Final = AddBinsByMutRate(Final)
	Final = AddCancerType(Final)
	if ByGroup == 'CancerType':
		Final = Final[~Final['type'].isna()] # Only want CNAs where we know the cancer type
		Final = Final.set_index(['Track','type','xbin','NumberOfTumorsInBin'])
		Final = AddNumberOfMutsInBins(Final)
		BreakpointFreq = sample(Final['fractional_overlap'], GetMeanOverlapsByType).CI().reset_index().rename(columns={
			'low':'breakpointFreq_low','true':'breakpointFreq_mean','high':'breakpointFreq_high'})
		FractionalOverlap = sample(Final['breakpoint_frequency'], GetMeanOverlapsByType).CI().reset_index().rename(columns={
			'low':'fractionalOverlap_low','true':'fractionalOverlap_mean','high':'fractionalOverlap_high'})
		Out = pd.concat([FractionalOverlap,BreakpointFreq[['breakpointFreq_low','breakpointFreq_mean','breakpointFreq_high']]], axis=1)
	elif ByGroup == 'CancerTypeBroad':
		Final = Final[~Final['type'].isna()] # Only want CNAs where we know the cancer type
		Final = Final.merge(GetCancerGroups(), left_on='type', right_on='type' , how='left')
		Final = Final.set_index(['Track','broad_type','xbin','NumberOfTumorsInBin'])
		Final = AddNumberOfMutsInBins(Final)
		BreakpointFreq = sample(Final['fractional_overlap'], GetMeanOverlapsByBroadType).CI().reset_index().rename(columns={
			'low':'breakpointFreq_low','true':'breakpointFreq_mean','high':'breakpointFreq_high'})
		FractionalOverlap = sample(Final['breakpoint_frequency'], GetMeanOverlapsByBroadType).CI().reset_index().rename(columns={
			'low':'fractionalOverlap_low','true':'fractionalOverlap_mean','high':'fractionalOverlap_high'})
		Out = pd.concat([FractionalOverlap,BreakpointFreq[['breakpointFreq_low','breakpointFreq_mean','breakpointFreq_high']]], axis=1)
	elif ByGroup == 'Length': # By CNA length
		Final = Final.set_index(['Track','Length_Category','xbin','NumberOfTumorsInBin'])
		Final = AddNumberOfMutsInBins(Final)
		BreakpointFreq = sample(Final['fractional_overlap'], GetMeanOverlapsByLength).CI().reset_index().rename(columns={
			'low':'breakpointFreq_low','true':'breakpointFreq_mean','high':'breakpointFreq_high'})
		FractionalOverlap = sample(Final['breakpoint_frequency'], GetMeanOverlapsByLength).CI().reset_index().rename(columns={
			'low':'fractionalOverlap_low','true':'fractionalOverlap_mean','high':'fractionalOverlap_high'})
		Out = pd.concat([FractionalOverlap,BreakpointFreq[['breakpointFreq_low','breakpointFreq_mean','breakpointFreq_high']]], axis=1)
	else: # Just group by TMB
		Final = Final.set_index(['Track','xbin','NumberOfTumorsInBin'])
		Final = AddNumberOfMutsInBins(Final)
		BreakpointFreq = sample(Final['fractional_overlap'], GetMeanOverlaps).CI().reset_index().rename(columns={
			'low':'breakpointFreq_low','true':'breakpointFreq_mean','high':'breakpointFreq_high'})
		FractionalOverlap = sample(Final['breakpoint_frequency'], GetMeanOverlaps).CI().reset_index().rename(columns={
			'low':'fractionalOverlap_low','true':'fractionalOverlap_mean','high':'fractionalOverlap_high'})
		Out = pd.concat([FractionalOverlap,BreakpointFreq[['breakpointFreq_low','breakpointFreq_mean','breakpointFreq_high']]], axis=1)
	return(Out)
