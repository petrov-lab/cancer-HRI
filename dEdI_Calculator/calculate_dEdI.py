'''
	 Functions used to calculate dE/dI via breakpoint frequency and fractional overlap of CNs
	 by patient sample and an optionally specified gene track.
'''

import pandas as pd
import numpy as np
import glob
import os
from gene_mapper import create_gene_map
from bootstrap import sample
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def ReadAndAnnotateCNAs(CNA_File):
	''' 
		Reads in a tab-separated file of CNAs.
		CNAs must have `Chromosome`, `Start`, `End`, and `Sample` in columns.
	'''
	CNAs = pd.read_csv(CNA_File, sep='\t')
	return(CNAs)

def GetGenomeTrack():
	'''
		Returns a dictionary of gene categories (keys) and the gene names/coordinates of those
		groups of genes (items). Dictionary contains one key ('All') of all mappabable genes.
	''' 
	chromosome_lengths = pd.read_csv(os.getcwd() + '/annotations/biomart_chromosomes.csv', index_col=['Chromosome'])['Length'].drop('MT') ### normal chromosome lengths
	all_genes = pd.read_pickle(os.getcwd() + '/annotations/biomart_gene_annotations.pkl.gz')
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
	gene_sets = pd.read_csv(gene_set_file_name, sep='\t')
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
				['Sample', 'Start', 'End'], append=True).drop('Chromosome', axis=1)*(1/normalization)
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

def GetMeanOverlaps(S):
	'''
		Calculates mean dE/dI of samples.
	'''
	return S.groupby(['Track','Sample']).mean()


def GetdEdI(CN_Path, out_file, gene_sets=""):
	'''
		Takes in a list of CNAs and outputs a table of dE/dI values
		(using fractional overlap and breakpoint frequencies).
	'''
	CNAs = ReadAndAnnotateCNAs(CN_Path)
	track_sets = GetGenomeTrack()
	all_tracks = track_sets['All']
	if len(gene_sets) != 0:
		track_sets = MapAdditionalGeneSetsToGenomeTrack(track_sets, gene_sets)
	Iter = {(track_name, Chr):track_subset for track_name, track in track_sets.items() for Chr, track_subset in track.groupby('Chromosome')}
	emptyList=[]
	for types in Iter.keys():
		if len(Iter.get(types)) != 0: 
			output = CalculateOverlaps(Iter.get(types), CNAs)
			output['Track'] = list(types)[0]
			output['Chromosome'] = list(types)[1]
			emptyList.append(output)
	final = pd.concat(emptyList).fillna(0).reset_index()
	final = final.set_index(['Track','Sample'])
	breakpointFreq = sample(final['fractional_overlap'], GetMeanOverlaps).CI().reset_index().rename(columns={'low':'breakpointFreq_low','true':'breakpointFreq_mean','high':'breakpointFreq_high'})
	fractionalOverlap = sample(final['breakpoint_frequency'], GetMeanOverlaps).CI().reset_index().rename(columns={'low':'fractionalOverlap_low','true':'fractionalOverlap_mean','high':'fractionalOverlap_high'})
	out = pd.concat([fractionalOverlap,breakpointFreq[['breakpointFreq_low','breakpointFreq_mean','breakpointFreq_high']]], axis=1)
	out.to_csv(out_file, sep='\t')




