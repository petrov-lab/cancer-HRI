'''
	This script contains all of the functions used to annotate different datasets used in the manuscript.
'''

#################
### Libraries ###
#################

import pandas as pd
import glob
import os
import numpy as np
import os.path

###################
### Directories ###
###################

RawDataDir="/labs/ccurtis2/tilk/scripts/hri-data/Raw/SNV" # Where the raw data for each dataset is stored
ProcessedDir="/labs/ccurtis2/tilk/scripts/hri-data/Processed/SNV/ProcessedWithGeneAnnotations/" # Output directory
AnnotationsDir="/labs/ccurtis2/tilk/scripts/hri-data/Annotations/" # Extra annotations for each dataset (e.g. tumor purity)

#################
### Functions ###
#################


def AddCancerType(muts):
	'''Reads in a dataframe of mutations in TCGA and annotates the cancer type of each patient barcode.'''
	SampleSourceCodes=pd.read_csv(AnnotationsDir + 'TCGA/TCGA_TSS_SourceSite.txt', sep='\t', dtype=str, keep_default_na=False)
	CancerTypeCodes=pd.read_csv(AnnotationsDir + 'TCGA/TCGA_CancerTypeNames.txt', sep='\t', header=None)
	MergedCodes = SampleSourceCodes.merge(CancerTypeCodes, left_on='Study Name', right_on=1, how='left')[['TSS Code',0]]
	MergedCodes['TSS Code'] = MergedCodes['TSS Code'].astype(str).str.zfill(2) ### these are not getting imported as 01, but as 1
	Barcodes=pd.DataFrame(muts['Tumor_Sample_Barcode'].str.split("-", 2).tolist(), columns = ['TCGA','sampleSource','restOfBarcode'])
	Barcodes['FullBarcode'] = muts['Tumor_Sample_Barcode']
	Barcodes = Barcodes.merge(MergedCodes, right_on='TSS Code', left_on='sampleSource', how='left').rename(columns={0:'subtype'})
	muts['type'] = Barcodes['subtype']
	return(muts)



def AddAgeOfPatients(muts):
	'''
	Reads in a dataframe of mutations in TCGA and annotates the age of each patient at diagnosis based on clinical info.
	Clinical info for each cancer type in TCGA was downloaded from cBioPortal.
	'''
	clinical = []
	ClinicalInfoDir = AnnotationsDir + 'TCGA/TCGAClinical'
	FilesWithinDir = glob.glob(os.path.join(ClinicalInfoDir, "*tsv")) 
	for f in FilesWithinDir:
		temp = pd.read_csv(f, sep='\t')
		temp = temp[['Sample ID','Diagnosis Age']]
		clinical.append(temp)
	AllClinical = pd.concat(clinical, ignore_index=True).dropna()
	Barcodes = muts['Tumor_Sample_Barcode'].str.slice(0,15).reset_index()
	Barcodes = Barcodes.merge(AllClinical, right_on='Sample ID', left_on='Tumor_Sample_Barcode', how='left')
	muts['Diagnosis_Age'] = Barcodes['Diagnosis Age']
	return(muts)


def AddStageOfTumor(muts):
	'''Reads in a dataframe of mutations in TCGA and annotates the stage of each tumor based on clinical info.'''
	TumorStage = pd.read_csv(AnnotationsDir + 'TCGA/TCGAClinicalOutcome.tsv.gz', sep='\t')[['tumor_stage','submitter_id']]
	Stage1 = dict.fromkeys(['stage i','stage ia', 'stage ib'], 1)  
	Stage2 = dict.fromkeys(['stage ii',  'stage iia', 'stage iib', 'stage iic',], 2)  
	Stage3 = dict.fromkeys(['stage iii','stage iiia', 'stage iiib', 'stage iiic'], 3)  
	Stage4 = dict.fromkeys(['stage iv','stage iva','stage ivb', 'stage ivc'], 4)  
	TumorStage['tumor_stage'] = TumorStage['tumor_stage'].replace(Stage1).replace(Stage2).replace(Stage3).replace(Stage4)
	TumorStage = TumorStage[(TumorStage['tumor_stage'] == 1) | (TumorStage['tumor_stage'] == 2) | (TumorStage['tumor_stage'] == 3) | (TumorStage['tumor_stage'] == 4)].drop_duplicates()
	muts['barcodeID'] = muts['Tumor_Sample_Barcode'].str.slice(0,12)
	muts = muts.merge(TumorStage, left_on=['barcodeID'], right_on=['submitter_id'], how='left' )
	return(muts)


def AnnotatePatientsByMutRate(muts):
	'''
	Adds a separate column that annotats total number of mutations a patient has, 
	as well as the number of intergenic and intronic mutations. 
	'''
	mutRate = muts.groupby(["Tumor_Sample_Barcode", "Variant_Classification"]).size().reset_index()
	intron = mutRate[mutRate['Variant_Classification'] == "Intron"][['Tumor_Sample_Barcode', 0]].rename(columns={0:'totalNumIntronMutations'})
	total = muts.groupby(["Tumor_Sample_Barcode"]).size().reset_index().rename(columns={0:'totalNumMutations'})
	muts = muts.merge(intron, left_on='Tumor_Sample_Barcode', right_on='Tumor_Sample_Barcode', how='outer').fillna(0)
	muts = muts.merge(total, left_on='Tumor_Sample_Barcode', right_on='Tumor_Sample_Barcode', how='outer').fillna(0)
	return(muts)

def AddGeneSets(muts):
	'''
	Reads in many gene set files and adds a boolean column of whether or not the mutation overlaps with a particular gene set.
	'''
	GeneSetDirectory="/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/"
	PancancerDriver_COSMIC_genes = pd.read_csv(GeneSetDirectory + 'Driver_Census_2016_10_18.csv')['Gene Symbol'] # From COSMIC's `Driver Gene Census`; downloaded Oct, 18 2016
	PancancerDriver_Bailey_genes = pd.read_csv(GeneSetDirectory + 'drivers_Bailey2018.txt', sep='\t', header=None)[0]
	PancancerDriver_Intogen = pd.read_csv('/labs/ccurtis2/tilk/02_DNDS/separateDatabaseFiles/intogen_cancer_drivers-2014.12/Drivers_type_role.tsv', sep='\t', comment='#')
	PancancerDriver_Intogen = PancancerDriver_Intogen[PancancerDriver_Intogen['Driver_type'] == 'MUTATION']['geneHGNCsymbol']
	housekeeping_genes = pd.read_csv(GeneSetDirectory + "housekeeping_genes.txt.gz", delim_whitespace=True)['HUGO']
	essential_genes = pd.read_csv(GeneSetDirectory + 'essential_genes.txt.gz', header=None)[0]
	ts=pd.read_csv(GeneSetDirectory + 'tumorSuppressorsAndOncogenes.txt', sep='\t')
	onco=ts[ts['Type'] == 'ONC']['Gene']
	tsg=ts[ts['Type'] == 'TSG']['Gene']
	Chr_Segregation = pd.read_csv(GeneSetDirectory + 'ChromosomeSegregation_GO_Term_ID_0007059.clean.txtwithUpperCase.txt', sep='\t')['Symbol']
	normalGTEX_bottom50Percent = pd.read_csv(GeneSetDirectory + 'TotalMeanExpression_Bottom_19503', sep=',')['GeneName']
	normalGTEX_top50Percent = pd.read_csv(GeneSetDirectory +  'TotalMeanExpression_Top_19503', sep=',')['GeneName']
	transcription_regulation= pd.read_csv(GeneSetDirectory + 'TranscriptionRegulation_GO_Term_ID_0140110.clean.txtwithUpperCase.txt', sep='\t')['Symbol']
	translation_regulation= pd.read_csv(GeneSetDirectory + 'TranslationRegulation_GO_Term_ID_0045182.clean.txtwithUpperCase.txt', sep='\t')['Symbol']
	all_interacting_proteins = pd.read_csv(GeneSetDirectory + 'interaction_counts.csv', sep=',')
	### Add annotations ###
	annotations = pd.DataFrame({'pancancerDriver_COSMIC': muts['Hugo_Symbol'].isin(PancancerDriver_COSMIC_genes),
								'pancancerDriver_Bailey': muts['Hugo_Symbol'].isin(PancancerDriver_Bailey_genes), 
								'pancancerDriver_Intogen': muts['Hugo_Symbol'].isin(PancancerDriver_Intogen), 
								'Housekeeping':muts['Hugo_Symbol'].isin(housekeeping_genes),
								'Essential':muts['Hugo_Symbol'].isin(essential_genes),
								'TranscriptionRegulation':muts['Hugo_Symbol'].isin(transcription_regulation),
								'TranslationRegulation': muts['Hugo_Symbol'].isin(translation_regulation), 
								'All_Interacting_Proteins': muts['Hugo_Symbol'].isin(all_interacting_proteins['gene']), 
								'Oncogene':muts['Hugo_Symbol'].isin(onco), 
								'Chr_Segregation': muts['Hugo_Symbol'].isin(Chr_Segregation),
								'Tumor_Suppressor': muts['Hugo_Symbol'].isin(tsg),
								'normalGTEX_bottom50Percent':muts['Gene'].isin(normalGTEX_bottom50Percent), 
								'normalGTEX_top50Percent': muts['Gene'].isin(normalGTEX_top50Percent)})
	SpecificDriverGenes = pd.read_csv(GeneSetDirectory + 'TCGA_CancerSpecificDrivers')
	muts['SpecificDriverGenes'] = False
	for CancerType in SpecificDriverGenes['1'].unique():
		muts['SpecificDriverGenes']  = np.where(((muts['type'] == CancerType) & (muts['Hugo_Symbol'].isin(SpecificDriverGenes['0']))), 
			True, muts['SpecificDriverGenes'])
	muts = pd.concat([muts, annotations], axis=1)
	return(muts)



def SplitByCancerType(muts, OutDir, ByCancerType):
	''' Files are split by cancer type or patient ID to allow for easier parallel processing. '''
	if ByCancerType:
		Groups=muts['type'].unique().tolist()
	else: # Split by patient
		Groups=muts['Tumor_Sample_Barcode'].unique().tolist()
	for SubGroup in Groups:
		if ByCancerType:
			temp = muts[muts['type'] == SubGroup]
		else:
			temp = muts[muts['Tumor_Sample_Barcode'] == SubGroup]
		temp.to_csv(OutDir + SubGroup + '.txt.gz', sep='\t', compression='gzip')
	print('finished splitting data!')


def AnnotateTumorPurity(df):
    ''' Reads in a dataframe of mutations and adds a column of tumor purity calls for each patient.'''
    AnnotationsDir="/labs/ccurtis2/tilk/scripts/hri-data/Annotations/"
    purity = pd.read_csv(AnnotationsDir + 'TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.txt', sep='\t')
    purity['Tumor_Sample_Barcode'] = purity['array']
    purity = purity[['Tumor_Sample_Barcode', 'purity']].rename(columns={'Tumor_Sample_Barcode':'Sample'})
    df['Sample'] = df['Tumor_Sample_Barcode'].str[0:15]
    if 'purity' in df.columns:
        df = df.drop(columns={'purity'})
    df = df.merge(purity, left_on=['Sample'], right_on=['Sample'], how='left')
    df['VAF'] = df['t_alt_count']/(df['t_depth'])
    df['PurityAdjusted_VAF'] = df['VAF']/df['purity']
    return(df)
    

def AnnotateMC3Data():
	muts = pd.read_csv(RawDataDir + '/mc3.v0.2.8.PUBLIC.maf.gz',sep='\t')
	muts = AddCancerType(muts)
	muts = AnnotatePatientsByMutRate(muts)
	muts = AddAgeOfPatients(muts)
	muts = AddStageOfTumor(muts)
	muts = AddGeneSets(muts)
	muts = AnnotateTumorPurity(muts)
	muts['SIFT'] = muts['SIFT'].str.split('(',expand=True)[0]
	muts['PolyPhen'] = muts['PolyPhen'].str.split('(',expand=True)[0]
	SplitByCancerType(muts, ProcessedDir + 'TCGA/', True)
	return(muts)


AnnotateMC3Data()

