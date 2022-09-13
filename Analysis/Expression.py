

'''
	This script calculates gene expression levels of different gene families in cancer samples.
'''

#################
### Libraries ###
#################

from scipy import stats
from scipy.stats import linregress
import statsmodels.api as sm
from dNdS import *
from dEdI import *

#################
### Functions ###
#################


def GetExpressionData():
	''' Read in raw expression data from COSMIC and annotates TCGA sample by cancer type.'''
	ExpressionDir="/labs/ccurtis2/tilk/02_DNDS/RNA-SeqAnalysis/COSMIC/"
	Expression = pd.read_csv(ExpressionDir + 'CosmicCompleteGeneExpression.tsv.gz', sep='\t')
	CancerType = GetData(Dataset='TCGA')[['Tumor_Sample_Barcode','type']].drop_duplicates()
	Expression = Expression.merge(CancerType.assign(SAMPLE_NAME= CancerType['Tumor_Sample_Barcode'].str[0:15]), 
		left_on='SAMPLE_NAME', right_on='SAMPLE_NAME', how='left')
	return(Expression)


def GetGeneNamesWithinGeneFamily(GeneFamilyName):
	'''
	Returns a series of the gene names (Hugo Symbols) associated with a specific genefamily.
	Options for @geneFamilyName are: 'Chaperonins', 'All', 'Proteasome' and 'HSP90'.
	'''
	ExpressionDir="/labs/ccurtis2/tilk/02_DNDS/RNA-SeqAnalysis/COSMIC/"		
	GeneSetDir="/labs/ccurtis2/tilk/scripts/hri-data/Annotations/GeneSets/"
	if GeneFamilyName == 'All': 
		return(pd.read_csv(ExpressionDir + 'CosmicCompleteGeneExpression.tsv.gz', sep='\t')['GENE_NAME'].drop_duplicates())
	elif GeneFamilyName == "HSP90":
		AllGenes = pd.read_csv(GeneSetDir + 'heatshockChaperoneGenes', sep='\t')
		AllGenes = AllGenes[AllGenes['Gene family description'] == 'Heat shock 90kDa proteins']
		return(AllGenes['Approved Symbol'])
	elif GeneFamilyName == "Chaperonins":
		AllGenes = pd.read_csv(GeneSetDir + 'heatshockChaperoneGenes', sep='\t')
		AllGenes = AllGenes[AllGenes['Gene family description'] == str(GeneFamilyName)]
		return(AllGenes['Approved Symbol'])
	elif GeneFamilyName == "Proteasome":
		return(np.array(['PSMC1','PSMC2','PSMC3','PSMC4','PSMC5','PSMC6','PSMD1','PSMD2','PSMD3']))


def GetExpressionDataForGeneFamiliesOfInterest():
	'''
	Returns gene expression levels of different pathways.
	'''
	OutputDir="/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/01_AnnotatedData/" # Save output
	ProteinFamilies = ['All','Chaperonins','HSP90','Proteasome']
	FinalDF = pd.DataFrame(columns = ['SAMPLE_NAME','SAMPLE_ID','type'])
	Expression = GetExpressionData()
	for Fam in ProteinFamilies:
		GeneSet = GetGeneNamesWithinGeneFamily(Fam)
		df = Expression[Expression['GENE_NAME'].isin(GeneSet)].groupby(
			['SAMPLE_NAME','SAMPLE_ID','type'])['Z_SCORE'].median().reset_index()
		FinalDF = FinalDF.merge(df, left_on= ['SAMPLE_NAME','SAMPLE_ID','type'], 
			right_on= ['SAMPLE_NAME','SAMPLE_ID','type'], how= 'right')
		FinalDF = FinalDF.rename(columns={'Z_SCORE': Fam})
	FinalDF.to_csv(OutputDir + "/ExpressionOfGeneSets_COSMIC.txt", sep='\t')
	return(FinalDF)


def GetSummarizedExpressionTable(AnnotatedExpression):
	''' Generates a pandas Dataframe of mean expression values with CI by gene groups of interest and adds metadata of counts in bins.'''
	Table = pd.concat([
				sample(AnnotatedExpression.set_index(['xbin'])['Chaperonins'], MeanByMutationBin).CI().assign(geneGroup='Chaperonins'),
				sample(AnnotatedExpression.set_index(['xbin'])['Proteasome'], MeanByMutationBin).CI().assign(geneGroup='Proteasome'),
				sample(AnnotatedExpression.set_index(['xbin'])['HSP90'], MeanByMutationBin).CI().assign(geneGroup='HSP90'),
				sample(AnnotatedExpression.set_index(['xbin'])['All'], MeanByMutationBin).CI().assign(geneGroup='All')]).reset_index()
	MetaData = pd.DataFrame({"NumMutsInBin" : AnnotatedExpression.groupby(['xbin']).size(),
                            "NumPatientsInBin" : AnnotatedExpression.groupby(['xbin'])['SAMPLE_NAME'].unique().apply(len)}).reset_index()
	return(Table.merge(MetaData, left_on='xbin',right_on='xbin'))


def MeanByMutationBin(df):
    ''' Calculate the mean for every bin of samples.'''
    return(df.groupby(['xbin']).mean())



def GetExpressionOfGeneSetsByMutationBin(MutType='SNV', ByCancerType=False):
	'''
	Returns expression of heat-shock and protein degradatiion pathways for SNVs and CNAs. 
	@mutType specifies whether to stratify tumors by string 'SNV' or 'CNA'.
	@byCancerType is a boolean that determines whether cancer-specific expression is calculated.
	'''
	Expression = GetExpressionDataForGeneFamiliesOfInterest() ### Generates the file "ExpressionOfGeneSets_COSMIC.txt"
	#Expression = pd.read_csv('/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/01_AnnotatedData/ExpressionOfGeneSets_COSMIC.txt', sep='\t')
	if MutType == 'SNV':
		Mut = GetData(Dataset='TCGA')
		if ByCancerType:
			Results = pd.DataFrame()
			for CancerSubgroup in Expression['type'].unique():
				MutsInGroup = Mut[Mut['type'] == CancerSubgroup]
				MutsInGroup = AnnotateMutationRate(MutsInGroup, MutType='kska', BinType='log').assign(SAMPLE_NAME = MutsInGroup['Tumor_Sample_Barcode'].str[0:15])
				AnnotatedExpression = Expression.merge(MutsInGroup[['SAMPLE_NAME','xbin']].drop_duplicates(), left_on = 'SAMPLE_NAME', right_on='SAMPLE_NAME')
				ResultsForOneCancerType = GetSummarizedExpressionTable(AnnotatedExpression).assign(type=CancerSubgroup)
				Results = pd.concat([Results, ResultsForOneCancerType])
			return(Results)
		else: # Look at pan-cancer expression by SNVs
			Mut = AnnotateMutationRate(Mut, MutType='kska', BinType='log')
			Mut['SAMPLE_NAME'] = Mut['Tumor_Sample_Barcode'].str[0:15]
			AnnotatedExpression = Expression.merge(Mut[['SAMPLE_NAME','xbin']].drop_duplicates(), left_on = 'SAMPLE_NAME', right_on='SAMPLE_NAME')
			return(GetSummarizedExpressionTable(AnnotatedExpression))
	elif MutType == 'CNV':
		CNAs = ReadAndAnnotateCNAs()
		if ByCancerType:
			Results = pd.DataFrame()
			for CancerSubgroup in Expression['type'].unique():
				CNAsInGroup = CNAs[CNAs['Sample'].isin(Expression[Expression['type'] == CancerSubgroup]['SAMPLE_ID'])]
				CNAsInGroup = AddBinsByMutRate(CNAsInGroup)
				AnnotatedExpression = Expression.merge(CNAsInGroup[['Sample','xbin']].drop_duplicates(), left_on = 'SAMPLE_ID', right_on='Sample')
				ResultsForOneCancerType = GetSummarizedExpressionTable(AnnotatedExpression).assign(type=CancerSubgroup)
				Results = pd.concat([Results, ResultsForOneCancerType])
			return(Results)
		else: # Look at pan-cancer expression by CNVs
			CNAs = AddBinsByMutRate(CNAs)
			AnnotatedExpression = Expression.merge(CNAs[['Sample','xbin']].drop_duplicates(), left_on = 'SAMPLE_ID', right_on='Sample')
			return(GetSummarizedExpressionTable(AnnotatedExpression))


def R2(df):
	X = df['Bin']
	Y = df['Z_SCORE']
	slope, intercept, r_value, p_value, std_err = linregress(X, Y)
	return(pd.DataFrame({'r_value':[r_value],'p_value':[p_value]}))


def GetExpressionByLoadForAllIndividualGenes():
	'''
	Calculates correlation coefficients of median expression by each mutational burden bin
	for all genes in the genome. Used for supplemental figure 22B.
	'''
	Expression = GetExpressionData()
	Mut = AnnotateMutationRate(GetData(Dataset='TCGA'), MutType='kska', BinType='log')
	Mut = Mut[['xbin','mutRate','Tumor_Sample_Barcode']].drop_duplicates()
	Mut['SAMPLE_NAME'] = Mut['Tumor_Sample_Barcode'].str[0:15]
	Expression = Expression.merge(Mut, left_on='SAMPLE_NAME', right_on='SAMPLE_NAME', how='left')
	MedianExpression = Expression.groupby(['GENE_NAME','xbin']).median().reset_index().dropna()
	NumTumorsInBin = MedianExpression.groupby(['GENE_NAME','xbin']).size().reset_index().rename(columns={0:"NumTumorsInBin"})
	MedianExpression = MedianExpression.merge(NumTumorsInBin, left_on=['GENE_NAME','xbin'], right_on = ['GENE_NAME','xbin'])
	MedianExpression['Bin'] = MedianExpression['xbin'].str.split('-', expand=True)[0]
	MedianExpression[['NumTumorsInBin','Bin','Z_SCORE']] = MedianExpression[['NumTumorsInBin','Bin','Z_SCORE']].astype(float)
	return(MedianExpression.groupby(['GENE_NAME']).apply(R2).reset_index())

