'''
    This contains the set of functions used in the manuscript to calculate dN/dS and variosu other statistics for the supplemental figures.
'''

from os import replace
from dNdS import *
from dEdI import *
from MutationalBias import * 
import subprocess
from GetNullCNVs import *
from Expression import GetGeneNamesWithinGeneFamily
import random


def GetIndependentlySampledMutationsForAxes(dNdScv=False):
    '''
        Randomly subsamples mutations to be independent for dN/dS and mutation rate.
    '''
    NumReps=10
    OutDF = pd.DataFrame()
    df = GetPermutedCounts().reset_index()
    for RepNum in np.arange(0,NumReps):
        MutationRateMuts = df.sample(frac=0.5) # Randomly sample half of mutations for mutation rate
        dNdSMuts = df[~df.index.isin(MutationRateMuts.index)]
        MutationRate = AnnotateMutationRate(MutationRateMuts, MutType='kska', BinType='log')[['mutRate','xbin','Tumor_Sample_Barcode']].drop_duplicates()
        dNdSMuts = dNdSMuts.merge(MutationRate, left_on = 'Tumor_Sample_Barcode', right_on = 'Tumor_Sample_Barcode')
        Sets = dict(
                Passengers = dNdSMuts[dNdSMuts['pancancerDriver_Bailey'] == False],
                Drivers_Bailey = dNdSMuts[dNdSMuts['pancancerDriver_Bailey'] == True]
        )
        Result = GetDnDsTable(Sets, dNdScv).assign(rep=RepNum)
        OutDF = OutDF.append(Result)
    return(OutDF)


def GetLowlyExpressedGenes():
    ''' Defined as having TPM less than 1 across all samples in GTEX.'''
    OutDir='/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/'
    samp=pd.read_csv(OutDir + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', skiprows=2, sep='\t')
    samp['Name'] = samp['Name'].str.split('.', expand=True)[0]
    psuedo = samp.set_index(['Description','Name'])
    psuedo.loc[(psuedo < 1).all(axis=1)].reset_index()[['Description','Name']].to_csv(
        OutDir + 'LessThan1_UnexpressedAcrossAllSamplesGTEX.txt', sep='\t')


def GetdNdSInLowlyExpressedGenes(dNdScv=False):
    ''' Calculates dN/dS by mutational burden in lowly expressed/pseudo genes.'''
    #GetLowlyExpressedGenes() # Generates the file LessThan1_UnexpressedAcrossAllSamplesGTEX.txt
    GeneSetDir='/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/'
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType='kska', BinType='log')
    Pseudo=pd.read_csv(GeneSetDir + 'LessThan1_UnexpressedAcrossAllSamplesGTEX.txt', sep='\t')
    df['Pseudo'] = df['Gene'].isin(Pseudo['Name'])
    Sets = dict(
        Pseudo=df[df['Pseudo'] == True],
        PseudoNoEssential=df[(df['Pseudo'] == True) & (df['Housekeeping'] == False) & (df['Essential'] == False)], 
        PseudoNoDrivers=df[(df['Pseudo'] == True) & (df['pancancerDriver_Bailey'] == False)] # Extra sanity checks
    )
    return(GetDnDsTable(Sets, dNdScv))


def GetdNdSAfterShortGenesRemoved():
    ''' Removes short genes that have fewer than 10 permutations and calculates dN/dS. '''
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType='kska', BinType='log')
    ShortGenes = df[(df['pKa'] < 10 ) & (df['pKs'] < 10)]['Gene'].unique()
    df = df[~df['Gene'].isin(ShortGenes)] # Remove short genes
    Sets = dict(
            Passengers = df[df['pancancerDriver_Bailey'] == False],
            Drivers_Bailey = df[df['pancancerDriver_Bailey'] == True]
    )
    return(GetDnDsTable(Sets, dNdScv=False))


def GetFunctionalGeneSets(dNdScv=False):
    '''
        Calculates dN/dS for various functional gene sets by low TMB (<10) vs high TMB (> 10)
    '''
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType='kska', BinType='log')
    # Re-map bins into low and high; dash is added for downstream function @GetDnDsTable
    df.loc[(df['xbin'] == '1-3') | (df['xbin'] == '3-10'), 'xbin2'] = 'low-low1'
    df.loc[(df['xbin'] != '1-3') & (df['xbin'] != '3-10'), 'xbin2'] = 'high-high1'
    df['xbin'] = df['xbin2']
    Sets = dict(
        All = df[(df['pancancerDriver_Bailey'] == False)],
        TranslationRegulation = df[(df['pancancerDriver_Bailey'] == False) & (df['TranslationRegulation'] == True)],
        TranscriptionRegulation = df[(df['pancancerDriver_Bailey'] == False) & (df['TranscriptionRegulation'] == True)],
        Interacting_Protein = df[(df['pancancerDriver_Bailey'] == False) & (df['All_Interacting_Proteins'] == True)],
        Chr_Segregation = df[(df['pancancerDriver_Bailey'] == False) & (df['Chr_Segregation'] == True)],
        Passengers_Essential = df[(df['pancancerDriver_Bailey'] == False) & (df['Essential'] == True)],
        Passengers_Housekeeping = df[(df['pancancerDriver_Bailey'] == False) & (df['Housekeeping'] == True)],
        normalGTEX_top50Percent_exclDrivers = df[(df['pancancerDriver_Bailey'] == False) & (df['normalGTEX_top50Percent'] == True)]
    )
    return(GetDnDsTable(Sets, dNdScv))


def GetSubsampledFromHighMutRate(NumReps, dNdScv=False):
    '''
        Randomly subsamples passengers from high mutation rate tumors in the same proportions of  
        muts as are binned in Fig. 2A, and then calculates dN/dS.
    '''
    list = [] # List where data gets appended to 
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType='kska', BinType='log')
    BinSizes = df.groupby(['xbin']).size().reset_index().rename(columns={0:'NumMutInBin'})
    # Sample mutations from high TMB tumors
    BarcodesToSampleFrom = df[~df['xbin'].isin(['1-3','3-10','10000-30000'])]['Tumor_Sample_Barcode'] 
    MutsToSampleFrom = df[(df['Tumor_Sample_Barcode'].isin(BarcodesToSampleFrom)) & (df['pancancerDriver_Bailey'] == False)]
    for rep in np.arange(1, NumReps, 1):
        for bin in BinSizes.index: 
            print(BinSizes['NumMutInBin'][bin])
            MutsInSamePropAsBin = MutsToSampleFrom.assign(xbin=1).set_index('xbin').sample(n = BinSizes['NumMutInBin'][bin], replace=True)
            if dNdScv:
                dNdSInBin = GetdNdScv(MutsInSamePropAsBin.reset_index()).rename(columns={'mle': 'true', 'cilow': 'low', 'cihigh':'high'}).assign(
                    type=BinSizes['xbin'][bin], rep=rep, BinSize = BinSizes['NumMutInBin'][bin])
            else:
                dNdSInBin = MutsInSamePropAsBin.groupby('xbin').sum().eval('(Ka/pKa)/(Ks/pKs)').reset_index().assign(
                    type=BinSizes['xbin'][bin], rep=rep, BinSize = BinSizes['NumMutInBin'][bin])
            list.append(dNdSInBin.reset_index())
    return(pd.concat(list))
    


def GetdNdSAfterRemovingByTumorPurity(Dataset, dNdScv=False):
    ''' Removes tumors based on purity thresholds and calculates dN/dS. '''
    OutList = [] # Empty DF where results are appended to
    df = GetPermutedCounts()
    if 'purity' not in df.columns:
        df = AnnotateTumorPurity(df)
    df = AnnotateMutationRate(df, MutType='kska', BinType='log')
    PurityThresholds = np.arange(0, 1, 0.1)
    for threshold in PurityThresholds:
        ImpureRemoved = df[df['purity'] >= threshold]
        Sets = dict(
            Passengers = ImpureRemoved[ImpureRemoved['pancancerDriver_Bailey'] == False],
            Drivers_Bailey = ImpureRemoved[ImpureRemoved['pancancerDriver_Bailey'] == True]
        )
        OutList.append(GetDnDsTable(Sets, dNdScv).assign(PurityThreshold = threshold))
    return(pd.concat(OutList))



def GetMutationalBias(NumReps):
    '''
        Subsets observed mutations according to specific mutation biases (e.g., Signature 1 vs Signature 30).
        Randomly removes mutations (simulating purifying selection) to specified true dN/dS values.
    '''
    MutationSet = GetPermutedCounts()
    dNdS_Values = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    SignatureNumbers = [1,2,3,4,5,6,8,9]
    ListOfProps = [1000, 10000, 100000]
    out = pd.DataFrame()
    for rep in np.arange(0,NumReps,1):
        for SignatureNum in SignatureNumbers:
            for ProportionOfSites in ListOfProps:
                print(SignatureNum, ProportionOfSites, rep)
                temp = GetMutationalBiasResults(dNdS_Values, SignatureNum, ProportionOfSites, MutationSet)
                out = out.append(temp)      
    return(out)




def GetCorrelationOfTumorPurityByMutationRate(Dataset):
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType='kska', BinType='log')
    if 'purity' not in df.columns:
        df = AnnotateTumorPurity(df)
    return(df[['Tumor_Sample_Barcode','xbin','purity','mutRate']].drop_duplicates())


def RunAnnovarToGetOverlapWithCommonVariants(Dataset, MAF):
    '''
    Runs Annovar to annotate mutations by their frequency in the 1000 Genomes Project.
	Required columns to run Annovar are: Chromosome, Start, Stop, Ref, Alt.
    Extra columns of patient barcode and mutation are added to merge files back later.
    '''
    AnnotationsDir="/labs/ccurtis2/tilk/scripts/hri-data/Annotations/"
    AnnovarDir="/labs/ccurtis2/tilk/software/annovar/"
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType = 'kska', BinType = 'log')
    df['ID'] = list(range(0, len(df)))
    AnnovarInput = df[['Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Tumor_Sample_Barcode','ID','Ka','Ks','xbin','pancancerDriver_Bailey']]
    AnnovarInput.to_csv(AnnotationsDir + Dataset + '/AnnovarInput' , sep= ' ', index=False, header=None)
    if MAF == '1':
        subprocess.call(['bash', '-c', 'perl ' + AnnovarDir + 'annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg19 -out ' +
            AnnotationsDir + Dataset + '/AnnovarOutput' + str(MAF) + ' ' + AnnotationsDir + Dataset + '/AnnovarInput ' + AnnovarDir + 'humandb'])
    else:
        subprocess.call(['bash', '-c', 'perl ' + AnnovarDir + 'annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg19 -out ' +
            AnnotationsDir + Dataset + '/AnnovarOutput' + str(MAF) + ' ' + AnnotationsDir + Dataset + '/AnnovarInput ' + AnnovarDir + 'humandb' +
            ' -maf ' + str(MAF)])
    


def GetOverlapWithCommonPolymoprhisms(PolymorphismFreq, Dataset='TCGA'):
    '''
        Takes in a value that determines what fequencies of polymorphisms to look at.
        Returns the number of polymorphisms with the specifided frequency that overlaps
        with variants in the 1000 Genomes project. Only available for WES data.
    '''
    AnnovarOutputDir = '/labs/ccurtis2/tilk/scripts/hri-data/Annotations/'
    #RunAnnovarToGetOverlapWithCommonVariants(Dataset, PolymorphismFreq) # Generates the Dropped and Filtered files
    Dropped = pd.read_csv(AnnovarOutputDir + Dataset + '/AnnovarOutput' + PolymorphismFreq + '.hg19_ALL.sites.2015_08_dropped',
        sep=' ', header=None)[[5,6,7,8,9]].assign(IsCommonPolymorphism = True) ### Has overlap with 1000 genomes
    Filtered = pd.read_csv(AnnovarOutputDir + Dataset + '/AnnovarOutput' + PolymorphismFreq + '.hg19_ALL.sites.2015_08_filtered',
        sep=' ', header=None)[[5,6,7,8,9]].assign(IsCommonPolymorphism = False) ### Has NO overlap with 1000 genomes
    merged = Dropped.append(Filtered)
    merged.columns = ['Tumor_Sample_Barcode','ID','Ka','Ks','xbin','IsCommonPolymorphism']
    df = merged.groupby(['xbin','IsCommonPolymorphism']).sum().reset_index()
    df[['BeginBin','EndBin']] = df['xbin'].str.split('-', expand=True)
    df = df[['BeginBin','Ka','Ks','IsCommonPolymorphism']].pivot(index='BeginBin',
        columns='IsCommonPolymorphism').reset_index().replace(np.nan, 0)
    df.columns = ['BeginBin','KaNoOverlap','KaOverlap','KsNoOverlap','KsOverlap']
    df['Fraction_Ka'] = df['KaOverlap']/df['KaNoOverlap']
    df['Fraction_Ks'] = df['KsOverlap']/df['KsNoOverlap']
    return(pd.melt(df[['BeginBin' , 'Fraction_Ka','Fraction_Ks']], id_vars=['BeginBin']))


def RecapitulateFigure(Dataset):
    '''
        Re-plot Martincorena et al 2017 (Fig 5) by calculating dN/dS in passengersa
        and drivers using dNdScv with the their binning scheme. Differences in their binning is
        - Remove samples below 600 substitutions and 3 tumor types.
        - Group "samples in 20 equal-sized bins according to mutation burden".
    '''
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType = 'kska', BinType = 'equal')
    df['xbin'] = df['xbin'].astype(str).str.strip('(]').str.replace(', ','-')
    df = df[~df['type'].isin(['CHOL', 'DLBC', 'UVM'])] # Remove cancer types not used in their analysis
    df = df[df['mutRate'] < 600]
    Sets = dict(
        All=df,
        Passengers = df[df['pancancerDriver_Bailey'] == False],
        Drivers_Bailey = df[df['pancancerDriver_Bailey'] == True]
    )
    return(GetDnDsTable(Sets, dNdScv=True))


def GetdEdIForNullCNA():
    NullCNAs = pd.concat([
        GetdEdIByMutationalBin(PermuteNullCNAs(RemoveFocal=True), ByGroup = '').assign(RemoveFocal=True),
        GetdEdIByMutationalBin(PermuteNullCNAs(RemoveFocal=False), ByGroup = '').assign(RemoveFocal=False)])
    return(NullCNAs)


def GetNullDistributionOfCorrelationWithMutBurdenForAllGenes():
    random.seed(10)
    r2 = pd.read_csv('/labs/ccurtis2/tilk/scripts/hri/Analysis/eLife_Tables/S22B_R2MutationalBinByExpressionForAllGenes')
    r2['Chaperonins'] = r2['GENE_NAME'].isin(GetGeneNamesWithinGeneFamily('Chaperonins'))
    r2['Proteasome'] = r2['GENE_NAME'].isin(GetGeneNamesWithinGeneFamily('Proteasome'))
    r2['HSP90'] = r2['GENE_NAME'].isin(GetGeneNamesWithinGeneFamily('HSP90'))
    SizeOfGeneSet = r2.groupby(['Chaperonins','Proteasome','HSP90']).size().reset_index().query(
        'Chaperonins == True or Proteasome == True or HSP90 == True')[0].sum()
    ### Randomly sample the same number of genes to calculate a null distribution of median R2 values
    NumberOfTimesToSample = 1000000
    OutList = []
    for Num in list(range(0, NumberOfTimesToSample)):
        MedianOfRandomSample = r2['0'].sample(n=SizeOfGeneSet).median()
        OutList.append(MedianOfRandomSample)
    return(pd.DataFrame(OutList))


def RemovePutativeGermlineVariants(Dataset, MAF, dNdScv=True):
    ''' Removes mutations that overlap common polymorphisms defined by @MAF and calculates dN/dS.'''
    AnnovarOutputDir = '/labs/ccurtis2/tilk/scripts/hri-data/Annotations/'
    #RunAnnovarToGetOverlapWithCommonVariants(Dataset, PolymorphismFreq) # Generates the Dropped and Filtered files
    Dropped = pd.read_csv(AnnovarOutputDir + Dataset + '/AnnovarOutput' + MAF + '.hg19_ALL.sites.2015_08_dropped',
        sep=' ', header=None).assign(IsCommonPolymorphism = True) ### Has overlap with 1000 genomes
    Filtered = pd.read_csv(AnnovarOutputDir + Dataset + '/AnnovarOutput' + MAF + '.hg19_ALL.sites.2015_08_filtered',
        sep=' ', header=None).assign(IsCommonPolymorphism = False) ### Has NO overlap with 1000 genomes
    merged = Dropped.append(Filtered).rename(columns={0: 'Chromosome', 1:'Start_Position', 2:'End_Position', 3:'Reference_Allele',
        4:'Tumor_Seq_Allele2', 5:'Tumor_Sample_Barcode', 6:'ID', 7:'Ka', 8:'Ks', 9:'xbin', 10:'pancancerDriver_Bailey'})
    df = merged[merged['IsCommonPolymorphism'] == False] # Remove putative common polymorphisms
    Sets = dict(
        All=df,
        Passengers = df[df['pancancerDriver_Bailey'] == False],
        Drivers_Bailey = df[df['pancancerDriver_Bailey'] == True]
    )
    return(GetDnDsTable(Sets, dNdScv))


