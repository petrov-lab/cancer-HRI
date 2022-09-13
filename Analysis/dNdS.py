'''
    This contains the set of functions used in the manuscript to calculate dN/dS for the main figures.
'''

#################
### Libraries ###
#################

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import os 
import glob
import pandas as pd
import numpy as np
from bootstrap import sample
import rpy2.robjects as ro 
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


#################
### Functions ###
#################


def GetPermutedCounts(Dataset='TCGA'):
    '''
        Takes a list of paths to concatenate together into one large table of Ka/Ks tallies.
        Five possible datasets saved are: 'ICCGandTCGA', 'TCGA' and 'Mutect2'.
        Returns an pandas dataframe of annotated mutations with permutation counts.
    '''
    List = []
    CountsDir="/labs/ccurtis2/tilk/scripts/hri-data/Processed/SNV/ProcessedWithPermutations"
    if Dataset == 'TCGA':
        ListOfFiles = glob.glob(os.path.join(CountsDir + '/TCGA/', "*ENST.counts.gz"))
    elif Dataset == 'ICCGandTCGA':
        ListOfFiles = list(glob.glob(os.path.join(CountsDir + '/ICGC/', "*ENST.counts.gz")) + 
                        glob.glob(os.path.join(CountsDir + '/TCGA/', "*ENST.counts.gz")))
    elif Dataset == 'Mutect2':
        ListOfFiles = glob.glob(os.path.join(CountsDir + '/TCGA/', "*ENST.counts.gz"))
    ### Merge datasets by common columns
    for FileName in ListOfFiles:
        df = pd.read_csv(FileName, sep='\t')
        List.append(df)
    Tallies = pd.concat(List).fillna(0)
    if Dataset == 'ICCGandTCGA':
        Tallies = Tallies[GetCommonColumnsInDataset(Dataset)]
    return(Tallies)


def GetCommonColumnsInDataset(Dataset):
    '''
    Common columns of interest that datasets (TCGA/ICGC) are merged on.
    '''
    PolyPhenColumn = ['PolyPhen'] # PolyPhen only available for TCGA
    PositionColumns = ['Chromosome','Start_Position','End_Position','Reference_Allele',
                        'Tumor_Seq_Allele2','Strand','Variant_Classification']
    GeneColumns = ['Gene','Transcript_ID','Oncogene','Tumor_Suppressor','Housekeeping', 'Essential',
                    'pancancerDriver_COSMIC','pancancerDriver_Bailey', 'pancancerDriver_Intogen', 
                    'normalGTEX_bottom50Percent','normalGTEX_top50Percent','All_Interacting_Proteins','Hugo_Symbol',
                    'TranslationRegulation','TranscriptionRegulation', 'SpecificDriverGenes','Chr_Segregation']
    KaKsColumns = ['Old_DNA_context','New_DNA_context','Loc','begin','end','HGVSc','Ka', 'Ks', 'pKa', 'pKs']
    MutColumns = [ 'type','Tumor_Sample_Barcode','totalNumIntronMutations', 'totalNumMutations', 'Diagnosis_Age']
    QCColumns = ['t_depth', 't_alt_count']
    if 'ICGC' in Dataset :
        return( PositionColumns + GeneColumns + KaKsColumns + MutColumns + QCColumns + ['totalIntergenic'])
    elif 'TCGA' in Dataset:
        return(PolyPhenColumn  + PositionColumns + GeneColumns + KaKsColumns + MutColumns + QCColumns)
    else: # COSMIC/Validation/Mutect2
        return(GeneColumns + KaKsColumns + ['Tumor_Sample_Barcode'])


def AnnotateMutationRate(df, MutType, BinType):
    '''
    Annotates patients by the mutation rate of tumors and bins them based on these values.
    '''
    if MutType == 'intron':
        mutRate = df[['Tumor_Sample_Barcode','totalNumIntronMutations']].drop_duplicates().reset_index().rename(columns={'totalNumIntronMutations':'mutRate'})
    elif MutType == 'intergenic': 
        mutRate = df[['Tumor_Sample_Barcode','totalIntergenic']].drop_duplicates().reset_index().rename(columns={'totalIntergenic':'mutRate'})
    elif MutType == 'kska':
        mutRate = df.groupby('Tumor_Sample_Barcode')[['Ka','Ks']].sum().sum(axis=1).reset_index().rename(columns={0:'mutRate'})
    elif MutType == 'all':
        mutRate = df[['Tumor_Sample_Barcode','totalNumMutations']].drop_duplicates().reset_index().rename(columns={'totalNumMutations':'mutRate'})
    elif MutType == 'age':
        mutRate = df[['Tumor_Sample_Barcode','Diagnosis_Age']].drop_duplicates().reset_index().rename(columns={'Diagnosis_Age':'mutRate'})
    if BinType == 'log':
        mutRate['xbin'] = pd.cut(mutRate['mutRate'], bins=[0,3,10,30,100,300,1000,3000,10000,30000], 
                        labels=['1-3','3-10','10-30','30-100','100-300','300-1000','1000-3000','3000-10000','10000-30000'])
    elif BinType == 'age':
        mutRate['xbin'] = pd.cut(mutRate['mutRate'], bins=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100], 
                        labels=['10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100'])
    elif BinType == 'equal':
        mutRate['xbin'] = pd.qcut(mutRate['mutRate'], 20, retbins=False) ### used for recapitulating Fig 5 in Martincorena et al 2017
    df = df.merge(mutRate , left_on = 'Tumor_Sample_Barcode', right_on = 'Tumor_Sample_Barcode')
    return(df)


def KaKs(df):
    """Calculate Ka/Ks for every bin of samples."""
    return df.groupby(level='xbin').sum().eval('(Ka/pKa)/(Ks/pKs)')


def PercentPathogenicPolyPhen(df):
    ''' Calculate the fraction of pathogenic mutations for every bin of samples.'''
    return(df.groupby(['xbin','pancancerDriver_Bailey']).sum().eval('(possibly_damaging + probably_damaging)/(benign + possibly_damaging + probably_damaging)'))


def BootstrapCI(df, NumBoot=10000):
    ''' 
        Calculates 95% confidence intervals by boostrapping mutations in each mutational load bin (i.e. by the column `xbin`).
        Outputs the true, lower and upper bounded confidence intervals for dN/dS (aka (Ka/pKa)/(Ks/pKs))
    '''
    RandomlySampledSet = [KaKs(df[['xbin','Ka','Ks','pKa','pKs']].groupby(['xbin']).apply(lambda x: x.sample(frac=1, replace=True))) for _ in range (NumBoot)]
    RandomlySampledSet = pd.concat(RandomlySampledSet).replace([np.inf, -np.inf], np.nan).dropna() # Remove infinities and nas
    mean = RandomlySampledSet.groupby(['xbin']).mean()
    std = RandomlySampledSet.groupby(['xbin']).std()
    z = 1.9
    err = z*std
    return(pd.DataFrame({'high' : (mean + err), 'true' :KaKs(df[['xbin','Ka','Ks','pKa','pKs']].set_index('xbin')), 'low' : mean - err}))

 
def GetdNdScv(df, RefGenome='hg19'):
    '''
    Uses the rpy2 package to run dNdScv in R. Takes an input of a pandas dataframe,
    converts it to a dataframe in R and returns a pandas dataframe.
    '''
    dndscv= importr('dndscv') # Load required packages to run dNdScv
    plyr= importr('plyr')
    ro.r.source(os.getcwd() + '/dNdScv.R') # Source R script that runs dNdScv
    df = df[['Tumor_Sample_Barcode','Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','xbin']].astype(str)
    with localconverter(ro.default_converter + pandas2ri.converter): 
        dfInR = ro.conversion.py2rpy(df) # Convert pandas dataframe to R 
    if RefGenome == 'hg19':
        dNdScv=ro.r.ApplydNdScvToHg19(dfInR) # Run script
    elif RefGenome == 'hg38':
        dNdScv=ro.r.ApplydNdScvToHg38(dfInR) # Run script
    with localconverter(ro.default_converter + pandas2ri.converter):
        dfInPd = ro.conversion.rpy2py(dNdScv) # Convert R dataframe back to Pandas
    return(dfInPd)



def GetDnDsTable(Sets, dNdScv=False, RefGenome='hg19'):
    Results = pd.DataFrame()
    GroupNames=[k  for  k in  Sets.keys()]
    for Group in GroupNames:
        if len(Sets[Group]) > 0: # Make sure mutations exist in each set
            MetaData = pd.DataFrame({"NumMutsInBin" : Sets[Group].groupby(['xbin']).size(),
                                    "NumPatientsInBin" : Sets[Group].groupby(['xbin'])['Tumor_Sample_Barcode'].unique().apply(len)})
            if dNdScv:
                if RefGenome == 'hg19':
                    CI = GetdNdScv(Sets[Group]).set_index('xbin').rename(columns={'mle': 'true', 'cilow': 'low', 'cihigh':'high'})
                elif RefGenome == 'hg38':
                    CI = GetdNdScv(Sets[Group],'hg38').set_index('xbin').rename(columns={'mle': 'true', 'cilow': 'low', 'cihigh':'high'})
            else:
                CI = BootstrapCI(Sets[Group])
                #CI=sample(Sets[Group].set_index('xbin'), KaKs).CI() # Calculates dN/dS for the group set and bootstraps to get CIs
            dNdS = pd.concat([MetaData, CI], axis=1).reset_index()
            dNdS[['BeginBin','EndBin']] = dNdS['xbin'].str.split('-', expand=True)
            dNdS_Long = pd.concat([ 
                dNdS[['xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
                dNdS[['xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','EndBin']].rename(columns={'EndBin':'Bin'})
            ]).assign(Group = Group)
            Results = Results.append(dNdS_Long)
    return(Results)


def GetDriversWithOncogenesAndTumorSupressors(Dataset, dNdScv=False, MutType='kska', BinType='log'):
    '''
        Calculates dN/dS for pan-cancer drivers, oncogenes, and tumor supressors.
    '''
    df = GetPermutedCounts(Dataset)
    df = AnnotateMutationRate(df, MutType, BinType)
    Sets = dict(
        All=df,
        Passengers = df[df['pancancerDriver_Bailey'] == False],
        Drivers_COSMIC = df[df['pancancerDriver_COSMIC'] == True],
        Drivers_Bailey = df[df['pancancerDriver_Bailey'] == True],
        Drivers_Intogen = df[df['pancancerDriver_Intogen'] == True],
        Oncogene = df[df['Oncogene'] == True],
        Tumor_Suppressor = df[df['Tumor_Suppressor'] == True],
        SpecificDriverGenes = df[df['SpecificDriverGenes'] == True]
    )
    return(GetDnDsTable(Sets, dNdScv))


def GetPathogenicPolyPhen(MutType='kska', BinType='log'):
    '''
        Calculates fraction of pathogenic mutations for drivers and passengers.
    '''
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType, BinType)
    PolyPhen = df.groupby(['xbin','pancancerDriver_Bailey','PolyPhen','Tumor_Sample_Barcode']).size().unstack(
                level='PolyPhen').replace(np.nan,0).reset_index()
    MetaData = pd.DataFrame(
                {"NumMutsInBin" : PolyPhen.groupby(['xbin','pancancerDriver_Bailey']).sum(
                    )[['benign','possibly_damaging','probably_damaging']].sum(axis=1),
                "NumPatientsInBin" : PolyPhen.groupby(['xbin','pancancerDriver_Bailey','Tumor_Sample_Barcode']).sum(
                    ).sum(axis=1).where(lambda x : x!=0).dropna().groupby(['xbin','pancancerDriver_Bailey']).size()
    }).reset_index()
    CI = sample(PolyPhen.set_index(['xbin','pancancerDriver_Bailey']), PercentPathogenicPolyPhen).CI().reset_index()   
    return(MetaData.merge(CI, left_on=['xbin','pancancerDriver_Bailey'], right_on=['xbin','pancancerDriver_Bailey']))


def GetClonalAndSubclonal(dNdScv=False, MutType='kska', BinType='log'):
    '''
        Calculates dN/dS for clonal and subclonal pan-cancer drivers.
    '''
    OutList = [] # Empty DF where results are appended to
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType, BinType)
    RangeOfAFs = np.arange(0.1, 1, 0.1)
    for AF in RangeOfAFs:
        Sets = dict(
            Passengers = df[df['pancancerDriver_Bailey'] == False],
            Drivers = df[df['pancancerDriver_Bailey'] == True],
            Clonal_Drivers = df[(df['pancancerDriver_Bailey'] == True) & (df['PurityAdjusted_VAF'] > AF)],
            Subclonal_Drivers = df[(df['pancancerDriver_Bailey'] == True) & (df['PurityAdjusted_VAF'] <= AF)],
            Clonal_Passengers = df[(df['pancancerDriver_Bailey'] == False) & (df['PurityAdjusted_VAF'] > AF)],
            Subclonal_Passengers = df[(df['pancancerDriver_Bailey'] == False) & (df['PurityAdjusted_VAF'] <= AF)]
        )
        OutList.append(GetDnDsTable(Sets, dNdScv).assign(AF=AF))
    return(pd.concat(OutList))



def GetCancerTypeBroad(dNdScv=False, MutType='kska', BinType='log'):
    '''
        Calculates dN/dS by broad cancer categories, e.g. 'Neuronal' or 'Circulatory'.
    '''
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
    OutList = [] # Empty DF where results are appended to
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType, BinType)
    for Group in CancerGroups.unique():
        MutsInGroup = df[df['type'].isin(CancerGroups[CancerGroups == Group].index.tolist())]
        Sets = dict(
            All=MutsInGroup,
            Passengers = MutsInGroup[MutsInGroup['pancancerDriver_Bailey'] == False],
            Drivers = MutsInGroup[MutsInGroup['pancancerDriver_Bailey'] == True]
        )
        OutList.append(GetDnDsTable(Sets, dNdScv).assign(CancerType=Group))
    return(pd.concat(OutList))


def GetCancerTypeSpecific(dNdScv=False, MutType='kska', BinType='log'):
    '''
        Calculates dN/dS by specific cancer-subtypes, e.g. 'BRCA' or 'LAML'.
    '''
    OutList = [] # Empty DF where results are appended to
    df = GetPermutedCounts()
    df = AnnotateMutationRate(df, MutType, BinType)
    for Subtype in df['type'].unique():
        MutsInGroup = df[df['type'] == Subtype]
        Sets = dict(
            All=MutsInGroup,
            Passengers = MutsInGroup[MutsInGroup['pancancerDriver_Bailey'] == False],
            Drivers = MutsInGroup[MutsInGroup['pancancerDriver_Bailey'] == True]
        )
        OutList.append(GetDnDsTable(Sets, dNdScv).assign(Subtype=Subtype))
    return(pd.concat(OutList))


    

