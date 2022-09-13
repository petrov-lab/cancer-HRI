'''
This script uses the rpy2 package to plot all main figures in R from the script PlotMainFigures.R
'''

import pandas as pd
import os
import numpy as np
from pandas.core.frame import DataFrame
from GetRCodeIntoPython import *
from dNdS import *
from bootstrap import sample
from Expression import *


def AvgdNdSInBin(df):
    """Calculate average dN/dS for randomly sampled sets of mutations."""
    return(df.groupby(level='xbin')['true'].mean())

def SortBins(df):
    ''' Sorts the dataframe so that the bins are visualized as staggered.'''
    BinOrder = pd.DataFrame({'BinID': list(range(0,9)), 'xbin' : 
        ['1-3','3-10','10-30','30-100','100-300','300-1000','1000-3000','3000-10000','10000-30000']})
    df = df.merge(BinOrder, left_on='xbin', right_on='xbin', how='left')
    return(df.sort_values(['dNdS','Group','Bin','BinID']))

def GetHighestAndLowestBins(df):
    df['LowBin'] = df.groupby(['CancerType','Group'])['Bin'].transform('min')
    df['HighBin'] = df.groupby(['CancerType','Group'])['Bin'].transform('max')
    df = df[(df['HighBin'] == df['Bin']) | (df['LowBin'] == df['Bin'])]
    df['MutationRateGroup'] = (df['Bin'] == df['LowBin']).map({True: 'Low Mutation Rate', False: 'High Mutation Rate'})
    return(df)

def GetDataForFigure2(FigNum):
    '''
        Input: @FigNum = a string of a figure panel for Fig.2 in the manuscript
        Output = Returns the raw input data used for plotting each sub-panel in Fig.2
        Commented out functions are used to process raw data and take while to run. 
    '''
    ResultsDir=os.getcwd() + '/Analysis/DataTables/'
    SetUpPlottingPackages()
    if FigNum == '2A': # dN/dS by drivers and passengers
        #GetDriversWithOncogenesAndTumorSupressors().to_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t')
        #GetDriversWithOncogenesAndTumorSupressors(dNdScv=True).to_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t')
        dnds = pd.concat([pd.read_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdScv'),
                            pd.read_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdS-permutation')])
        dnds = dnds[dnds['Group'].isin(['Passengers','Drivers_Bailey'])]
        dnds = dnds.replace(np.inf, np.nan).dropna() # Remove any noisy estimates where CI can't be inferred
        dnds = dnds[~((dnds['Group'] == 'Drivers_Bailey') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv
        dnds = SortBins(dnds)
        return(ConvertPandasDFtoR(dnds))
    elif FigNum == '2B': # Percent pathogenic via PolyPhen
        #GetPathogenicPolyPhen().to_csv(ResultsDir + 'PolyPhenDriverAndPassengersTCGA', sep='\t')
        poly = pd.read_csv(ResultsDir + 'PolyPhenDriverAndPassengersTCGA', sep='\t')
        poly[['BeginBin','EndBin']] = poly['xbin'].str.split('-', expand=True)
        poly['Group'] = poly['pancancerDriver_Bailey']
        Poly_Long = pd.concat([ 
            poly[['Group','xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
            poly[['Group','xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','EndBin']].rename(columns={'EndBin':'Bin'})
        ])
        Poly_Long = Poly_Long[~((Poly_Long['Group'] == True) & (Poly_Long['xbin'] == '1-3'))] # Remove noisy bin that reaches entire CI
        Poly = SortBins(Poly_Long.assign(dNdS=0).astype(str))
        return(ConvertPandasDFtoR(Poly))
    elif FigNum == '2C': # dE/dI
        #GetdEdIByMutationalBin(ByCancerType=False).to_csv( ResultsDir + 'dEdITCGA', sep='\t')
        dedi = pd.read_csv(ResultsDir + 'dEdITCGA', sep='\t')
        dedi = dedi[dedi['Track'].isin(['Drivers','Passengers'])]
        dedi[['BeginBin','EndBin']] = dedi['xbin'].str.split('-', expand=True)
        dEdI_Long = pd.concat([ 
            dedi[['Track','Length_Category','xbin','NumMutsInBin','NumberOfTumorsInBin','breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
            dedi[['Track','Length_Category','xbin','NumMutsInBin','NumberOfTumorsInBin','breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','EndBin']].rename(columns={'EndBin':'Bin'})])
        dEdI_Long['Group'] = dEdI_Long['Track']
        dEdI_Long['dNdS'] = dEdI_Long['Length_Category']
        dEdI_Long = SortBins(dEdI_Long)
        dEdI_Long = dEdI_Long.replace(np.inf, np.nan).dropna()
        return(ConvertPandasDFtoR(dEdI_Long))
    elif FigNum == '2D': # Clonal vs subclonal
        #GetClonalAndSubclonal('TCGA',True).to_csv(ResultsDir + 'dNdScvClonalSubclonalDriverAndPassengersTCGA', sep='\t')
        #GetClonalAndSubclonal('TCGA',False).to_csv(ResultsDir + 'dNdSClonalSubclonalDriverAndPassengersTCGA', sep='\t')
        clon = pd.concat([pd.read_csv(ResultsDir + 'dNdScvClonalSubclonalDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdScv'),
                        pd.read_csv(ResultsDir + 'dNdSClonalSubclonalDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdS-permutation')])
        clon = clon[~clon['Group'].isin(['All'])].drop(columns={'Unnamed: 0.1', 'Unnamed: 0' })
        clon['AF'] = round(clon['AF'], 1).astype(float)
        clon = clon[clon['AF'] == 0.2]
        clon = clon.replace(np.inf, np.nan).dropna() # Remove any noisy estimates where CI can't be inferred
        clon = clon[~((clon['Group'] == 'Subclonal_Drivers') & (clon['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv
        #clon = clon[~((clon['Group'] == 'Drivers') & (clon['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv
        clon = clon[~((clon['Group'] == 'Clonal_Drivers') & (clon['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv
        clon[['Bin','NumMutsInBin','NumPatientsInBin']] = clon[['Bin','NumMutsInBin','NumPatientsInBin']].astype(int)
        return(ConvertPandasDFtoR(SortBins(clon)))
    elif FigNum == '2E': # Broad cancer type
        #GetCancerTypeBroad(Dataset='TCGA', dNdScv=False).to_csv(ResultsDir + 'dNdSCancerTypeBroadDriverAndPassengersTCGA', sep='\t')
        broad = pd.read_csv(ResultsDir + 'dNdSCancerTypeBroadDriverAndPassengersTCGA', sep='\t')
        broad = broad[broad['NumMutsInBin'] > 5]
        broad = broad.replace(np.inf, np.nan).dropna() # Remove any noisy estimates where CI can't be inferred
        NumberSamples = broad.query('Group == "All"')[['CancerType','xbin','NumPatientsInBin']].drop_duplicates().groupby(
            'CancerType')['NumPatientsInBin'].sum().reset_index()
        broad = broad[broad['Group'].isin(['Passengers','Drivers'])]
        broad = GetHighestAndLowestBins(broad)
        broad[['NumPatientsInBin','NumMutsInBin']] = broad[['NumPatientsInBin','NumMutsInBin']].astype(int)
        NumberSamples['Label'] = NumberSamples['CancerType'] + ' (n=' + NumberSamples['NumPatientsInBin'].astype(int).astype(str) + ')'
        broad = broad.merge(NumberSamples[['Label','CancerType']], left_on='CancerType', right_on='CancerType')
        broad['dNdS'] = broad['Label']
        broad = SortBins(broad)
        return(ConvertPandasDFtoR(broad)) 
    elif FigNum == '2F-1': # Specific cancer type
        #GetCancerTypeSpecific(Dataset='TCGA').to_csv(ResultsDir + 'dNdSCancerTypeSpecificDriverAndPassengersTCGA', sep='\t')
        spec = pd.read_csv(ResultsDir + 'dNdSCancerTypeSpecificDriverAndPassengersTCGA', sep='\t').rename(columns={
            'Subtype' : 'CancerType' })
        spec = spec[spec['Group'].isin(['Passengers','Drivers'])]
        spec = spec.replace(np.inf, np.nan).dropna() # Remove any noisy estimates where CI can't be inferred
        spec = GetHighestAndLowestBins(spec)
        spec = spec[spec['CancerType'].isin(['BRCA','PRAD','LGG','LUAD','UCEC','HNSC'])]
        return(ConvertPandasDFtoR(spec))
    elif FigNum == '2F-2': # dE/dI by specific cancer types
        #GetdEdIByMutationalBin(ByCancerType=True).to_csv(ResultsDir + 'dEdICancerTypeSpecificDriverAndPassengersTCGA', sep='\t')
        dedi = pd.read_csv(ResultsDir + 'dEdICancerTypeSpecificDriverAndPassengersTCGA', sep='\t')
        dedi = dedi[dedi['Track'].isin(['Drivers','Passengers'])]
        dedi[['BeginBin','EndBin']] = dedi['xbin'].str.split('-', expand=True)
        dedi = dedi.replace(np.inf, np.nan).dropna() # Remove any noisy estimates where CI can't be inferred
        dEdI_Long = pd.concat([ 
            dedi[['Track','type','xbin', 'breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
            dedi[['Track','type','xbin', 'breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','EndBin']].rename(columns={'EndBin':'Bin'}) ])
        dEdI_Long = GetHighestAndLowestBins(dEdI_Long.rename(columns={'type':'CancerType','Track':'Group'}))
        dEdI_Long = dEdI_Long[dEdI_Long['CancerType'].isin(['BRCA','PRAD','LGG','LUAD','UCEC','HNSC'])]
        return(ConvertPandasDFtoR(dEdI_Long))
    elif FigNum == '2G': # Expression    
        # GetExpressionOfGeneSetsByMutationBin().to_csv(ResultsDir + 'SNVExpressionTCGA', sep='\t')
        exp = pd.read_csv(ResultsDir + 'SNVExpressionTCGA', sep='\t').rename(columns={'geneGroup':'Group'})
        exp[['BeginBin','EndBin']] = exp['xbin'].str.split('-', expand=True)
        Exp_Long = pd.concat([ 
            exp[['Group','xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
            exp[['Group','xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','EndBin']].rename(columns={'EndBin':'Bin'})])
        Exp_Long = SortBins(Exp_Long.assign(dNdS=0))
        return(ConvertPandasDFtoR(Exp_Long))



def GetDataForSupplementalFigure(FigNum):
    '''
        Input: @FigNum = a string of a figure panel for the supplemental section of the manuscript
        Output = Returns the raw input data used for plotting each sub-panel/panel
        Commented out functions are used to process raw data and take while to run. 
    '''
    ResultsDir=os.getcwd() + '/Analysis/DataTables/'
    SetUpPlottingPackages()
    if FigNum == '2A': # Correcting bias from simulations of mutations by mutational signatures
        #GetMutationalBias(NumReps=100).to_csv(ResultsDir + 'S2A_NullModelBySignatures100RepsWithRandomSamplingofPkaPksTCGA')
        bias = pd.read_csv(ResultsDir + 'S2A_NullModelBySignatures100RepsWithRandomSamplingofPkaPksTCGA')
        return(ConvertPandasDFtoR(bias))
    elif FigNum == '2B': # dN/dS subsampled from high TMB tumors in same prop as Fig.2A (dNdS-permutation)
        #GetSubsampledFromHighMutRate(Dataset='TCGA', NumReps=1000, dNdScv=False).to_csv(ResultsDir + 'S2B_dNdSSubsampledFromHighMutationRateTCGA1000Reps')
        RandSamp = pd.read_csv(ResultsDir + 'S2B_dNdSSubsampledFromHighMutationRateTCGA1000Reps').drop(
                        columns={'xbin'}).rename(columns={'type':'xbin','0':'true'})
        RandSamp = sample(RandSamp.set_index(['xbin']), AvgdNdSInBin).CI().reset_index().assign(Type='RandSamp')
        RandSamp['Bin'] = RandSamp['xbin'].str.split('-', expand=True)[0].astype(int)
        RandSamp = SortBins(RandSamp.assign(dNdS='', Group=''))
        TruedNdS = SortBins(pd.read_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t').query('Group == "Passengers"').assign(
                    Type='TruedNdS', dNdS=''))
        Combined = pd.concat([RandSamp, TruedNdS[['xbin','low','true','high','Type','Bin','dNdS','Group']]])
        return(ConvertPandasDFtoR(Combined.astype(str).reset_index()))
    elif FigNum == '2C': # dN/dS only in lowly expressed genes (dNdScv and dNdS-permutation)
        #GetdNdSInLowlyExpressedGenes(Dataset='TCGA').to_csv(ResultsDir + 'S2C_dNdSLowlyExpressedGenesTCGA')
        low = pd.read_csv(ResultsDir + 'S2C_dNdSLowlyExpressedGenesTCGA').assign(dNdS='dNdS-permutation')
        low = low.replace(np.inf, np.nan).dropna() # Remove any noisy estimates where CI can't be inferred
        low = low[low['Group'] == 'Pseudo']
        low = SortBins(low)
        return(ConvertPandasDFtoR(low))
    elif FigNum == '3': # Oncogenes and tumor suppressors 
        #GetDriversWithOncogenesAndTumorSupressors().to_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t')
        #GetDriversWithOncogenesAndTumorSupressors(dNdScv=True).to_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t') 
        drivers = pd.concat([pd.read_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t').assign(dNdS = 'dNdScv'),
                            pd.read_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t').assign(dNdS = 'dNdS-permutation')])
        drivers = drivers[drivers['Group'].isin(['Passengers', 'Drivers_Bailey', 'Oncogene', 'Tumor_Suppressor'])]
        drivers = drivers[~((drivers['Group'] == 'Drivers_Bailey') & (drivers['xbin'] == '1-3'))] # Remove same noisy bin as original data
        drivers = drivers.replace(np.inf, np.nan).dropna() 
        drivers = drivers[~((drivers['Group'] == 'Tumor_Suppressor') & (drivers['xbin'] == '3-10'))] 
        drivers = SortBins(drivers)
        return(ConvertPandasDFtoR(drivers))
    elif FigNum == '4A': # Overlap with common polymorphisms
        # RunAnnovarToGetOverlapWithCommonVariants(Dataset='TCGA', MAF='1')
        # GetOverlapWithCommonPolymoprhisms('1', Dataset='TCGA').to_csv(ResultsDir + 'S7A_1_OverlapWithCommonPolymorphismsTCGA')
        overlap = pd.read_csv(ResultsDir + 'S7A_1_OverlapWithCommonPolymorphismsTCGA')
        return(ConvertPandasDFtoR(overlap))
    elif FigNum == '4B':
        # RunAnnovarToGetOverlapWithCommonVariants(Dataset='TCGA', MAF='0.05')
        # GetOverlapWithCommonPolymoprhisms('0.05', Dataset='TCGA').to_csv(ResultsDir + 'S7B_0.05_OverlapWithCommonPolymorphismsTCGA')
        overlap = pd.read_csv(ResultsDir + 'S7B_0.05_OverlapWithCommonPolymorphismsTCGA')
        return(ConvertPandasDFtoR(overlap))
    elif FigNum == '4C':
        # RunAnnovarToGetOverlapWithCommonVariants(Dataset='TCGA', MAF='0.01')
        # GetOverlapWithCommonPolymoprhisms('0.01', Dataset='TCGA').to_csv(ResultsDir + 'S7C_0.01_OverlapWithCommonPolymorphismsTCGA')
        overlap = pd.read_csv(ResultsDir + 'S7C_0.01_OverlapWithCommonPolymorphismsTCGA')
        return(ConvertPandasDFtoR(overlap))
    elif FigNum == '4D':
        # RunAnnovarToGetOverlapWithCommonVariants(Dataset='TCGA', MAF='0.005')
        # GetOverlapWithCommonPolymoprhisms('0.005', Dataset='TCGA').to_csv(ResultsDir + 'S7D_0.005_OverlapWithCommonPolymorphismsTCGA')
        overlap = pd.read_csv(ResultsDir + 'S7D_0.005_OverlapWithCommonPolymorphismsTCGA')
        return(ConvertPandasDFtoR(overlap))
    elif FigNum == '5': # Cancer specific drivers
        #GetDriversWithOncogenesAndTumorSupressors().to_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t')
        #GetDriversWithOncogenesAndTumorSupressors(dNdScv=True).to_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t') 
        dnds = pd.concat([pd.read_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t').assign(dNdS = 'dNdScv'),
                            pd.read_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t').assign(dNdS = 'dNdS-permutation')])
        dnds = dnds[~((dnds['Group'] == 'Drivers_COSMIC') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin as original data
        dnds = dnds[~((dnds['Group'] == 'Drivers_Bailey') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin as original data
        dnds = dnds[~((dnds['Group'] == 'SpecificDriverGenes') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin as original data
        dnds = dnds[dnds['Group'].isin(['Passengers', 'Drivers_COSMIC', 'Drivers_Bailey','Drivers_Intogen','SpecificDriverGenes'])]
        dnds = SortBins(dnds)
        dnds = dnds.replace(np.inf, np.nan).dropna() 
        dnds = dnds[~((dnds['Group'] == 'Drivers_Intogen') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin as original data
        return(ConvertPandasDFtoR(dnds))
    elif FigNum == '6A': # Raw correlation between purity and mutations 
        #GetCorrelationOfTumorPurityByMutationRate(Dataset='TCGA').to_csv( ResultsDir + 'CorrelationBetweenPurityAndMutationRateTCGA')
        pur = pd.read_csv(ResultsDir + 'CorrelationBetweenPurityAndMutationRateTCGA')
        return(ConvertPandasDFtoR(pur.astype(str)))
    elif FigNum == '6B': # dN/dS by the bottom most mutatinal burden bins vs rest
        pur = pd.read_csv(ResultsDir + 'CorrelationBetweenPurityAndMutationRateTCGA')
        pur['xbin2'] = np.where( ((pur['xbin'] == '1-3') | (pur['xbin'] == '3-10')),  pur['xbin'], '10-10000')
        return(ConvertPandasDFtoR(pur.astype(str)))
    elif FigNum == '6C': # dN/dS after removing tumor by purity thresholds
        # GetdNdSAfterRemovingByTumorPurity(Dataset='TCGA', dNdScv=False).to_csv(ResultsDir + 'S9_dNdSRemovingTumorsByPurityTCGA')
        dnds = pd.read_csv(ResultsDir + 'S9_dNdSRemovingTumorsByPurityTCGA')
        dnds = SortBins(dnds.assign(dNdS=dnds['PurityThreshold']))
        dnds = dnds.replace(np.inf, np.nan).dropna() 
        dnds = dnds[~((dnds['Group'] == 'Drivers_Bailey') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin as original data
        return(ConvertPandasDFtoR(dnds))
    elif FigNum == '6D': # Num of mutations in each bin after removing tumor by purity thresholds and calculating dN/dS
        # GetdNdSAfterRemovingByTumorPurity(Dataset='TCGA', dNdScv=False).to_csv(ResultsDir + 'S9_dNdSRemovingTumorsByPurityTCGA')
        dnds = pd.read_csv(ResultsDir + 'S9_dNdSRemovingTumorsByPurityTCGA')
        dnds = SortBins(dnds.assign(dNdS=dnds['PurityThreshold']))
        dnds = dnds.replace(np.inf, np.nan).dropna() 
        dnds = dnds[~((dnds['Group'] == 'Drivers_Bailey') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin as original data
        NumberMuts = dnds[['PurityThreshold','xbin','NumMutsInBin']].drop_duplicates().groupby([
            'xbin','PurityThreshold'])['NumMutsInBin'].sum().reset_index().rename(columns={'NumMutsInBin':'SummedMutsInBins'})
        NumberMuts['Bin'] = NumberMuts['xbin'].str.split('-', expand=True)[0]
        return(ConvertPandasDFtoR(NumberMuts))
    elif FigNum == '7A': # recapitulating Martincorene et. al.
        #GetDriversWithOncogenesAndTumorSupressors().to_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t')
        #GetDriversWithOncogenesAndTumorSupressors(dNdScv=True).to_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t') 
        dnds = pd.concat([pd.read_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdScv'),
                     pd.read_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdS-permutation')])
        dnds = dnds.replace(np.inf, np.nan).dropna() # Remove any noisy estimates where CI can't be inferred
        dnds = dnds[dnds['Group'].isin(['Passengers', 'Drivers_Bailey','All'])]
        dnds = dnds[~((dnds['Group'] == 'Drivers_Bailey') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv
        dnds = SortBins(dnds)
        return(ConvertPandasDFtoR(dnds))
    elif FigNum == '7B': # recapitulating Martincorene et. al.
        #RecapitulateFigure(Dataset='TCGA').to_csv(ResultsDir + 'S10B_dNdScvBinnedSameAsFig5OfMartincorenaEtAl')
        dnds = pd.read_csv(ResultsDir + 'S10B_dNdScvBinnedSameAsFig5OfMartincorenaEtAl')
        dnds['Bin'] = round(dnds['Bin'].astype(float)).astype(int)
        dnds = dnds.replace(np.inf, np.nan).dropna() 
        dnds = dnds.sort_values(['Bin','xbin','Group'])
        return(ConvertPandasDFtoR(dnds))
    elif FigNum == '8': # Null dE/dI
        #GetdEdIForNullCNA().to_csv(ResultsDir + 'S23_COSMICNullCNAsSinglePermutationOfCNAs')
        dedi = pd.read_csv(ResultsDir + 'S23_COSMICNullCNAsSinglePermutationOfCNAs')
        dedi = dedi.replace(np.inf, np.nan).dropna() 
        dedi = dedi[dedi['xbin'] != '10000-30000'] # Don't look at bins over >10,000
        dedi = dedi[dedi['Track'].isin(['All'])]
        dedi[['BeginBin','EndBin']] = dedi['xbin'].str.split('-', expand=True)
        dEdI_Long = pd.concat([ 
            dedi[['xbin','RemoveFocal','fractionalOverlap_low', 'fractionalOverlap_mean','fractionalOverlap_high','BeginBin']].rename(
                columns={'BeginBin':'Bin','fractionalOverlap_low':'low', 'fractionalOverlap_mean':'mean', 'fractionalOverlap_high':'high'}).assign(metric='Fractional Overlap'),
            dedi[['xbin','RemoveFocal','fractionalOverlap_low', 'fractionalOverlap_mean','fractionalOverlap_high','EndBin']].rename(
                columns={'EndBin':'Bin','fractionalOverlap_low':'low', 'fractionalOverlap_mean':'mean', 'fractionalOverlap_high':'high'}).assign(metric='Fractional Overlap'),
            dedi[['xbin','RemoveFocal','breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','BeginBin']].rename(
                columns={'BeginBin':'Bin','breakpointFreq_low':'low', 'breakpointFreq_mean':'mean', 'breakpointFreq_high':'high'}).assign(metric='Breakpoint Frequency'),
            dedi[['xbin','RemoveFocal','breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','EndBin']].rename(
                columns={'EndBin':'Bin','breakpointFreq_low':'low', 'breakpointFreq_mean':'mean', 'breakpointFreq_high':'high'}).assign(metric='Breakpoint Frequency')
        ]).rename(columns={'RemoveFocal':'Group'}).replace({"Group": {True:'FocalKept', False : "FocalRemoved"} })
        dEdI_Long = SortBins(dEdI_Long.assign(dNdS=dEdI_Long['metric'])).sort_values(['dNdS','Group','BinID'])
        return(ConvertPandasDFtoR(dEdI_Long))
    elif FigNum == '9': # dE/dI using fractional overlap 
        #GetdEdIByMutationalBin(ByCancerType=False).to_csv(ResultsDir + 'dEdITCGA', sep='\t')
        dedi = pd.read_csv(ResultsDir + 'dEdITCGA', sep='\t')
        dedi = dedi[dedi['Track'].isin(['Drivers','Passengers'])]
        dedi[['BeginBin','EndBin']] = dedi['xbin'].str.split('-', expand=True)
        dEdI_Long = pd.concat([ 
            dedi[['Track','Length_Category','xbin','NumberOfTumorsInBin', 'fractionalOverlap_low', 'fractionalOverlap_mean','fractionalOverlap_high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
            dedi[['Track','Length_Category','xbin','NumberOfTumorsInBin', 'fractionalOverlap_low', 'fractionalOverlap_mean','fractionalOverlap_high','EndBin']].rename(columns={'EndBin':'Bin'})
        ])
        dEdI_Long['Group'] = dEdI_Long['Track']
        dEdI_Long['dNdS'] = dEdI_Long['Length_Category']
        dEdI_Long = SortBins(dEdI_Long)
        dEdI_Long = dEdI_Long.replace(np.inf, np.nan).dropna()
        return(ConvertPandasDFtoR(dEdI_Long))
    elif FigNum == '10': # dN/dS aftering different threshold of clonal vs subclonal
        #GetClonalAndSubclonal('TCGA',False).to_csv(ResultsDir + 'dNdSClonalSubclonalDriverAndPassengersTCGA', sep='\t') 
        dnds = pd.read_csv(ResultsDir + 'dNdSClonalSubclonalDriverAndPassengersTCGA', sep='\t')
        dnds = dnds.replace(np.inf, np.nan).dropna() # Remove any noisy estimates where CI can't be inferred
        dnds = dnds[~((dnds['Group'] == 'Subclonal_Drivers') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv
        dnds = dnds[~((dnds['Group'] == 'Drivers') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv
        dnds = dnds[~((dnds['Group'] == 'Clonal_Drivers') & (dnds['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv
        dnds = dnds.replace(np.inf, np.nan).dropna() 
        dnds = SortBins(dnds.assign(dNdS=0))
        return(ConvertPandasDFtoR(dnds))
    elif FigNum == '11': # dN/dS in functional gene sets
        #GetFunctionalGeneSets(dNdScv=False).to_csv(ResultsDir + 'dNdSFunctionalGeneSetsTCGA', sep='\t')
        dnds = pd.read_csv(ResultsDir + 'dNdSFunctionalGeneSetsTCGA', sep='\t')
        dnds = dnds.replace(np.inf, np.nan).dropna() 
        dnds = SortBins(dnds.assign(dNdS=0))
        dnds["Group"].replace({"Chr_Segregation": "Chromosome Segregtion", "Interacting_Protein": "InteractingProtein",
        'Passengers_Essential': 'Essential Genes', 'Passengers_Housekeeping' : 'Housekeeping Genes', 'TranscriptionRegulation':
        'Transcription Regulation', 'TranslationRegulation': 'Translation Regulation', 'normalGTEX_top50Percent_exclDrivers':
        'Highly-Expressed Passengers'}, inplace=True)
        dnds = GetHighestAndLowestBins(dnds.assign(CancerType=0))
        return(ConvertPandasDFtoR(dnds.astype(str)))
    elif FigNum == '12A': # dN/dS in broad cancer groups
        #GetCancerTypeBroad().to_csv(ResultsDir + 'dNdSCancerTypeBroadDriverAndPassengersTCGA', sep='\t')
        dnds = pd.read_csv(ResultsDir + 'dNdSCancerTypeBroadDriverAndPassengersTCGA', sep='\t')
        dnds = dnds[dnds['NumMutsInBin'] > 5]
        dnds = dnds[dnds['Group'].isin(['Passengers','Drivers'])]
        dnds = SortBins(dnds.assign(dNdS=dnds['CancerType']))
        dnds = dnds.replace(np.inf, np.nan).dropna() 
        return(ConvertPandasDFtoR(dnds))
    elif FigNum == '12B': # dN/dS by cancer types for all 
        #GetCancerTypeSpecific(Dataset='TCGA').to_csv(ResultsDir + 'dNdSCancerTypeSpecificDriverAndPassengersTCGA', sep='\t')
        specific = pd.read_csv(ResultsDir + 'dNdSCancerTypeSpecificDriverAndPassengersTCGA', sep='\t')
        NumbersAcrossAllBins = specific[['NumPatientsInBin','Group','Subtype','xbin']].drop_duplicates().groupby(
                    ['Group', 'Subtype'])['NumPatientsInBin'].sum().reset_index().query('Group == "All"').rename(columns={
                        'NumPatientsInBin' : 'TotalPatientsInAllBins' })
        specific = specific.merge(NumbersAcrossAllBins[['Subtype','TotalPatientsInAllBins']], left_on='Subtype', right_on='Subtype')
        specific = specific[specific['Group'].isin(['Passengers','Drivers'])]
        specific = specific[specific['true'] != 0] # Choose next bin if dN/dS can't be estimated and is 0
        specific = specific.replace(np.inf, np.nan).dropna() # Remove any noisy estimates where CI can't be inferred
        specific = GetHighestAndLowestBins(specific.rename(columns={'Subtype':'CancerType'}))
        return(ConvertPandasDFtoR(specific))
    elif FigNum == '13A': # dE/dI for specific cancer types with large # of samples (fractional overlap)
        #GetdEdIByMutationalBin(ByCancerType=True).to_csv(ResultsDir + 'dEdICancerTypeSpecificDriverAndPassengersTCGA', sep='\t')
        dedi = pd.read_csv(ResultsDir + 'dEdICancerTypeSpecificDriverAndPassengersTCGA', sep='\t')
        dedi = dedi[dedi['type'].isin(['UCEC','PRAD','LGG','HNSC','LUAD','BRCA'])]
        dedi = dedi[dedi['Track'].isin(['Drivers','Passengers'])]
        dedi = dedi.replace(np.inf, np.nan).dropna() 
        dedi[['BeginBin','EndBin']] = dedi['xbin'].str.split('-', expand=True)
        dEdI_Long = pd.concat([ 
            dedi[['Track','type','xbin', 'fractionalOverlap_low', 'fractionalOverlap_mean','fractionalOverlap_high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
            dedi[['Track','type','xbin', 'fractionalOverlap_low', 'fractionalOverlap_mean','fractionalOverlap_high','EndBin']].rename(columns={'EndBin':'Bin'})
        ])
        dEdI_Long['Group'] = dEdI_Long['Track']
        dEdI_Long['dNdS'] = dEdI_Long['type']
        dEdI_Long = SortBins(dEdI_Long)
        return(ConvertPandasDFtoR(dEdI_Long))
    elif FigNum == '13B': # dE/dI by broad cancer type groups (fractional overlap)
        #GetdEdIByMutationalBin(CNAs=ReadAndAnnotateCNAs(), ByGroup = 'CancerTypeBroad').to_csv(ResultsDir + 'dEdICancerTypeBroadriverAndPassengersTCGA', sep='\t')
        dedi = pd.read_csv(ResultsDir + 'dEdICancerTypeBroadriverAndPassengersTCGA', sep='\t')
        dedi = dedi[dedi['Track'].isin(['Drivers','Passengers'])]
        dedi = dedi.replace(np.inf, np.nan).dropna() 
        dedi[['BeginBin','EndBin']] = dedi['xbin'].str.split('-', expand=True)
        dEdI_Long = pd.concat([ 
                    dedi[['Track','broad_type','xbin','fractionalOverlap_low', 'fractionalOverlap_mean','fractionalOverlap_high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
                    dedi[['Track','broad_type','xbin','fractionalOverlap_low', 'fractionalOverlap_mean','fractionalOverlap_high','EndBin']].rename(columns={'EndBin':'Bin'})
        ])
        dEdI_Long['Group'] = dEdI_Long['Track']
        dEdI_Long['dNdS'] = dEdI_Long['broad_type']
        dEdI_Long = SortBins(dEdI_Long)
        return(ConvertPandasDFtoR(dEdI_Long))
    elif FigNum == '13C': # dE/dI for specific cancer types with large # of samples (breakpoint freq)
        #GetdEdIByMutationalBin(ByCancerType=True).to_csv(ResultsDir + 'dEdICancerTypeSpecificDriverAndPassengersTCGA', sep='\t') 
        dedi = pd.read_csv(ResultsDir + 'dEdICancerTypeSpecificDriverAndPassengersTCGA', sep='\t')
        dedi = dedi[dedi['type'].isin(['UCEC','PRAD','LGG','HNSC','LUAD','BRCA'])]
        dedi = dedi[dedi['Track'].isin(['Drivers','Passengers'])]
        dedi = dedi.replace(np.inf, np.nan).dropna() 
        dedi[['BeginBin','EndBin']] = dedi['xbin'].str.split('-', expand=True)
        dEdI_Long = pd.concat([ 
                    dedi[['Track','type','xbin','breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
                    dedi[['Track','type','xbin','breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','EndBin']].rename(columns={'EndBin':'Bin'})
        ])
        dEdI_Long['Group'] = dEdI_Long['Track']
        dEdI_Long['dNdS'] = dEdI_Long['type']
        dEdI_Long = SortBins(dEdI_Long)
        return(ConvertPandasDFtoR(dEdI_Long))
    elif FigNum == '13D': # dE/dI by broad cancer type groups (breakpoint freq)
        #GetdEdIByMutationalBin(CNAs=ReadAndAnnotateCNAs(), ByGroup = 'CancerTypeBroad').to_csv(ResultsDir + 'dEdICancerTypeBroadriverAndPassengersTCGA', sep='\t')
        dedi = pd.read_csv(ResultsDir + 'dEdICancerTypeBroadriverAndPassengersTCGA', sep='\t')
        dedi = dedi[dedi['Track'].isin(['Drivers','Passengers'])]
        dedi = dedi.replace(np.inf, np.nan).dropna() 
        dedi[['BeginBin','EndBin']] = dedi['xbin'].str.split('-', expand=True)
        dEdI_Long = pd.concat([ 
                    dedi[['Track','broad_type','xbin', 'breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
                    dedi[['Track','broad_type','xbin', 'breakpointFreq_low', 'breakpointFreq_mean','breakpointFreq_high','EndBin']].rename(columns={'EndBin':'Bin'})
        ])
        dEdI_Long['Group'] = dEdI_Long['Track']
        dEdI_Long['dNdS'] = dEdI_Long['broad_type']
        dEdI_Long = SortBins(dEdI_Long)
        return(ConvertPandasDFtoR(dEdI_Long))
    elif FigNum == '14A': # Expression by numbers of CNVs
        #GetExpressionOfGeneSetsByMutationBin('CNV').to_csv(ResultsDir + 'CNVExpressionTCGA', sep='\t') 
        exp = pd.read_csv(ResultsDir + 'CNVExpressionTCGA', sep='\t').rename(columns={'geneGroup':'Group'})
        exp = exp[exp['NumPatientsInBin'] > 5]
        exp[['BeginBin','EndBin']] = exp['xbin'].str.split('-', expand=True)
        Exp_Long = pd.concat([ 
                    exp[['Group','xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','BeginBin']].rename(columns={'BeginBin':'Bin'}),
                    exp[['Group','xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','EndBin']].rename(columns={'EndBin':'Bin'})
        ])
        Exp_Long = SortBins(Exp_Long.assign(dNdS=0))
        return(ConvertPandasDFtoR(Exp_Long))
    elif FigNum == '14B':
        #GetExpressionByLoadForAllIndividualGenes().to_csv(ResultsDir + 'S22B_R2MutationalBinByExpressionForAllGenes')
        r2 = pd.read_csv(ResultsDir + 'S22B_R2MutationalBinByExpressionForAllGenes')
        r2['Chaperonins'] = r2['GENE_NAME'].isin(GetGeneNamesWithinGeneFamily('Chaperonins'))
        r2['Proteasome'] = r2['GENE_NAME'].isin(GetGeneNamesWithinGeneFamily('Proteasome'))
        r2['HSP90'] = r2['GENE_NAME'].isin(GetGeneNamesWithinGeneFamily('HSP90'))
        r2= r2.rename(columns={'0':'r2'})
        return(ConvertPandasDFtoR(r2.astype(str)))
    elif FigNum == '14C':
        #GetNullDistributionOfCorrelationWithMutBurdenForAllGenes().to_csv(ResultsDir + 'S22C_NullDistributionOfR2ForSameSizeGeneSet', sep='\t')
        null = pd.read_csv(ResultsDir + 'S22C_NullDistributionOfR2ForSamrSizeGeneSet', sep='\t').rename(columns={'0':'null'})
        return(ConvertPandasDFtoR(null))
    elif FigNum == '14D':
        #GetExpressionOfGeneSetsByMutationBin('CNV', True).to_csv(ResultsDir + 'S22D_CNVExpressionByCancerTypeTCGA', sep='\t')
        exp = pd.read_csv(ResultsDir + 'S22D_CNVExpressionByCancerTypeTCGA', sep='\t')
        exp = exp[exp['xbin'] != '10000-30000'].rename(columns={'type':'CancerType','geneGroup':'Group'})
        exp[['BeginBin','EndBin']] = exp['xbin'].str.split('-', expand=True)
        exp = exp[['Group','CancerType','xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','BeginBin']].rename(columns={'BeginBin':'Bin'})
        exp = GetHighestAndLowestBins(exp)
        exp = exp.pivot_table(index=['Group','CancerType'], columns='MutationRateGroup', values=['true']).reset_index()
        exp.columns = ['Group','CancerType','High','Low']
        return(ConvertPandasDFtoR(exp))
    elif FigNum == '14E':
        #GetExpressionOfGeneSetsByMutationBin('SNV', True).to_csv(ResultsDir + 'S22E_SNVExpressionByCancerTypeTCGA', sep='\t')
        exp = pd.read_csv(ResultsDir + 'S22E_SNVExpressionByCancerTypeTCGA', sep='\t')
        exp = exp[exp['xbin'] != '10000-30000'].rename(columns={'type':'CancerType','geneGroup':'Group'})
        exp[['BeginBin','EndBin']] = exp['xbin'].str.split('-', expand=True)
        exp = exp[['Group','CancerType','xbin','NumMutsInBin','NumPatientsInBin', 'low', 'true','high','BeginBin']].rename(columns={'BeginBin':'Bin'})
        exp = GetHighestAndLowestBins(exp)
        exp = exp.pivot_table(index=['Group','CancerType'], columns='MutationRateGroup', values=['true']).reset_index()
        exp.columns = ['Group','CancerType','High','Low']
        return(ConvertPandasDFtoR(exp))
    elif FigNum == '21': # dN/dS of different mutation calls for TCGA
        #GetDriversWithOncogenesAndTumorSupressors(Dataset='Mutect2', dNdScv=True, MutType='kska', BinType='log').to_csv(ResultsDir + 'S3B_dNdScvMutect2')
        #GetDriversWithOncogenesAndTumorSupressors(Dataset='Mutect2', dNdScv=False, MutType='kska', BinType='log').to_csv(ResultsDir + 'S3B_dNdSMutect2')
        mutect = pd.concat([pd.read_csv(ResultsDir + 'S3B_dNdScvMutect2').assign(dNdS='dNdScv'), 
                pd.read_csv(ResultsDir + 'S3B_dNdSMutect2').assign(dNdS='dNdS-permutation')])
        #GetDriversWithOncogenesAndTumorSupressors().to_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t')
        #GetDriversWithOncogenesAndTumorSupressors(dNdScv=True).to_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t')
        mc3 = pd.concat([pd.read_csv(ResultsDir + 'dNdScvDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdScv'), 
                pd.read_csv(ResultsDir + 'dNdSDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdS-permutation')])
        df = pd.concat([SortBins(mc3).assign(Dataset='MC3 SNP Calls'), SortBins(mutect).assign(Dataset='Mutect2 SNP Calls')])
        df = df[df['Group'].isin(['Drivers_Bailey','Passengers'])]
        df = df.replace(np.inf, np.nan).dropna()
        df = df[~((df['Group'] == 'Drivers_Bailey') & (df['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv  
        return(ConvertPandasDFtoR(df))


    



def PlotReviewerResponses(Figure):
    SetUpPlottingPackages()
    ResultsDir=os.getcwd() + '/Analysis/DataTables/'
    if Figure == 'TCGA_and_ICGC':
        df = pd.concat([pd.read_csv(DataDir + 'S3C_dNdScvICGCandTCGA').assign(dNdS='dNdScv'), 
            pd.read_csv(DataDir + 'dNdSDriverAndPassengersICGCandTCGA', sep='\t').assign(dNdS='dNdS-permutation')])
        df = df[df['Group'].isin(['Drivers_Bailey','Passengers'])]
        df = df.replace(np.inf, np.nan).dropna()
        df = df[~((df['Group'] == 'Drivers_Bailey') & (df['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv
        df = SortBins(df)
        ro.r.PlotReviewerResponse_1(ConvertPandasDFtoR(df)) # ReviewerResponse_TCGAandICGC_PostFiltering.pdf
    elif Figure == 'TCGA_MutationCalls':
        mutect = pd.concat([pd.read_csv(DataDir + 'S3B_dNdScvMutect2').assign(dNdS='dNdScv'), 
                pd.read_csv(DataDir + 'S3B_dNdSMutect2').assign(dNdS='dNdS-permutation')])
        mc3 = pd.concat([pd.read_csv(DataDir + 'dNdScvDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdScv'), 
                pd.read_csv(DataDir + 'dNdSDriverAndPassengersTCGA', sep='\t').assign(dNdS='dNdS-permutation')])
        df = pd.concat([SortBins(mc3).assign(Dataset='MC3 SNP Calls'), SortBins(mutect).assign(Dataset='Mutect2 SNP Calls')])
        df = df[df['Group'].isin(['Drivers_Bailey','Passengers'])]
        df = df.replace(np.inf, np.nan).dropna()
        df = df[~((df['Group'] == 'Drivers_Bailey') & (df['xbin'] == '1-3'))] # Remove same noisy bin in dNdScv  
        ro.r.PlotReviewerResponse_2(ConvertPandasDFtoR(df)) # ReviewerResponse_TCGA_MutationCalls.pdf



def PlotFigure(FigNum):
    SetUpPlottingPackages()
    if FigNum == 'Fig2': #dN/dS by burden
        ro.r.PlotCombinedPanels(FigNum, GetDataForFigure2('2A'), GetDataForFigure2('2B'), GetDataForFigure2('2C'),
                GetDataForFigure2('2D'), GetDataForFigure2('2E'), GetDataForFigure2('2F-1'), GetDataForFigure2('2F-2'),
                GetDataForFigure2('2G'))
    elif FigNum == 'FigS2': # Correcting for mut biases
        ro.r.PlotCombinedPanels(FigNum, GetDataForSupplementalFigure(FigNum='2A'), GetDataForSupplementalFigure(FigNum='2B'),
            GetDataForSupplementalFigure(FigNum='2C'))
    elif FigNum == 'FigS3': # Oncogenes and tumor suppressors
        ro.r.PlotSupFig3(GetDataForSupplementalFigure(FigNum='3'))
    elif FigNum == 'FigS4': # Overlap with germline polymorphisms
        ro.r.PlotCombinedPanels(FigNum, GetDataForSupplementalFigure(FigNum='4A'), GetDataForSupplementalFigure(FigNum='4B'), 
                GetDataForSupplementalFigure(FigNum='4C'), GetDataForSupplementalFigure(FigNum='4D'))
    elif FigNum == 'FigS5': # Cancer specific drivers
        ro.r.PlotSupFig5(GetDataForSupplementalFigure(FigNum='5'))
    elif FigNum == 'FigS6': # dN/dS and purity
        ro.r.PlotCombinedPanels(FigNum, GetDataForSupplementalFigure(FigNum='6A'), GetDataForSupplementalFigure(FigNum='6B'), 
            GetDataForSupplementalFigure(FigNum='6C'),GetDataForSupplementalFigure(FigNum='6D') )
    elif FigNum == 'FigS7': # Figure comparison
        ro.r.PlotCombinedPanels(FigNum, GetDataForSupplementalFigure('7A'), GetDataForSupplementalFigure(FigNum='7B'))
    elif FigNum == 'FigS8': # Null dE/dI
        ro.r.PlotSupFig8(GetDataForSupplementalFigure(FigNum='8'))
    elif FigNum == 'FigS9': # Fractional overlap
        ro.r.PlotSupFig9(GetDataForSupplementalFigure(FigNum='9'))
    elif FigNum == 'FigS10': # Clonal/subclonal
        ro.r.PlotSupFig10(GetDataForSupplementalFigure(FigNum='10'))
    elif FigNum == 'FigS11': # Gene sets
        ro.r.PlotSupFig11(GetDataForSupplementalFigure(FigNum='11'))
    elif FigNum == 'FigS12': # dN/dS in cancer groups
        ro.r.PlotCombinedPanels(FigNum, GetDataForSupplementalFigure(FigNum='12A'), GetDataForSupplementalFigure(FigNum='12B'))
    elif FigNum == 'FigS13': # dE/dI in cancer groups
        ro.r.PlotCombinedPanels(FigNum, GetDataForSupplementalFigure(FigNum='13A'), GetDataForSupplementalFigure(FigNum='13B'),
                GetDataForSupplementalFigure(FigNum='13C'), GetDataForSupplementalFigure(FigNum='13D'))
    elif FigNum == 'FigS14': # Expression 
        ro.r.PlotCombinedPanels(FigNum, GetDataForSupplementalFigure(FigNum='14A'), GetDataForSupplementalFigure(FigNum='14B'),
            GetDataForSupplementalFigure(FigNum='14C'), GetDataForSupplementalFigure(FigNum='14D'), 
            GetDataForSupplementalFigure(FigNum='14E'))
    elif FigNum == 'FigS21': # dN/dS by mutation quality calls 
        ro.r.PlotSupFig21(GetDataForSupplementalFigure(FigNum='21'))
    elif FigNum == 'FigS22': # Counts of mutations in bin
        ro.r.PlotCombinedPanels(FigNum, GetDataForFigure2('2A'), GetDataForFigure2('2B'), GetDataForFigure2('2C'),
                GetDataForFigure2('2D'))
    else:
        print('Figure number not found.')
