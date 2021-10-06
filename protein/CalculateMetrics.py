
import pandas as pd
import os
import numpy as np
from pandas.core.frame import DataFrame
import rpy2.robjects as ro 
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from GetData import * 

def SetUpPlottingPackages():
    ggplot = importr('ggplot2')
    #scales = importr('scales')
    #ggpmisc = importr('ggpmisc')
    #cowplot = importr('cowplot')
    lme4 = importr('lme4')
    glmnet = importr('glmnet')
    data_table = importr('data.table')
    ro.r.source('/labs/ccurtis2/tilk/scripts/protein/GetRegressionStats.R') # Source R script for plotting

def ConvertPandasDFtoR(df):
    with localconverter(ro.default_converter + pandas2ri.converter): 
        dfInR = ro.conversion.py2rpy(df) # Convert pandas dataframe to R 
    return(dfInR)

def ConvertRDataframetoPandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        dfInPd = ro.conversion.rpy2py(df) # Convert R dataframe back to Pandas
    return(dfInPd)

def GetRegressionStatsInput(Dataset, DataType):
    ''' 
    Calculates mixed effects regression (controlling for cancer type) coefficients for
    the effect of expression/drug response/rnai response to mutational load.
    @Dataset = `CCLE`, `GTEX` or `TCGA`
    @DataType = `Expression`, `RNAi` or `Drug`
    '''
    if DataType == 'RNAi':
        Values = GetRNAiResponseData()
    elif DataType == 'Drug':
        Values = GetDrugResponseData(Screen='primary')
    elif DataType == 'Expression':
        Values = GetExpressionData(Dataset)
    Muts = AnnotateMutationalLoad(GetMutations(Dataset), MutType='SNV')
    Purity = GetTumorPurity()
    Tissue = GetTissueType(Dataset)
    Stats = Values.merge(Muts, left_on='Barcode', right_on='Barcode')
    Stats['LogScore'] = np.log10(Stats['MutLoad'])
    if Dataset == 'TCGA':
        Covariates = Tissue.merge(Purity, left_on='Barcode', right_on='Barcode')
    else:
        Covariates = Tissue
    return(Stats.merge(Covariates, left_on='Barcode', right_on='Barcode'))

def GetRegressionEstimates(Dataset, DataType):
    SetUpPlottingPackages()
    df = GetRegressionStatsInput(Dataset, DataType)
    R = ConvertPandasDFtoR(df)
    return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'Mixed-Effect')))


GetRegressionEstimates(Dataset='TCGA', DataType='Expression')