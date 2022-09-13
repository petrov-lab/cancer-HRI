'''
    This script contains the set of functions required to run R code in python via the rpy2 package.
'''


import os
import rpy2.robjects as ro 
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


def SetUpPlottingPackages():
    '''Uses rpy2 package to be able to run R code in python. Code in R is used for plotting and statistics.'''
    ggplot = importr('ggplot2')
    scales = importr('scales')
    ggpmisc = importr('ggpmisc')
    cowplot = importr('cowplot')
    data_table = importr('data.table')
    ro.r.source('/labs/ccurtis2/tilk/scripts/hri/Analysis/AllPlottingFunctions.R') # Source R script for plotting

def ConvertPandasDFtoR(df):
    with localconverter(ro.default_converter + pandas2ri.converter): 
        dfInR = ro.conversion.py2rpy(df) # Convert pandas dataframe to R 
    return(dfInR)

def ConvertRDataframetoPandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        dfInPd = ro.conversion.rpy2py(df) # Convert R dataframe back to Pandas
    return(dfInPd)