'''
    This script generated 'null CNAs' by randomly permuting the start and stop position of observed CNA.
'''

import numpy as np
import pandas as pd
import random

def AnnotateFocal(patient): 
    '''
    Focal begins and ends within a chromosome arm, and non-focal begins and/or ends at a telomere and/or centromere.   
    '''
    patient['Focal'] = np.nan
    for row in np.arange(1, (len(patient) + 1)):
        chr = patient['Chromosome'][row-1:row].values[0]
        if (chr != "X") and (chr != "Y"):
            chr = patient['Chromosome'][row-1:row].astype(int).values[0]
            start = patient['Start'][row - 1:row].astype(int).values[0]
            stop = patient['End'][row-1:row].astype(int).values[0]
            lengthOfCN = patient['End'][row - 1:row].astype(int).values[0] - patient['Start'][row-1:row].astype(int).values[0]
            ChrCoord = GetChromosomalCoordinates(chr)   
            ### checks beginning of the chromosome up to centromere or end of the chromosome up to centromere   
            if ((start > ChrCoord[1]) and (stop < ChrCoord[2])) or ((start > ChrCoord[3]) and (stop < ChrCoord[0])):
                patient['Focal'][row-1:row] = True
            else:
                patient['Focal'][row-1:row] = False
        else:
            patient['Focal'][row-1:row] = False
    return (patient)


def GetChromosomalCoordinates(chr):
    '''
        Returns 4 values for each chromosome in the following order that corresponds to:
        1) End of the chromosome, not including telomere
        2) Start of the chromosome, not including telomere
        3) Position where the centrosome starts
        4) Position where the centrosome ends
    '''
    GeneSetDir="/labs/ccurtis2/tilk/scripts/hri-data/Annotations/GeneSets/"
    CentromereAndTelomere = pd.read_csv(GeneSetDir + 'CentromeresAndTelomeresUsedToMask_vCOSMIC_CNAs.csv', sep='\t')
    if chr == 23: ### X chromosome
        out = CentromereAndTelomere[CentromereAndTelomere['Chromosome'] == "X"]
    elif chr == 24: ### Y chromosome
        out = CentromereAndTelomere[CentromereAndTelomere['Chromosome'] == "Y"]
    else:
        out = CentromereAndTelomere[CentromereAndTelomere['Chromosome'] == str(chr)]
    return( out['end'].values[0], out['start'].values[0], out['cen_start'].values[0], out['cen_end'].values[0])



def GetNullCNPerPatient(df, NumPermutations, IncludeFocal):
    '''
    Generates random, "null" CN alterations per patient. Takes in a list of CN alterations and permutes its locations.
    Ensures that the permuted CN never falls between a centromere or telomere.
    '''
    numCN = len(df)
    chrOutList = []; startOutList = []; endOutList = []; barcodeList = []; gainOrLossList = []; studyList = [] ; nameList = []### initialize lists
    for row in np.arange(1, numCN + 1):
        chr = df['Chromosome'][row-1:row].astype(int).values[0]
        end = df['End'][row-1:row].astype(int).values[0] 
        start = df['Start'][row-1:row].astype(int).values[0]
        barcode = df['ID_SAMPLE'][row-1:row].values[0]
        study = df['ID_STUDY'][row-1:row].values[0]
        sample_name = df['SAMPLE_NAME'][row-1:row].values[0]
        gainOrLoss = df['MUT_TYPE'][row-1:row].values[0]
        lengthOfCN = end - start
        chromosomeValues = GetChromosomalCoordinates(chr)
        End = chromosomeValues[0]; Start = chromosomeValues[1]; CentromereStart = chromosomeValues[2]; CentromereEnd = chromosomeValues[3]
        for perm in np.arange(1, NumPermutations + 1):
            if IncludeFocal: ### no filtering of CN alterations to make sure that they don't overlap centromeres  
                if (End - lengthOfCN) < Start:
                    print ("Could not generate randomly permuted CN that doesn't overlap a telomere: " + str(chr) + " " + str(start) + " " + str(end))
                    continue
                else:
                    startOut = random.randint( Start, (End - lengthOfCN)) ### pick a copy number alteration that doesn't overlap telomere
            else: ### filter out CN alteration that overlap centromeres
                if ((CentromereStart - Start) - lengthOfCN) < 0: ### length of CN is larger than start of chr to centromere beginning
                    startOut=random.randint(CentromereEnd,(End - lengthOfCN))     ### pick a value from end of centromere to end of chr
                elif ((End - CentromereEnd) - lengthOfCN) < 0:  ### length of CN is larger than end of centromere to end of chr
                    startOut=random.randint(Start, (CentromereStart - lengthOfCN))
                elif (((End - CentromereEnd) - lengthOfCN) > 0) and (((CentromereStart - Start) - lengthOfCN) > 0):
                    bin = random.randint(1,2)
                    if bin == 1:  ### pick a number in the beginning of the chromosome up to centromere
                        startOut=random.randint(Start,(CentromereStart - lengthOfCN))
                    else: ### pick a number in the end of the chromosome up to centromere
                        startOut=random.randint(CentromereEnd,( End - lengthOfCN)) 
                else: 
                    print ("Could not generate randomly permuted CN for: " + str(chr) + " " + str(start) + " " + str(end))
                    continue
            endOut=startOut + lengthOfCN
            chrOutList.insert(0, chr); nameList.insert(0, sample_name); startOutList.insert(0, startOut); endOutList.insert(0, endOut); barcodeList.insert(0, barcode); studyList.insert(0, study); gainOrLossList.insert(0, gainOrLoss)
    DataframeOut = pd.DataFrame({'Chromosome': chrOutList,'Start': startOutList, 'End': endOutList, 'ID_SAMPLE': barcodeList, 'ID_STUDY': studyList, 'MUT_TYPE': gainOrLossList, 'SAMPLE_NAME' : nameList})
    return(DataframeOut)



def PermuteNullCNAs(RemoveFocal):
    '''
    Permutes the start and stop of each CNA, while maintaing its length to generate null CNAs where 
    selection is not expected to persist.
    @RemoveFocal = boolean that filters out CN alterations that overlap centromeres.
    '''
    CN_Dir = "/labs/ccurtis2/tilk/scripts/hri-data/Raw/CNV/"
    CNs = pd.read_csv(CN_Dir + 'CosmicCompleteCNA.tsv.gz', sep='\t')
    expanded = CNs['Chromosome:G_Start..G_Stop'].str.extract("(?P<Chromosome>[0-9]+):(?P<Start>\d+)..(?P<End>\d+)") # weird Chr:Start..Stop annotation to extract
    CNs[['Chromosome','Start','End']] = expanded
    AllPatients = CNs['ID_SAMPLE'].unique()
    FinalDF = pd.DataFrame()
    for PatientBarcode in AllPatients:
        patient = CNs.query('ID_SAMPLE == @PatientBarcode')
        if RemoveFocal:
            patient = AnnotateFocal(patient) ### annotate CNs by focal status
            patient = patient.loc[patient['Focal'] == True] ### select for only focal CN
            NullPatient = GetNullCNPerPatient(patient, NumPermutations=1, IncludeFocal = False)
        else:
            NullPatient = GetNullCNPerPatient(patient, NumPermutations=1, IncludeFocal = True)
        NullPatient = NullPatient.rename(columns={'ID_SAMPLE':'Sample','ID_STUDY':'Study'})
        NullPatient['MUT_TYPE'] = NullPatient['MUT_TYPE'].map({'LOSS':-1, 'GAIN':+1, np.nan:0}).astype(int)
        FinalDF = FinalDF.append(NullPatient.rename(columns={'MUT_TYPE':'CN'}))
    return(FinalDF)

