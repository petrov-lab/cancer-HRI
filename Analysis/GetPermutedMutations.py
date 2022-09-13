'''
    This script contains the set of functions used to generate permuted tallies forsynonymous and non-synonymous mutations
    based on 3-nucleotide sequence context for every observed mutation in TCGA to calculate dNdS-perm in the manuscript. 
'''

#################
### Libraries ###
#################

import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import pandas as pd
import os 
import glob
import re
from GetGeneSequences import *
from bootstrap import sample


#################
### Functions ###
#################

def ReadAnnotatedMutations(CancerType):
    ''' Reads in annnotated mutations from @GetAnnotatedMutations.py'''
    ProcessedDir="/labs/ccurtis2/tilk/scripts/hri-data/Processed/SNV/ProcessedWithGeneAnnotations/" # Input directory
    muts = pd.read_csv(ProcessedDir + "TCGA/" + CancerType + ".txt.gz", sep='\t') # Mutations from @AnnotateMC3Data()
    return(muts)


def GetGeneticCode():
    '''
        Returns a dictionary of codons to one letter amino acids.
    '''
    genetic_code = dict(
    ttt='F',tct='S',tat='Y',tgt='C',
    ttc='F',tcc='S',tac='Y',tgc='C',
    tta='L',tca='S',taa='*',tga='*',
    ttg='L',tcg='S',tag='*',tgg='W',
    ctt='L',cct='P',cat='H',cgt='R',
    ctc='L',ccc='P',cac='H',cgc='R',
    cta='L',cca='P',caa='Q',cga='R',
    ctg='L',ccg='P',cag='Q',cgg='R',
    att='I',act='T',aat='N',agt='S',
    atc='I',acc='T',aac='N',agc='S',
    ata='I',aca='T',aaa='K',aga='R',
    atg='M',acg='T',aag='K',agg='R',
    gtt='V',gct='A',gat='D',ggt='G',
    gtc='V',gcc='A',gac='D',ggc='G',
    gta='V',gca='A',gaa='E',gga='G',
    gtg='V',gcg='A',gag='E',ggg='G')
    return(genetic_code)


def AnnotateMutationsAsKaOrKs(muts):
    '''
        Adds annotations of 0 or 1 when a mutation is non-synonymous (Ka) or synonymous (Ks).
        Returns protein-coding mutations.
    '''
    MutationClasses = dict(Ka={'Missense_Mutation', 'Nonsense_Mutation'},  Ks={'Silent'})
    KaOrKs = {k:muts['Variant_Classification'].isin(frozenset(v)) for k, v in MutationClasses.items()}
    muts = pd.concat([muts, pd.DataFrame(KaOrKs) * 1], axis=1)
	muts = muts[~((muts['Ka'] == 0) & (muts['Ks'] == 0))] ### select only protein coding mutations
    return(muts)


def ParseCodingSequence(CodingSequenceColumn):
    '''
        Takes in a pandas Series of HGVSc annotations, e.g. "c.1270G>A", a G to A mutation
        at position 1270 in the open reading frame. Splits column into 'begin','end',and 'Loc'
        of the mutation.
    '''
    df = CodingSequenceColumn.str.extract('(?P<begin>[ACGT]+)?>(?P<end>[ACGT]+)?', expand=True)
    df['begin'] = df['begin'].str.lower()
    df['end'] = df['end'].str.lower()
    df['Loc'] = CodingSequenceColumn.str.extract('(\d+)').fillna(0).astype(int)
    df['SNV_or_MNV'] = CodingSequenceColumn.str.contains('\+|-|dup|del|ins|delins|inv').map(pd.Series({False:'SNV', True:'MNV'}))
    return(df)


def GetCodingSequenceOfSingleGeneMapped(ENST, ENSG, HUGO, all_genes, GeneType):
    '''
        @all_genes is a series with indicies that correspond to differente gene identifiers
        (ENST,ENSG, and Hugo Symbol), and values that are the coding sequence of the gene.
        This function iterates through different gene identifiers that are provided and
        returns the matching coding sequence of the gene. Only required gene identifier is
        the parameter (), the Hugo Symbol. All others gene identifiers help further map genes
        to coding sequences. 
    '''
    if GeneType == 'ENST':
        try: 
            return all_genes.loc[ENST,:,:].values[0].lower()
        except:
            return ''
    elif GeneType == 'ENSG':
        try: 
            return all_genes.loc[:,ENSG,:].values[0].lower()
        except:
            return ''
    elif GeneType == 'HUGO':
        try: 
            return all_genes.loc[:,:,HUGO].values[0].lower()
        except:
            return ''

def GetCodingSequence(row, RefGenome, AllMappedGenes, GeneType='ENST'):
    ENST = row.loc['Transcript_ID'] if 'Transcript_ID' in (row.index) else ""
    ENSG = row.loc['Gene'] if 'Gene' in (row.index) else ""
    HUGO = row.loc['Hugo_Symbol'] if 'Hugo_Symbol' in (row.index) else ""
    seq = GetCodingSequenceOfSingleGeneMapped(ENST, ENSG, HUGO, AllMappedGenes, GeneType)
    return(seq.replace(" ",""))



def PermuteMutation(row, RefGenome, AllMappedGenes, GeneType='ENST'):
    '''
        Permutes an observed mutation in all ways possible that still maintain the same gene & 3-nt context.
        Returns a pandas.Series with 3 entries: 
        pKa ~ # of permutations that lead to a non-synonymous change
        pKs ~ # of permutations  "        "        synonymous change
        pKn ~ # of permutations  "        "        nonsense mutation
    '''
    output = pd.Series(dict(pKs=0, pKa=0, pKn=0))
    genetic_code=GetGeneticCode()
    seq = GetCodingSequence(row, RefGenome, AllMappedGenes, GeneType)
    Loc = row['Loc']
    if type(Loc) == float or np.isnan(Loc):
        # Location is NaN
        return output
    Loc = int(row['Loc'])
    if seq is None or Loc + 1 >= len(seq):
        # If we can't find the DNA sequence or if the mutation's location is at the end of the gene return nothing
        return output
    DNA_context = seq[Loc - 2: Loc + 1]
    if 'n' in DNA_context:
        # Can't permute trinucleotide with N in context
        return output
    output['Old_DNA_context'] = DNA_context
    new_context = DNA_context[0] + row['end'] + DNA_context[2] 
    output['New_DNA_context'] = new_context
    new_seq = seq.replace(DNA_context, new_context)
    for mut in re.finditer(DNA_context, seq):
        position = mut.start() + 1
        codon_start = position  - (position % 3)
        old_codon = seq[codon_start:codon_start+3]
        new_codon = new_seq[codon_start:codon_start+3]
        if old_codon in genetic_code and new_codon in genetic_code:
            old_AA = genetic_code[old_codon]
            new_AA = genetic_code[new_codon]
            if old_AA == new_AA:
                output['pKs'] += 1
            else:
                output['pKa'] += 1
                if new_AA == '*':
                    output['pKn'] += 1
    return output



def GetPermutationCounts(CancerType, GeneType='ENST'):
    '''
    Generates a new set of files that contain pKa/pKs columns for each mutationn.
    @CancerType = a string used for reading in a specific file within each dataset
    TCGA data is split by cancer type, e.g. @CancerType would be "BRCA" for TCGA
    WGS data is split by patient, e.g. @CancerType would be "DO50409" for ICGC
    '''
    print(CancerType)
    OutDir="/labs/ccurtis2/tilk/scripts/hri-data/Processed/SNV/ProcessedWithPermutations/" 
    muts = ReadAnnotatedMutations(CancerType)
    muts = AnnotateMutationsAsKaOrKs(muts) ### will also only select for non-synonymous and synonymous mutations
    muts = pd.concat([muts, ParseCodingSequence(muts['HGVSc']) ], axis=1)
    muts = muts[(muts['SNV_or_MNV'] == "SNV") & (muts['Loc'] > 2) ] ### remove MNV variants or mutations that land the very beginnin of CDS
    permutations = muts.apply(PermuteMutation, args=("GRCH37", GetAllGenesMappedToGRCH37(), GeneType), axis=1)
    out = pd.concat([muts, permutations], axis=1)
    out = out[~out['Old_DNA_context'].isnull()] ### remove permutations where mutation landed at end of gene
    out.to_csv(OutDir + '/TCGA/' + CancerType + "." + GeneType + ".counts.gz", sep='\t', compression="gzip")  # Save file


GetPermutationCounts(CancerType = sys.argv[1], GeneType=sys.argv[2])