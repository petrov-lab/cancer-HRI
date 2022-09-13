
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

 
def GetAllGenesMappedToGRCH38(filename=os.getcwd() + "/annotations/Homo_sapiens.GRCh38.cds.all.fa.gz"):
	''' 
		Returns pandas.Series containing DNA sequence of every gene mapped to reference genome GRCH38.
		The datasets that are mapped to GRCh38 are SNPS called by mutect2.
	'''
	S = pd.read_table(filename, header=None, lineterminator='>')[0]
	df = S.str.partition('\n')
	seqs = df[2].str.replace('\n', '')
	S2 = pd.Series(seqs.values, 
	index=pd.MultiIndex.from_tuples(df[0].str.replace('.', ' ').str.split(), names=['ENST','one','ENSG','one2','Hugo']), name='Sequence').sort_index()
	S2.index = S2.index.droplevel(['one','one2']) 
	return(S2)

def GetAllGenesMappedToGRCH37(filename=os.getcwd() + "/annotations/Homo_sapiens.GRCh37.cds.all.fa.gz"):
	'''
		Returns pandas.Series containing DNA sequence of every gene mapped to reference genome GRCH37.
		The datasets that are mapped to GRCh38 are SNPS from MC3 dataset.
	'''
	S = pd.read_table(filename, header=None, lineterminator='>')[0]
	df = S.str.partition('\n')
	seqs = df[2].str.replace('\n', '').str.strip()
	S2 = pd.Series(seqs.values, 
	index=pd.MultiIndex.from_tuples(df[0].str.split(), names=['ENST','ENSG','Hugo']), name='Sequence').sort_index()
	return(S2)


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
