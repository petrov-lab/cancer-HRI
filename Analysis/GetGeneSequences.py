'''
	This script contains the functions used to read in the DNA sequences of every gene mapped
	to COSMIC, GRCh38 and GRCh37 reference genomes.
'''

import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

 
def GetAllGenesMappedToGRCH38(filename=os.getcwd() + "/GeneSequences/Homo_sapiens.GRCh38.cds.all.fa.gz"):
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

def GetAllGenesMappedToGRCH37(filename=os.getcwd() + "/GeneSequences/Homo_sapiens.GRCh37.cds.all.fa.gz"):
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



def GetAllGenesMappedToCOSMIC(filename=os.getcwd() + "/GeneSequences/All_COSMIC_Genes.fasta.gz"):
	'''
		Returns pandas.Series containing DNA sequence of every gene used in the COSMIC database. 
		Contains dummy column ENSG.
	'''
	S = pd.read_table(filename, header=None, lineterminator='>')[0]
	df = S.str.partition('\n')
	seqs = df[2].str.replace('\n', '')
	df[0] += [" None"]
	S2 = pd.Series(seqs.values, index=pd.MultiIndex.from_tuples(df[0].str.split(), names=['Hugo', 'ENST','ENSG']), name='Sequence').sort_index()
	S2.index = S2.index.set_levels(S2.index.levels[1].str.replace(r'_ENST.*$', ''), level=1)
	S2 = S2.reorder_levels(['ENST','ENSG','Hugo']).sort_index()
	return(S2)


