'''
	Subsets observed mutations according to specific mutation biases (i.e. Signature 1 vs Signature 30).
	Randomly removes mutations (simulating purifying selection) to specified true dN/dS values.
	Generates table used in Supplemental Figure 2 of the manuscript.
'''


#################
### Libraries ###
#################

import numpy as np
import pandas as pd
from numpy.random import choice


#################
### Functions ###
#################

def GetMutationSignatures():
	'''
		Raw dataframe of the proportion of each nucleotide context (e.g. ATA T>G) within each COSMIC signature. 
	'''
	MutationBias = pd.read_csv('/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/00_RawData/COSMICSignatures/sigProfiler_exome_SBS_signatures.csv', sep=',')
	MutationBias[['begin','end']] = MutationBias['Type'].str.lower().str.split('>', expand=True)
	MutationBias['SubType'] = MutationBias['SubType'].str.lower()
	return(MutationBias)



def ReadPermutedCountsByMutationSignature(MutationSet):
	'''
		Reads in MC3 SNP calls, that have annotations of the 3-nucleotide context of each mutation. 
		Randomly subsets 1000 mutations of each context that contain the pKa and pKs counts of the gene in which the mutation arose.
		Outputs a table of "neutral" pKa and pKs counts by each nucleotide context and the proportion of that context in each COSMIC signature. 
	'''
	NumToSample=10000
	sig = GetMutationSignatures()
	NeutralPermutedCounts = MutationSet.groupby(['begin','end','Old_DNA_context'])[['pKa','pKs']].apply(
		pd.DataFrame.sample, n=NumToSample, replace=True).reset_index().groupby(['begin','end','Old_DNA_context']).sum().reset_index()
	SigMerged = NeutralPermutedCounts.merge(sig, left_on= ['Old_DNA_context','begin', 'end'], 
		right_on=['SubType','begin','end'], indicator=True, how='left').dropna()
	return(SigMerged)


def GetMutBiasBySignature(SignatureNum, MutationSet):
	'''
		Outputs a dataframe of proportions of each trinucleotide context according to a specified specific COSMIC Signatures.
	'''
	MutationBias = ReadPermutedCountsByMutationSignature(MutationSet)
	SignatureColumn = 'SBS' + str(SignatureNum)
	MutationBias[['begin','end']] = MutationBias['Type'].str.lower().str.split('>', expand=True)
	MutationBias = MutationBias[["SubType",SignatureColumn,'begin','end', 'pKs','pKa']]
	MutationBias['SubType'] = MutationBias['SubType'].str.lower()
	return(MutationBias)


def SimulateMutations(bias, ProportionOfSites):
	'''
		Returns a dataframe of biased mutations according to a specified signature.
		Randomly subsets each mutational context to contain a specific proportion of mutations,
		(i.e. 1000 total possible mutations * bias of occuring at 2%).
		Changes observed mutations to non-synonymous or synonymous according to 
		pKA/(pKA + pKs).
	'''
	SignatureColumnName = [col for col in bias.columns if 'SBS' in col] ### i.e. SBS1 or COSMIC signature 1
	bias['probKa'] = bias['pKa']/ (bias['pKa'] + bias['pKs'])
	bias['probKs'] = 1 - bias['probKa'] 
	bias['ProportionOfSites'] = round((bias[SignatureColumnName]*ProportionOfSites))
	BiasGrouped = bias.groupby([SignatureColumnName[0], 'SubType', 'begin', 'end', 'pKa', 'pKs'])[['probKa','probKs','ProportionOfSites']].sum() ### keep index info
	KaKsCounts = BiasGrouped.apply(KaOrKs, axis=1) ### apply it to every row
	return(pd.concat([BiasGrouped, KaKsCounts], axis=1))


def KaOrKs(weights):
	'''
		Generates Ka or Ks tallies based on pKs or pKa tallies.
	'''
	#print(weights)
	output = pd.Series(dict(Ks=0, Ka=0))
	MutationClass = ['Ka','Ks']
	NumSites = int(weights['ProportionOfSites'])
	if NumSites > 0: ### if there are any sites to simulate
		for i in range(1, NumSites + 1):
			KaOrKs_Choice = choice(MutationClass, p=[weights['probKa'],weights['probKs']])
			if KaOrKs_Choice == 'Ka': ### non-syn mutation
				output['Ka'] += 1
			else: ### syn mutation
				output['Ks'] += 1
	return(output)

def SimulatePurifyingSelection(counts, dNdS_Value):
	''' 
		Randomly removes non-synonymous mutations at specified proportions (dNdS values).
	'''
	SummedCounts = counts[['Ks','Ka','pKa','pKs']].sum()
	SummedCounts['neutral_Ka'] = SummedCounts['Ka']
	NumNonSynMuts = int(SummedCounts['Ka'])
	FinalCounts = 0 
	for i in range(1, NumNonSynMuts + 1):
		count = choice([1,0], p=[dNdS_Value, 1 - dNdS_Value])
		FinalCounts += count
	SummedCounts['Ka'] = FinalCounts
	return(SummedCounts)

def AvgdNdS(df):
	'''
		Calculates the average dN/dS.
	'''
	return(df.groupby(['SignatureNumber','ProportionOfSitesWithEachContext','True_Ka_Ks']).mean())


def GetMutationalBiasResults(dNdS_Values, SignatureNum, ProportionOfSites, MutationSet):
	''' 
		Generates a final dataframe of dNdS values according to specified parameters. 
	'''
	output = pd.DataFrame()
	bias = GetMutBiasBySignature(SignatureNum, MutationSet) ### biased set of mutations according to mutation signatures
	counts = SimulateMutations(bias, ProportionOfSites).reset_index()
	counts['pKa'] = counts['pKa'] * counts['ProportionOfSites'] # Don't add permutation counts to mutations that weren't simulated
	counts['pKs'] = counts['pKs'] * counts['ProportionOfSites']
	print((counts.sum()['Ka']/counts.sum()['Ks'])/(counts.sum()['pKa']/counts.sum()['pKs'])) # Expected neutral dN/dS from random sampling
	for value in dNdS_Values:
		df = SimulatePurifyingSelection(counts, value)
		true_dNdS = value
		raw_Ka_Ks = df['Ka']/(df['Ks'])
		neutral_dNdS = (df['neutral_Ka']/df['pKa'])/(df['Ks']/df['pKs'])
		Ka_pKa_Ks_pKs = (df['Ka']/df['pKa'])/(df['Ks']/df['pKs'])
		param = pd.DataFrame({'SignatureNumber': SignatureNum, 
							'ProportionOfSitesWithEachContext' : ProportionOfSites, 
							'True_Ka_Ks' : true_dNdS, 
							'neutral_dNdS' : neutral_dNdS,
							'raw_Ka_Ks' : raw_Ka_Ks, 
							'Ka_pKa_Ks_pKs' : Ka_pKa_Ks_pKs}, 
							index=[0])
		output=output.append(param)
	return(output)



	


