import pandas as pd
import numpy as np
import os

# Key attributes of nomenclatures
prefixes = dict(ensembl_gene_id='ENSG', 
                refseq_accession='NM_', 
                ccds_id='CCDS',
                ucsd_id='uc',
                vega_id='OTTHUMG')
numeric_ids = {'hgnc_id', 'entrez_id', 'mgd_id', 'rgd_id', 'omim_id'}
symbolic = {'symbol', 'alias_symbol', 'prev_symbol', 'uniprot_ids', 'ena', 'cosmic'}
non_unique = {'refseq_accession', 'pubmed_id', 'ena', 'alias_symbol', 'prev_symbol'}
#########


def get_hgnc_database():
    return(pd.read_csv('/labs/ccurtis2/tilk/scripts/hri-data/Annotations/GeneSets/hgnc_complete_set.txt', sep='\t', low_memory=False))

def invert_non_unique_map(S):
    out = pd.Series({v_i:k for k, v in S.dropna().str.split('|').items() for v_i in v}, name=S.index.names[0])
    out.index.names = [S.name]
    return out

def drop_numeric_prefixes(hgnc_df):
    for col in numeric_ids - {'entrez_id'}:
        hgnc_df.loc[:, col] = hgnc_df[col].str.rpartition(':')[2].values
    return hgnc_df

backups = dict( symbol=['prev_symbol', 'alias_symbol', 'uniprot_ids'],
                cosmic=['symbol', 'prev_symbol', 'alias_symbol', 'uniprot_ids'],
                ensembl_gene_id=['refseq_accession'])

def create_gene_map(genes_in, output_type, 
    input_type='infer', hgnc_database=None, match_threshold=0.95,
    symbolic_to_infer=symbolic-non_unique,
    backups=backups, no_output='None', case_sensitive=True, verbose=True):
    
#    if not verbose:
#        print = lambda _: None

    df = hgnc_df.copy() if hgnc_database is None else hgnc_database

    if type(genes_in) != pd.Series:
        genes_in = pd.Series(genes_in)
    assert output_type in df.columns, "{:} is not a nomenclature within the hgnc database.".format(output_type)  # Break early
#    assert output_type not in non_unique, """{:} is a nomenclature where several gene names map to the same HGNC ID. 
#    Therefore, I cannot provide a unique mapping to it.""".format(output_type)

    if input_type == 'infer':
        print('Inferring the nomenclature based on gene prefixes...')
        matches = {in_type:(genes_in.str.slice(0, len(prefix)) == prefix).mean() for in_type, prefix in prefixes.items()}
        best_match = max(matches, key=matches.get)
        if matches[best_match] >= match_threshold:
            print("Based on gene prefixes, I conclude that the input are", best_match, '{:.2%} of input genes have this prefix.'.format(matches[best_match]))
            input_type = best_match
        else:
            print("Input genes have no known prefix.")

    if input_type == 'infer':
        print("Inferring the nomenclature based on overlap of names...")
        overlaps = pd.Series({in_type:genes_in.isin(df[in_type]).mean() for in_type in symbolic_to_infer})
        acceptable_overlaps = overlaps[overlaps >= match_threshold]
        if len(acceptable_overlaps) >= 1:
            input_type = acceptable_overlaps.idxmax()
            print('Based on set overlap, I concluded that the input are', input_type, '{:.2%} of the input genes fit this type.'.format(acceptable_overlaps.loc[input_type]))
            if len(acceptable_overlaps) > 1:
                print('These other nomenclatures had an overlap >= {:.1%}'.format(match_threshold))
                print(acceptable_overlaps.drop(input_type).to_string())
        else:
            print("Couldn't infer a nomenclature.")
            print(overlaps.idxmax(), 'is my best guess ({:.2%} overlap with the input genes).'.format(overlaps.max()))
            raise RuntimeError
    

    df[output_type] = df[output_type].fillna(no_output)

    def robust_map(input_type, output_type):
        return df.set_index(input_type, drop=False)[output_type] if not input_type in non_unique else invert_non_unique_map(hgnc_df.set_index(output_type)[input_type])

    Map = robust_map(input_type, output_type) 
    clean_map = Map.loc[~Map.index.isnull()]
    if clean_map.index.duplicated().any():
        clean_map = clean_map.groupby(level=0).agg(lambda s: s.iloc[0] if len(s.unique()) == 1 else np.nan)

    reduced = clean_map.reindex(genes_in.values)
    if input_type in backups:
        alternates = backups[input_type]
        if not case_sensitive and input_type in symbolic:
            alternates.insert(0, input_type)        # The first map is always case-sensitive, only cleanup measures are case-insensitive. 
        for alternate in alternates:
            unmapped_genes = reduced.loc[reduced.isnull()]
            if len(unmapped_genes) == 0:
                break
            try:
                extra_map = robust_map(alternate, output_type)
            except KeyError:
                continue
            keep = extra_map.index.isin(unmapped_genes.index) if case_sensitive else extra_map.index.str.upper().isin(unmapped_genes.index.str.upper())
            extras = extra_map.loc[keep]
            if len(extras) > 0:
                print('Your input nomenclature are {:}, but {:} ({:.1%}) of the un-mapped symbols are {:}. Adding these alternatives to the map.'.format(
                  input_type, len(extras), len(extras)/len(unmapped_genes), alternate))
                reduced.loc[extras.index] = extras.values

    failed = reduced == no_output
    print('{:.1%} of genes mapped to no output'.format(failed.mean()))
    failed_genes = reduced[failed]
    if len(failed_genes) < 50:
       print('The following genes were not mapped:')
       print(failed_genes.to_string())
    return reduced 

