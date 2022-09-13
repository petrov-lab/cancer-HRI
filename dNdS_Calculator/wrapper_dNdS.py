import argparse
from calculate_dNdS import *

### Description message ###
parser = argparse.ArgumentParser(description="Calculates dN/dS using a permutation-based null model of mutagenesis per sample.", epilog="Please cite: Tilk, S., Tkachenko, S., Curtis, C., Petrov, D. & McFarland, C. D. Most cancers carry a substantial deleterious load due to Hill-Robertson interference. eLife. (2022) https://doi.org/10.7554/eLife.67790")

### Required arguments ###
requiredArgs = parser.add_argument_group('required arguments')
requiredArgs.add_argument("-input", help = "Tab-separated file of SNPs that are in MAF format. Must contain the following columns: 'Variant_Classification', 'HGVSc', 'Tumor_Sample_Barcode', and a column for gene names (either 'Hugo_Symbol', 'Gene' or 'Transcript_ID'.)", required=True)
requiredArgs.add_argument("-output", help = "Path to output file of dN/dS values.", default='stdout', required=True)
requiredArgs.add_argument("-ref", help ="Reference genome that SNPs are mapped to, either 'GRCH37' or 'GRCH38'.", default='GRCH37', required=True)

### Parse arguments ###
args = parser.parse_args()

### Calculate dN/dS ###
GetdNdS(args.input, args.output)

