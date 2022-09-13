import argparse
from calculate_dEdI import *

### Description message ###
parser = argparse.ArgumentParser(description="Calculates dE/dI for all CNAs per sample.",epilog="Please cite: Tilk, S., Tkachenko, S., Curtis, C., Petrov, D. & McFarland, C. D. Most cancers carry a substantial deleterious load due to Hill-Robertson interference. eLife. (2022) https://doi.org/10.7554/eLife.67790")

### Required arguments ###
requiredArgs = parser.add_argument_group('required arguments')
requiredArgs.add_argument("-input", help = "Tab-separated file of CNAs that contains the following columns:'Sample','Start','End'.", required=True )
requiredArgs.add_argument("-output", help = "Path to output file of dE/dI values.", default='stdout', required=True)

### Optional arguments ###
parser.add_argument("-AddGeneSet", default=False, help = "Tab-separated file with a column of gene names ('GeneName') and a column of strings ('Type') that identifies the gene set. This will calculate dE/dI per sample for all inputted CNAs, and separately for these sets of genes.", required=False)

### Parse arguments ###
args = parser.parse_args()

### Calculate dE/dI ###
GetdEdI(args.input, args.output)


