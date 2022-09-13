# Hill-Robertson Interference (HRI) in cancer

## Paper

> Tilk, S., Tkachenko, S., Curtis, C., Petrov, D. & McFarland, C. D. Most cancers carry a substantial deleterious load due to Hill-Robertson interference. eLife. (2022) https://doi.org/10.7554/eLife.67790. https://elifesciences.org/articles/67790

Scripts used to generate and plot data in this manuscript can be found in the Analysis directory. User-friendly code to calculate dE/dI and dN/dS (with a permutation-based null model of mutagenesis) can also be found below. 

## Environment
All dependencies required to run the code can be found in Environment directory. Create the virtual environment and load it using conda.

```
conda create -f Environment/environment.yaml
source activate cancer-HRI
```

## dN/dS Calculator
A permutation-based, nonparametric (parameter-free) estimation of dN/dS. In this approach, every observed mutation is permuted while preserving the gene, patient samples, specific base change (e.g. A>T) and tri-nucleotide context. The permutations are then tallied for both nonsynonymous and synonymous substitutions and used as proportional estimates of the observed number of nonsynonymous and synonymous mutations in the absence of selection. 

### Usage
```
usage: wrapper.py [-h] -input INPUT -output OUTPUT -ref REF

Calculates dN/dS using a permutation-based null model of mutagenesis per sample.

optional arguments:
  -h, --help      show this help message and exit

required arguments:
  -input INPUT    Tab-separated file of SNPs that are in MAF format. Must
                  contain the following columns: 'Variant_Classification',
                  'HGVSc', 'Tumor_Sample_Barcode', and a column for gene names
                  (either 'Hugo_Symbol', 'Gene' or 'Transcript_ID'.)
  -output OUTPUT  Path to output file of dN/dS values.
  -ref REF        Reference genome that SNPs are mapped to, either 'GRCH37' or
                  'GRCH38'.
```
### Examples
Wrapper script should be ran from within the dNdS_Calculator directory.

```
python wrapper_dNdS.py -input examples/dNdS_Input.txt -output examples/dNdS_Output.txt -ref "GRCH37"
```

## dE/dI Calculator
Quantifies selection in CNAs using two alternative measures: Breakpoint Frequency and Fractional Overlap.
For both measures, we compare the number of CNAs that either terminate (Breakpoint Frequency) within or partially overlap (Fractional Overlap) Exonic regions of the genome relative to non-coding (Intergenic and Intronic) regions (dE/dI). Like dN/dS, dE/dI is expected to be <1 in genomic regions experiencing negative selection, >1 in regions experiencing positive = selection (e.g. driver genes), and approximately 1 when selection is absent or inefficient. Confidence intervals (95%) are determined by bootstrapping.

### Usage

```
usage: wrapper_dEdI.py [-h] -input INPUT -output OUTPUT [-AddGeneSet ADDGENESET]

Calculates dE/dI for all CNAs per sample.

optional arguments:
  -h, --help            show this help message and exit
  -AddGeneSet ADDGENESET
                        Tab-separated file with a column of gene names
                        ('GeneName') and a column of strings ('Type') that
                        identifies the gene set. This will calculate dE/dI for
                        all inputted CNAs, and separately for these sets of
                        genes.

required arguments:
  -input INPUT          Tab-separated file of CNAs that contains the following
                        columns:'Sample','Start','End'.
  -output OUTPUT        Path to output file of dE/dI values.
```

### Examples
Wrapper script should be ran from within the dEdI_Calculator directory.
To calculate dE/dI of CNAs in a set of genes, e.g. drivers, add optional argument -AddGeneSet.

```
python wrapper_dEdI.py -input examples/exampleCNAInput.txt 
                  -output examples/CNAOutputWithGeneTracks.txt 
                  -AddGeneSet examples/drivers_Bailey2018.txt
```

