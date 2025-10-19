# Reproducing Figures from "Group Normalization for Genomic Data"

This repository contains R scripts to exactly reproduce Figures 1 and 3 from the PLOS One paper "Group Normalization for Genomic Data" by Ghandi et al. (2012).

## Required Data

To reproduce the exact figures from the paper, you need to download the original raw data:

1. **Raw CEL files from Lee et al. (2007)**:
   - Available from: http://chemogenomics.stanford.edu/supplements/03nuc/datasets.html
   - Dataset S3: "Raw Microarray Data" - Contains CEL files output from Affymetrix GeneChip operating system

2. **Processed Data (if raw data is unavailable)**:
   - The analyzed data described in the paper can be found at the same URL
   - Dataset S1: "Analyzed data"

## Setup

1. Create a directory called "data"
2. Download and extract the raw CEL files into the data directory
3. Update the file paths in the script to match your downloaded CEL files

## Data Processing

The script follows the exact methodology described in the paper:

1. **Reading CEL files**: Uses the `affy` package to read and process Affymetrix CEL files.
2. **TAS-based Normalization**: Implements normalization similar to Affymetrix Tiling Analysis Software (TAS) as specified in the paper: "All 3 nucleosomal and 3 control CEL files were imported into Affymetrix Tiling Analysis Software (TAS) as a 2-sample analysis."
3. **Group Normalization**: Implements the Group Normalization algorithm exactly as described in the paper.
4. **Visualization**: Creates Figures 1 and 3 exactly as shown in the paper.

## Figure 1: Genomic Hybridization Signals

This figure demonstrates that genomic hybridization signals from different replicates are highly reproducible (high correlation), but have significant variation in efficiency across the genome (non-uniform signal).

- **Part A**: Shows two replicate genomic hybridization signals along chromosome III
- **Part B**: Shows the correlation between two genomic hybridizations across the genome

## Figure 3: Group Normalization Method

This figure illustrates how the Group Normalization method works:

- Probes are sorted by their values in a reference condition (black line)
- For each probe, a reference set consisting of probes with similar response is identified (dashed boxes)
- High (red) and low (green) signals from experimental conditions are normalized based on these reference sets

## Requirements

The script requires the following R packages:
- BiocManager (for installing Bioconductor packages)
- affy (for processing Affymetrix CEL files)
- oligo (alternative for processing CEL files)
- affxparser (for low-level CEL file parsing if needed)
- GenomicRanges
- rtracklayer

## Citation

Ghandi M, Lee D, Mohammad-Noori M, Beer MA (2012) Group Normalization for Genomic Data. PLoS ONE 7(8): e38695. https://doi.org/10.1371/journal.pone.0038695

Lee W, Tillo D, Bray N, Morse RH, Davis RW, et al. (2007) A high-resolution atlas of nucleosome occupancy in yeast. Nat Genet 39: 1235–1244. https://doi.org/10.1038/ng2117
