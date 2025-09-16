# Technical Description of Affymetrix Tiling Arrays

## 1. Array Design: Probe Design, Number, and Positioning

### Core Design Principles

Affymetrix tiling arrays represent a specialized high-density microarray platform designed to provide comprehensive, unbiased coverage of genomic regions. Unlike traditional microarrays that target specific genes or exons, tiling arrays systematically interrogate contiguous genomic regions or entire genomes with regularly spaced oligonucleotide probes, regardless of known annotations.

### Probe Design

- **Oligonucleotide Length**: Typically 25-nucleotide (25-mer) DNA oligonucleotides
- **Probe Types**:
  - **Perfect Match (PM)**: Sequences exactly complementary to target genomic regions
  - **Mismatch (MM)**: Identical to PM probes except for a single nucleotide mismatch at the central position (13th base), serving as specificity controls
- **Selection Criteria**:
  - Avoidance of repetitive elements
  - Balancing of GC content
  - Minimizing cross-hybridization potential
  - Optimizing thermodynamic properties for consistent hybridization

### Probe Number and Density

- **Probe Quantity**: Modern Affymetrix tiling arrays contain millions of probe features
  - Example: Drosophila Genome Tiling Arrays contain approximately 6.4 million oligonucleotides, with approximately 3.2 million PM probes
  - Human genome tiling arrays may contain 40-50 million probes
- **Density/Resolution**: Described by the distance between the start positions of adjacent probes
  - Common resolutions:
    - High-resolution: 5-35 base pairs (bp) between probe starts
    - Medium-resolution: 35-100 bp between probe starts
    - Low-resolution: >100 bp between probe starts
  - For example, the Affymetrix S. cerevisiae Tiling 1.0R Array averages one probe every 5 bp

### Probe Positioning

- **Tiling Strategy**:
  - **End-to-end tiling**: Probes positioned directly adjacent to each other
  - **Overlapping tiling**: Probes overlap with adjacent probes
  - **Gapped tiling**: Probes spaced at regular intervals with gaps between them
- **Coverage Types**:
  - **Genome-wide**: Complete coverage of an entire genome
  - **Targeted**: Focused coverage of specific chromosomal regions or gene clusters
  - **Promoter-focused**: Concentrated coverage around transcription start sites
- **Strand Representation**:
  - **Single-strand**: Probes designed to hybridize to one DNA strand
  - **Double-strand**: Probes designed to hybridize to both DNA strands

## 2. Probe Data: Collection and Interpretation

### Data Collection

- **Sample Preparation**:
  - Target DNA/RNA is fragmented, labeled with biotin or fluorescent dyes
  - Labeled targets hybridize to complementary probes on the array
  - After washing to remove non-specifically bound material, arrays are scanned
- **Raw Data Acquisition**:
  - Scanner captures fluorescence intensity at each probe position
  - Intensity values correspond to the amount of target bound to each probe
  - Data is stored as CEL files containing raw intensity values for each probe

### Data Preprocessing

- **Background Correction**:
  - Removes non-specific hybridization signal
  - Methods include: MAS5.0, RMA, GC-RMA, and probe-specific models
- **Normalization**:
  - Removes technical variation between arrays
  - Common methods: quantile normalization, scaling, and Group Normalization (as described by Ghandi and Beer)
- **Probe-level Summarization**:
  - For arrays with PM/MM probe pairs, MM signal may be subtracted from PM signal
  - Neighboring probes often smoothed or averaged in sliding windows
  - Specialized algorithms like MAT (Model-based Analysis of Tiling arrays) model probe behavior based on sequence composition

### Data Interpretation

- **Signal Representation**:
  - Processed data typically visualized as intensity values along genomic coordinates
  - For ChIP-chip, data often represented as log2(ChIP/Input) ratios
  - For expression studies, intensities indicate transcriptional activity
- **Feature Detection**:
  - Peak calling algorithms identify regions of significant signal enrichment
  - Hidden Markov Models (HMMs) detect transitions between bound/unbound states
  - Sliding window approaches with statistical thresholds
- **Statistical Analysis**:
  - Multiple testing correction for millions of probes (FDR, Bonferroni)
  - Significance determined by comparing to background models
  - Permutation tests to establish empirical significance thresholds

## 3. Example Arrays and Associated Analyses

### Affymetrix S. cerevisiae Tiling 1.0R Array

- **Specifications**:
  - Covers entire 12.5 Mb S. cerevisiae genome
  - ~3.2 million perfect match probes
  - 5 bp resolution (average distance between adjacent probes)
- **Notable Applications**:
  - **Nucleosome Mapping**: Lee et al. (2007) created high-resolution maps of nucleosome occupancy
    - Methodology: Micrococcal nuclease digestion followed by IP of histone H3
    - Analysis: Sliding window approach to identify nucleosome positions
    - Findings: Revealed nucleosome-free regions upstream of genes and positioning relative to transcription start sites
  - **Transcriptome Analysis**: David et al. (2006) conducted comprehensive transcriptome analysis
    - Methodology: RNA hybridization under various conditions
    - Analysis: Segmentation algorithms to define transcriptional units
    - Findings: Identified hundreds of novel non-coding RNAs and revised gene boundaries

### Affymetrix Drosophila Genome Tiling Arrays

- **Specifications**:
  - ~3.2 million PM probes covering Drosophila genome
  - One probe per 35 bp of genomic sequence on average
- **Notable Applications**:
  - **Chromatin Modification Mapping**: By Kharchenko et al. (2011)
    - Methodology: ChIP-chip for various histone modifications
    - Analysis: Hidden Markov Models to classify chromatin states
    - Findings: Identified distinct chromatin signatures associated with functional elements
  - **Transcription Factor Binding**: By Li et al. (2008)
    - Methodology: ChIP-chip for developmental transcription factors
    - Analysis: MAT algorithm for peak detection
    - Findings: Revealed combinatorial binding patterns during development

### Affymetrix Human Tiling Arrays

- **Specifications**:
  - Multiple array sets covering human genome
  - 35-100 bp resolution depending on array version
- **Notable Applications**:
  - **ENCODE Project**: Comprehensive analysis of functional elements
    - Methodology: Multiple assays including ChIP-chip, DNase-seq
    - Analysis: Integrated analysis across multiple data types
    - Findings: Mapped transcription factor binding sites, chromatin modifications, and transcriptionally active regions
  - **DNA Methylation Analysis**: By Weber et al. (2005)
    - Methodology: Methylated DNA immunoprecipitation (MeDIP)
    - Analysis: Sliding window approach with significance testing
    - Findings: Identified methylation patterns in promoters and gene bodies

## 4. Technical Challenges and Considerations

### Experimental Challenges

- **Cross-hybridization**: Non-specific binding of targets to similar probe sequences
- **Probe Behavior**: Variable hybridization efficiency due to sequence composition
- **Dynamic Range**: Limited range of signal detection compared to sequencing approaches
- **Spatial Artifacts**: Technical biases in different regions of the physical array

### Analytical Challenges

- **Multiple Testing**: Statistical challenges due to millions of probes tested simultaneously
- **Normalization**: Accounting for probe-specific behavior and array biases
- **Spatial Resolution**: Limitations in precisely defining boundaries of features
- **Data Integration**: Combining tiling array data with other genomic data types

### Modern Context

While high-throughput sequencing technologies (ChIP-seq, RNA-seq) have largely superseded tiling arrays for many applications, tiling arrays remain valuable in specific contexts:

- Well-established analysis pipelines with known error characteristics
- Lower cost for certain applications
- Lower computational requirements for analysis
- Historical datasets valuable for meta-analyses and reanalysis with improved methods

## 5. Advanced Analysis Approaches

### Signal Enhancement Methods

- **Probe Sequence Normalization**: Models accounting for sequence-specific hybridization behaviors
- **Wavelet-based Denoising**: Removal of noise while preserving signal structures
- **Group Normalization**: Addressing both global variation and probe-specific biases
- **Cross Normalization**: Amplifying biologically relevant differences between datasets

### Integrative Analysis

- **Multi-platform Integration**: Combining tiling array data with sequencing data
- **Multi-factor Analysis**: Integrating multiple transcription factors or histone modifications
- **Time-course Analysis**: Studying dynamic changes in binding or expression
- **Network Inference**: Reconstructing regulatory networks from binding and expression data

---

*This technical description provides an overview of Affymetrix tiling array technology, focusing on design principles, data processing, and applications. While specific details may vary between array versions and experimental protocols, the fundamental concepts remain consistent across the platform.*