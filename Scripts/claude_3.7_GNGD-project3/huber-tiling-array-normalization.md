# Normalization in Huber et al. 2006 Tiling Array Paper

## Paper Reference
Huber, W., Toedling, J., and Steinmetz, L.M. (2006). Transcript mapping with high-density oligonucleotide tiling arrays. Bioinformatics 22, 1963-1970.

## The Biology
- High-density oligonucleotide tiling arrays used for unbiased transcript mapping
- Detects expressed regions without relying on existing gene annotations
- Can identify novel transcripts, non-coding RNAs, and alternative splice variants
- Used for mapping transcriptomes in organisms like yeast

## Normalization Challenges
- Sequence-dependent probe response (different oligonucleotides have different hybridization affinities)
- Need to accurately determine transcript boundaries along genomic coordinates
- Signal variability between different probes makes quantitative comparison difficult

## Normalization Approaches Considered

### 1. DNA Reference Normalization (Their Main Method)
- Uses genomic DNA hybridization as a reference
- Adjusts for sequence-specific probe response
- Employs empirical probe response parameters from reference hybridizations
- Accounts for intrinsic hybridization properties of each probe

### 2. Perfect Match (PM) and Mismatch (MM) Probes
- PM probes: match genomic sequence exactly
- MM probes: contain single nucleotide change to estimate non-specific binding
- Not applicable for all platforms (e.g., Nimblegen doesn't use MM probes)

### 3. Statistical Normalization Methods
- Alternatives like GC content normalization (dividing signal by median of probes with identical GC content)
- Quantile normalization extended to multivariate probe sequence space
- Naef and Magnasco model-based approaches

## Evaluation Methods
- Application to real biological data (yeast transcriptome)
- Comparison with biological expectations rather than just statistical criteria
- Used alongside segmentation algorithm (piecewise constant expression profile)
- Statistical metrics like AIC (Akaike Information Criterion) and BIC (Bayesian Information Criterion)
- Performance in identifying known and novel transcripts

## Effectiveness
- Successfully identified transcript boundaries, structure and expression levels
- Detected expected transcripts and novel features (operon-like transcripts, complex architectures)
- Mapped positions of 3' and 5' UTRs of coding genes
- Identified hundreds of antisense RNA transcripts
- For yeast data, corrections showed minimal improvements over raw data

## Key Findings on Control Features
- Control features showed high enrichment requiring correction
- Sequence normalization methods reduced enrichment of control features while retaining enrichment of target features
- Different effectiveness depending on organism and experimental context
- DNA reference normalization provided clean normalized signal for accurate transcript boundary identification

## Conclusion
The DNA reference normalization method combined with their segmentation algorithm was effective for transcript mapping in yeast, providing accurate transcript boundary identification and allowing detection of novel transcriptional features. The approach is adaptable to different organisms and experimental conditions.