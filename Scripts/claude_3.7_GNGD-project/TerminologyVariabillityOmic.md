# Comprehensive Terminology for Variability in Genomic Data Analysis

## The vocabulary of variability in genomic data analysis

The field of genomic data normalization has developed a rich vocabulary for describing sources of variation and methods for their correction. This terminology, largely formalized through the work of Terry Speed's group at UC Berkeley and WEHI, and recently unified by Gerard and Stephens, provides the conceptual framework for understanding group normalization approaches.

## Core Concepts of Variation

### Fundamental Variation Dichotomies

**Wanted versus Unwanted Variation** represents the central organizing principle in genomic normalization. **Wanted variation** encompasses the biological signals of scientific interest—treatment effects, disease states, genetic differences, or temporal changes that researchers aim to detect and quantify. In contrast, **unwanted variation** includes all sources of variability that confound the analysis but are not of primary scientific interest. This concept was formalized in Gagnon-Bartsch and Speed's seminal 2012 paper in *Biostatistics*, which introduced the RUV (Remove Unwanted Variation) framework. The distinction is context-dependent: batch effects that are unwanted in most analyses might become wanted variation when studying laboratory reproducibility.

**Technical versus Biological Variation** divides variability by its origin. **Technical variation** arises from experimental procedures rather than biological differences, including pipetting errors, RNA degradation, PCR amplification biases, sequencing depth differences, and platform-specific artifacts. **Biological variation** represents legitimate differences between samples due to actual biological phenomena, including both the variation of interest and biological confounders like age or sex effects. Yang et al. (2002) in *Nucleic Acids Research* established methods for separating these sources in microarray data.

**Observable versus Unobservable Factors** distinguishes between variation sources based on measurement capability. **Observable factors** (also called known covariates) are directly measured and recorded during experiments—batch identifiers, processing dates, technician identities, or sample collection times. **Unobservable factors** (hidden or latent factors) affect gene expression but weren't measured or recorded. Leek and Storey's 2007 surrogate variable analysis (SVA) paper in *PLoS Genetics* provided the statistical framework for estimating these hidden factors.

**Systematic versus Random Variation** categorizes effects by their predictability. **Systematic variation** follows reproducible patterns that can be modeled and corrected—dye biases in two-color arrays, spatial effects on array surfaces, or consistent batch effects. **Random variation** represents unpredictable fluctuations following statistical distributions, typically modeled as measurement error or biological noise that cannot be deterministically removed.

## The Architecture of Normalization

### Levels of Normalization

**Within-sample normalization** corrects for biases affecting different features within the same sample. For RNA-seq, this primarily involves adjusting for gene length and GC content biases, enabling fair comparison of expression levels across genes of different lengths. The RPKM, FPKM, and TPM units introduced by Mortazavi et al. (2008) and refined by Li and Dewey (2011) standardized this approach.

**Between-sample normalization** enables valid comparisons across different samples by correcting for global differences in library size, RNA composition, and technical processing. The field has developed sophisticated methods including **quantile normalization** (Bolstad et al., 2003), which forces all samples to have identical distributions, and **TMM normalization** (Robinson and Oshlack, 2010), which uses a trimmed mean of M-values to account for compositional biases.

**Between-group normalization** specifically addresses systematic differences between experimental groups that aren't of biological interest. This includes corrections for batch effects when treatment groups are processed separately, or adjustments for confounding variables that differ systematically between groups. The ComBat method (Johnson et al., 2007) exemplifies this approach, using empirical Bayes to remove known batch effects while preserving biological variation.

## Control Features and Reference Systems

### Types of Control Genes

The concept of **negative control genes** is fundamental to modern normalization methods. These are genes assumed to be unaffected by the conditions of interest, with expression variation reflecting only unwanted technical or biological factors. The RUV methods formalized their use, with RUVg specifically using negative controls for factor estimation. Selection can be based on biological knowledge (housekeeping genes), external spike-ins (ERCC controls), or empirical identification from the data itself.

**Positive control genes** have known, expected expression changes between conditions. While less commonly used in normalization, they provide validation that methods preserve true biological signals. **Empirical control genes** represent an important innovation where controls are identified computationally from the data, typically as the least differentially expressed genes in an initial analysis pass.

**Housekeeping genes** are endogenous genes with presumed stable expression across conditions, traditionally used for qPCR normalization. However, RNA-seq has revealed that many classical housekeeping genes show substantial variation, leading to more sophisticated control selection methods.

## Mathematical Framework and Factor Models

### The RUV Model Structure

The RUV framework represents gene expression as **Y = Xβ + Wα + ε**, where Y is the expression matrix, X captures factors of wanted variation with coefficients β, W represents factors of unwanted variation with loadings α, and ε is random error. This decomposition, introduced by Gagnon-Bartsch and Speed (2012), provides the mathematical foundation for separating wanted from unwanted effects.

**Factor loadings** (the α matrix) describe how each unwanted factor affects each gene's expression. These are typically estimated from control genes or replicate samples. The **number of unwanted factors** (k) is a crucial parameter, often determined by cross-validation or scree plots. Gerard and Stephens (2021) showed how different choices of k relate to different RUV variants.

### RUV Method Variants

The proliferation of RUV methods reflects different strategies for estimating unwanted factors:

**RUVg** uses negative control genes exclusively, performing factor analysis on this subset to estimate W. **RUVs** uses replicate samples or negative control samples where the factors of interest are constant, estimating unwanted factors from the variation within these groups. **RUVr** uses residuals from a first-pass regression of all genes on the factors of interest, essentially treating genes that show little differential expression as empirical controls.

Gerard and Stephens' **RUV\*** (2021) provided a unifying framework showing these methods as special cases of a general approach. Their RUVB extends this with Bayesian factor analysis, accounting for uncertainty in factor estimation—a limitation of earlier methods.

## Batch Effects and Confounding

### Understanding Batch Effects

**Batch effects** are systematic differences between groups of samples processed together. The term encompasses any non-biological variation correlating with experimental batches—different days, technicians, reagent lots, or equipment runs. Leek et al. (2010) demonstrated their pervasive impact, showing batch effects often exceed biological signals in magnitude.

**Confounding** occurs when unwanted variation correlates with factors of interest, making it impossible to separate technical from biological effects without additional information or assumptions. **Perfect confounding** represents the extreme case where batch completely aligns with treatment groups, making correction impossible without biological replicates spanning batches.

### Surrogate Variable Analysis

**Surrogate variables** are statistically constructed variables representing unmeasured sources of variation. Introduced by Leek and Storey (2007), SVA estimates these latent factors by identifying systematic patterns in expression data not explained by known variables. The method assumes these patterns reflect batch effects or other technical artifacts rather than interesting biology.

**Hidden factors** and **latent factors** are related concepts referring to unmeasured sources of systematic variation. While sometimes used interchangeably with surrogate variables, they can also refer to true biological factors that weren't measured but affect expression patterns.

## Statistical Properties and Distributions

### Count Data Characteristics

**Overdispersion** describes the common observation that RNA-seq count variance exceeds what the Poisson distribution predicts. This extra-Poisson variation arises from both technical sources (varying capture efficiency) and biological heterogeneity. The negative binomial distribution, with its additional dispersion parameter, has become the standard model for RNA-seq counts.

**Zero inflation** refers to an excess of zero counts beyond what count distributions predict. In single-cell RNA-seq, this manifests as **dropout events**—genes appearing unexpressed due to technical failures rather than biological absence. Methods like ZINB (zero-inflated negative binomial) explicitly model these excess zeros.

### Compositional Effects

**Compositional data** constraints arise because RNA-seq measures relative rather than absolute abundance. If one gene's expression increases dramatically, others appear relatively decreased even if their absolute expression is unchanged. This **compositional bias** motivated the development of TMM normalization and similar methods that assume most genes are unchanged.

## Modern Extensions and Unified Frameworks

### The Gerard-Stephens Contributions

David Gerard and Matthew Stephens expanded the RUV framework in several important directions. Their **MOUTHWASH** method (2020) combines empirical Bayes shrinkage with unwanted variation removal, providing better uncertainty quantification through **local false sign rates** rather than traditional FDR.

The **RUV\* framework** (2021) unified previously disparate RUV methods, showing them as special cases of a general matrix imputation problem. This theoretical advance enables development of new methods and clearer understanding of when different approaches are equivalent.

**RUVB** introduces full Bayesian treatment of unwanted variation, propagating uncertainty from factor estimation through to differential expression testing—addressing a key limitation where earlier methods treated estimated factors as known.

## Visualization and Quality Assessment

### Diagnostic Plots

**MA plots**, plotting log-ratios (M) against average intensities (A), reveal intensity-dependent biases and the effectiveness of normalization. Introduced for microarrays by Yang et al. (2002), they remain valuable for RNA-seq visualization.

**RLE plots** (Relative Log Expression) display the distribution of each gene's expression relative to its median across samples. Gandolfo and Speed (2018) showed how RLE plots effectively reveal unwanted variation patterns, with tight, centered distributions indicating successful normalization.

**Principal component analysis (PCA)** projects high-dimensional expression data onto major axes of variation. PCA plots reveal batch effects when samples cluster by batch rather than biological groups, and successful normalization should reduce batch-driven separation.

## Implementation Terminology

### Normalization Units

The evolution from **RPKM** to **FPKM** to **TPM** reflects growing understanding of compositional effects. While RPKM normalizes for sequencing depth then gene length, TPM reverses this order, ensuring the sum of all TPMs equals one million in every sample—providing true proportional abundances that are comparable across samples.

### Key Software Implementations

The terminology is embodied in software packages that have become standard tools: **limma** for linear modeling and batch correction, **sva** for surrogate variable analysis and ComBat, **RUVSeq** for applying RUV methods to RNA-seq, **edgeR** and **DESeq2** for differential expression with integrated normalization, and **vicar** for the Gerard-Stephens methods.

## Conclusion

This vocabulary provides the conceptual foundation for understanding group normalization in genomics. The terminology has evolved from simple scaling methods to sophisticated factor models that separate wanted from unwanted variation. The recent unification by Gerard and Stephens suggests the field is maturing toward a coherent theoretical framework, though challenges remain in handling increasingly complex experimental designs and single-cell data. Understanding this terminology is essential for both implementing normalization methods correctly and interpreting their results in biological context.