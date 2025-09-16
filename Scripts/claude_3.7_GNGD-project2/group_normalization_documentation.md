# Documentation: Group Normalization for Genomic Data

This document provides documentation for the R code implementation of the "Group Normalization for Genomic Data" method by Mahmoud Ghandi and Michael A. Beer (PLOS ONE, 2012).

## Overview

Group Normalization is a method for normalizing genomic data (like microarray or sequencing data) that addresses both global variations between experiments and local probe-specific effects in one integrated step. Unlike other normalization methods, it doesn't assume identical signal distributions across conditions and can efficiently correct for nonlinear probe effects.

The paper introduces two main algorithms:
1. **Group Normalization (GN)** - For normalizing probe values in experimental conditions
2. **Cross Normalization** - For amplifying biologically relevant differences between two conditions

## Implementation Files

Two R script artifacts are provided:

1. **Group Normalization Code** - Complete implementation with all algorithms and code to reproduce the paper's figures
2. **Group Normalization Example** - Simplified demonstration using simulated data

## Core Concepts

### The Problem

In genomic data analysis (e.g., tiling arrays for nucleosome positioning):
- Different probes have varying hybridization efficiencies (the "probe effect")
- Even in control conditions with uniform DNA concentration, the signal varies significantly across probes
- These probe-specific variations are highly reproducible between experiments
- Traditional normalization methods often make assumptions about signal distributions that don't hold for all experimental designs

### Group Normalization Solution

The key insight of Group Normalization:
- For each probe, find a set of "reference probes" that have similar dynamic response properties
- Use these reference probes to estimate how the probe responds to both high and low biological signals
- Use this estimated response to normalize the probe's value in experimental conditions

### Mathematical Model

The observed signal y<sub>i</sub> for a probe i is modeled as:

y<sub>i</sub> = A<sub>i</sub> × x<sub>i</sub> + B<sub>i</sub> + ε

Where:
- x<sub>i</sub> is the desired biological signal (normalized to 0-1)
- A<sub>i</sub> is a probe-specific efficiency factor
- B<sub>i</sub> is a probe-specific background signal
- ε is random noise

Group Normalization estimates A<sub>i</sub> and B<sub>i</sub> for each probe by finding reference probes with similar behavior.

## Algorithms

### Binary Group Normalization

This method estimates two signal levels (high and low) for each probe:

1. For each probe p<sub>i</sub>, find a set of reference probes Ref(p<sub>i</sub>) that have similar signals in a reference condition
2. Sort the reference probes by their values in the experimental condition
3. Calculate μ<sub>i,low</sub> (mean of the lowest 30% of reference probes) and μ<sub>i,high</sub> (mean of the highest 30% of reference probes)
4. Normalize the probe value using: x<sub>i</sub> = (y<sub>i</sub> - μ<sub>i,low</sub>) / (μ<sub>i,high</sub> - μ<sub>i,low</sub>)

### Quantile-based Group Normalization

An alternative approach that doesn't require defining high and low signal ranges:

1. For each probe, find its reference set as in the binary method
2. Apply quantile normalization to the reference set
3. Use the probe's rank in the reference set to determine its normalized value

### Cross Normalization

For comparing two experimental conditions:

1. Use condition A as the reference for normalizing condition B
2. Use condition B as the reference for normalizing condition A
3. The correspondence between these asymmetric normalizations highlights significant differences

### Reference Group Assignment

Three methods for finding reference probes:

1. **Single Reference Method** - Sort probes by value in one reference condition
2. **Sorted Average Method** - Sort probes by average value across multiple reference conditions
3. **Minimum Distance Method** - Find probes with minimum distance in a multidimensional reference space

## Implementation Details

### Key Functions

- `find_reference_probes()` - Identifies reference probes for a given probe
- `binary_group_normalize()` - Implements Binary Group Normalization
- `quantile_group_normalize()` - Implements Quantile-based Group Normalization
- `cross_normalize()` - Implements Cross Normalization
- `calculate_signal_quality()` - Calculates signal quality metric for evaluation

### Signal Quality Evaluation

The implementation includes the signal quality metric defined in the paper:
- Signal (S): Mean square change in signal for significantly changed probes between conditions
- Noise (N): Mean square change in signal between replicates
- Signal Quality in dB = 10 × log<sub>10</sub>(S/N)

## Usage

### Basic Usage

```r
# Load genomic data
reference_data <- read.table("genomic_reference.txt")
experiment_data <- read.table("nucleosome_enriched.txt")

# Apply Binary Group Normalization
normalized_data <- binary_group_normalize(
  reference_data, 
  experiment_data,
  n_refs = 1000,
  low_range = c(0.1, 0.4),
  high_range = c(0.6, 0.9)
)

# Apply Cross Normalization for two conditions
condition_a <- read.table("before_treatment.txt")
condition_b <- read.table("after_treatment.txt")
cross_norm_result <- cross_normalize(condition_a, condition_b)
```

### Handling Repetitive Regions

The implementation allows excluding repetitive regions:

```r
# Identify repetitive regions
repeat_indices <- which(is_repetitive_region)

# Apply normalization excluding repeats
normalized_data <- binary_group_normalize(
  reference_data, 
  experiment_data,
  exclude_repeats = TRUE,
  repeat_indices = repeat_indices
)
```

## Visualization Functions

The implementation includes functions to recreate the figures from the paper:

- `plot_joint_distribution()` - Visualizes probe distribution before and after normalization (Figure 5A)
- `plot_nucleosome_occupancy()` - Plots nucleosome occupancy at genomic regions (Figure 5B)
- `plot_cross_normalization()` - Visualizes differential signals from cross normalization (Figure 5C)
- `plot_signal_quality()` - Compares signal quality across methods (Figure 7)

## Advantages of Group Normalization

As demonstrated in the implementation:

1. It significantly improves signal quality compared to MAS5, MAT, and Quantile normalization
2. It does not require the assumption that treatment and control have identical signal distributions
3. It is flexible enough to correct for nonlinear and higher-order probe effects
4. Cross normalization efficiently amplifies biologically relevant differences
5. The method is computationally efficient and easy to implement

## Applications

The Group Normalization method can be applied to various genomic datasets:

1. Nucleosome positioning data (primary example in the paper)
2. ChIP-chip and ChIP-seq experiments
3. RNA-seq data for differential expression analysis
4. Any genomic dataset with probe-specific or sequence-specific effects

## Extension to Other Technologies

While the implementation focuses on tiling array data, the Group Normalization approach can be adapted to other high-throughput technologies that suffer from sequence-specific effects:

- Sequence-dependent shearing rates
- Endonuclease sequence cleavage preferences
- Sequence-specific priming efficiencies in massively parallel sequencing
- GC content biases in sequencing technologies

## Conclusion

The Group Normalization method provides a data-driven approach to normalization that estimates probe parameters from similar probes rather than using explicit models. This makes it robust to various types of probe effects and applicable to diverse genomic datasets.

The implementation provided in these R scripts demonstrates the effectiveness of the method using simulated data and provides tools to apply it to real genomic datasets.
