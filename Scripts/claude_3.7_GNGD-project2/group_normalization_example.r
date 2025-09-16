# Example Demonstration of Group Normalization with Simulated Data
# Based on "Group Normalization for Genomic Data" by Ghandi & Beer, 2012

# This script demonstrates the Group Normalization algorithm with simulated data
# We'll create a simple dataset and demonstrate the key concepts from the paper

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("preprocessCore", quietly = TRUE))
  BiocManager::install("preprocessCore")

library(preprocessCore)

############################################################
# Group Normalization Implementation
############################################################

#' Binary Group Normalization
#'
#' @param reference Probe values in reference condition
#' @param experiment Probe values in experiment condition
#' @param n_refs Number of reference probes to use
#' @param low_range Range for low signal (e.g., c(0.1, 0.4))
#' @param high_range Range for high signal (e.g., c(0.6, 0.9))
#' @return Normalized probe values
binary_group_normalize <- function(reference, experiment, n_refs = 1000,
                                  low_range = c(0.1, 0.4), high_range = c(0.6, 0.9)) {
  
  # Number of probes
  n_probes <- length(reference)
  
  # Initialize result
  normalized <- numeric(n_probes)
  
  # Sort probes by their value in the reference condition
  sorted_idx <- order(reference)
  
  # For each probe
  for (i in 1:n_probes) {
    # Find rank of probe i in reference
    rank_i <- which(sorted_idx == i)
    
    # Find reference set (probes with similar signal in reference)
    start_idx <- max(1, rank_i - floor(n_refs/2))
    end_idx <- min(n_probes, start_idx + n_refs - 1)
    
    # Adjust start if end is at the limit
    if (end_idx == n_probes) {
      start_idx <- max(1, end_idx - n_refs + 1)
    }
    
    # Get reference set indices
    ref_indices <- sorted_idx[start_idx:end_idx]
    
    # Get experiment values for reference probes
    ref_exp_values <- experiment[ref_indices]
    
    # Sort reference probes by their experiment values
    sorted_ref_exp <- sort(ref_exp_values)
    
    # Calculate low and high signal levels
    n_low_start <- ceiling(low_range[1] * length(sorted_ref_exp))
    n_low_end <- floor(low_range[2] * length(sorted_ref_exp))
    n_high_start <- ceiling(high_range[1] * length(sorted_ref_exp))
    n_high_end <- floor(high_range[2] * length(sorted_ref_exp))
    
    mu_low <- mean(sorted_ref_exp[n_low_start:n_low_end])
    mu_high <- mean(sorted_ref_exp[n_high_start:n_high_end])
    
    # Normalize probe signal
    normalized[i] <- (experiment[i] - mu_low) / (mu_high - mu_low)
  }
  
  return(normalized)
}

#' Cross Normalization
#'
#' @param condition_a Probe values in condition A
#' @param condition_b Probe values in condition B
#' @param n_refs Number of reference probes to use
#' @param low_range Range for low signal
#' @param high_range Range for high signal
#' @return List with normalized values for A vs B and B vs A
cross_normalize <- function(condition_a, condition_b, n_refs = 1000,
                           low_range = c(0.1, 0.4), high_range = c(0.6, 0.9)) {
  
  # Normalize A using B as reference
  a_norm <- binary_group_normalize(condition_b, condition_a, n_refs, low_range, high_range)
  
  # Normalize B using A as reference
  b_norm <- binary_group_normalize(condition_a, condition_b, n_refs, low_range, high_range)
  
  return(list(a_norm = a_norm, b_norm = b_norm))
}

############################################################
# Simulate Dataset for Demonstration
############################################################

set.seed(123) # For reproducibility

# Number of probes
n_probes <- 5000

# Simulate probe efficiency and background
probe_efficiency <- runif(n_probes, 0.5, 1.5)  # Each probe has different efficiency
probe_background <- runif(n_probes, 0, 0.5)    # Each probe has different background

# Create simulated biological conditions
# Genomic DNA reference - uniform concentration (=1) for all probes
genomic_reference <- probe_efficiency * 1 + probe_background + rnorm(n_probes, 0, 0.1)

# Create a nucleosome pattern (alternating bound/unbound regions)
nucleosome_signal <- rep(0, n_probes)
for (i in seq(1, n_probes, 200)) {
  start <- i
  end <- min(i + 149, n_probes)
  nucleosome_signal[start:end] <- 1  # Nucleosome bound regions
}

# Nucleosome-enriched condition
nucleosome_condition <- probe_efficiency * nucleosome_signal + probe_background + rnorm(n_probes, 0, 0.1)

# Create two experimental conditions with different nucleosome patterns
# Condition 1: before treatment
before_signal <- rep(0, n_probes)
for (i in seq(1, n_probes, 400)) {
  start <- i
  end <- min(i + 149, n_probes)
  before_signal[start:end] <- 1  # Nucleosome bound in specific regions
}
before_condition <- probe_efficiency * before_signal + probe_background + rnorm(n_probes, 0, 0.1)

# Condition 2: after treatment (different binding pattern)
after_signal <- rep(0, n_probes)
for (i in seq(201, n_probes, 400)) {
  start <- i
  end <- min(i + 149, n_probes)
  after_signal[start:end] <- 1  # Different nucleosome binding pattern
}
after_condition <- probe_efficiency * after_signal + probe_background + rnorm(n_probes, 0, 0.1)

############################################################
# Demonstrate Group Normalization
############################################################

# Apply Binary Group Normalization
cat("Applying Binary Group Normalization...\n")
normalized_nucleosome <- binary_group_normalize(
  genomic_reference, 
  nucleosome_condition,
  n_refs = 500,
  low_range = c(0.1, 0.4),
  high_range = c(0.6, 0.9)
)

# Apply Cross Normalization to before/after conditions
cat("Applying Cross Normalization...\n")
cross_norm_results <- cross_normalize(
  before_condition,
  after_condition,
  n_refs = 500
)

############################################################
# Visualize Results
############################################################

# Create plotting function
plot_signals <- function(region = 1:500) {
  par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
  
  # Plot 1: Raw signals showing probe effect
  plot(region, genomic_reference[region], type = "l", col = "black", lwd = 2,
       main = "Raw Signals: Probe Effect in Genomic Reference",
       xlab = "Probe Position", ylab = "Signal Intensity")
  lines(region, nucleosome_condition[region], col = "red", lwd = 2)
  legend("topright", c("Genomic Reference", "Nucleosome Enriched"), 
         col = c("black", "red"), lwd = 2, bty = "n")
  
  # Plot 2: Before vs. After normalized with Group Normalization
  plot(region, normalized_nucleosome[region], type = "l", col = "blue", lwd = 2,
       main = "Group Normalized Nucleosome Signal",
       xlab = "Probe Position", ylab = "Normalized Signal")
  abline(h = 0, lty = 2, col = "gray")
  abline(h = 1, lty = 2, col = "gray")
  text(min(region) + 10, 0.1, "Unbound", pos = 4)
  text(min(region) + 10, 0.9, "Bound", pos = 4)
  
  # Plot 3: Cross normalization of before/after conditions
  plot(region, cross_norm_results$a_norm[region], type = "l", col = "green", lwd = 2,
       main = "Cross Normalization: Differential Binding",
       xlab = "Probe Position", ylab = "Cross Normalized Signal")
  lines(region, cross_norm_results$b_norm[region], col = "purple", lwd = 2, lty = 2)
  abline(h = 0, lty = 2, col = "gray")
  legend("topright", c("Before vs After", "After vs Before"), 
         col = c("green", "purple"), lwd = 2, lty = c(1, 2), bty = "n")
  
  par(mfrow = c(1, 1))
}

# Plot the signals for a region
plot_signals(1:1000)

############################################################
# Comparison with Other Normalization Methods
############################################################

# Quantile Normalization
quantile_norm <- function(data_matrix) {
  return(normalize.quantiles(data_matrix))
}

# Simple MAS5-like Normalization
mas5_normalize <- function(data_matrix) {
  # Scale to a target median
  target_median <- 500
  scale_factors <- target_median / apply(data_matrix, 2, median, na.rm = TRUE)
  scaled <- sweep(data_matrix, 2, scale_factors, "*")
  return(scaled)
}

# Apply alternative normalizations
cat("Applying alternative normalization methods for comparison...\n")
data_matrix <- cbind(nucleosome_condition, genomic_reference)
quant_norm_result <- quantile_norm(data_matrix)[, 1]
mas5_norm_result <- mas5_normalize(data_matrix)[, 1]

# Compare normalization methods
compare_methods <- function(region = 1:500) {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  # Plot original data
  plot(region, genomic_reference[region], type = "l", col = "black", lwd = 2,
       main = "Raw Data",
       xlab = "Probe Position", ylab = "Signal Intensity")
  lines(region, nucleosome_condition[region], col = "red", lwd = 2)
  legend("topright", c("Genomic Reference", "Nucleosome Enriched"), 
         col = c("black", "red"), lwd = 2, bty = "n")
  
  # Plot Group Normalization
  plot(region, normalized_nucleosome[region], type = "l", col = "blue", lwd = 2,
       main = "Group Normalization",
       xlab = "Probe Position", ylab = "Normalized Signal")
  abline(h = c(0, 1), lty = 2, col = "gray")
  
  # Plot Quantile Normalization
  plot(region, quant_norm_result[region], type = "l", col = "green", lwd = 2,
       main = "Quantile Normalization",
       xlab = "Probe Position", ylab = "Normalized Signal")
  
  # Plot MAS5 Normalization
  plot(region, mas5_norm_result[region], type = "l", col = "purple", lwd = 2,
       main = "MAS5-like Normalization",
       xlab = "Probe Position", ylab = "Normalized Signal")
  
  par(mfrow = c(1, 1))
}

# Compare different normalization methods
compare_methods(1:1000)

############################################################
# Evaluate Signal Quality
############################################################

# Create a replicate of the nucleosome condition with noise
nucleosome_replicate <- probe_efficiency * nucleosome_signal + probe_background + rnorm(n_probes, 0, 0.1)

# Function to calculate signal quality
calculate_signal_quality <- function(signal_a, signal_b, replicate_b, window = 100, top_percent = 0.02) {
  # Calculate difference between conditions
  diff_ab <- signal_a - signal_b
  
  # Apply smoothing
  smoothed_diff <- stats::filter(diff_ab, rep(1/window, window), sides = 2)
  # Handle NAs at the edges
  na_idx <- which(is.na(smoothed_diff))
  if (length(na_idx) > 0) {
    smoothed_diff[na_idx] <- diff_ab[na_idx]
  }
  
  # Find top differentially bound probes
  n_top <- ceiling(top_percent * length(smoothed_diff))
  top_indices <- order(abs(smoothed_diff), decreasing = TRUE)[1:n_top]
  
  # Calculate signal power (S)
  signal_power <- mean((signal_a[top_indices] - signal_b[top_indices])^2)
  
  # Calculate noise power (N)
  noise_power <- mean((signal_b - replicate_b)^2)
  
  # Calculate Signal Quality in dB
  signal_quality_db <- 10 * log10(signal_power / noise_power)
  
  return(signal_quality_db)
}

# Calculate signal quality for different normalization methods
sq_gn <- calculate_signal_quality(normalized_nucleosome, nucleosome_condition, nucleosome_replicate)
sq_quant <- calculate_signal_quality(quant_norm_result, nucleosome_condition, nucleosome_replicate)
sq_mas5 <- calculate_signal_quality(mas5_norm_result, nucleosome_condition, nucleosome_replicate)

# Display signal quality results
cat("\nSignal Quality Comparison (dB):\n")
cat("Group Normalization:", round(sq_gn, 2), "dB\n")
cat("Quantile Normalization:", round(sq_quant, 2), "dB\n")
cat("MAS5-like Normalization:", round(sq_mas5, 2), "dB\n")

# Create barplot of signal quality
barplot(c(sq_mas5, sq_quant, sq_gn), 
        names.arg = c("MAS5", "Quantile", "Group Norm"),
        col = c("purple", "green", "blue"),
        main = "Signal Quality Comparison",
        ylab = "Signal Quality (dB)")

############################################################
# Visualize Joint Distribution (Figure 5A from paper)
############################################################

# Plot joint distribution before and after normalization
plot_joint_distribution <- function(sample_size = 1000) {
  # Sample points for plotting
  idx <- sample(1:n_probes, sample_size)
  
  # Create a 1x2 layout
  par(mfrow = c(1, 2))
  
  # Plot before normalization
  plot(genomic_reference[idx], nucleosome_condition[idx], 
       pch = 20, cex = 0.5, col = "darkblue", 
       xlab = "Genomic Reference", ylab = "Nucleosome Enriched",
       main = "Before Normalization")
  
  # Calculate correlation
  cor_before <- cor(genomic_reference, nucleosome_condition)
  text(min(genomic_reference) + 0.1 * diff(range(genomic_reference)), 
       max(nucleosome_condition) - 0.1 * diff(range(nucleosome_condition)),
       paste("Correlation =", round(cor_before, 3)),
       pos = 4)
  
  # Plot after normalization
  plot(genomic_reference[idx], normalized_nucleosome[idx], 
       pch = 20, cex = 0.5, col = "darkred", 
       xlab = "Genomic Reference", ylab = "Normalized Signal",
       main = "After Normalization")
  
  # Calculate correlation
  cor_after <- cor(genomic_reference, normalized_nucleosome)
  text(min(genomic_reference) + 0.1 * diff(range(genomic_reference)), 
       max(normalized_nucleosome) - 0.1 * diff(range(normalized_nucleosome)),
       paste("Correlation =", round(cor_after, 3)),
       pos = 4)
  
  # Reset layout
  par(mfrow = c(1, 1))
}

# Visualize joint distribution
plot_joint_distribution(2000)

cat("\nGroup Normalization demonstration complete.\n")

# The next section contains code to reproduce specific figures from the paper
# Using the real data would require access to the original datasets
# This example demonstrates the concepts and methods using simulated data
