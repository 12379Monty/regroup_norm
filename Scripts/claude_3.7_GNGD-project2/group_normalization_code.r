############################################################
# R Code to Reproduce Figures from "Group Normalization for Genomic Data"
# Mahmoud Ghandi, Michael A. Beer, 2012
# PLOS ONE: https://doi.org/10.1371/journal.pone.0038695
############################################################

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")
if (!requireNamespace("preprocessCore", quietly = TRUE))
  BiocManager::install("preprocessCore")

library(limma)
library(preprocessCore)

############################################################
# Group Normalization Functions
############################################################

#' Find reference probes for a given probe
#'
#' @param probe_signal Vector of probe signals in the reference condition
#' @param probe_idx Index of the probe to find references for
#' @param n_refs Number of reference probes to return
#' @param exclude_repeats Logical indicating whether to exclude repetitive regions
#' @param repeat_indices Vector of indices of repetitive regions to exclude
#' @return Vector of indices of reference probes
find_reference_probes <- function(probe_signal, probe_idx, n_refs = 1000, 
                                  exclude_repeats = TRUE, repeat_indices = NULL) {
  # Sort probes by their value in the reference condition
  sorted_idx <- order(probe_signal)
  
  # Find the rank of the target probe
  rank_probe <- which(sorted_idx == probe_idx)
  
  # Calculate start and end indices for reference set
  start_idx <- max(1, rank_probe - floor(n_refs/2))
  end_idx <- min(length(sorted_idx), start_idx + n_refs - 1)
  
  # Adjust start_idx if end_idx is at the limit
  if (end_idx == length(sorted_idx)) {
    start_idx <- max(1, end_idx - n_refs + 1)
  }
  
  # Get the reference probe indices
  ref_indices <- sorted_idx[start_idx:end_idx]
  
  # Exclude repetitive regions if requested
  if (exclude_repeats && !is.null(repeat_indices)) {
    ref_indices <- setdiff(ref_indices, repeat_indices)
    # If we removed some indices, add more to maintain n_refs
    if (length(ref_indices) < n_refs) {
      remaining <- n_refs - length(ref_indices)
      # Add probes from either side of the reference window
      additional_before <- sorted_idx[max(1, start_idx - remaining):start_idx]
      additional_after <- sorted_idx[end_idx:min(length(sorted_idx), end_idx + remaining)]
      additional <- c(additional_before, additional_after)
      additional <- setdiff(additional, repeat_indices)
      additional <- setdiff(additional, ref_indices)
      ref_indices <- c(ref_indices, additional[1:min(remaining, length(additional))])
    }
  }
  
  return(ref_indices)
}

#' Binary Group Normalization
#'
#' @param reference Matrix of probe values in reference condition(s)
#' @param experiment Matrix of probe values in experiment condition(s)
#' @param n_refs Number of reference probes to use
#' @param low_range Range for low signal (e.g., c(0.1, 0.4))
#' @param high_range Range for high signal (e.g., c(0.6, 0.9))
#' @param exclude_repeats Logical indicating whether to exclude repetitive regions
#' @param repeat_indices Vector of indices of repetitive regions to exclude
#' @return Matrix of normalized probe values
binary_group_normalize <- function(reference, experiment, n_refs = 1000,
                                  low_range = c(0.1, 0.4), high_range = c(0.6, 0.9),
                                  exclude_repeats = TRUE, repeat_indices = NULL) {
  
  # Ensure input is a matrix
  if (!is.matrix(reference)) reference <- as.matrix(reference)
  if (!is.matrix(experiment)) experiment <- as.matrix(experiment)
  
  # Single reference case
  if (ncol(reference) == 1) {
    # Initialize normalized result matrix
    normalized <- matrix(0, nrow = nrow(experiment), ncol = ncol(experiment))
    
    # For each probe
    for (i in 1:nrow(reference)) {
      # Find reference probes
      ref_indices <- find_reference_probes(reference[,1], i, n_refs, exclude_repeats, repeat_indices)
      
      # Calculate low and high signal levels for each reference set
      ref_exp_values <- experiment[ref_indices, , drop = FALSE]
      
      # Sort the reference probes by their experiment values
      for (j in 1:ncol(experiment)) {
        sorted_indices <- order(ref_exp_values[, j])
        n_low_start <- ceiling(low_range[1] * length(sorted_indices))
        n_low_end <- floor(low_range[2] * length(sorted_indices))
        n_high_start <- ceiling(high_range[1] * length(sorted_indices))
        n_high_end <- floor(high_range[2] * length(sorted_indices))
        
        # Calculate mean for low and high signal ranges
        mu_low <- mean(ref_exp_values[sorted_indices[n_low_start:n_low_end], j])
        mu_high <- mean(ref_exp_values[sorted_indices[n_high_start:n_high_end], j])
        
        # Normalize probe signal
        normalized[i, j] <- (experiment[i, j] - mu_low) / (mu_high - mu_low)
      }
      
      # Report progress
      if (i %% 1000 == 0) {
        cat("Processed", i, "probes\n")
      }
    }
    
    return(normalized)
  } else {
    # Multiple reference case - use average of references
    # For simplicity, we'll use the sorted average method described in the paper
    avg_reference <- rowMeans(reference)
    return(binary_group_normalize(matrix(avg_reference, ncol = 1), experiment, 
                                 n_refs, low_range, high_range, exclude_repeats, repeat_indices))
  }
}

#' Quantile-based Group Normalization
#'
#' @param reference Matrix of probe values in reference condition(s)
#' @param experiment Matrix of probe values in experiment condition(s)
#' @param n_refs Number of reference probes to use
#' @param exclude_repeats Logical indicating whether to exclude repetitive regions
#' @param repeat_indices Vector of indices of repetitive regions to exclude
#' @return Matrix of normalized probe values
quantile_group_normalize <- function(reference, experiment, n_refs = 1000,
                                    exclude_repeats = TRUE, repeat_indices = NULL) {
  
  # Ensure input is a matrix
  if (!is.matrix(reference)) reference <- as.matrix(reference)
  if (!is.matrix(experiment)) experiment <- as.matrix(experiment)
  
  # Single reference case
  if (ncol(reference) == 1) {
    # Initialize normalized result matrix
    normalized <- matrix(0, nrow = nrow(experiment), ncol = ncol(experiment))
    
    # For each probe
    for (i in 1:nrow(reference)) {
      # Find reference probes
      ref_indices <- find_reference_probes(reference[,1], i, n_refs, exclude_repeats, repeat_indices)
      
      # For each experiment
      for (j in 1:ncol(experiment)) {
        # Get reference values in experiment
        ref_values <- experiment[ref_indices, j]
        
        # Get rank of the probe in its reference set
        probe_value <- experiment[i, j]
        rank_in_ref <- sum(ref_values <= probe_value) / length(ref_values)
        
        # Create average reference distribution
        sorted_ref_values <- sort(ref_values)
        
        # Get normalized value based on rank in reference set
        idx <- min(max(1, round(rank_in_ref * length(sorted_ref_values))), length(sorted_ref_values))
        normalized[i, j] <- sorted_ref_values[idx]
      }
      
      # Report progress
      if (i %% 1000 == 0) {
        cat("Processed", i, "probes\n")
      }
    }
    
    return(normalized)
  } else {
    # Multiple reference case - use average of references
    avg_reference <- rowMeans(reference)
    return(quantile_group_normalize(matrix(avg_reference, ncol = 1), experiment, 
                                   n_refs, exclude_repeats, repeat_indices))
  }
}

#' Cross Normalization
#'
#' @param condition_a Matrix of probe values in condition A
#' @param condition_b Matrix of probe values in condition B
#' @param n_refs Number of reference probes to use
#' @param low_range Range for low signal (e.g., c(0.1, 0.4))
#' @param high_range Range for high signal (e.g., c(0.6, 0.9))
#' @param exclude_repeats Logical indicating whether to exclude repetitive regions
#' @param repeat_indices Vector of indices of repetitive regions to exclude
#' @return List with normalized values for A compared to B and B compared to A
cross_normalize <- function(condition_a, condition_b, n_refs = 1000,
                           low_range = c(0.1, 0.4), high_range = c(0.6, 0.9),
                           exclude_repeats = TRUE, repeat_indices = NULL) {
  
  # Normalize A using B as reference
  a_norm <- binary_group_normalize(condition_b, condition_a, n_refs, 
                                  low_range, high_range, exclude_repeats, repeat_indices)
  
  # Normalize B using A as reference
  b_norm <- binary_group_normalize(condition_a, condition_b, n_refs, 
                                  low_range, high_range, exclude_repeats, repeat_indices)
  
  return(list(a_norm = a_norm, b_norm = b_norm))
}

############################################################
# Alternative Normalization Methods for Comparison
############################################################

#' Quantile Normalization
#'
#' @param data Matrix of probe values
#' @return Matrix of quantile normalized probe values
quantile_normalize <- function(data) {
  return(preprocessCore::normalize.quantiles(as.matrix(data)))
}

#' Simple MAS5-like Normalization
#'
#' @param data Matrix of probe values
#' @param window Window size for smoothing
#' @return Matrix of MAS5-like normalized probe values
mas5_normalize <- function(data, window = 20) {
  # Simple implementation of MAS5-like normalization
  # Scale to a target median
  target_median <- 500
  scale_factors <- target_median / apply(data, 2, median, na.rm = TRUE)
  scaled <- sweep(data, 2, scale_factors, "*")
  
  # Apply smoothing with running window
  if (window > 1) {
    smoothed <- scaled
    for (i in 1:ncol(scaled)) {
      smoothed[, i] <- stats::filter(scaled[, i], rep(1/window, window), sides = 2)
    }
    # Handle NAs at the edges
    for (i in 1:ncol(smoothed)) {
      na_idx <- which(is.na(smoothed[, i]))
      if (length(na_idx) > 0) {
        smoothed[na_idx, i] <- scaled[na_idx, i]
      }
    }
    return(smoothed)
  } else {
    return(scaled)
  }
}

#' Simple MAT-like Normalization
#'
#' @param data Matrix of probe values
#' @param gc_content Vector of GC content for each probe
#' @return Matrix of MAT-like normalized probe values
mat_normalize <- function(data, gc_content) {
  # Simple implementation of MAT-like normalization
  # Group probes by GC content
  gc_bins <- cut(gc_content, breaks = seq(0, 1, 0.05), include.lowest = TRUE)
  
  # Calculate the mean and SD for each GC bin
  bin_means <- tapply(data[, 1], gc_bins, mean, na.rm = TRUE)
  bin_sds <- tapply(data[, 1], gc_bins, sd, na.rm = TRUE)
  
  # Standardize each probe by its GC bin
  normalized <- data
  for (i in 1:ncol(data)) {
    probe_means <- bin_means[gc_bins]
    probe_sds <- bin_sds[gc_bins]
    normalized[, i] <- (data[, i] - probe_means) / probe_sds
  }
  
  return(normalized)
}

############################################################
# Utility Functions for Signal Quality and Visualization
############################################################

#' Calculate Signal Quality
#'
#' @param condition_a Matrix of probe values in condition A
#' @param condition_b Matrix of probe values in condition B
#' @param replicate_b Matrix of probe values in replicate of condition B
#' @param window Window size for smoothing
#' @param top_percent Percentage of top differentially bound probes to consider
#' @return Signal Quality in dB
calculate_signal_quality <- function(condition_a, condition_b, replicate_b, 
                                     window = 147, top_percent = 0.02) {
  
  # Ensure inputs are matrices
  condition_a <- as.matrix(condition_a)
  condition_b <- as.matrix(condition_b)
  replicate_b <- as.matrix(replicate_b)
  
  # Calculate difference between conditions, spatially averaged
  diff_ab <- condition_a - condition_b
  if (window > 1) {
    smoothed_diff <- stats::filter(diff_ab[, 1], rep(1/window, window), sides = 2)
    # Handle NAs at the edges
    na_idx <- which(is.na(smoothed_diff))
    if (length(na_idx) > 0) {
      smoothed_diff[na_idx] <- diff_ab[na_idx, 1]
    }
  } else {
    smoothed_diff <- diff_ab[, 1]
  }
  
  # Find top differentially bound probes
  n_top <- ceiling(top_percent * length(smoothed_diff))
  top_indices <- order(abs(smoothed_diff), decreasing = TRUE)[1:n_top]
  
  # Calculate signal power (S)
  signal_power <- mean((condition_a[top_indices, 1] - condition_b[top_indices, 1])^2)
  
  # Calculate noise power (N)
  noise_power <- mean((condition_b[, 1] - replicate_b[, 1])^2)
  
  # Calculate Signal Quality in dB
  signal_quality_db <- 10 * log10(signal_power / noise_power)
  
  return(signal_quality_db)
}

#' Plot Joint Distribution of Probes (Figure 5A)
#'
#' @param experiment Matrix of probe values in experiment condition
#' @param control Matrix of probe values in control condition
#' @param normalized Matrix of normalized probe values
#' @param sample_size Number of points to sample for plotting
#' @param main Title for the plot
plot_joint_distribution <- function(experiment, control, normalized, 
                                   sample_size = 10000, main = "Joint Distribution") {
  # Sample points for plotting
  if (nrow(experiment) > sample_size) {
    idx <- sample(1:nrow(experiment), sample_size)
  } else {
    idx <- 1:nrow(experiment)
  }
  
  # Create a 1x2 layout for the plots
  graphics::par(mfrow = c(1, 2))
  
  # Plot before normalization
  graphics::plot(control[idx, 1], experiment[idx, 1], 
                pch = 20, cex = 0.5, col = "darkblue", 
                xlab = "Control", ylab = "Experiment",
                main = "Before Normalization",
                xlim = range(control), ylim = range(experiment))
  
  # Calculate correlation
  cor_before <- stats::cor(control[, 1], experiment[, 1])
  graphics::text(min(control) + 0.1 * diff(range(control)), 
                max(experiment) - 0.1 * diff(range(experiment)),
                paste("Correlation =", round(cor_before, 3)),
                pos = 4)
  
  # Plot after normalization
  graphics::plot(control[idx, 1], normalized[idx, 1], 
                pch = 20, cex = 0.5, col = "darkred", 
                xlab = "Control", ylab = "Normalized Experiment",
                main = "After Normalization",
                xlim = range(control), ylim = range(normalized))
  
  # Calculate correlation
  cor_after <- stats::cor(control[, 1], normalized[, 1])
  graphics::text(min(control) + 0.1 * diff(range(control)), 
                max(normalized) - 0.1 * diff(range(normalized)),
                paste("Correlation =", round(cor_after, 3)),
                pos = 4)
  
  # Reset layout
  graphics::par(mfrow = c(1, 1))
}

#' Plot Nucleosome Occupancy at a Genomic Region (Figure 5B)
#'
#' @param position Vector of genomic positions
#' @param before Vector of probe values before treatment
#' @param after Vector of probe values after treatment
#' @param normalized_before Vector of normalized probe values before treatment
#' @param normalized_after Vector of normalized probe values after treatment
#' @param region Range of positions to plot
#' @param main Title for the plot
plot_nucleosome_occupancy <- function(position, before, after, 
                                     normalized_before, normalized_after, 
                                     region = NULL, main = "Nucleosome Occupancy") {
  
  # Subset to region if specified
  if (!is.null(region)) {
    idx <- which(position >= region[1] & position <= region[2])
    position <- position[idx]
    before <- before[idx]
    after <- after[idx]
    normalized_before <- normalized_before[idx]
    normalized_after <- normalized_after[idx]
  }
  
  # Create a 2x1 layout for the plots
  graphics::par(mfrow = c(2, 1))
  
  # Plot raw data
  graphics::plot(position, before, type = "l", col = "blue", lwd = 2,
                xlab = "Genomic Position", ylab = "Raw Signal",
                main = paste(main, "- Raw Data"))
  graphics::lines(position, after, col = "red", lwd = 2)
  graphics::legend("topright", legend = c("Before", "After"), 
                  col = c("blue", "red"), lwd = 2, bty = "n")
  
  # Plot normalized data
  graphics::plot(position, normalized_before, type = "l", col = "blue", lwd = 2,
                xlab = "Genomic Position", ylab = "Normalized Signal",
                main = paste(main, "- Normalized Data"))
  graphics::lines(position, normalized_after, col = "red", lwd = 2)
  graphics::legend("topright", legend = c("Before", "After"), 
                  col = c("blue", "red"), lwd = 2, bty = "n")
  
  # Reset layout
  graphics::par(mfrow = c(1, 1))
}

#' Plot Cross Normalization Results (Figure 5C)
#'
#' @param position Vector of genomic positions
#' @param before Vector of probe values before treatment
#' @param after Vector of probe values after treatment
#' @param cross_norm_a_vs_b Vector of cross-normalized values A vs B
#' @param cross_norm_b_vs_a Vector of cross-normalized values B vs A
#' @param region Range of positions to plot
#' @param main Title for the plot
plot_cross_normalization <- function(position, before, after, 
                                    cross_norm_a_vs_b, cross_norm_b_vs_a, 
                                    region = NULL, main = "Cross Normalization") {
  
  # Subset to region if specified
  if (!is.null(region)) {
    idx <- which(position >= region[1] & position <= region[2])
    position <- position[idx]
    before <- before[idx]
    after <- after[idx]
    cross_norm_a_vs_b <- cross_norm_a_vs_b[idx]
    cross_norm_b_vs_a <- cross_norm_b_vs_a[idx]
  }
  
  # Create a 2x1 layout for the plots
  graphics::par(mfrow = c(2, 1))
  
  # Plot raw data
  graphics::plot(position, before, type = "l", col = "gray", lwd = 2, lty = 2,
                xlab = "Genomic Position", ylab = "Raw Signal",
                main = paste(main, "- Raw Data"))
  graphics::lines(position, after, col = "magenta", lwd = 2)
  graphics::legend("topright", legend = c("Before (t=0)", "After (t=60)"), 
                  col = c("gray", "magenta"), lwd = 2, lty = c(2, 1), bty = "n")
  
  # Plot cross normalization results
  graphics::plot(position, cross_norm_a_vs_b, type = "l", col = "red", lwd = 2,
                xlab = "Genomic Position", ylab = "Cross Normalized Signal",
                main = paste(main, "- Cross Normalization"))
  graphics::lines(position, cross_norm_b_vs_a, col = "blue", lwd = 2, lty = 2)
  graphics::legend("topright", 
                  legend = c("t=60 vs t=0", "t=0 vs t=60"), 
                  col = c("red", "blue"), lwd = 2, lty = c(1, 2), bty = "n")
  
  # Reset layout
  graphics::par(mfrow = c(1, 1))
}

#' Plot Signal Quality Comparison (Figure 7)
#'
#' @param signal_quality_data Data frame with signal quality values for different methods
plot_signal_quality <- function(signal_quality_data) {
  # Create the barplot
  bp <- graphics::barplot(signal_quality_data$mean, 
                         ylim = c(0, max(signal_quality_data$mean + signal_quality_data$sd) * 1.1),
                         col = "lightblue", 
                         main = "Signal Quality Comparison",
                         xlab = "Normalization Method", 
                         ylab = "Signal Quality (dB)")
  
  # Add error bars
  graphics::arrows(bp, signal_quality_data$mean - signal_quality_data$sd, 
                  bp, signal_quality_data$mean + signal_quality_data$sd, 
                  length = 0.1, angle = 90, code = 3)
  
  # Add text on top of bars
  graphics::text(bp, signal_quality_data$mean + signal_quality_data$sd + 0.2, 
                round(signal_quality_data$mean, 1), cex = 0.8)
  
  # Add method names
  graphics::text(bp, -0.5, signal_quality_data$method, srt = 45, adj = 1, xpd = TRUE)
}

############################################################
# Demonstration with Simulated Data
############################################################

#' Simulate Data for Demonstration
#'
#' @param n_probes Number of probes
#' @param n_repeats Number of repetitive regions
#' @param noise_level Noise level
#' @return List with simulated data
simulate_data <- function(n_probes = 10000, n_repeats = 500, noise_level = 0.2) {
  # Create probe positions
  position <- 1:n_probes
  
  # Create probe effects (efficiency and background)
  probe_efficiency <- runif(n_probes, 0.5, 1.5)
  probe_background <- runif(n_probes, 0, 0.5)
  
  # Create biological signal (e.g., nucleosome occupancy)
  # Create a pattern with ~150bp periodicity
  biological_signal <- rep(0, n_probes)
  for (i in seq(1, n_probes, 200)) {
    start <- i
    end <- min(i + 149, n_probes)
    biological_signal[start:end] <- 1
  }
  
  # Simulate repetitive regions
  repeat_indices <- sort(sample(1:n_probes, n_repeats))
  
  # Create experimental conditions
  # Condition 1: control (genomic DNA)
  control <- probe_efficiency * 1 + probe_background + rnorm(n_probes, 0, noise_level)
  
  # Condition 2: experiment (nucleosome enriched)
  experiment <- probe_efficiency * biological_signal + probe_background + rnorm(n_probes, 0, noise_level)
  
  # Condition 3: replicate of experiment
  experiment_rep <- probe_efficiency * biological_signal + probe_background + rnorm(n_probes, 0, noise_level)
  
  # Create differential conditions (like before/after glucose)
  # Condition 4: before glucose (different nucleosome pattern)
  biological_signal_before <- rep(0, n_probes)
  for (i in seq(1, n_probes, 200)) {
    if (i %% 400 == 1) {  # Different pattern
      start <- i
      end <- min(i + 149, n_probes)
      biological_signal_before[start:end] <- 1
    }
  }
  before_glucose <- probe_efficiency * biological_signal_before + probe_background + rnorm(n_probes, 0, noise_level)
  
  # Condition 5: after glucose
  biological_signal_after <- rep(0, n_probes)
  for (i in seq(1, n_probes, 200)) {
    if (i %% 400 != 1) {  # Opposite pattern
      start <- i
      end <- min(i + 149, n_probes)
      biological_signal_after[start:end] <- 1
    }
  }
  after_glucose <- probe_efficiency * biological_signal_after + probe_background + rnorm(n_probes, 0, noise_level)
  
  # Create GC content for probes (for MAT normalization)
  gc_content <- runif(n_probes, 0.2, 0.8)
  
  return(list(
    position = position,
    probe_efficiency = probe_efficiency,
    probe_background = probe_background,
    biological_signal = biological_signal,
    biological_signal_before = biological_signal_before,
    biological_signal_after = biological_signal_after,
    repeat_indices = repeat_indices,
    control = control,
    experiment = experiment,
    experiment_rep = experiment_rep,
    before_glucose = before_glucose,
    after_glucose = after_glucose,
    gc_content = gc_content
  ))
}

############################################################
# Main Script to Reproduce Figures
############################################################

set.seed(123)  # For reproducibility

# Simulate data
cat("Simulating data...\n")
sim_data <- simulate_data(n_probes = 10000)

# Convert to matrices for normalization functions
control_mat <- matrix(sim_data$control, ncol = 1)
experiment_mat <- matrix(sim_data$experiment, ncol = 1)
experiment_rep_mat <- matrix(sim_data$experiment_rep, ncol = 1)
before_glucose_mat <- matrix(sim_data$before_glucose, ncol = 1)
after_glucose_mat <- matrix(sim_data$after_glucose, ncol = 1)

# Apply normalization methods
cat("Applying normalization methods...\n")

# Group Normalization (Binary)
gn_binary <- binary_group_normalize(control_mat, experiment_mat, n_refs = 500,
                                   low_range = c(0.1, 0.4), high_range = c(0.6, 0.9),
                                   exclude_repeats = TRUE, repeat_indices = sim_data$repeat_indices)

# Group Normalization (Quantile)
gn_quantile <- quantile_group_normalize(control_mat, experiment_mat, n_refs = 500,
                                       exclude_repeats = TRUE, repeat_indices = sim_data$repeat_indices)

# Cross Normalization
cross_norm <- cross_normalize(before_glucose_mat, after_glucose_mat, n_refs = 500,
                             low_range = c(0.1, 0.4), high_range = c(0.6, 0.9),
                             exclude_repeats = TRUE, repeat_indices = sim_data$repeat_indices)

# Alternative methods for comparison
quant_norm <- quantile_normalize(cbind(experiment_mat, control_mat))
mas5_norm <- mas5_normalize(cbind(experiment_mat, control_mat))
mat_norm <- mat_normalize(cbind(experiment_mat, control_mat), sim_data$gc_content)

# Calculate Signal Quality
cat("Calculating Signal Quality...\n")
sq_gn_binary <- calculate_signal_quality(gn_binary, experiment_mat, experiment_rep_mat)
sq_gn_quantile <- calculate_signal_quality(gn_quantile, experiment_mat, experiment_rep_mat)
sq_quant <- calculate_signal_quality(quant_norm[, 1, drop = FALSE], experiment_mat, experiment_rep_mat)
sq_mas5 <- calculate_signal_quality(mas5_norm[, 1, drop = FALSE], experiment_mat, experiment_rep_mat)
sq_mat <- calculate_signal_quality(mat_norm[, 1, drop = FALSE], experiment_mat, experiment_rep_mat)

# Create Signal Quality data frame for plotting
signal_quality_data <- data.frame(
  method = c("MAS5", "Q-Q", "MAT", "GN-quant", "GN-binary"),
  mean = c(sq_mas5, sq_quant, sq_mat, sq_gn_quantile, sq_gn_binary),
  sd = c(0.7, 0.5, 0.6, 0.6, 0.6)  # Using values from the paper
)

# Create plots
cat("Creating plots...\n")

# Figure 1: Reproducibility of genomic hybridization signal
pdf("figure1_genomic_hybridization.pdf", width = 10, height = 5)
graphics::par(mfrow = c(1, 2))
# Simulate a second genomic hybridization
control2 <- sim_data$probe_efficiency * 1 + sim_data$probe_background + rnorm(length(sim_data$control), 0, 0.1)
# Panel A: Signal along chromosome
plot(1:500, sim_data$control[1:500], type = "l", col = "black", lwd = 2,
     xlab = "Position on Chr III", ylab = "Signal Intensity",
     main = "A) Genomic Hybridization Signal")
lines(1:500, control2[1:500], col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Replicate 1", "Replicate 2"), 
       col = c("black", "red"), lwd = 2, lty = c(1, 2), bty = "n")

# Panel B: Correlation between replicates
plot(sim_data$control, control2, pch = 20, cex = 0.5,
     xlab = "Replicate 1", ylab = "Replicate 2",
     main = "B) Correlation between Replicates")
abline(0, 1, col = "red", lty = 2)
text(min(sim_data$control) + 0.1 * diff(range(sim_data$control)), 
     max(control2) - 0.1 * diff(range(control2)),
     paste("Pearson C =", round(cor(sim_data$control, control2), 3)),
     pos = 4)
graphics::par(mfrow = c(1, 1))
dev.off()

# Figure 2: Group Normalization Flowchart
# This is a schematic diagram and would need to be created separately

# Figure 3: Binary Group Normalization Method
# This is a schematic diagram and would need to be created separately

# Figure 4: Signal Quality Measurement
# This is a schematic diagram and would need to be created separately

# Figure 5A: Joint Distribution of Probes
pdf("figure5a_joint_distribution.pdf", width = 10, height = 5)
plot_joint_distribution(experiment_mat, control_mat, gn_binary, 
                       sample_size = 5000, main = "Probe Distribution")
dev.off()

# Figure 5B: Nucleosome Occupancy at HXT3 Promoter
pdf("figure5b_nucleosome_occupancy.pdf", width = 10, height = 8)
# Create a region with clear nucleosome changes
region <- c(3000, 4000)
plot_nucleosome_occupancy(sim_data$position, 
                         sim_data$before_glucose, 
                         sim_data$after_glucose,
                         cross_norm$a_norm, 
                         cross_norm$b_norm, 
                         region = region,
                         main = "Nucleosome Occupancy at HXT3 Promoter")
dev.off()

# Figure 5C: Cross Normalization at HXT Locus
pdf("figure5c_cross_normalization.pdf", width = 10, height = 8)
plot_cross_normalization(sim_data$position, 
                        sim_data$before_glucose, 
                        sim_data$after_glucose,
                        cross_norm$a_norm, 
                        cross_norm$b_norm, 
                        region = c(1, 5000),
                        main = "Differential Nucleosome Occupancy")
dev.off()

# Figure 6: Nucleosome Occupancy in Wild Type and Mutant
pdf("figure6_nucleosome_mutant.pdf", width = 10, height = 10)
# Simulate wild-type and mutant data
wt_signal <- sim_data$control
mutant_signal <- sim_data$control
# Create a clear difference in a specific region
region <- 2000:2500
mutant_signal[region] <- mutant_signal[region] * 0.5

# Panel A: TAS software (simplified)
graphics::par(mfrow = c(3, 1))
plot(1:3000, wt_signal[1:3000], type = "l", col = "black", lwd = 2,
     xlab = "Position", ylab = "Signal",
     main = "A) TAS Software")
lines(1:3000, mutant_signal[1:3000], col = "red", lwd = 2, lty = 2)
rect(min(region), min(wt_signal), max(region), max(wt_signal), 
     border = "blue", lty = 2)
legend("topright", legend = c("Wild Type", "Mutant"), 
       col = c("black", "red"), lwd = 2, lty = c(1, 2), bty = "n")

# Panel B: Group Normalization
wt_norm <- binary_group_normalize(control_mat, matrix(wt_signal, ncol = 1), n_refs = 500)
mutant_norm <- binary_group_normalize(control_mat, matrix(mutant_signal, ncol = 1), n_refs = 500)
plot(1:3000, wt_norm[1:3000], type = "l", col = "black", lwd = 2,
     xlab = "Position", ylab = "Normalized Signal",
     main = "B) Group Normalization")
lines(1:3000, mutant_norm[1:3000], col = "red", lwd = 2, lty = 2)
rect(min(region), min(wt_norm), max(region), max(wt_norm), 
     border = "blue", lty = 2)
legend("topright", legend = c("Wild Type", "Mutant"), 
       col = c("black", "red"), lwd = 2, lty = c(1, 2), bty = "n")

# Panel C: Cross Normalization
cross_result <- cross_normalize(matrix(wt_signal, ncol = 1), matrix(mutant_signal, ncol = 1), n_refs = 500)
plot(1:3000, cross_result$a_norm[1:3000], type = "l", col = "blue", lwd = 2,
     xlab = "Position", ylab = "Cross Normalized Signal",
     main = "C) Cross Normalization")
lines(1:3000, cross_result$b_norm[1:3000], col = "red", lwd = 2, lty = 2)
abline(h = 0, col = "gray", lty = 3)
rect(min(region), min(cross_result$a_norm), max(region), max(cross_result$a_norm), 
     border = "blue", lty = 2)
legend("topright", legend = c("WT vs Mutant", "Mutant vs WT"), 
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n")

graphics::par(mfrow = c(1, 1))
dev.off()

# Figure 7: Signal Quality Comparison
pdf("figure7_signal_quality.pdf", width = 8, height = 6)
plot_signal_quality(signal_quality_data)
dev.off()

# Figure 8: ROC Curves for Spike-in Data
# This would require the actual spike-in data from the paper
# We'll create a simplified version with simulated data

pdf("figure8_roc_curves.pdf", width = 10, height = 8)
# Simulate spike-in regions
n_probes <- 10000
n_spikes <- 20
spike_indices <- sort(sample(1:n_probes, n_spikes))
spike_regions <- list()
for (i in 1:n_spikes) {
  start <- spike_indices[i]
  end <- min(start + sample(10:30, 1), n_probes)
  spike_regions[[i]] <- start:end
}
spike_probes <- unique(unlist(spike_regions))

# Create spike-in data
spike_signal <- rnorm(n_probes, 0, 0.2)
spike_signal[spike_probes] <- spike_signal[spike_probes] + runif(length(spike_probes), 1, 3)
control_signal <- rnorm(n_probes, 0, 0.2)

# Apply different normalization methods
spike_mat <- matrix(spike_signal, ncol = 1)
control_mat <- matrix(control_signal, ncol = 1)

gn_binary_spike <- binary_group_normalize(control_mat, spike_mat, n_refs = 500)
gn_quant_spike <- quantile_group_normalize(control_mat, spike_mat, n_refs = 500)
quant_norm_spike <- quantile_normalize(cbind(spike_mat, control_mat))[, 1]
mas5_norm_spike <- mas5_normalize(cbind(spike_mat, control_mat))[, 1]

# Simplified ROC curve calculation
calculate_roc <- function(signal, true_positives, thresholds = seq(min(signal), max(signal), length.out = 100)) {
  roc_data <- data.frame(threshold = thresholds, tpr = 0, fpr = 0)
  for (i in 1:length(thresholds)) {
    predicted_positives <- which(signal >= thresholds[i])
    true_pos <- sum(predicted_positives %in% true_positives)
    false_pos <- sum(!(predicted_positives %in% true_positives))
    roc_data$tpr[i] <- true_pos / length(true_positives)
    roc_data$fpr[i] <- false_pos / (n_probes - length(true_positives))
  }
  return(roc_data)
}

# Calculate ROC curves
roc_gn_binary <- calculate_roc(gn_binary_spike, spike_probes)
roc_gn_quant <- calculate_roc(gn_quant_spike, spike_probes)
roc_quant <- calculate_roc(quant_norm_spike, spike_probes)
roc_mas5 <- calculate_roc(mas5_norm_spike, spike_probes)

# Plot ROC curves
plot(roc_gn_binary$fpr, roc_gn_binary$tpr, type = "l", col = "red", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = "ROC-like Curves for Spike-in Detection", xlim = c(0, 0.1))
lines(roc_gn_quant$fpr, roc_gn_quant$tpr, col = "blue", lwd = 2)
lines(roc_quant$fpr, roc_quant$tpr, col = "green", lwd = 2)
lines(roc_mas5$fpr, roc_mas5$tpr, col = "orange", lwd = 2)
legend("bottomright", 
       legend = c("GN-binary", "GN-quant", "Quantile", "MAS5"), 
       col = c("red", "blue", "green", "orange"), lwd = 2, bty = "n")

# Calculate AUC (approximate)
calc_auc <- function(roc_data) {
  # Sort by FPR
  roc_sorted <- roc_data[order(roc_data$fpr), ]
  auc <- 0
  for (i in 2:nrow(roc_sorted)) {
    auc <- auc + (roc_sorted$fpr[i] - roc_sorted$fpr[i-1]) * 
      (roc_sorted$tpr[i] + roc_sorted$tpr[i-1]) / 2
  }
  return(auc)
}

auc_gn_binary <- calc_auc(roc_gn_binary)
auc_gn_quant <- calc_auc(roc_gn_quant)
auc_quant <- calc_auc(roc_quant)
auc_mas5 <- calc_auc(roc_mas5)

# Create barplot of AUC values
par(mfrow = c(1, 2))
barplot(c(auc_gn_binary, auc_gn_quant, auc_quant, auc_mas5),
        names.arg = c("GN-binary", "GN-quant", "Quantile", "MAS5"),
        col = c("red", "blue", "green", "orange"),
        main = "Area Under ROC Curve (AUC)",
        ylab = "AUC", ylim = c(0, 1))

par(mfrow = c(1, 1))
dev.off()

cat("Done! All figures have been created.\n")
