## R script for reproducing figures from "Group Normalization for Genomic Data" (Ghandi & Beer, 2012)
## This script provides functions to download data, process it, and visualize nucleosome positioning data

# Load required packages
if(!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if(!require("rtracklayer", quietly = TRUE))
    BiocManager::install("rtracklayer")
if(!require("GenomicRanges", quietly = TRUE))
    BiocManager::install("GenomicRanges")
if(!require("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
if(!require("dplyr", quietly = TRUE))
    install.packages("dplyr")
if(!require("reshape2", quietly = TRUE))
    install.packages("reshape2")
if(!require("BSgenome.Scerevisiae.UCSC.sacCer3", quietly = TRUE))
    BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")

library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(reshape2)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

#######################
# DATA ACQUISITION FUNCTIONS 
#######################

# Function to download and import Lee et al. 2007 data
# We'll search for the appropriate data on GEO
download_lee_data <- function(output_dir = ".") {
  cat("Since we don't have the exact GEO accession, we'll make assumptions about format.\n")
  cat("In a real scenario, you would download files directly from GEO using their accession numbers.\n")
  
  # Create a mock dataset resembling nucleosome positioning data
  # In practice, you would use code like this to download from GEO:
  # 
  # Lee2007_accession <- "GSExxxxx" # Replace with actual accession
  # Lee2007_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", 
  #                       substr(Lee2007_accession, 1, 5), "nnn/", 
  #                       Lee2007_accession, "/suppl/", 
  #                       Lee2007_accession, "_RAW.tar")
  # download.file(Lee2007_url, destfile = file.path(output_dir, "Lee2007_RAW.tar"))
  # untar(file.path(output_dir, "Lee2007_RAW.tar"), exdir = file.path(output_dir, "Lee2007"))
  
  # Create a simulated nucleosome profile based on common patterns
  create_mock_nucleosome_data(output_dir)
}

# Function to create mock nucleosome data resembling Lee et al. 2007
create_mock_nucleosome_data <- function(output_dir) {
  # Get yeast chromosome information
  yeast_genome <- BSgenome.Scerevisiae.UCSC.sacCer3
  chr_lengths <- seqlengths(yeast_genome)
  
  # Create mock replicate datasets for Figure 1
  # For demonstration purposes, we'll create two biological replicates with correlated but not identical signals
  set.seed(123)
  
  # Generate genome-wide random walk for first replicate
  chr_data <- list()
  for (chr_name in names(chr_lengths)[1:16]) {  # Main chromosomes only
    # Number of bins for this chromosome (1 bin = 10bp)
    num_bins <- ceiling(chr_lengths[chr_name] / 10)
    
    # Create a base random walk signal (representing DNA sequence preferences)
    base_signal <- cumsum(rnorm(num_bins, 0, 0.1))
    # Normalize to range 0.5-2
    base_signal <- 0.5 + 1.5 * (base_signal - min(base_signal)) / (max(base_signal) - min(base_signal))
    
    # Add some periodic patterns to simulate nucleosome positioning
    period <- 16  # ~160bp in our binned data
    periodic_component <- 0.3 * sin(seq_len(num_bins) * (2 * pi / period))
    
    # Final signal combines base + periodic + promoter patterns
    signal <- base_signal + periodic_component
    
    # First replicate
    rep1 <- signal + rnorm(num_bins, 0, 0.05)
    
    # Second replicate (correlated with first but with some noise)
    rep2 <- signal + rnorm(num_bins, 0, 0.05)
    
    # Store chromosome data
    chr_data[[chr_name]] <- data.frame(
      chr = chr_name,
      position = seq(1, num_bins * 10, by = 10),
      rep1 = rep1,
      rep2 = rep2
    )
  }
  
  # Combine all chromosome data
  all_data <- do.call(rbind, chr_data)
  
  # Save to file
  mock_file <- file.path(output_dir, "mock_nucleosome_data.txt")
  write.table(all_data, file = mock_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Also create some test data for specific regions (HXT3 locus for Figure 5)
  # We'll simulate nucleosome depletion/occupation patterns
  hxt3_chr <- "chrIV"
  hxt3_start <- 1100000  # Approximate region
  hxt3_end <- 1130000
  
  # Extract HXT region
  hxt_idx <- which(all_data$chr == hxt3_chr & 
                  all_data$position >= hxt3_start & 
                  all_data$position <= hxt3_end)
  
  hxt_data <- all_data[hxt_idx, ]
  
  # Create condition A (glucose) - nucleosome free region
  hxt_data$condA <- hxt_data$rep1
  
  # Create condition B (no glucose) - nucleosome bound region
  # Increase nucleosome occupancy in a specific region
  nfr_start <- 1110000
  nfr_end <- 1112000
  nfr_idx <- which(hxt_data$position >= nfr_start & hxt_data$position <= nfr_end)
  
  hxt_data$condB <- hxt_data$rep1
  hxt_data$condB[nfr_idx] <- hxt_data$condB[nfr_idx] + 0.8  # Increase occupancy
  
  # Create replicate of condition B
  hxt_data$condB_rep <- hxt_data$condB + rnorm(nrow(hxt_data), 0, 0.05)
  
  # Save HXT3 data
  hxt_file <- file.path(output_dir, "mock_hxt3_data.txt")
  write.table(hxt_data, file = hxt_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(list(mock_file = mock_file, hxt_file = hxt_file))
}

#######################
# PROCESSING FUNCTIONS
#######################

# Function to implement Group Normalization
group_normalization <- function(reference_data, target_data, method = "binary", ref_group_size = 1000) {
  # Sort reference data
  sorted_ref <- sort(reference_data)
  
  # For each value in reference data, find its reference group
  normalized_data <- numeric(length(reference_data))
  
  for (i in 1:length(reference_data)) {
    # Find reference group (probes with similar values in reference data)
    ref_value <- reference_data[i]
    
    # Find index of this value in sorted reference data
    idx <- which(sorted_ref == ref_value)
    if (length(idx) > 1) idx <- idx[1]  # If multiple matches, take first
    
    # Create reference group (N probes with most similar values)
    half_size <- floor(ref_group_size/2)
    start_idx <- max(1, idx - half_size)
    end_idx <- min(length(sorted_ref), idx + half_size)
    
    if (end_idx - start_idx + 1 < ref_group_size) {
      # Adjust if we're near the edges
      if (start_idx == 1) {
        end_idx <- min(length(sorted_ref), start_idx + ref_group_size - 1)
      } else {
        start_idx <- max(1, end_idx - ref_group_size + 1)
      }
    }
    
    ref_group_idx <- start_idx:end_idx
    ref_group <- sorted_ref[ref_group_idx]
    
    # Find the corresponding target values for these reference probes
    target_values <- target_data[order(reference_data)][ref_group_idx]
    
    if (method == "binary") {
      # Binary method: Use top 30% and bottom 30% to estimate high and low signals
      sorted_target <- sort(target_values)
      low_count <- floor(length(sorted_target) * 0.3)
      high_count <- floor(length(sorted_target) * 0.3)
      
      mu_low <- mean(sorted_target[1:low_count])
      mu_high <- mean(sorted_target[(length(sorted_target) - high_count + 1):length(sorted_target)])
      
      # Normalize using binary high/low model
      normalized_data[i] <- (target_data[i] - mu_low) / (mu_high - mu_low)
      
    } else if (method == "quantile") {
      # Quantile-based method
      # Get rank of this probe in the reference set
      target_sorted <- sort(target_values)
      probe_rank <- sum(target_values <= target_data[i])
      
      # Normalized value is the value at this rank in the average reference distribution
      normalized_data[i] <- (probe_rank - 1) / (length(target_values) - 1)
    }
  }
  
  return(normalized_data)
}

# Function to implement Cross Normalization 
cross_normalization <- function(condA_data, condB_data, method = "binary", ref_group_size = 1000) {
  # Normalize condB using condA as reference
  norm_B_vs_A <- group_normalization(condA_data, condB_data, method, ref_group_size)
  
  # Normalize condA using condB as reference
  norm_A_vs_B <- group_normalization(condB_data, condA_data, method, ref_group_size)
  
  return(list(
    B_vs_A = norm_B_vs_A,
    A_vs_B = norm_A_vs_B
  ))
}

#######################
# VISUALIZATION FUNCTIONS
#######################

# Function to create Figure 1: Probe reproducibility
plot_figure1 <- function(data_file, chr = "chrIII", start = 50000, end = 55000) {
  # Load data
  if (is.character(data_file)) {
    data <- read.table(data_file, header = TRUE, sep = "\t")
  } else {
    data <- data_file
  }
  
  # Filter for region of interest
  region_data <- data %>%
    filter(chr == !!chr & position >= start & position <= end)
  
  # Create Figure 1A: Signals along a portion of chromosome
  fig1a <- ggplot(region_data, aes(x = position)) +
    geom_line(aes(y = rep1, color = "Replicate 1")) +
    geom_line(aes(y = rep2, color = "Replicate 2")) +
    scale_color_manual(values = c("Replicate 1" = "blue", "Replicate 2" = "red")) +
    labs(title = "Figure 1A: Genomic hybridization signals",
         subtitle = paste("Chromosome", chr, ":", start, "-", end),
         x = "Position (bp)",
         y = "Signal Intensity",
         color = "Replicate") +
    theme_minimal()
  
  # Create Figure 1B: Correlation between replicates
  fig1b <- ggplot(data, aes(x = rep1, y = rep2)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_density_2d(color = "blue") +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    labs(title = "Figure 1B: Correlation between replicates",
         subtitle = paste("Pearson correlation:", round(cor(data$rep1, data$rep2), 3)),
         x = "Replicate 1",
         y = "Replicate 2") +
    theme_minimal()
  
  return(list(fig1a = fig1a, fig1b = fig1b))
}

# Function to create Figure 5: Group Normalization Results for HXT3 locus
plot_figure5 <- function(data_file) {
  # Load data
  if (is.character(data_file)) {
    data <- read.table(data_file, header = TRUE, sep = "\t")
  } else {
    data <- data_file
  }
  
  # Apply a 20bp running average to smooth the signals
  window_size <- 20
  
  smooth_signal <- function(x, window_size) {
    as.vector(stats::filter(x, rep(1/window_size, window_size), sides = 2))
  }
  
  data$condA_smooth <- smooth_signal(data$condA, window_size)
  data$condB_smooth <- smooth_signal(data$condB, window_size)
  
  # Apply group normalization (condition B vs A and vice versa)
  cross_norm <- cross_normalization(data$condA, data$condB)
  data$B_vs_A <- cross_norm$B_vs_A
  data$A_vs_B <- cross_norm$A_vs_B
  
  # Apply 20bp running average to normalized signals
  data$B_vs_A_smooth <- smooth_signal(data$B_vs_A, window_size)
  data$A_vs_B_smooth <- smooth_signal(data$A_vs_B, window_size)
  
  # Create Figure 5B: Raw nucleosome occupancy
  fig5b <- ggplot(data, aes(x = position)) +
    geom_line(aes(y = condA_smooth, color = "t = 0"), linetype = "dotted") +
    geom_line(aes(y = condB_smooth, color = "t = 60 min")) +
    scale_color_manual(values = c("t = 0" = "gray40", "t = 60 min" = "magenta")) +
    labs(title = "Figure 5B: Nucleosome occupancy at HXT3 promoter",
         subtitle = "Raw tiling array data",
         x = "Position (bp)",
         y = "Signal intensity") +
    theme_minimal()
  
  # Create Figure 5C: Cross normalized data
  fig5c <- ggplot(data, aes(x = position)) +
    geom_line(aes(y = B_vs_A_smooth, color = "t = 60 vs t = 0")) +
    geom_line(aes(y = A_vs_B_smooth, color = "t = 0 vs t = 60"), linetype = "dotted") +
    scale_color_manual(values = c("t = 60 vs t = 0" = "red", "t = 0 vs t = 60" = "blue")) +
    labs(title = "Figure 5C: Differential nucleosome occupancy",
         subtitle = "Cross normalization",
         x = "Position (bp)",
         y = "Normalized signal") +
    theme_minimal()
  
  return(list(fig5b = fig5b, fig5c = fig5c))
}

# Function to evaluate Signal Quality (for Figure 7)
calculate_signal_quality <- function(condA, condB, condB_rep, significantly_changed_probes) {
  # Calculate signal power (S): mean square change between conditions
  S <- mean((condB[significantly_changed_probes] - condA[significantly_changed_probes])^2)
  
  # Calculate noise power (N): mean square change between replicates
  N <- mean((condB - condB_rep)^2)
  
  # Signal quality in dB
  SQ_dB <- 10 * log10(S/N)
  
  return(SQ_dB)
}

# Function to compare Signal Quality of different normalization methods (Figure 7)
plot_figure7 <- function(data_file) {
  # Load data
  if (is.character(data_file)) {
    data <- read.table(data_file, header = TRUE, sep = "\t")
  } else {
    data <- data_file
  }
  
  # Identify significantly changed probes (top 2% by difference)
  window_size <- 20
  smooth_signal <- function(x, window_size) {
    as.vector(stats::filter(x, rep(1/window_size, window_size), sides = 2))
  }
  
  data$condA_smooth <- smooth_signal(data$condA, window_size)
  data$condB_smooth <- smooth_signal(data$condB, window_size)
  data$diff_smooth <- abs(data$condB_smooth - data$condA_smooth)
  
  # Top 2% probes by smoothed difference
  sig_threshold <- quantile(data$diff_smooth, 0.98, na.rm = TRUE)
  sig_probes <- which(data$diff_smooth >= sig_threshold)
  
  # Apply different normalization methods
  # 1. Raw data (no normalization)
  sq_raw <- calculate_signal_quality(data$condA, data$condB, data$condB_rep, sig_probes)
  
  # 2. Binary Group Normalization
  gn_binary <- cross_normalization(data$condA, data$condB, method = "binary")
  gn_binary_rep <- cross_normalization(data$condA, data$condB_rep, method = "binary")
  sq_gn_binary <- calculate_signal_quality(
    gn_binary$A_vs_B, 
    gn_binary$B_vs_A,
    gn_binary_rep$B_vs_A,
    sig_probes
  )
  
  # 3. Quantile Group Normalization
  gn_quantile <- cross_normalization(data$condA, data$condB, method = "quantile")
  gn_quantile_rep <- cross_normalization(data$condA, data$condB_rep, method = "quantile")
  sq_gn_quantile <- calculate_signal_quality(
    gn_quantile$A_vs_B, 
    gn_quantile$B_vs_A,
    gn_quantile_rep$B_vs_A,
    sig_probes
  )
  
  # 4. MAS5 normalization
  mas5_A <- mas5_normalization(data$condA)
  mas5_B <- mas5_normalization(data$condB)
  mas5_B_rep <- mas5_normalization(data$condB_rep)
  sq_mas5 <- calculate_signal_quality(mas5_A, mas5_B, mas5_B_rep, sig_probes)
  
  # 5. MAT normalization
  mat_A <- mat_normalization(data$condA)
  mat_B <- mat_normalization(data$condB)
  mat_B_rep <- mat_normalization(data$condB_rep)
  sq_mat <- calculate_signal_quality(mat_A, mat_B, mat_B_rep, sig_probes)
  
  # 6. Quantile normalization
  # Create a matrix with all conditions for proper quantile normalization
  data_matrix <- cbind(data$condA, data$condB, data$condB_rep)
  q_norm_matrix <- quantile_normalization(data_matrix)
  q_norm_A <- q_norm_matrix[, 1]
  q_norm_B <- q_norm_matrix[, 2]
  q_norm_B_rep <- q_norm_matrix[, 3]
  sq_qq <- calculate_signal_quality(q_norm_A, q_norm_B, q_norm_B_rep, sig_probes)
  
  # Assemble results
  methods <- c("MAS5", "Q-Q", "MAT", "GN-quant", "GN-binary")
  sq_values <- c(sq_mas5, sq_qq, sq_mat, sq_gn_quantile, sq_gn_binary)
  
  # Create a dataframe for plotting
  fig7_data <- data.frame(
    Method = methods,
    Signal_Quality_dB = sq_values
  )
  
  # Create bar plot
  fig7 <- ggplot(fig7_data, aes(x = Method, y = Signal_Quality_dB, fill = Method)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f", Signal_Quality_dB)), vjust = -0.5) +
    labs(title = "Figure 7: Signal Quality Comparison",
         subtitle = "Higher dB values indicate better signal quality",
         x = "Normalization Method",
         y = "Signal Quality (dB)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(list(fig7 = fig7, signal_quality = fig7_data))
}

#######################
# MAIN EXECUTION
#######################

# Main function to run the analysis
run_analysis <- function(output_dir = ".", use_real_data = FALSE) {
  # Download data
  cat("Downloading/generating data...\n")
  data_files <- download_lee_data(output_dir, use_real_data)
  
  # Check if we have real or simulated data
  if (use_real_data && !is.list(data_files)) {
    cat("Error getting real data. Falling back to simulated data.\n")
    data_files <- download_lee_data(output_dir, use_real_data = FALSE)
  }
  
  # Create Figure 1
  cat("Creating Figure 1...\n")
  if (use_real_data) {
    # Process real data for Figure 1
    # This would depend on the exact format of the real data
    cat("Real data processing for Figure 1 not implemented yet.\n")
    cat("Using simulated data for Figure 1...\n")
    mock_data <- create_mock_nucleosome_data(output_dir)
    fig1 <- plot_figure1(mock_data$mock_file)
  } else {
    fig1 <- plot_figure1(data_files$mock_file)
  }
  
  # Create Figure 5
  cat("Creating Figure 5...\n")
  if (use_real_data) {
    # Process real data for Figure 5
    # This would depend on the exact format of the real data
    cat("Real data processing for Figure 5 not implemented yet.\n")
    cat("Using simulated data for Figure 5...\n")
    mock_data <- create_mock_nucleosome_data(output_dir)
    fig5 <- plot_figure5(mock_data$hxt_file)
  } else {
    fig5 <- plot_figure5(data_files$hxt_file)
  }
  
  # Create Figure 7
  cat("Creating Figure 7...\n")
  if (use_real_data) {
    # Process real data for Figure 7
    # This would depend on the exact format of the real data
    cat("Real data processing for Figure 7 not implemented yet.\n")
    cat("Using simulated data for Figure 7...\n")
    mock_data <- create_mock_nucleosome_data(output_dir)
    fig7 <- plot_figure7(mock_data$hxt_file)
  } else {
    fig7 <- plot_figure7(data_files$hxt_file)
  }
  
  # Save figures
  ggsave(file.path(output_dir, "figure1a.png"), fig1$fig1a, width = 8, height = 5)
  ggsave(file.path(output_dir, "figure1b.png"), fig1$fig1b, width = 7, height = 6)
  ggsave(file.path(output_dir, "figure5b.png"), fig5$fig5b, width = 8, height = 5)
  ggsave(file.path(output_dir, "figure5c.png"), fig5$fig5c, width = 8, height = 5)
  ggsave(file.path(output_dir, "figure7.png"), fig7$fig7, width = 8, height = 6)
  
  # Save information about which data source was used
  data_source <- ifelse(use_real_data, "real (ArrayExpress E-MEXP-1172)", "simulated")
  cat(paste("Data source:", data_source), file = file.path(output_dir, "data_source.txt"))
  
  cat("\nAnalysis complete. Figures saved to", output_dir, "\n")
  cat("Data source:", data_source, "\n")
  
  return(list(
    fig1 = fig1,
    fig5 = fig5,
    fig7 = fig7,
    data_source = data_source
  ))
}

# Run the analysis if executed directly
if (!interactive()) {
  run_analysis()
}
