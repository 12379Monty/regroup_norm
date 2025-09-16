# Group Normalization for Genomic Data in Base R
# Based on Ghandi & Beer 2012, reproducing Figure 3
# Using only base R graphics (no ggplot2)

#' Import Lee 2007 dataset from ArrayExpress
#' 
#' This function downloads and processes the Lee 2007 nucleosome dataset (E-MEXP-1172)
#' from the ArrayExpress database.
#' @param save_path Path to save the downloaded data
#' @return A list containing the processed reference and experimental data
import_lee_2007_data <- function(save_path = "lee_2007_data") {
  # Create directory if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path)
  }
  
  # URLs for the data (assuming these are the correct URLs)
  # Note: In reality, you would need to find the actual URLs for the data files
  # from ArrayExpress using the accession code E-MEXP-1172
  reference_url <- "https://www.ebi.ac.uk/arrayexpress/files/E-MEXP-1172/E-MEXP-1172.raw.1.zip"
  experimental_url <- "https://www.ebi.ac.uk/arrayexpress/files/E-MEXP-1172/E-MEXP-1172.raw.2.zip"
  
  # Define local file paths
  reference_zip <- file.path(save_path, "reference.zip")
  experimental_zip <- file.path(save_path, "experimental.zip")
  
  # Download the files (commented out for safety)
  # download.file(reference_url, reference_zip, mode = "wb")
  # download.file(experimental_url, experimental_zip, mode = "wb")
  
  # Since direct download may not work due to access restrictions,
  # alternative approach is to use ArrayExpress package or manual download
  
  cat("Note: ArrayExpress data often requires manual download.\n")
  cat("Please go to https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MEXP-1172\n")
  cat("and download the data files manually.\n\n")
  
  # Check if files are already downloaded
  reference_file <- file.path(save_path, "genomic_reference.txt")
  experimental_file <- file.path(save_path, "nucleosome_experimental.txt")
  
  # Parse the data if files exist (this is a placeholder for the actual parsing code)
  if (file.exists(reference_file) && file.exists(experimental_file)) {
    reference_data <- read.delim(reference_file, sep = "\t")
    experimental_data <- read.delim(experimental_file, sep = "\t")
    
    # Process data
    # ... (specific processing steps would depend on the actual file format)
    
    return(list(
      reference = reference_data,
      experimental = experimental_data
    ))
  } else {
    cat("Data files not found. Using simulated data instead.\n")
    return(NULL)
  }
}

#' Generate simulated data mimicking Lee 2007 dataset
#' 
#' This function generates synthetic data with properties similar to the 
#' Lee 2007 nucleosome positioning data.
#' @param n_probes Number of probes to simulate
#' @param log_intensity_range Range of log intensity values
#' @return A data frame with simulated probe data
simulate_lee_2007_data <- function(n_probes = 5000, 
                                  log_intensity_range = c(2.0, 4.5)) {
  # Set seed for reproducibility
  set.seed(42)
  
  # Create probe indices
  probe_indices <- 1:n_probes
  
  # Generate reference condition (genomic hybridization) with log-normal distribution
  # Scale to the specified log intensity range
  min_val <- log_intensity_range[1]
  max_val <- log_intensity_range[2]
  
  # Generate log-normal data and transform to log scale
  reference_log <- min_val + (max_val - min_val) * 
    (pnorm(rnorm(n_probes, mean = 0, sd = 1)))
  
  # Sort by log intensity
  sorted_idx <- order(reference_log)
  reference_log <- reference_log[sorted_idx]
  
  # Generate experimental condition (nucleosome-enriched)
  # The experimental data has a relationship with reference data but with variations
  noise <- rnorm(n_probes, mean = 0, sd = 0.15)
  experimental_log <- reference_log + noise
  
  # Ensure we stay within the specified range
  experimental_log <- pmin(pmax(experimental_log, min_val), max_val)
  
  # Create a data frame
  data <- data.frame(
    probe_index = probe_indices,
    reference_log = reference_log,
    experimental_log = experimental_log
  )
  
  return(data)
}

#' Find reference set for a given probe
#' 
#' This function finds the n_ref probes with the closest reference log intensity
#' to the specified probe.
#' @param probe_idx Index of the probe
#' @param data Data frame containing probe data
#' @param n_ref Number of probes in the reference set
#' @return A data frame containing the reference set
find_reference_set <- function(probe_idx, data, n_ref = 1000) {
  # Get reference log intensity for the specified probe
  ref_val <- data$reference_log[probe_idx]
  
  # Calculate distance from each probe to the reference value
  dist <- abs(data$reference_log - ref_val)
  
  # Sort by distance and take the first n_ref
  ref_set <- data[order(dist)[1:n_ref], ]
  
  return(ref_set)
}

#' Create Figure 3 from Ghandi & Beer 2012 paper using base R graphics
#' 
#' This function recreates Figure 3 and the additional plot from the paper 
#' "Group Normalization for Genomic Data" by Ghandi and Beer (2012).
#' @param data Data frame containing probe data
#' @param probe1_idx Index of the first highlighted probe
#' @param probe2_idx Index of the second highlighted probe
#' @param n_ref Number of probes in each reference set
#' @param percentile_low Percentile to use for low signal
#' @param percentile_high Percentile to use for high signal
#' @param create_pdf Whether to create PDF files of the plots
create_figure3_base <- function(data, 
                              probe1_idx = NULL, 
                              probe2_idx = NULL,
                              n_ref = 1000,
                              percentile_low = 0.3,
                              percentile_high = 0.7,
                              create_pdf = TRUE) {
  
  # If probe indices are not provided, select them automatically
  if (is.null(probe1_idx)) {
    probe1_idx <- round(nrow(data) * 0.3)  # 30% point
  }
  if (is.null(probe2_idx)) {
    probe2_idx <- round(nrow(data) * 0.7)  # 70% point
  }
  
  # Find reference sets
  ref_set1 <- find_reference_set(probe1_idx, data, n_ref)
  ref_set2 <- find_reference_set(probe2_idx, data, n_ref)
  
  # Calculate low and high signal levels
  # For the first reference set
  ref_set1_sorted <- ref_set1[order(ref_set1$experimental_log), ]
  n1 <- nrow(ref_set1_sorted)
  low_cutoff1 <- ceiling(n1 * percentile_low)
  high_cutoff1 <- ceiling(n1 * percentile_high)
  
  low_signal1 <- mean(ref_set1_sorted$experimental_log[1:low_cutoff1])
  high_signal1 <- mean(ref_set1_sorted$experimental_log[high_cutoff1:n1])
  
  # For the second reference set
  ref_set2_sorted <- ref_set2[order(ref_set2$experimental_log), ]
  n2 <- nrow(ref_set2_sorted)
  low_cutoff2 <- ceiling(n2 * percentile_low)
  high_cutoff2 <- ceiling(n2 * percentile_high)
  
  low_signal2 <- mean(ref_set2_sorted$experimental_log[1:low_cutoff2])
  high_signal2 <- mean(ref_set2_sorted$experimental_log[high_cutoff2:n2])
  
  # Get probe indices for reference sets
  ref_set1_min_idx <- min(ref_set1$probe_index)
  ref_set1_max_idx <- max(ref_set1$probe_index)
  ref_set2_min_idx <- min(ref_set2$probe_index)
  ref_set2_max_idx <- max(ref_set2$probe_index)
  
  # Define plot margins and layout
  old_par <- par(no.readonly = TRUE)  # Save original parameters
  
  # Create PDF if requested
  if (create_pdf) {
    pdf("figure3_base_r.pdf", width = 10, height = 12)
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), oma = c(1, 1, 2, 1))
  } else {
    # Just create the plots in the current device
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), oma = c(1, 1, 2, 1))
  }
  
  # Create first plot (similar to Figure 3)
  # ----------------------------------------
  
  # Set up the plot region
  x_range <- range(data$probe_index)
  y_range <- range(c(data$reference_log, data$experimental_log))
  
  # Create empty plot
  plot(NA, type = "n", 
       xlim = x_range, 
       ylim = y_range, 
       xlab = "Probe Index (sorted by reference condition)",
       ylab = "Log Intensity",
       main = "Figure 3 Recreation: Group Normalization Principle",
       cex.main = 1.2,
       cex.lab = 1.1)
  
  # Add reference set rectangles
  rect(ref_set1_min_idx, y_range[1], 
       ref_set1_max_idx, y_range[2], 
       col = rgb(0, 0, 1, 0.1), border = NA)
  
  rect(ref_set2_min_idx, y_range[1], 
       ref_set2_max_idx, y_range[2], 
       col = rgb(0, 0, 1, 0.1), border = NA)
  
  # Add data points - use a subset for efficiency
  subset_factor <- max(1, floor(nrow(data) / 2000))  # Show at most 2000 points
  subset_idx <- seq(1, nrow(data), by = subset_factor)
  
  # Add reference and experimental data points
  points(data$probe_index[subset_idx], data$reference_log[subset_idx], 
         pch = 20, col = rgb(0, 0, 0, 0.5), cex = 0.5)
  points(data$probe_index[subset_idx], data$experimental_log[subset_idx], 
         pch = 20, col = rgb(0.5, 0.5, 0.5, 0.3), cex = 0.5)
  
  # Highlight the probes of interest
  points(data$probe_index[c(probe1_idx, probe2_idx)], 
         data$reference_log[c(probe1_idx, probe2_idx)], 
         pch = 19, col = "blue", cex = 1.5)
  
  # Add horizontal lines for signal levels
  segments(ref_set1_min_idx, low_signal1, ref_set1_max_idx, low_signal1, 
           col = "green3", lwd = 2, lty = 2)
  segments(ref_set1_min_idx, high_signal1, ref_set1_max_idx, high_signal1, 
           col = "red", lwd = 2, lty = 2)
  segments(ref_set2_min_idx, low_signal2, ref_set2_max_idx, low_signal2, 
           col = "green3", lwd = 2, lty = 2)
  segments(ref_set2_min_idx, high_signal2, ref_set2_max_idx, high_signal2, 
           col = "red", lwd = 2, lty = 2)
  
  # Add legend
  legend("topright", 
         legend = c("Reference condition", "Experimental condition", 
                    "Reference set", "Low signal level", "High signal level"),
         col = c("black", "gray", rgb(0, 0, 1, 0.5), "green3", "red"),
         pch = c(20, 20, NA, NA, NA),
         lty = c(NA, NA, 1, 2, 2),
         lwd = c(NA, NA, 10, 2, 2),
         cex = 0.8,
         bg = "white")
  
  # Create second plot (additional plot)
  # -----------------------------------
  
  # Set up the plot region for the second plot
  x_range <- range(data$reference_log)
  y_range <- range(data$experimental_log)
  
  # Create empty plot
  plot(NA, type = "n", 
       xlim = x_range, 
       ylim = y_range, 
       xlab = "Log Reference Probe Intensity",
       ylab = "Log Experimental Intensity",
       main = "Additional Plot: Reference Log Intensity vs Experimental Signal",
       cex.main = 1.2,
       cex.lab = 1.1)
  
  # Add reference set rectangles in log space
  rect(min(ref_set1$reference_log), y_range[1], 
       max(ref_set1$reference_log), y_range[2], 
       col = rgb(0, 0, 1, 0.1), border = NA)
  
  rect(min(ref_set2$reference_log), y_range[1], 
       max(ref_set2$reference_log), y_range[2], 
       col = rgb(0, 0, 1, 0.1), border = NA)
  
  # Add data points
  points(data$reference_log[subset_idx], data$experimental_log[subset_idx], 
         pch = 20, col = rgb(0.5, 0.5, 0.5, 0.3), cex = 0.5)
  
  # Highlight the probes of interest
  points(data$reference_log[c(probe1_idx, probe2_idx)], 
         data$experimental_log[c(probe1_idx, probe2_idx)], 
         pch = 19, col = "blue", cex = 1.5)
  
  # Add horizontal lines for signal levels
  segments(min(ref_set1$reference_log), low_signal1, 
           max(ref_set1$reference_log), low_signal1, 
           col = "green3", lwd = 2, lty = 2)
  segments(min(ref_set1$reference_log), high_signal1, 
           max(ref_set1$reference_log), high_signal1, 
           col = "red", lwd = 2, lty = 2)
  segments(min(ref_set2$reference_log), low_signal2, 
           max(ref_set2$reference_log), low_signal2, 
           col = "green3", lwd = 2, lty = 2)
  segments(min(ref_set2$reference_log), high_signal2, 
           max(ref_set2$reference_log), high_signal2, 
           col = "red", lwd = 2, lty = 2)
  
  # Add legend
  legend("topleft", 
         legend = c("Probe signal", "Reference set", 
                    "Low signal level", "High signal level"),
         col = c("gray", rgb(0, 0, 1, 0.5), "green3", "red"),
         pch = c(20, NA, NA, NA),
         lty = c(NA, 1, 2, 2),
         lwd = c(NA, 10, 2, 2),
         cex = 0.8,
         bg = "white")
  
  # Add overall title
  mtext("Group Normalization for Genomic Data (Ghandi & Beer 2012)", 
        outer = TRUE, cex = 1.5)
  
  # Close PDF if created
  if (create_pdf) {
    dev.off()
    cat("Plots saved to figure3_base_r.pdf\n")
  }
  
  # Reset graphics parameters
  par(old_par)
  
  # Return important calculated values
  return(list(
    probe1_idx = probe1_idx,
    probe2_idx = probe2_idx,
    ref_set1 = ref_set1,
    ref_set2 = ref_set2,
    low_signal1 = low_signal1,
    high_signal1 = high_signal1,
    low_signal2 = low_signal2,
    high_signal2 = high_signal2
  ))
}

#' Group Normalization implementation
#' 
#' This function implements the Group Normalization method from Ghandi & Beer 2012.
#' @param reference Reference condition data (log intensity)
#' @param experimental Experimental condition data (log intensity)
#' @param n_ref Number of probes in each reference set
#' @param percentile_low Percentile to use for low signal
#' @param percentile_high Percentile to use for high signal
#' @return A vector of normalized values
group_normalize <- function(reference, experimental, 
                           n_ref = 1000,
                           percentile_low = 0.3,
                           percentile_high = 0.7) {
  
  n_probes <- length(reference)
  normalized <- numeric(n_probes)
  
  # Create a data frame for sorting
  data <- data.frame(
    index = 1:n_probes,
    reference = reference,
    experimental = experimental
  )
  
  # Sort by reference value
  data <- data[order(data$reference), ]
  
  # For each probe, find its reference set and normalize
  for (i in 1:n_probes) {
    # Find the reference set (n_ref probes with closest reference values)
    start_idx <- max(1, i - n_ref/2)
    end_idx <- min(n_probes, i + n_ref/2 - 1)
    ref_set <- data[start_idx:end_idx, ]
    
    # Sort reference set by experimental value
    ref_set <- ref_set[order(ref_set$experimental), ]
    
    # Calculate low and high signal levels
    n_ref_actual <- nrow(ref_set)
    low_cutoff <- ceiling(n_ref_actual * percentile_low)
    high_cutoff <- ceiling(n_ref_actual * percentile_high)
    
    mu_low <- mean(ref_set$experimental[1:low_cutoff])
    mu_high <- mean(ref_set$experimental[high_cutoff:n_ref_actual])
    
    # Normalize the experimental value
    if (mu_high > mu_low) {
      normalized[data$index[i]] <- (data$experimental[i] - mu_low) / (mu_high - mu_low)
    } else {
      normalized[data$index[i]] <- 0
    }
  }
  
  return(normalized)
}

#' Plot normalized data using base R
#' 
#' @param data Data frame containing original and normalized data
#' @param create_pdf Whether to create a PDF file
plot_normalized_data <- function(data, create_pdf = TRUE) {
  # Create PDF if requested
  if (create_pdf) {
    pdf("normalized_data_base_r.pdf", width = 8, height = 6)
  }
  
  # Set up the plot
  par(mar = c(4, 4, 3, 1))
  
  # Create the plot
  plot(data$probe_index, data$normalized, 
       type = "p", pch = 20, cex = 0.5, col = rgb(0, 0, 1, 0.5),
       xlab = "Probe Index (sorted by reference condition)",
       ylab = "Normalized Signal",
       main = "Group Normalized Signal",
       cex.main = 1.2,
       cex.lab = 1.1)
  
  # Add a smooth trend line
  smooth_line <- smooth.spline(data$probe_index, data$normalized, spar = 0.5)
  lines(smooth_line, col = "red", lwd = 2)
  
  # Add legend
  legend("topright", 
         legend = c("Normalized data", "Smoothed trend"),
         col = c(rgb(0, 0, 1, 0.5), "red"),
         pch = c(20, NA),
         lty = c(NA, 1),
         lwd = c(NA, 2),
         cex = 0.8,
         bg = "white")
  
  # Close PDF if created
  if (create_pdf) {
    dev.off()
    cat("Normalized data plot saved to normalized_data_base_r.pdf\n")
  }
}

# Main execution

# Try to import real data
lee_data <- import_lee_2007_data()

# If real data is not available, use simulated data
if (is.null(lee_data)) {
  cat("Using simulated data...\n")
  sim_data <- simulate_lee_2007_data(n_probes = 5000, 
                                    log_intensity_range = c(2.0, 4.5))
  
  # Select probe indices that will have experimental log intensities 
  # close to 3.1 and 3.75
  # First, find probes with experimental values close to the targets
  target1 <- 3.1
  target2 <- 3.75
  dist1 <- abs(sim_data$experimental_log - target1)
  dist2 <- abs(sim_data$experimental_log - target2)
  
  probe1_idx <- which.min(dist1)
  probe2_idx <- which.min(dist2)
  
  cat("Selected probe indices:\n")
  cat("Probe 1 (target ~3.1):", probe1_idx, "with value", 
      sim_data$experimental_log[probe1_idx], "\n")
  cat("Probe 2 (target ~3.75):", probe2_idx, "with value", 
      sim_data$experimental_log[probe2_idx], "\n\n")
  
  # Create the figures using base R
  results <- create_figure3_base(sim_data, probe1_idx, probe2_idx)
  
  # Create an interactive plot as well
  x11(width = 10, height = 12)  # Opens a new interactive graphics window
  create_figure3_base(sim_data, probe1_idx, probe2_idx, create_pdf = FALSE)
  
  # Perform Group Normalization on the simulated data
  normalized_data <- group_normalize(sim_data$reference_log, 
                                    sim_data$experimental_log)
  
  # Add normalized data to the data frame
  sim_data$normalized <- normalized_data
  
  # Plot the normalized data
  plot_normalized_data(sim_data)
  
  # Also create an interactive plot of the normalized data
  x11(width = 8, height = 6)
  plot_normalized_data(sim_data, create_pdf = FALSE)
  
  cat("All plots created successfully!\n")
} else {
  # Process real data
  # ...specific code for processing the actual Lee 2007 data would go here
  cat("Processing real Lee 2007 data...\n")
}
