import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
import matplotlib.patches as patches

# Set random seed for reproducibility
np.random.seed(42)

# Create mock data to simulate Figure 3 from Ghandi & Beer 2012
# This simulates the Lee 2007 data mentioned in the paper

# Number of probes
n_probes = 5000

# Generate reference condition data (genomic hybridization) with log-normal distribution
reference_data = np.exp(np.random.normal(0, 0.5, n_probes))
reference_data = np.sort(reference_data)  # Sort by intensity

# Generate experimental condition data (nucleosome-enriched) 
# The experimental data has a relationship with reference data but with variations
experimental_data = reference_data * (0.7 + 0.6 * np.random.random(n_probes)) 

# Create a DataFrame to store the data
data = pd.DataFrame({
    'probe_index': np.arange(n_probes),
    'reference': reference_data,
    'experimental': experimental_data
})

# Sort by reference values
data = data.sort_values('reference').reset_index(drop=True)

# Function to find the 1000 nearest neighbors for a specific probe
def find_reference_set(probe_idx, data, n_ref=1000):
    ref_val = data.loc[probe_idx, 'reference']
    data['dist'] = np.abs(data['reference'] - ref_val)
    ref_set = data.sort_values('dist').head(n_ref)
    return ref_set

# Choose two probes to highlight (one with high signal, one with low signal)
# These are at different positions in the sorted reference data
low_signal_probe_idx = 1500
high_signal_probe_idx = 3500

# Get reference sets
low_signal_ref_set = find_reference_set(low_signal_probe_idx, data)
high_signal_ref_set = find_reference_set(high_signal_probe_idx, data)

# Create Figure 3 from the paper
plt.figure(figsize=(10, 6))

# Plot reference condition (black)
plt.scatter(data['probe_index'], data['reference'], s=1, color='black', alpha=0.7, label='Reference condition')

# Plot experimental condition (gray)
plt.scatter(data['probe_index'], data['experimental'], s=1, color='gray', alpha=0.5, label='Experimental condition')

# Highlight the two probe reference sets with rectangles
rect_height = max(data['experimental'].max(), data['reference'].max()) * 1.1
rect_width = 1000  # Number of probes in reference set

# Low signal probe reference set (blue rectangle)
low_signal_rect = patches.Rectangle(
    (low_signal_probe_idx - rect_width/2, 0), 
    rect_width, rect_height, 
    linewidth=1, edgecolor='blue', facecolor='blue', alpha=0.1
)
plt.gca().add_patch(low_signal_rect)

# High signal probe reference set (blue rectangle)
high_signal_rect = patches.Rectangle(
    (high_signal_probe_idx - rect_width/2, 0), 
    rect_width, rect_height, 
    linewidth=1, edgecolor='blue', facecolor='blue', alpha=0.1
)
plt.gca().add_patch(high_signal_rect)

# Mark the probes of interest
plt.scatter([low_signal_probe_idx], [data.loc[low_signal_probe_idx, 'reference']], 
            s=50, color='blue', edgecolor='black', zorder=10)
plt.scatter([high_signal_probe_idx], [data.loc[high_signal_probe_idx, 'reference']], 
            s=50, color='blue', edgecolor='black', zorder=10)

# Show low and high signal levels in the experimental condition
# For the low signal reference set
low_ref_set_sorted = low_signal_ref_set.sort_values('experimental')
low_signal_low = low_ref_set_sorted['experimental'].iloc[:300].mean()  # Bottom 30%
low_signal_high = low_ref_set_sorted['experimental'].iloc[-300:].mean()  # Top 30%

# For the high signal reference set
high_ref_set_sorted = high_signal_ref_set.sort_values('experimental')
high_signal_low = high_ref_set_sorted['experimental'].iloc[:300].mean()  # Bottom 30%
high_signal_high = high_ref_set_sorted['experimental'].iloc[-300:].mean()  # Top 30%

# Add horizontal lines for signal levels
plt.hlines(low_signal_low, low_signal_probe_idx - rect_width/2, low_signal_probe_idx + rect_width/2, 
           colors='green', linestyles='dashed', linewidth=2, label='Low signal level')
plt.hlines(low_signal_high, low_signal_probe_idx - rect_width/2, low_signal_probe_idx + rect_width/2, 
           colors='red', linestyles='dashed', linewidth=2, label='High signal level')

plt.hlines(high_signal_low, high_signal_probe_idx - rect_width/2, high_signal_probe_idx + rect_width/2, 
           colors='green', linestyles='dashed', linewidth=2)
plt.hlines(high_signal_high, high_signal_probe_idx - rect_width/2, high_signal_probe_idx + rect_width/2, 
           colors='red', linestyles='dashed', linewidth=2)

plt.xlabel('Probe Index (sorted by reference condition)')
plt.ylabel('Signal Intensity')
plt.title('Figure 3 Recreation: Group Normalization Principle')
plt.legend()
plt.tight_layout()

# Create the additional plot with reference log probe intensity on x-axis
plt.figure(figsize=(10, 6))

# Plot log of reference vs experimental
plt.scatter(np.log(data['reference']), data['experimental'], s=1, color='gray', alpha=0.5)

# Highlight the same two probes and their reference sets
log_low_ref = np.log(data.loc[low_signal_probe_idx, 'reference'])
log_high_ref = np.log(data.loc[high_signal_probe_idx, 'reference'])

# Calculate the boundaries of the vertical rectangles in log-reference space
log_low_min = np.log(low_signal_ref_set['reference'].min())
log_low_max = np.log(low_signal_ref_set['reference'].max())
log_high_min = np.log(high_signal_ref_set['reference'].min())
log_high_max = np.log(high_signal_ref_set['reference'].max())

# Add vertical rectangles
low_rect_x_width = log_low_max - log_low_min
high_rect_x_width = log_high_max - log_high_min

low_signal_rect_log = patches.Rectangle(
    (log_low_min, 0), 
    low_rect_x_width, rect_height, 
    linewidth=1, edgecolor='blue', facecolor='blue', alpha=0.1
)
plt.gca().add_patch(low_signal_rect_log)

high_signal_rect_log = patches.Rectangle(
    (log_high_min, 0), 
    high_rect_x_width, rect_height, 
    linewidth=1, edgecolor='blue', facecolor='blue', alpha=0.1
)
plt.gca().add_patch(high_signal_rect_log)

# Mark the probes of interest
plt.scatter([log_low_ref], [data.loc[low_signal_probe_idx, 'experimental']], 
            s=50, color='blue', edgecolor='black', zorder=10)
plt.scatter([log_high_ref], [data.loc[high_signal_probe_idx, 'experimental']], 
            s=50, color='blue', edgecolor='black', zorder=10)

# Add horizontal lines for signal levels
plt.hlines(low_signal_low, log_low_min, log_low_max, 
           colors='green', linestyles='dashed', linewidth=2, label='Low signal level')
plt.hlines(low_signal_high, log_low_min, log_low_max, 
           colors='red', linestyles='dashed', linewidth=2, label='High signal level')

plt.hlines(high_signal_low, log_high_min, log_high_max, 
           colors='green', linestyles='dashed', linewidth=2)
plt.hlines(high_signal_high, log_high_min, log_high_max, 
           colors='red', linestyles='dashed', linewidth=2)

plt.xlabel('Log Reference Probe Intensity')
plt.ylabel('Experimental Signal Intensity')
plt.title('Additional Plot: Reference Log Probe Intensity vs Experimental Signal')
plt.legend()
plt.tight_layout()

plt.show()

# Let's also create a more realistic simulation of the Lee 2007 data with actual probe intensity distributions

def simulate_lee_2007_data(n_probes=5000, n_highlight=100):
    """
    Simulate data similar to Lee 2007 nucleosome positioning data with more realistic distributions
    """
    # Create probe indices
    probe_indices = np.arange(n_probes)
    
    # Generate reference condition (genomic hybridization)
    # Use a mixture of log-normal distributions to create a more realistic shape
    reference_mix1 = np.exp(np.random.normal(0, 0.3, n_probes))
    reference_mix2 = np.exp(np.random.normal(0.5, 0.2, n_probes))
    reference_data = 0.7 * reference_mix1 + 0.3 * reference_mix2
    
    # Add a smoothly varying component to simulate probe effects
    x = np.linspace(0, 2*np.pi, n_probes)
    smooth_variation = 0.5 + 0.2 * np.sin(x) + 0.1 * np.sin(3*x) + 0.05 * np.sin(7*x)
    reference_data *= smooth_variation
    
    # Sort by reference intensity
    sorted_indices = np.argsort(reference_data)
    reference_data = reference_data[sorted_indices]
    
    # Generate experimental condition (nucleosome-enriched)
    # Create regions of high and low nucleosome occupancy
    experimental_data = np.zeros(n_probes)
    
    # Base level relation to reference data (probe effect)
    experimental_data = 0.5 * reference_data + 0.5 * np.random.random(n_probes)
    
    # Add nucleosome peaks (high occupancy regions)
    for i in range(10):
        center = np.random.randint(500, n_probes-500)
        width = np.random.randint(100, 200)  # ~150bp for nucleosomes
        height = 1.5 + 0.5 * np.random.random()
        
        # Create a Gaussian peak
        peak = height * np.exp(-0.5 * ((probe_indices - center) / (width/3))**2)
        experimental_data += peak
    
    # Select some highlighted probes that are part of a nucleosome peak
    # and their 100 nearest neighbors
    highlight_centers = []
    for i in range(2):
        center = np.random.randint(1000, n_probes-1000)
        highlight_centers.append(center)
    
    return probe_indices, reference_data, experimental_data, highlight_centers

# Generate more realistic data
probe_indices, reference_data, experimental_data, highlight_centers = simulate_lee_2007_data()

# Create DataFrame
realistic_data = pd.DataFrame({
    'probe_index': probe_indices,
    'reference': reference_data,
    'experimental': experimental_data
})

# Find reference sets for highlighted probes
highlight1_ref_set = find_reference_set(highlight_centers[0], realistic_data, n_ref=100)
highlight2_ref_set = find_reference_set(highlight_centers[1], realistic_data, n_ref=100)

# Create a more realistic plot
plt.figure(figsize=(12, 8))

# Plot the data
plt.subplot(2, 1, 1)
plt.plot(realistic_data['probe_index'], realistic_data['reference'], '-', color='black', alpha=0.7, linewidth=1, label='Reference (genomic)')
plt.plot(realistic_data['probe_index'], realistic_data['experimental'], '-', color='red', alpha=0.5, linewidth=1, label='Experimental (nucleosome)')

# Highlight the reference sets
highlight1_min = highlight1_ref_set['probe_index'].min()
highlight1_max = highlight1_ref_set['probe_index'].max()
highlight2_min = highlight2_ref_set['probe_index'].min()
highlight2_max = highlight2_ref_set['probe_index'].max()

# Add vertical rectangles
rect_height = max(realistic_data['experimental'].max(), realistic_data['reference'].max()) * 1.1

plt.axvspan(highlight1_min, highlight1_max, color='blue', alpha=0.2)
plt.axvspan(highlight2_min, highlight2_max, color='green', alpha=0.2)

plt.xlabel('Probe Index (sorted by reference condition)')
plt.ylabel('Signal Intensity')
plt.title('Realistic Simulation of Lee 2007 Data')
plt.legend()

# Plot log reference vs experimental (the additional plot requested)
plt.subplot(2, 1, 2)
plt.scatter(np.log(realistic_data['reference']), realistic_data['experimental'], s=1, color='gray', alpha=0.5)

# Highlight the same regions in log space
log_highlight1_min = np.log(highlight1_ref_set['reference'].min())
log_highlight1_max = np.log(highlight1_ref_set['reference'].max())
log_highlight2_min = np.log(highlight2_ref_set['reference'].min())
log_highlight2_max = np.log(highlight2_ref_set['reference'].max())

plt.axvspan(log_highlight1_min, log_highlight1_max, color='blue', alpha=0.2)
plt.axvspan(log_highlight2_min, log_highlight2_max, color='green', alpha=0.2)

plt.xlabel('Log Reference Probe Intensity')
plt.ylabel('Experimental Signal Intensity')
plt.title('Additional Plot: Log Reference vs Experimental Signal')
plt.tight_layout()

plt.show()
