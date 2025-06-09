#%% md
# # Jupyter notebook sample
#%%
import numpy as np
# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import re
#%%
TARGET_PTM = ('Citrullination', r'UniMod:7', 'R')

PTM_DATA = '/Users/cgu3/Documents/diann/data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_citrullination_report-lib.parquet'
#%%
metadata = pd.read_csv("data/metadata.csv")
metadata
#%%
dat_raw = pd.read_parquet(PTM_DATA)
#%%
dat_raw
#%%
# left join the metadata with the citrullination data
dat = dat_raw.merge(metadata, on='ID', how='left')
#%%
dat.head()
#%%
# peptide-level sum intensity for each run and color by cancer type, only citrullination PTM
#%%
# check total number of row where the Modified.Sequence contain "(*)" and * means any number of any character
dat['has_target_PTM'] = dat['Modified.Sequence'].str.contains(TARGET_PTM[1])
dat['has_target_site'] = dat['Stripped.Sequence'].str.contains(TARGET_PTM[2])
dat['num_PTM'] = dat['Modified.Sequence'].str.count(TARGET_PTM[1])
#%%
total_intensity = dat.groupby(['Run', 'ID', 'assay', 'has_target_PTM']).agg(total_intensity=('Precursor.Quantity', 'sum')).reset_index()
total_intensity['total_intensity'] = np.log10(total_intensity['total_intensity'])
total_intensity = total_intensity.merge(metadata, on='ID')
#%%
total_intensity
#%%
import matplotlib.pyplot as plt

# Drop rows with missing 'group' values
total_intensity_labeled = total_intensity.dropna(subset=['group'])

# Aggregate the data by summing the intensities, grouped by 'group', 'Run', and 'assay'
aggregated_intensity = total_intensity_labeled.groupby(['group', 'ID', 'assay']).agg(
    total_intensity=('total_intensity', 'sum')
).reset_index()

# Define colors for 'group'
case_colors = ['blue', 'orange']

# Plot histograms in the same plot with different colors
plt.figure(figsize=(10, 6))
for assay in aggregated_intensity['assay'].unique():
    plt.figure(figsize=(10, 6))
    for case_label, case_color in zip(aggregated_intensity['group'].unique(), case_colors):
        subset = aggregated_intensity[(aggregated_intensity['group'] == case_label) & (aggregated_intensity['assay'] == assay)]
        subset['total_intensity'].hist(bins=100, alpha=0.5, label=f'{case_label}', color=case_color)

    plt.xlabel('Total intensity (log10)')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of total intensity by Case/Control for assay {assay}')
    plt.legend()
    plt.show()
#%%
import matplotlib.pyplot as plt

# Drop rows with missing 'group' values
total_intensity_labeled = total_intensity.dropna(subset=['group'])

# Define colors for each combination of 'group' and 'has_target_PTM'
color_mapping = {
    ('Case', True): 'green',
    ('Case', False): 'orange',
    ('Control', True): 'red',
    ('Control', False): 'yellow'
}

# Plot histograms in the same plot with different colors
plt.figure(figsize=(10, 6))
for (case_label, mod_label), color in color_mapping.items():
    subset = total_intensity_labeled[(total_intensity_labeled['group'] == case_label) &
                                     (total_intensity_labeled['has_target_PTM'] == mod_label)]
    label = f'{case_label}, {"Modified peptide" if mod_label else "Unmodified peptide"}'
    subset['total_intensity'].hist(bins=100, alpha=0.5, label=label, color=color)

plt.xlabel('Total intensity (log10)')
plt.ylabel('Frequency')
plt.title(f'Distribution of total intensity of peptide with {TARGET_PTM[0]} on site {TARGET_PTM[2]}')
plt.legend()
plt.show()
#%%
import matplotlib.pyplot as plt

# Drop rows with missing 'Cancer Type' values
total_intensity_labeled = total_intensity[(total_intensity['has_target_PTM']) & (total_intensity['assay'] == 'IgB')].dropna(subset=['Cancer Type'])

# Plot histograms in the same plot with different colors
plt.figure(figsize=(10, 6))
colors = plt.cm.get_cmap('tab10', len(total_intensity_labeled['Cancer Type'].unique()))

for i, (label, color) in enumerate(zip(total_intensity_labeled['Cancer Type'].unique(), colors.colors)):
    subset = total_intensity_labeled[total_intensity_labeled['Cancer Type'] == label]
    subset['total_intensity'].hist(bins=30, alpha=0.3, label=label, color=color)

plt.xlabel('Total intensity (log10)')
plt.ylabel('Frequency')
plt.title(f'Distribution of total intensity of peptide with {TARGET_PTM[0]} on site {TARGET_PTM[2]} by Cancer Type')
plt.legend()
plt.show()
#%%
import matplotlib.pyplot as plt

# Drop rows with missing 'Cancer Type' values
total_intensity_labeled = total_intensity[total_intensity['has_target_PTM']].dropna(subset=['Cancer Type'])

# Filter for specific cancer types and controls
desired_types = ['Breast', 'LEAP Control', 'MERIT Control']
# desired_types = ['Breast', 'Lung']
total_intensity_labeled = total_intensity_labeled[total_intensity_labeled['Cancer Type'].isin(desired_types)]

# Loop through each unique assay value
for assay in total_intensity_labeled['assay'].unique():
    subset_assay = total_intensity_labeled[total_intensity_labeled['assay'] == assay]

    # Plot histograms in the same plot with different colors
    plt.figure(figsize=(10, 6))
    colors = plt.cm.get_cmap('viridis', len(subset_assay['Cancer Type'].unique()))

    for i, (label, color) in enumerate(zip(subset_assay['Cancer Type'].unique(), colors.colors)):
        subset = subset_assay[subset_assay['Cancer Type'] == label]
        subset['total_intensity'].hist(bins=30, alpha=0.4, label=label, color=color)

    plt.xlabel('Total intensity (log10)')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of total intensity of modified peptide with {TARGET_PTM[0]} on site {TARGET_PTM[2]} by Cancer Type for assay {assay}')
    plt.legend()
    plt.show()
#%%
peptide_level_aggregation = dat[(dat['has_target_site'] ==True)].groupby(['ID', 'assay', 'Stripped.Sequence']).agg({'has_target_PTM': 'any'}).reset_index()

protein_level_aggregation = dat[(dat['has_target_site'] ==True)].groupby(['ID', 'assay',  'Protein.Group']).agg({'has_target_PTM': 'any'}).reset_index()
#%%
peptide_level_aggregation['has_target_PTM'].groupby(peptide_level_aggregation['assay']).mean()
#%%
protein_level_aggregation['has_target_PTM'].groupby(protein_level_aggregation['assay']).mean()
#%%
import matplotlib.pyplot as plt

# Merge the data
dat_labeled = protein_level_aggregation.merge(metadata, on='ID')

# Drop rows with missing 'group' values
dat_labeled = dat_labeled.dropna(subset=['group'])

# Calculate the mean of 'has_PTM' for each 'Run', 'group', and 'assay'
mean_has_PTM = dat_labeled.groupby(['ID', 'assay', 'group'])['has_target_PTM'].mean().reset_index()

# Get unique assays
unique_assays = mean_has_PTM['assay'].unique()

# Create subplots
fig, axes = plt.subplots(nrows=len(unique_assays), ncols=1, figsize=(10, 6 * len(unique_assays)))

# Ensure axes is always iterable
if len(unique_assays) == 1:
    axes = [axes]

# Loop through each unique assay and plot histograms
for ax, assay in zip(axes, unique_assays):
    subset_assay = mean_has_PTM[mean_has_PTM['assay'] == assay]
    for label, color in zip(subset_assay['group'].unique(), ['blue', 'orange']):
        subset = subset_assay[subset_assay['group'] == label]
        subset['has_target_PTM'].hist(bins=50, alpha=0.5, label=label, color=color, ax=ax)

    ax.set_xlabel('Proportion of modified peptides')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Distribution of proportion of modified protein with {TARGET_PTM[0]} on site {TARGET_PTM[2]} for assay {assay}')
    ax.legend()

plt.tight_layout()
plt.show()
#%%
import matplotlib.pyplot as plt

# Merge the data
dat_labeled = peptide_level_aggregation.merge(metadata, on='ID')

# Drop rows with missing 'group' values
dat_labeled = dat_labeled.dropna(subset=['group'])

# Calculate the total number of peptides with PTM for each 'Run', 'group', and 'assay'
total_peptides_with_PTM = dat_labeled.groupby(['ID', 'group', 'assay'])['has_target_PTM'].sum().reset_index()

# Get unique assays
unique_assays = total_peptides_with_PTM['assay'].unique()

# Create subplots
fig, axes = plt.subplots(nrows=len(unique_assays), ncols=1, figsize=(10, 6 * len(unique_assays)))

# Ensure axes is always iterable
if len(unique_assays) == 1:
    axes = [axes]

# Loop through each unique assay and plot histograms
for ax, assay in zip(axes, unique_assays):
    subset_assay = total_peptides_with_PTM[total_peptides_with_PTM['assay'] == assay]
    for label, color in zip(subset_assay['group'].unique(), ['blue', 'orange']):
        subset = subset_assay[subset_assay['group'] == label]
        subset['has_target_PTM'].hist(bins=50, alpha=0.5, label=label, color=color, ax=ax)

    ax.set_xlabel('Total number of peptides with PTM')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Histogram of Total Number of Peptides with {TARGET_PTM[0]} on site {TARGET_PTM[2]} by group for assay {assay}')
    ax.legend()

plt.tight_layout()
plt.show()
