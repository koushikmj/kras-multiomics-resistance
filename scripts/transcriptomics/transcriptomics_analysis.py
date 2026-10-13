#!/usr/bin/env python3
# ============================================================================
# Transcriptomics Analysis: KRAS A146 vs Control  (SYNTHETIC DATA)
# BINF6420 Project - Resistance Mechanisms in KRAS-mutant Cancers
#
# Inputs : data/transcriptomics/counts_kras_a146.csv, meta_kras_a146.csv
# Outputs: results/transcriptomics/*.csv, figures/transcriptomics/*.png
#
# Method : per-gene two-sided Student t-test on raw counts, log2FC with +1
#          pseudocount, |log2FC| > 1 and p < 0.05. This is a teaching-level
#          approach; see README "Limitations" (no library-size normalisation,
#          no multiple-testing correction; DESeq2 is the proper tool).
# ============================================================================

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "transcriptomics"
RESULTS_DIR = REPO_ROOT / "results" / "transcriptomics"
FIGURES_DIR = REPO_ROOT / "figures" / "transcriptomics"
COUNTS_FILE = DATA_DIR / "counts_kras_a146.csv"
META_FILE = DATA_DIR / "meta_kras_a146.csv"

for f in (COUNTS_FILE, META_FILE):
    if not f.exists():
        sys.exit(f"Missing input: {f}\nSee data/transcriptomics/README.md for how to add the count matrix.")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# ============================================================================
# STEP 1: Load RNA-seq count data
# ============================================================================

print("Loading data...")
counts = pd.read_csv(COUNTS_FILE)
meta = pd.read_csv(META_FILE)

# Set gene names as index
counts = counts.rename(columns={'Unnamed: 0': 'gene'})
counts = counts.set_index('gene')

print(f"✓ Data loaded: {counts.shape[0]} genes, {counts.shape[1]} samples")

# Identify sample groups
ctrl_cols = [col for col in counts.columns if 'CTRL' in col]
kras_cols = [col for col in counts.columns if 'A146' in col]

print(f"  Control samples (n={len(ctrl_cols)}): {ctrl_cols}")
print(f"  KRAS A146 samples (n={len(kras_cols)}): {kras_cols}")

# ============================================================================
# STEP 2: Differential Expression Analysis
# ============================================================================

print("\nPerforming differential expression analysis...")

de_results = pd.DataFrame()
de_results['gene'] = counts.index

# Calculate mean expression for each group
de_results['control_mean'] = counts[ctrl_cols].mean(axis=1).values
de_results['kras_mean'] = counts[kras_cols].mean(axis=1).values

# Calculate fold change (add pseudocount to avoid division by zero)
de_results['fold_change'] = (de_results['kras_mean'] + 1) / (de_results['control_mean'] + 1)
de_results['log2FC'] = np.log2(de_results['fold_change'])

# Calculate p-values using independent t-test
pvalues = []
for gene in counts.index:
    ctrl_vals = counts.loc[gene, ctrl_cols].values
    kras_vals = counts.loc[gene, kras_cols].values

    try:
        t_stat, p_val = stats.ttest_ind(ctrl_vals, kras_vals)
        pvalues.append(p_val)
    except:
        pvalues.append(1.0)

de_results['pvalue'] = pvalues
de_results['-log10_pvalue'] = -np.log10(de_results['pvalue'] + 1e-300)

# Remove infinite values
de_results = de_results.replace([np.inf, -np.inf], np.nan)
de_results = de_results.dropna()

# Define significance thresholds
FC_THRESHOLD = 1.0  # log2 fold change (2-fold)
PVAL_THRESHOLD = 0.05

# Categorize genes
de_results['category'] = 'Not Significant'
de_results.loc[(de_results['log2FC'] > FC_THRESHOLD) & 
               (de_results['pvalue'] < PVAL_THRESHOLD), 'category'] = 'Upregulated'
de_results.loc[(de_results['log2FC'] < -FC_THRESHOLD) & 
               (de_results['pvalue'] < PVAL_THRESHOLD), 'category'] = 'Downregulated'

# Summary statistics
n_up = len(de_results[de_results['category'] == 'Upregulated'])
n_down = len(de_results[de_results['category'] == 'Downregulated'])
n_ns = len(de_results[de_results['category'] == 'Not Significant'])

print(f"✓ Analysis complete:")
print(f"  Upregulated genes: {n_up}")
print(f"  Downregulated genes: {n_down}")
print(f"  Not significant: {n_ns}")

# Save results
de_results.to_csv(RESULTS_DIR / 'differential_expression_results.csv', index=False)
print("  Saved: differential_expression_results.csv")

# ============================================================================
# STEP 3: Create Volcano Plot
# ============================================================================

print("\nCreating volcano plot...")

fig, ax = plt.subplots(figsize=(12, 8))

# Define colors
colors = {'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'gray'}

# Plot each category
for category in ['Not Significant', 'Downregulated', 'Upregulated']:
    subset = de_results[de_results['category'] == category]
    ax.scatter(subset['log2FC'], subset['-log10_pvalue'], 
               c=colors[category], alpha=0.6, s=30, 
               label=f'{category} (n={len(subset)})', edgecolors='none')

# Add threshold lines
pval_log_threshold = -np.log10(PVAL_THRESHOLD)
ax.axhline(y=pval_log_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.axvline(x=FC_THRESHOLD, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.axvline(x=-FC_THRESHOLD, color='black', linestyle='--', linewidth=1, alpha=0.5)

# Label top significant genes
significant_genes = de_results[de_results['category'] != 'Not Significant']
top_genes = significant_genes.nlargest(8, '-log10_pvalue')

for idx, row in top_genes.iterrows():
    ax.annotate(row['gene'], xy=(row['log2FC'], row['-log10_pvalue']),
               xytext=(5, 5), textcoords='offset points', fontsize=8, alpha=0.8,
               bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3),
               arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=0.5))

# Formatting
ax.set_xlabel('log₂(Fold Change)', fontsize=14, fontweight='bold')
ax.set_ylabel('-log₁₀(p-value)', fontsize=14, fontweight='bold')
ax.set_title('Volcano Plot: KRAS A146 vs Control', fontsize=16, fontweight='bold', pad=20)
ax.legend(loc='upper right', frameon=True, fontsize=10)
ax.grid(True, alpha=0.3, linestyle=':')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'volcano_plot_kras.png', dpi=300, bbox_inches='tight')
print("  Saved: volcano_plot_kras.png")
plt.close()

# ============================================================================
# STEP 4: Principal Component Analysis (PCA)
# ============================================================================

print("\nPerforming PCA analysis...")

# Transpose data (samples as rows, genes as columns)
counts_t = counts.T

# Log2 transformation
counts_log = np.log2(counts_t + 1)

# Standardize
scaler = StandardScaler()
counts_scaled = scaler.fit_transform(counts_log)

# Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(counts_scaled)

# Create PCA dataframe
pca_df = pd.DataFrame({
    'sample': counts.columns,
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'group': ['Control' if 'CTRL' in c else 'KRAS A146' for c in counts.columns]
})

var1 = pca.explained_variance_ratio_[0] * 100
var2 = pca.explained_variance_ratio_[1] * 100

print(f"✓ PCA complete:")
print(f"  PC1 explains {var1:.1f}% of variance")
print(f"  PC2 explains {var2:.1f}% of variance")

# Save PCA results
pca_df.to_csv(RESULTS_DIR / 'pca_results.csv', index=False)
print("  Saved: pca_results.csv")

# Create PCA plot
fig, ax = plt.subplots(figsize=(10, 8))

# Plot points
for group in ['Control', 'KRAS A146']:
    subset = pca_df[pca_df['group'] == group]
    color = 'blue' if group == 'Control' else 'red'
    ax.scatter(subset['PC1'], subset['PC2'], c=color, s=150, 
               alpha=0.7, edgecolors='black', linewidth=1.5, label=group)

    # Add sample labels
    for idx, row in subset.iterrows():
        ax.annotate(row['sample'], (row['PC1'], row['PC2']),
                   xytext=(8, 8), textcoords='offset points', fontsize=9)

ax.set_xlabel(f'PC1 ({var1:.1f}% variance)', fontsize=12, fontweight='bold')
ax.set_ylabel(f'PC2 ({var2:.1f}% variance)', fontsize=12, fontweight='bold')
ax.set_title('PCA: KRAS A146 vs Control Samples', fontsize=14, fontweight='bold')
ax.legend(loc='best', frameon=True, fontsize=11)
ax.grid(True, alpha=0.3, linestyle=':')
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
ax.axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'pca_plot_kras.png', dpi=300, bbox_inches='tight')
print("  Saved: pca_plot_kras.png")
plt.close()

# ============================================================================
# STEP 5: Heatmap of Top Differentially Expressed Genes
# ============================================================================

print("\nCreating heatmap of top 50 genes...")

# Get top 50 most significant genes
significant_genes = de_results[de_results['category'] != 'Not Significant']
top50_genes = significant_genes.nsmallest(50, 'pvalue')

print(f"  Selected 50 genes:")
print(f"    Upregulated: {len(top50_genes[top50_genes['category'] == 'Upregulated'])}")
print(f"    Downregulated: {len(top50_genes[top50_genes['category'] == 'Downregulated'])}")

# Extract counts for these genes
heatmap_data = counts.loc[top50_genes['gene'], ctrl_cols + kras_cols].copy()

# Log2 transform
heatmap_data_log = np.log2(heatmap_data + 1)

# Z-score normalization (per gene across samples)
heatmap_data_zscore = heatmap_data_log.sub(heatmap_data_log.mean(axis=1), axis=0)
heatmap_data_zscore = heatmap_data_zscore.div(heatmap_data_log.std(axis=1), axis=0)

# Create heatmap
fig, ax = plt.subplots(figsize=(10, 16))

# Plot heatmap
im = ax.imshow(heatmap_data_zscore, aspect='auto', cmap='RdBu_r', 
               vmin=-2.5, vmax=2.5, interpolation='nearest')

# Set ticks and labels
ax.set_xticks(range(len(heatmap_data.columns)))
ax.set_xticklabels(heatmap_data.columns, rotation=45, ha='right', fontsize=10)
ax.set_yticks(range(len(top50_genes)))
ax.set_yticklabels(top50_genes['gene'].values, fontsize=8)

# Add column separator between Control and KRAS
ax.axvline(x=len(ctrl_cols) - 0.5, color='black', linewidth=2)

# Labels
ax.set_xlabel('Samples', fontsize=12, fontweight='bold')
ax.set_ylabel('Genes', fontsize=12, fontweight='bold')
ax.set_title('Top 50 Differentially Expressed Genes: KRAS A146 vs Control', 
             fontsize=14, fontweight='bold', pad=20)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Z-score', fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'heatmap_top50_genes.png', dpi=300, bbox_inches='tight')
print("  Saved: heatmap_top50_genes.png")
plt.close()

# Save heatmap data
heatmap_data_zscore.to_csv(RESULTS_DIR / 'heatmap_data.csv')
print("  Saved: heatmap_data.csv")

# ============================================================================
# STEP 6: Summary Report
# ============================================================================

print("\n" + "="*80)
print("ANALYSIS SUMMARY")
print("="*80)

print(f"\nTotal genes analyzed: {len(de_results)}")
print(f"Significantly altered: {n_up + n_down}")
print(f"  - Upregulated: {n_up}")
print(f"  - Downregulated: {n_down}")

print("\nTop 10 Upregulated Genes:")
top_up = de_results[de_results['category'] == 'Upregulated'].nlargest(10, 'log2FC')
for idx, row in top_up.iterrows():
    fc = 2**row['log2FC']
    print(f"  {row['gene']:15s}: {fc:6.1f}x higher (log2FC={row['log2FC']:.2f}, p={row['pvalue']:.2e})")

print("\nTop 10 Downregulated Genes:")
top_down = de_results[de_results['category'] == 'Downregulated'].nsmallest(10, 'log2FC')
for idx, row in top_down.iterrows():
    fc = 2**row['log2FC']
    print(f"  {row['gene']:15s}: {fc:6.1f}x lower (log2FC={row['log2FC']:.2f}, p={row['pvalue']:.2e})")

print("\n" + "="*80)
print("FILES GENERATED")
print("="*80)
print("  1. differential_expression_results.csv - Full DE analysis results")
print("  2. pca_results.csv - PCA coordinates")
print("  3. heatmap_data.csv - Z-score normalized expression")
print("  4. volcano_plot_kras.png - Volcano plot visualization")
print("  5. pca_plot_kras.png - PCA plot visualization")
print("  6. heatmap_top50_genes.png - Heatmap visualization")
print("="*80)
print("\nAnalysis complete! ✓")