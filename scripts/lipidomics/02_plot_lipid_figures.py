#!/usr/bin/env python3
"""Step 2: plot lipidomics Figures 1-3 and the PS-species volcano.

Reads the tables written by 01_summarize_lipid_classes.py and writes
PNG + PDF figures to figures/lipidomics/.
"""

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lipid_io import FIGURES_DIR, RESULTS_DIR
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Set publication-quality style
from matplotlib import font_manager
_available = {f.name for f in font_manager.fontManager.ttflist}
plt.rcParams['font.family'] = next((f for f in ['Arial', 'Helvetica'] if f in _available), 'DejaVu Sans')
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['figure.dpi'] = 300

# Color palette
color_wt = '#3498db'  # Blue
color_kras = '#e74c3c'  # Red
color_treatment = {
    'Control': ['#3498db', '#e74c3c'],
    'UGCGi': ['#2874a6', '#c0392b'],
    '2DG': ['#85c1e9', '#f1948a'],
    'Glucose Starvation': ['#5dade2', '#ec7063']
}

print("=" * 70)
print("KRAS LIPIDOME VISUALIZATION SUITE")
print("=" * 70)

# ============================================================================
# FIGURE 1: Lipid Class Abundance - WT vs KRAS G12V (Control)
# ============================================================================

print("\nGenerating Figure 1: Baseline Lipid Class Comparison...")

fig1_data = pd.read_csv(RESULTS_DIR / 'figure1_control_comparison.csv')

fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(fig1_data))
width = 0.35

# Create bars
bars1 = ax.bar(x - width/2, fig1_data['WT_mean'], width, 
               yerr=fig1_data['WT_std'], 
               label='WT Control', 
               color=color_wt,
               capsize=5,
               edgecolor='black',
               linewidth=1.2)

bars2 = ax.bar(x + width/2, fig1_data['KRAS_mean'], width,
               yerr=fig1_data['KRAS_std'], 
               label='KRAS G12V Control',
               color=color_kras,
               capsize=5,
               edgecolor='black',
               linewidth=1.2)

# Customize axes
ax.set_xlabel('Lipid Class', fontsize=12, fontweight='bold')
ax.set_ylabel('Total Abundance (AU)', fontsize=12, fontweight='bold')
ax.set_title('Baseline Lipid Class Abundance: WT vs KRAS G12V', 
             fontsize=14, fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(fig1_data['Lipid_Class'], fontsize=11)
ax.legend(fontsize=10, frameon=True, shadow=True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(0, max(fig1_data['WT_mean'].max(), fig1_data['KRAS_mean'].max()) * 1.2)

# Add significance stars
for i, row in fig1_data.iterrows():
    if row['significance'] != 'ns':
        y_pos = max(row['WT_mean'] + row['WT_std'], 
                    row['KRAS_mean'] + row['KRAS_std']) + 1
        ax.text(i, y_pos, row['significance'], 
                ha='center', va='bottom', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'Figure1_Lipid_Class_Comparison.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / 'Figure1_Lipid_Class_Comparison.pdf', bbox_inches='tight')
print("✓ Saved: Figure1_Lipid_Class_Comparison.png/pdf")
plt.close()

# ============================================================================
# FIGURE 2: UGCGi Treatment Effects on Glycosphingolipids
# ============================================================================

print("\nGenerating Figure 2: UGCGi Treatment Effects...")

fig2_data = pd.read_csv(RESULTS_DIR / 'figure2_ugcgi_effects.csv')

fig, ax = plt.subplots(figsize=(12, 6))

# Prepare data for grouped bar chart
lipid_classes = fig2_data['Lipid_Class'].unique()
genotypes = fig2_data['Genotype'].unique()
x = np.arange(len(lipid_classes))
width = 0.2

# Plot bars for each genotype and condition
for i, genotype in enumerate(genotypes):
    data_subset = fig2_data[fig2_data['Genotype'] == genotype]

    # Control bars
    control_pos = x + (i * 2 - 1) * width
    ax.bar(control_pos, data_subset['Control_mean'], width,
           yerr=data_subset['Control_std'],
           label=f'{genotype} Control',
           color=color_treatment['Control'][i],
           capsize=4,
           edgecolor='black',
           linewidth=1)

    # UGCGi bars
    ugcgi_pos = x + (i * 2) * width
    ax.bar(ugcgi_pos, data_subset['UGCGi_mean'], width,
           yerr=data_subset['UGCGi_std'],
           label=f'{genotype} UGCGi',
           color=color_treatment['UGCGi'][i],
           capsize=4,
           edgecolor='black',
           linewidth=1,
           hatch='//')

# Customize axes
ax.set_xlabel('Glycosphingolipid Class', fontsize=12, fontweight='bold')
ax.set_ylabel('Total Abundance (AU)', fontsize=12, fontweight='bold')
ax.set_title('Effect of UGCG Inhibition on Glycosphingolipid Abundance', 
             fontsize=14, fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(lipid_classes, fontsize=11)
ax.legend(fontsize=9, frameon=True, shadow=True, ncol=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add significance indicators
for i, lipid_class in enumerate(lipid_classes):
    for j, genotype in enumerate(genotypes):
        row = fig2_data[(fig2_data['Lipid_Class'] == lipid_class) & 
                        (fig2_data['Genotype'] == genotype)].iloc[0]
        if row['significance'] != 'ns':
            x_pos = x[i] + (j * 2 - 0.5) * width
            y_pos = max(row['Control_mean'] + row['Control_std'],
                       row['UGCGi_mean'] + row['UGCGi_std']) + 0.5
            ax.text(x_pos, y_pos, row['significance'],
                   ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'Figure2_UGCGi_Effects.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / 'Figure2_UGCGi_Effects.pdf', bbox_inches='tight')
print("✓ Saved: Figure2_UGCGi_Effects.png/pdf")
plt.close()

# ============================================================================
# FIGURE 3: PS Species Distribution (Horizontal Bar Chart)
# ============================================================================

print("\nGenerating Figure 3: PS Species Distribution...")

fig3_data = pd.read_csv(RESULTS_DIR / 'figure3_ps_species.csv')

# Filter for top species by absolute Log2FC and those with p < 0.05
fig3_filtered = fig3_data[(fig3_data['p_value'] < 0.05) | 
                          (abs(fig3_data['Log2_FC']) > 0.5)]
fig3_filtered = fig3_filtered.sort_values('Log2_FC', ascending=True).tail(15)

fig, ax = plt.subplots(figsize=(10, 8))

# Create horizontal bar chart
colors_bars = [color_kras if x > 0 else color_wt for x in fig3_filtered['Log2_FC']]
bars = ax.barh(range(len(fig3_filtered)), fig3_filtered['Log2_FC'], 
               color=colors_bars, edgecolor='black', linewidth=1)

# Customize axes
ax.set_yticks(range(len(fig3_filtered)))
ax.set_yticklabels(fig3_filtered['PS_Species'], fontsize=9)
ax.set_xlabel('Log₂ Fold Change (KRAS G12V / WT)', fontsize=12, fontweight='bold')
ax.set_title('Phosphatidylserine Species Distribution', 
             fontsize=14, fontweight='bold', pad=20)
ax.axvline(x=0, color='black', linestyle='-', linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add significance stars
for i, (idx, row) in enumerate(fig3_filtered.iterrows()):
    x_pos = row['Log2_FC']
    offset = 0.1 if x_pos > 0 else -0.1
    ax.text(x_pos + offset, i, row['significance'],
           ha='left' if x_pos > 0 else 'right',
           va='center', fontsize=9, fontweight='bold')

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=color_kras, label='Enriched in KRAS G12V'),
                   Patch(facecolor=color_wt, label='Enriched in WT')]
ax.legend(handles=legend_elements, loc='lower right', fontsize=10)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'Figure3_PS_Species.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / 'Figure3_PS_Species.pdf', bbox_inches='tight')
print("✓ Saved: Figure3_PS_Species.png/pdf")
plt.close()



# ============================================================================
# BONUS: Volcano Plot for WT vs KRAS G12V
# ============================================================================

print("\nGenerating Bonus Figure: Volcano Plot...")

fig3_data_full = pd.read_csv(RESULTS_DIR / 'figure3_ps_species.csv')

# Calculate -log10(p-value)
fig3_data_full['-log10_pval'] = -np.log10(fig3_data_full['p_value'].replace(0, 1e-300))

# Create volcano plot
fig, ax = plt.subplots(figsize=(10, 8))

# Color points based on significance and fold change
colors_volcano = []
for _, row in fig3_data_full.iterrows():
    if row['p_value'] < 0.05 and row['Log2_FC'] > 0.5:
        colors_volcano.append(color_kras)  # Enriched in KRAS
    elif row['p_value'] < 0.05 and row['Log2_FC'] < -0.5:
        colors_volcano.append(color_wt)  # Enriched in WT
    else:
        colors_volcano.append('gray')  # Not significant

ax.scatter(fig3_data_full['Log2_FC'], 
          fig3_data_full['-log10_pval'],
          c=colors_volcano,
          s=50,
          alpha=0.6,
          edgecolors='black',
          linewidth=0.5)

# Add threshold lines
ax.axvline(x=0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.axvline(x=-0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=1, alpha=0.5)

# Customize
ax.set_xlabel('Log₂ Fold Change (KRAS G12V / WT)', fontsize=12, fontweight='bold')
ax.set_ylabel('-Log₁₀ (p-value)', fontsize=12, fontweight='bold')
ax.set_title('Phosphatidylserine Species: WT vs KRAS G12V', 
             fontsize=14, fontweight='bold', pad=20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=color_kras, label='Enriched in KRAS G12V (p<0.05, FC>1.4)'),
    Patch(facecolor=color_wt, label='Enriched in WT (p<0.05, FC<0.7)'),
    Patch(facecolor='gray', label='Not significant')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'Figure_Bonus_Volcano.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / 'Figure_Bonus_Volcano.pdf', bbox_inches='tight')
print("✓ Saved: Figure_Bonus_Volcano.png/pdf")
plt.close()

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("ALL FIGURES GENERATED SUCCESSFULLY!")
print("=" * 70)
print("\nGenerated files:")
print("  • Figure1_Lipid_Class_Comparison.png/pdf")
print("  • Figure2_UGCGi_Effects.png/pdf")
print("  • Figure3_PS_Species.png/pdf")
print("  • Figure_Bonus_Volcano.png/pdf")
print("\nAll figures are publication-ready at 300 DPI.")
print("=" * 70)