"""
Step 3: Signature scoring analysis using ssGSEA (Python implementation)
Correlate METRN expression with 7 pathway signature scores in GSE124272
"""

import numpy as np
import pandas as pd
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# Try gsva package; if not working, use manual ssGSEA
out_dir = "D:/MedResearch/METRN_IVDD"

# ---- Load gene-level expression matrix from R output ----
# We have the FDR CSV files from Step 1 which have probe-level data
# But we need gene-level. Let's reconstruct from the series matrix.

# Since R is unstable, use GEOparse in Python
import GEOparse

gse1 = GEOparse.get_GEO(geo="GSE124272", destdir=os.path.join(out_dir, "geo_cache"), silent=True)

# Extract expression matrix
expr1 = gse1.pivot_samples('VALUE').T  # samples as rows, probes as columns
# Actually GEOparse gives probes as index, samples as columns
expr1_raw = gse1.pivot_samples('VALUE')
print(f"Expression matrix: {expr1_raw.shape}")

# Map probes to gene symbols
probe_to_gene = {}
for gsm_name, gsm in gse1.gsms.items():
    pass  # Not needed, use table_data

# Use the platform annotation
gpl = gse1.gpls[gse1.metadata['platform_id'][0]]
print(f"Platform: {gpl.metadata['platform_id']}")
print(f"GPL table columns: {list(gpl.table.columns)}")

# Create probe-to-gene mapping
if 'GENE_SYMBOL' in gpl.table.columns:
    probe_gene_map = dict(zip(gpl.table['ID'], gpl.table['GENE_SYMBOL']))
elif 'Gene Symbol' in gpl.table.columns:
    probe_gene_map = dict(zip(gpl.table['ID'], gpl.table['Gene Symbol']))
else:
    # Search for gene symbol column
    for col in gpl.table.columns:
        if 'gene' in col.lower() and 'symbol' in col.lower():
            probe_gene_map = dict(zip(gpl.table['ID'], gpl.table[col]))
            print(f"Using column: {col}")
            break

# Filter out empty gene symbols
probe_gene_map = {k: v for k, v in probe_gene_map.items() if v and v.strip() != ''}

# Collapse probes to genes (median of multiple probes)
print(f"Mapped probes: {len(probe_gene_map)}")

# Create gene-level expression
gene_expr_dict = {}
for probe, gene in probe_gene_map.items():
    if probe in expr1_raw.index:
        vals = expr1_raw.loc[probe].values
        if gene in gene_expr_dict:
            # Take median of multiple probes
            existing = gene_expr_dict[gene]
            gene_expr_dict[gene] = np.median([existing, vals], axis=0)
        else:
            gene_expr_dict[gene] = vals

# Build gene expression DataFrame
sample_names = expr1_raw.columns
gene_expr = pd.DataFrame(gene_expr_dict, index=sample_names).T
print(f"Gene-level matrix: {gene_expr.shape}")

# ---- Define gene sets ----
gene_sets = {
    'Mechanotransduction': ['PIEZO1', 'PIEZO2', 'TRPV4', 'YAP1', 'WWTR1', 
        'CTGF', 'CYR61', 'ITGA1', 'ITGA2', 'ITGA3', 'ITGA5', 'ITGA7',
        'ITGB1', 'ITGB3', 'ITGB5', 'ITGB6',
        'PTK2', 'RHOA', 'ROCK1', 'ROCK2', 'LIMK1', 'CFL1', 
        'VCL', 'TLN1', 'FOSB', 'CREB3L2', 'TMEM102'],
    
    'ECM_Homeostasis': ['COL2A1', 'ACAN', 'SOX9', 'COL1A1', 'FN1', 'COMP', 
        'MMP1', 'MMP3', 'MMP13', 'ADAMTS4', 'ADAMTS5', 'COL10A1', 'COL3A1'],
    
    'Autophagy': ['ATG5', 'ATG7', 'ATG12', 'ULK1', 'BECN1', 
        'MAP1LC3B', 'SQSTM1', 'ATG3', 'ATG4B', 'ATG16L1'],
    
    'Senescence': ['CDKN1A', 'CDKN2A', 'TP53', 'IL6', 'CXCL8', 
        'MMP1', 'MMP3', 'MMP13', 'SERPINE1', 'RB1', 'CCND1'],
    
    'ER_Stress': ['CALR', 'HSPA5', 'ATF4', 'DDIT3', 'XBP1', 
        'ERN1', 'CREB3L2', 'HSP90B1', 'PDIA3', 'EIF2AK3'],
    
    'Inflammation': ['IL1B', 'IL6', 'TNF', 'NFKB1', 'CXCL8', 
        'CCL2', 'TLR2', 'TLR4', 'IL10', 'TNFAIP3'],
    
    'Ferroptosis': ['GPX4', 'SLC7A11', 'ACSL4', 'FTH1', 'TFRC', 
        'NCOA4', 'ALOX15', 'POR', 'FDFT1', 'HSPB1']
}

# Filter gene sets for available genes
available_genes = set(gene_expr.index)
filtered_sets = {}
for name, genes in gene_sets.items():
    found = [g for g in genes if g in available_genes]
    missing = [g for g in genes if g not in available_genes]
    print(f"{name}: found {len(found)}/{len(genes)} genes, missing: {missing}")
    if len(found) >= 2:
        filtered_sets[name] = found

print(f"\nUsable gene sets: {list(filtered_sets.keys())}")

# ---- Manual ssGSEA implementation ----
def ssgsea_score(expr_matrix, gene_set):
    """
    Simplified ssGSEA: rank-based enrichment score
    For each sample, rank genes by expression, then compute 
    running sum of gene set membership.
    """
    scores = {}
    for sample in expr_matrix.columns:
        # Rank genes by expression (descending)
        sample_vals = expr_matrix[sample].sort_values(ascending=False)
        ranks = pd.Series(range(1, len(sample_vals) + 1), index=sample_vals.index)
        
        # Compute enrichment score
        N = len(ranks)
        Nh = len([g for g in gene_set if g in ranks.index])
        if Nh == 0:
            scores[sample] = 0
            continue
        
        # Running sum
        hit_indices = ranks.index.isin(gene_set)
        running_sum = 0
        max_sum = 0
        
        for i, gene in enumerate(ranks.index):
            if gene in gene_set:
                running_sum += ranks[gene] / sum(ranks[gene] for g in gene_set if g in ranks.index)
            else:
                running_sum -= 1.0 / (N - Nh)
            if abs(running_sum) > abs(max_sum):
                max_sum = running_sum
        
        scores[sample] = max_sum
    
    return pd.Series(scores)

# Alternative: simple mean-based scoring (more robust for small samples)
def mean_signature_score(expr_matrix, gene_set):
    """Mean z-score of gene set members"""
    available = [g for g in gene_set if g in expr_matrix.index]
    if len(available) < 2:
        return pd.Series(0, index=expr_matrix.columns)
    
    subset = expr_matrix.loc[available]
    # Z-score each gene across samples, then average
    z_scores = (subset - subset.mean(axis=1).values.reshape(-1, 1)) / subset.std(axis=1).values.reshape(-1, 1)
    # Replace NaN (constant genes) with 0
    z_scores = z_scores.fillna(0)
    return z_scores.mean(axis=0)

# ---- Compute signature scores ----
print("\n=== Computing signature scores ===")
signature_scores = pd.DataFrame(index=gene_expr.columns)
for name, genes in filtered_sets.items():
    # Use mean z-score method (more robust)
    scores = mean_signature_score(gene_expr, genes)
    signature_scores[name] = scores
    print(f"{name}: computed")

print(f"Signature scores shape: {signature_scores.shape}")

# ---- Correlate METRN with signature scores ----
print("\n=== METRN-Signature correlations ===")
metrn_vals = gene_expr.loc['METRN']

cor_results = []
for sig_name in signature_scores.columns:
    sig_vals = signature_scores[sig_name]
    r, p = stats.spearmanr(metrn_vals.values, sig_vals.values)
    cor_results.append({
        'Dataset': 'GSE124272',
        'Signature': sig_name,
        'Spearman_r': r,
        'P_value': p,
        'n_genes': len(filtered_sets[sig_name])
    })
    print(f"{sig_name}: r={r:.3f}, P={p:.4f}")

cor_df = pd.DataFrame(cor_results)

# FDR correction
from scipy.stats import false_discovery_control
cor_df['FDR'] = false_discovery_control(cor_df['P_value'].values, method='bh') if len(cor_df) > 0 else [1.0]*len(cor_df)

print("\nWith FDR correction:")
print(cor_df.to_string())

# ---- METRN-PIEZO1 correlation ----
print("\n=== METRN-PIEZO1 correlation ===")
piezo1_vals = gene_expr.loc['PIEZO1']
r_all, p_all = stats.spearmanr(metrn_vals.values, piezo1_vals.values)
r_pearson, p_pearson = stats.pearsonr(metrn_vals.values, piezo1_vals.values)
print(f"All samples: Spearman r={r_all:.3f}, P={p_all:.4f}")
print(f"All samples: Pearson r={r_pearson:.3f}, P={p_pearson:.4f}")

# ---- Top co-expressed genes ----
print("\n=== Top METRN co-expressed genes ===")
gene_cors = []
for gene in gene_expr.index:
    if gene == 'METRN':
        continue
    vals = gene_expr.loc[gene].values
    if np.std(vals) == 0:
        continue
    r, p = stats.spearmanr(metrn_vals.values, vals)
    gene_cors.append({'Gene': gene, 'Spearman_r': r, 'P_value': p})

gene_cor_df = pd.DataFrame(gene_cors).sort_values('Spearman_r', key=abs, ascending=False)

# FDR for all gene correlations
if len(gene_cor_df) > 0:
    gene_cor_df['FDR'] = false_discovery_control(gene_cor_df['P_value'].values, method='bh')

# Top 20
top20 = gene_cor_df.head(20)
print("Top 20 co-expressed genes:")
print(top20.to_string())

# Key genes
key_genes = ['PIEZO1', 'CALR', 'LIMK1', 'FOSB', 'CREB3L2', 'TMEM102', 'MESP1', 'ATG12', 'ULK1', 'COL2A1', 'MMP13', 'SOX9']
print("\nKey gene correlations:")
for g in key_genes:
    row = gene_cor_df[gene_cor_df['Gene'] == g]
    if len(row) > 0:
        print(f"  {g}: r={row['Spearman_r'].values[0]:.3f}, P={row['P_value'].values[0]:.4f}, FDR={row['FDR'].values[0]:.4f}")
    else:
        print(f"  {g}: not found")

# ---- Save all results ----
print("\n=== Saving results ===")
cor_df.to_csv(os.path.join(out_dir, 'signature_correlation_results.csv'), index=False)
signature_scores.to_csv(os.path.join(out_dir, 'GSVA_scores_GSE124272.csv'))
gene_cor_df.to_csv(os.path.join(out_dir, 'METRN_coexpression_all_genes_GSE124272.csv'))
top20.to_csv(os.path.join(out_dir, 'METRN_top20_coexpressed_GSE124272.csv'), index=False)

# ---- Generate plots ----
print("\n=== Generating plots ===")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Bar plot of METRN-signature correlations
fig, ax = plt.subplots(figsize=(8, 5))
sig_names = cor_df['Signature'].values
r_vals = cor_df['Spearman_r'].values
fdr_vals = cor_df['FDR'].values

colors = ['#F44336' if f < 0.05 else '#9E9E9E' for f in fdr_vals]
bars = ax.bar(range(len(sig_names)), r_vals, color=colors)
ax.set_xticks(range(len(sig_names)))
ax.set_xticklabels(sig_names, rotation=45, ha='right')
ax.axhline(y=0, color='gray', linestyle='--')
ax.set_ylabel('Spearman r')
ax.set_title('METRN-Signature Score Correlations (GSE124272)')
# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#F44336', label='FDR < 0.05'),
                   Patch(facecolor='#9E9E9E', label='FDR >= 0.05')]
ax.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'METRN_signature_correlation_GSE124272.png'), dpi=150)
plt.savefig(os.path.join(out_dir, 'METRN_signature_correlation_GSE124272.pdf'))
print("Bar plot saved")

# METRN-PIEZO1 scatter
fig, ax = plt.subplots(figsize=(7, 5))
# Get labels
pdata = gse1.phenotype_data
labels = []
for gsm in metrn_vals.index:
    title = pdata.loc[gsm, 'title'] if gsm in pdata.index else ''
    if 'patient' in title.lower() or 'IDD' in title:
        labels.append('IDD')
    elif 'volunteer' in title.lower() or 'healthy' in title.lower() or 'control' in title.lower():
        labels.append('Control')
    else:
        labels.append('Unknown')

idd_mask = [l == 'IDD' for l in labels]
ctrl_mask = [l == 'Control' for l in labels]

ax.scatter(metrn_vals.values[idd_mask], piezo1_vals.values[idd_mask], c='#F44336', s=60, label='IDD', zorder=3)
ax.scatter(metrn_vals.values[ctrl_mask], piezo1_vals.values[ctrl_mask], c='#4CAF50', s=60, label='Control', zorder=3)

# Add regression line
from numpy.polynomial.polynomial import polyfit
b, m = polyfit(metrn_vals.values, piezo1_vals.values, 1)
x_line = np.linspace(metrn_vals.values.min(), metrn_vals.values.max(), 100)
ax.plot(x_line, b + m * x_line, 'k-', linewidth=1)

ax.set_xlabel('METRN Expression')
ax.set_ylabel('PIEZO1 Expression')
ax.set_title(f'METRN vs PIEZO1 (GSE124272)\nSpearman r={r_all:.3f}, P={p_all:.4f}')
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'METRN_PIEZO1_scatter_GSE124272.png'), dpi=150)
plt.savefig(os.path.join(out_dir, 'METRN_PIEZO1_scatter_GSE124272.pdf'))
print("Scatter plot saved")

print("\n=== Step 3 Complete ===")
print(f"All output saved to: {out_dir}")
