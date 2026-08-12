import os
import GEOparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ── Output directories ─────────────────────────────────────────────────────
os.makedirs("data/geo", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# ── Step 1: Load GEO dataset ───────────────────────────────────────────────
print("=" * 60)
print("Step 1: Loading GSE10072...")
print("=" * 60)

gse = GEOparse.get_GEO(geo="GSE10072", destdir="data/geo")
pivot = gse.pivot_samples('VALUE')

print(f"Samples: {len(gse.gsms)}")
print(f"Expression matrix shape: {pivot.shape}")

# Build metadata table
metadata_rows = []
for gsm_name, gsm in gse.gsms.items():
    title = gsm.metadata['title'][0]
    source = gsm.metadata.get('source_name_ch1', ['N/A'])[0]
    metadata_rows.append({'sample': gsm_name, 'title': title, 'source': source})

metadata_df = pd.DataFrame(metadata_rows)
print(f"\nGroups found: {metadata_df['source'].unique()}")

# Split into tumor and normal
tumor = pivot[metadata_df[metadata_df['source'] == 'Adenocarcinoma of the Lung']['sample'].tolist()]
normal = pivot[metadata_df[metadata_df['source'] == 'Normal Lung Tissue']['sample'].tolist()]
print(f"Tumor: {tumor.shape[1]} samples | Normal: {normal.shape[1]} samples\n")

# ── Step 2: Differential expression — Welch's t-test + BH-FDR ────────────
print("=" * 60)
print("Step 2: Running differential expression analysis...")
print("=" * 60)

results = []
for probe in pivot.index:
    fc = tumor.loc[probe].mean() - normal.loc[probe].mean()
    _, p = stats.ttest_ind(tumor.loc[probe], normal.loc[probe], equal_var=False)
    results.append([probe, fc, p])

deg_df = pd.DataFrame(results, columns=['probe', 'log2FC', 'pvalue'])
_, deg_df['padj'], _, _ = multipletests(deg_df['pvalue'], method='fdr_bh')

sig = deg_df[(abs(deg_df['log2FC']) >= 1) & (deg_df['padj'] < 0.05)].sort_values('padj')
print(f"Significant DEGs found: {len(sig)}")

sig.to_csv("results/DEGs.csv", index=False)
print("Saved: results/DEGs.csv\n")

# ── Step 3: Probe-to-gene symbol mapping (GPL96) ──────────────────────────
print("=" * 60)
print("Step 3: Mapping probe IDs to gene symbols...")
print("=" * 60)

gpl = gse.gpls['GPL96']
annot = gpl.table[['ID', 'Gene Symbol']].copy()
annot.columns = ['probe', 'gene_symbol']
annot = annot[annot['gene_symbol'].notna()]
annot = annot[annot['gene_symbol'] != '']

sig_named = sig.merge(annot, on='probe', how='left')
sig_named = sig_named.dropna(subset=['gene_symbol'])
sig_named = sig_named.drop_duplicates(subset='probe')

print(f"DEGs with gene names: {len(sig_named)}")
print("\nTop 10 DEGs:")
print(sig_named[['gene_symbol', 'log2FC', 'padj']].head(10))

sig_named.to_csv("results/DEGs_named.csv", index=False)
print("Saved: results/DEGs_named.csv\n")

# ── Step 4: Pathway enrichment (GSEApy Enrichr ORA) ───────────────────────
print("=" * 60)
print("Step 4: Running pathway enrichment analysis...")
print("=" * 60)

up_genes = sig_named[sig_named['log2FC'] > 0]['gene_symbol'].tolist()
down_genes = sig_named[sig_named['log2FC'] < 0]['gene_symbol'].tolist()

print(f"Upregulated genes: {len(up_genes)}")
print(f"Downregulated genes: {len(down_genes)}\n")

enr_up = gp.enrichr(
    gene_list=up_genes,
    gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023'],
    organism='human',
    outdir="results/enrichment_up"
)

enr_down = gp.enrichr(
    gene_list=down_genes,
    gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023'],
    organism='human',
    outdir="results/enrichment_down"
)

up_sig = enr_up.results[enr_up.results['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
down_sig = enr_down.results[enr_down.results['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')

print("=== TOP UPREGULATED PATHWAYS ===")
print(up_sig[['Gene_set', 'Term', 'Adjusted P-value']].head(5).to_string())
print("\n=== TOP DOWNREGULATED PATHWAYS ===")
print(down_sig[['Gene_set', 'Term', 'Adjusted P-value']].head(5).to_string())

up_sig.to_csv("results/pathways_up.csv", index=False)
down_sig.to_csv("results/pathways_down.csv", index=False)
print("\nSaved: results/pathways_up.csv, results/pathways_down.csv\n")

# ── Step 5a: Volcano plot ──────────────────────────────────────────────────
print("=" * 60)
print("Step 5: Generating figures...")
print("=" * 60)

colors = []
for _, row in deg_df.iterrows():
    if row['padj'] < 0.05 and row['log2FC'] >= 1:
        colors.append('red')
    elif row['padj'] < 0.05 and row['log2FC'] <= -1:
        colors.append('blue')
    else:
        colors.append('grey')

plt.figure(figsize=(10, 7))
plt.scatter(deg_df['log2FC'], -np.log10(deg_df['pvalue']),
            c=colors, alpha=0.4, s=8)
plt.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
plt.axvline(1, color='black', linestyle='--', linewidth=0.8)
plt.axvline(-1, color='black', linestyle='--', linewidth=0.8)

top_genes = sig_named.nlargest(5, 'log2FC')[['gene_symbol', 'log2FC', 'pvalue']]
bot_genes = sig_named.nsmallest(5, 'log2FC')[['gene_symbol', 'log2FC', 'pvalue']]
for _, row in pd.concat([top_genes, bot_genes]).iterrows():
    plt.annotate(row['gene_symbol'],
                 (row['log2FC'], -np.log10(row['pvalue'])),
                 fontsize=8, fontweight='bold')

plt.xlabel('log2 Fold Change', fontsize=12)
plt.ylabel('-log10(p-value)', fontsize=12)
plt.title('Volcano Plot — Lung Adenocarcinoma vs Normal', fontsize=14)
plt.tight_layout()
plt.savefig("figures/volcano_plot.png", dpi=150)
plt.close()
print("Saved: figures/volcano_plot.png")

# ── Step 5b: Pathway bar chart ─────────────────────────────────────────────
up_plot = up_sig.head(8)[['Term', 'Adjusted P-value']].copy()
up_plot['Direction'] = 'Upregulated'
down_plot = down_sig.head(8)[['Term', 'Adjusted P-value']].copy()
down_plot['Direction'] = 'Downregulated'

combined = pd.concat([up_plot, down_plot])
combined['-log10(padj)'] = -np.log10(combined['Adjusted P-value'])
combined['Term'] = combined['Term'].str.split('(').str[0].str.strip().str[:55]

plt.figure(figsize=(13, 9))
colors_map = {'Upregulated': '#d62728', 'Downregulated': '#1f77b4'}
for direction, group in combined.groupby('Direction'):
    plt.barh(group['Term'], group['-log10(padj)'],
             color=colors_map[direction], alpha=0.8, label=direction)

plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
plt.title('Top Enriched Pathways — Lung Adenocarcinoma vs Normal', fontsize=14)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("figures/pathway_plot.png", dpi=150)
plt.close()
print("Saved: figures/pathway_plot.png")

# ── Step 5c: Heatmap ───────────────────────────────────────────────────────
top_up = sig_named.nlargest(12, 'log2FC').drop_duplicates('gene_symbol')
top_down = sig_named.nsmallest(13, 'log2FC').drop_duplicates('gene_symbol')
top_probes = pd.concat([top_up, top_down])

heat_data = pivot.loc[top_probes['probe'].values]
probe_to_gene = dict(zip(top_probes['probe'], top_probes['gene_symbol']))
heat_data.index = heat_data.index.map(probe_to_gene)

tumor_cols = metadata_df[metadata_df['source'] == 'Adenocarcinoma of the Lung']['sample'].tolist()
normal_cols = metadata_df[metadata_df['source'] == 'Normal Lung Tissue']['sample'].tolist()
heat_data = heat_data[tumor_cols + normal_cols]

heat_data_z = heat_data.subtract(heat_data.mean(axis=1), axis=0).divide(heat_data.std(axis=1), axis=0)

plt.figure(figsize=(16, 9))
sns.heatmap(heat_data_z, cmap='RdBu_r', center=0,
            xticklabels=False, yticklabels=True,
            linewidths=0, cbar_kws={'label': 'Z-score'})

plt.axvline(x=len(tumor_cols), color='black', linewidth=2)
plt.text(len(tumor_cols)/2, -1.2, 'Tumor (n=58)',
         ha='center', fontsize=11, fontweight='bold', color='darkred')
plt.text(len(tumor_cols) + len(normal_cols)/2, -1.2, 'Normal (n=49)',
         ha='center', fontsize=11, fontweight='bold', color='steelblue')

plt.title('Top 25 DEGs — Lung Adenocarcinoma vs Normal', fontsize=14)
plt.ylabel('Gene', fontsize=12)
plt.tight_layout()
plt.savefig("figures/heatmap.png", dpi=150)
plt.close()
print("Saved: figures/heatmap.png")

print("\nAll steps complete. Run 03_survival_analysis.py next.")
