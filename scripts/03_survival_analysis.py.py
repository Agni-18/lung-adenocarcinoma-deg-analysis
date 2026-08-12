import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# ── Output directories ─────────────────────────────────────────────────────
os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# ── Step 1: Load clinical survival data ───────────────────────────────────
print("=" * 60)
print("Step 1: Loading TCGA-LUAD clinical data...")
print("=" * 60)

clinical = pd.read_csv("data/tcga/LUAD_clinical.tsv", sep='\t')
surv = clinical[['sampleID', 'vital_status', 'days_to_death', 'days_to_last_followup']].copy()
surv['time'] = surv['days_to_death'].fillna(surv['days_to_last_followup'])
surv['event'] = (surv['vital_status'] == 'DECEASED').astype(int)
surv = surv.dropna(subset=['time'])
surv = surv[surv['time'] > 0]

print(f"Patients with survival data: {len(surv)}")
print(f"Deaths recorded: {surv['event'].sum()}\n")

# ── Step 2: Load expression data ───────────────────────────────────────────
print("=" * 60)
print("Step 2: Loading TCGA-LUAD expression data...")
print("=" * 60)

expr = pd.read_csv("data/tcga/LUAD_expr.txt", sep='\t', index_col=0)
expr = expr.drop(columns=['Entrez_Gene_Id'], errors='ignore')
print(f"Expression matrix shape: {expr.shape}\n")

# ── Step 3: Match samples ─────────────────────────────────────────────────
print("=" * 60)
print("Step 3: Matching samples across datasets...")
print("=" * 60)

# Trim TCGA sample IDs to match clinical format (first 4 fields)
expr.columns = ['-'.join(c.split('-')[:4]) for c in expr.columns]

common = list(set(surv['sampleID']) & set(expr.columns))
print(f"Matched samples: {len(common)}\n")

surv = surv[surv['sampleID'].isin(common)].set_index('sampleID')
expr = expr[common]

# ── Step 4 & 5: Kaplan-Meier survival analysis ────────────────────────────
print("=" * 60)
print("Step 4: Running Kaplan-Meier survival analysis...")
print("=" * 60)

genes_to_test = ['SPP1', 'MMP1', 'MMP12', 'TOP2A', 'AGER', 'SFTPC', 'COL11A1', 'GREM1']
available = [g for g in genes_to_test if g in expr.index]
print(f"Genes found in expression data: {available}\n")

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()
results = []

for i, gene in enumerate(available):
    gene_expr = expr.loc[gene]
    median_val = gene_expr.median()

    high_samples = gene_expr[gene_expr >= median_val].index.tolist()
    low_samples = gene_expr[gene_expr < median_val].index.tolist()

    high_surv = surv.loc[[s for s in high_samples if s in surv.index]]
    low_surv = surv.loc[[s for s in low_samples if s in surv.index]]

    if len(high_surv) < 5 or len(low_surv) < 5:
        print(f"  {gene}: not enough samples, skipping")
        continue

    result = logrank_test(
        high_surv['time'], low_surv['time'],
        event_observed_A=high_surv['event'],
        event_observed_B=low_surv['event']
    )
    p_val = result.p_value
    results.append({'gene': gene, 'p_value': round(p_val, 4)})

    ax = axes[i]
    kmf = KaplanMeierFitter()
    kmf.fit(high_surv['time'], high_surv['event'], label=f'High {gene}')
    kmf.plot_survival_function(ax=ax, color='red', ci_show=False)
    kmf.fit(low_surv['time'], low_surv['event'], label=f'Low {gene}')
    kmf.plot_survival_function(ax=ax, color='blue', ci_show=False)

    star = '★' if p_val < 0.05 else ''
    ax.set_title(f'{gene}  p={p_val:.4f} {star}', fontsize=11)
    ax.set_xlabel('Days')
    ax.set_ylabel('Survival Probability')

plt.suptitle('Kaplan-Meier Survival Curves — TCGA-LUAD', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig("figures/KM_curves.png", dpi=150)
plt.close()
print("Saved: figures/KM_curves.png\n")

# ── Summary ────────────────────────────────────────────────────────────────
results_df = pd.DataFrame(results).sort_values('p_value')

print("=== SURVIVAL RESULTS ===")
print(results_df.to_string(index=False))

print("\nSignificant genes (p < 0.05):")
print(results_df[results_df['p_value'] < 0.05].to_string(index=False))

results_df.to_csv("results/survival_results.csv", index=False)
print("\nSaved: results/survival_results.csv")
print("\nAnalysis complete.")