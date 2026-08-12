import os
import urllib.request
import GEOparse

# ── Create output directories ──────────────────────────────────────────────
os.makedirs("data/geo", exist_ok=True)
os.makedirs("data/tcga", exist_ok=True)

# ── 1. Download GSE10072 from NCBI GEO ────────────────────────────────────
print("=" * 60)
print("Step 1: Downloading GSE10072 from NCBI GEO...")
print("=" * 60)

gse = GEOparse.get_GEO(geo="GSE10072", destdir="data/geo")
print(f"Samples: {len(gse.gsms)}")
print(f"Platforms: {len(gse.gpls)}")
print("GSE10072 download complete.\n")

# ── 2. Download TCGA-LUAD clinical data ───────────────────────────────────
print("=" * 60)
print("Step 2: Downloading TCGA-LUAD clinical data...")
print("=" * 60)

clinical_url = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.LUAD.sampleMap%2FLUAD_clinicalMatrix"
)
clinical_path = "data/tcga/LUAD_clinical.tsv"

urllib.request.urlretrieve(clinical_url, clinical_path)
print(f"Clinical data saved to: {clinical_path}\n")

# ── 3. Download TCGA-LUAD expression data ─────────────────────────────────
print("=" * 60)
print("Step 3: Downloading TCGA-LUAD expression data (~74MB)...")
print("This may take 2-3 minutes depending on your connection.")
print("=" * 60)

expr_url = (
    "https://media.githubusercontent.com/media/cBioPortal/datahub/master/"
    "public/luad_tcga_pan_can_atlas_2018/"
    "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
)
expr_path = "data/tcga/LUAD_expr.txt"

urllib.request.urlretrieve(expr_url, expr_path)
print(f"Expression data saved to: {expr_path}\n")

# ── Verify downloads ───────────────────────────────────────────────────────
print("=" * 60)
print("Verifying downloads...")
print("=" * 60)

files = [
    "data/geo/GSE10072_family.soft.gz",
    "data/tcga/LUAD_clinical.tsv",
    "data/tcga/LUAD_expr.txt"
]

for f in files:
    if os.path.exists(f):
        size_mb = os.path.getsize(f) / (1024 * 1024)
        print(f"  ✓ {f}  ({size_mb:.1f} MB)")
    else:
        print(f"  ✗ {f}  NOT FOUND")

print("\nAll downloads complete. Run 02_deg_analysis.py next.")
