"""Single-cell QC, normalisation and BBKNN integration."""

import scanpy as sc
import numpy as np

input_file = "data/single_cell_raw.h5ad"
output_file = "results/single_cell_BBKNN.h5ad"
batch_key = "dataset"
library_key = "library_id"

adata = sc.read_h5ad(input_file)
adata.layers["counts"] = adata.X.copy()

# Doublet detection
sc.external.pp.scrublet(
    adata,
    batch_key=library_key,
    expected_doublet_rate=0.05,
    threshold=0.20,
    random_state=42,
)

# Cell-level QC
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
keep = (
    (adata.obs["n_genes_by_counts"] > 200)
    & (adata.obs["n_genes_by_counts"] < 10000)
    & (adata.obs["total_counts"] > 500)
    & (adata.obs["total_counts"] < 50000)
    & (adata.obs["pct_counts_mt"] < 15)
    & (adata.obs["doublet_score"] <= 0.20)
)
adata = adata[keep].copy()

# Retain protein-coding genes detected in at least 30 cells
detected = np.asarray((adata.layers["counts"] > 0).sum(axis=0)).ravel()
keep_gene = (adata.var["gene_biotype"] == "protein_coding") & (detected >= 30)
adata = adata[:, keep_gene].copy()

# BBKNN integration
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key=batch_key)
sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=True, random_state=42)
sc.external.pp.bbknn(adata, batch_key=batch_key)
sc.tl.umap(adata, random_state=42)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden", random_state=42)

adata.write_h5ad(output_file)
