"""Joint malignant epithelial-cell identification with inferCNVpy and scMalignantFinder."""

import numpy as np
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
from scipy import sparse
from scMalignantFinder import classifier

input_file = "data/single_cell_annotated.h5ad"
gene_position_file = "data/gene_coordinates.tsv"
training_file = "data/scMalignantFinder_training.h5ad"
feature_file = "data/scMalignantFinder_features.txt"
output_file = "results/epithelial_malignancy.tsv"

adata = sc.read_h5ad(input_file)
immune_types = [
    "CD4T", "CD8T", "Treg", "NKT", "NK", "B", "Plasma",
    "Monocyte", "Macrophage", "DC", "Neutrophil", "Mast",
]
adata = adata[adata.obs["celltype"].isin(["Epithelial"] + immune_types)].copy()

# Add genomic coordinates
pos = pd.read_csv(gene_position_file, sep="\t", index_col="gene")
pos = pos.reindex(adata.var_names)
adata.var["chromosome"] = pos["chromosome"]
adata.var["start"] = pos["start"]
adata.var["end"] = pos["end"]
keep = adata.var[["chromosome", "start", "end"]].notna().all(axis=1)
keep &= ~adata.var_names.str.startswith("HLA-")
adata = adata[:, keep].copy()

# inferCNVpy
cnv.tl.infercnv(
    adata,
    reference_key="celltype",
    reference_cat=immune_types,
    window_size=100,
    n_jobs=1,
)
xcnv = adata.obsm["X_cnv"]
cnv_score = (
    np.asarray(xcnv.multiply(xcnv).mean(axis=1)).ravel()
    if sparse.issparse(xcnv)
    else np.mean(xcnv**2, axis=1)
)
adata.obs["cnv_score"] = cnv_score

# Determine CNV thresholds separately for each sample
adata.obs["cnv_malignant"] = False
for sample in adata.obs["sample_id"].unique():
    sample_cells = adata.obs["sample_id"] == sample
    epi = sample_cells & (adata.obs["celltype"] == "Epithelial")
    ref = sample_cells & adata.obs["celltype"].isin(immune_types)
    if ref.sum() > 20:
        cutoff = adata.obs.loc[ref, "cnv_score"].mean() + 2 * adata.obs.loc[ref, "cnv_score"].std()
        adata.obs.loc[epi, "cnv_malignant"] = adata.obs.loc[epi, "cnv_score"] > cutoff
    else:
        cutoff = adata.obs.loc[epi, "cnv_score"].median()
        adata.obs.loc[epi, "cnv_malignant"] = adata.obs.loc[epi, "cnv_score"] >= cutoff

# scMalignantFinder logistic regression
epi = adata[adata.obs["celltype"] == "Epithelial"].copy()
epi.X = epi.layers["counts"].copy()
epi.raw = None
temp_file = "results/epithelial_for_scMalignantFinder.h5ad"
epi.write_h5ad(temp_file)

model = classifier.scMalignantFinder(
    test_input=temp_file,
    train_h5ad_path=training_file,
    feature_path=feature_file,
    model_method="LogisticRegression",
    norm_type=True,
    n_thread=1,
)
model.load()
prediction = model.predict()

result = adata.obs.loc[epi.obs_names, ["sample_id", "cnv_score", "cnv_malignant"]].copy()
result["scMF_prediction"] = prediction.obs.loc[result.index, "scMalignantFinder_prediction"]
result["malignancy"] = np.where(
    result["cnv_malignant"] & (result["scMF_prediction"] == "Malignant"),
    "Malignant",
    "Operational_normal",
)
result.to_csv(output_file, sep="\t", index_label="cell_id")
