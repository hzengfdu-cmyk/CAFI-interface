"""MFP microenvironment classification and bulk TIDE."""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier

platform = "RNA-seq"  # Alternatively, use "Microarray"
expression_file = f"results/bulk_ComBat_{platform}.tsv"
mfp_root = Path("data/MFP")
tide_root = Path("data/TIDEpy")

sys.path.insert(0, str(mfp_root))
from portraits.utils import read_gene_sets, ssgsea_formula

expression = pd.read_csv(expression_file, sep="\t", index_col=0)
training = pd.read_csv(mfp_root / "Cohorts/Pan_TCGA/signatures.tsv", sep="\t", index_col=0).T
annotation = pd.read_csv(mfp_root / "Cohorts/Pan_TCGA/annotation.tsv", sep="\t", index_col=0)
gene_sets = read_gene_sets(str(mfp_root / "signatures/gene_signatures.gmt"))

scores = ssgsea_formula(expression.T, gene_sets)
mad = scores.sub(scores.mean(axis=0), axis=1).abs().mean(axis=0)
scaled = scores.sub(scores.median(axis=0), axis=1).div(mad, axis=1)
model = KNeighborsClassifier(n_neighbors=35)
model.fit(training.clip(-2, 2), annotation.loc[training.index, "MFP"])
mfp = model.predict(scaled.loc[:, training.columns].clip(-2, 2))
pd.DataFrame({"sample_id": scaled.index, "MFP": mfp}).to_csv(
    f"results/MFP_{platform}.tsv", sep="\t", index=False
)

# Disable repeated internal TIDE normalisation for ComBat-corrected matrices
sys.path.insert(0, str(tide_root))
from tidepy.pred import TIDE

tide = TIDE(expression, cancer="Other", ignore_norm=True, force_normalize=False)
tide.to_csv(f"results/TIDE_{platform}.tsv", sep="\t")
