# Single-cell differential expression and GSEA
suppressPackageStartupMessages({library(Seurat); library(clusterProfiler); library(readr)})
set.seed(42)

sce <- readRDS("results/single_cell_states.rds")
gene_sets <- read.delim("data/gsea_term2gene.tsv")

epi <- subset(sce, cells = colnames(sce)[sce$state %in% c("Hi_inva_Epithelial", "Low_inva_Epithelial")])
Idents(epi) <- "state"
deg <- FindMarkers(
  epi,
  ident.1 = "Hi_inva_Epithelial",
  ident.2 = "Low_inva_Epithelial",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0,
  only.pos = FALSE
)
deg$gene <- rownames(deg)
write.table(deg, "results/DEG_Hi_vs_Low.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Invadopodiahi malignant-cell signature
signature <- subset(deg, avg_log2FC > 1 & pct.1 > 0.2 & p_val_adj < 0.05)$gene
writeLines(signature, "results/InvadopodiaHi_signature.txt")

# Rank all tested genes for GSEA
ranked <- sort(setNames(deg$avg_log2FC, deg$gene), decreasing = TRUE)
gsea <- GSEA(
  ranked,
  TERM2GENE = gene_sets[, c("term", "gene")],
  minGSSize = 5,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  seed = TRUE
)
write.table(as.data.frame(gsea), "results/GSEA_Hi_vs_Low.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
