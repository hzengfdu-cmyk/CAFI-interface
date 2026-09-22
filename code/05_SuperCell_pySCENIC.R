# SuperCell construction and pySCENIC commands
suppressPackageStartupMessages({library(Seurat); library(SuperCell); library(Matrix); library(data.table)})
set.seed(42)

sce <- readRDS("results/single_cell_states.rds")
epi <- subset(sce, cells = colnames(sce)[sce$state %in% c(
  "Hi_inva_Epithelial", "Low_inva_Epithelial", "Normal_Epithelial"
)])
epi <- NormalizeData(epi) |> FindVariableFeatures(nfeatures = 2000)

sc <- SCimplify(
  GetAssayData(epi, assay = "RNA", layer = "data"),
  gamma = 100,
  k.knn = 20,
  n.pc = 30,
  genes.use = VariableFeatures(epi),
  cell.annotation = as.character(epi$state),
  cell.split.condition = epi$state,
  do.approx = FALSE,
  seed = 42
)

raw_counts <- GetAssayData(epi, assay = "RNA", layer = "counts")
super_counts <- supercell_GE(raw_counts, groups = sc$membership, mode = "sum")
saveRDS(list(SuperCell = sc, counts = super_counts), "results/SuperCell.rds")
fwrite(as.data.table(t(as.matrix(super_counts)), keep.rownames = "Cell"),
       "results/SuperCell_counts.csv")

# Example pySCENIC commands
pyscenic_commands <- c(
  "pyscenic grn results/SuperCell.loom data/allTFs_hg38.txt --method grnboost2 --seed 42 --num_workers 40 -o results/adjacencies.tsv",
  "pyscenic ctx results/adjacencies.tsv data/hg38_10kb_v10.rankings.feather --annotations_fname data/motifs-v10nr.tbl --expression_mtx_fname results/SuperCell.loom --mask_dropouts --num_workers 40 -o results/regulons.csv",
  "pyscenic aucell results/SuperCell.loom results/regulons.csv --seed 42 --num_workers 40 -o results/SCENIC.loom"
)
writeLines(pyscenic_commands, "results/run_pySCENIC.sh")
