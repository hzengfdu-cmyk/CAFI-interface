# Invadopodia scoring, epithelial states and sample groups
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(dplyr)})
set.seed(42)

input_file <- "data/single_cell_annotated.rds"
malignancy_file <- "results/epithelial_malignancy.tsv"
output_file <- "results/single_cell_states.rds"

sce <- readRDS(input_file)
mal <- read.delim(malignancy_file, row.names = 1)
sce$malignancy <- "Not_epithelial"
sce$malignancy[match(rownames(mal), colnames(sce))] <- mal$malignancy

signature23 <- c(
  "ACTB", "ACTG1", "ACTR2", "ACTR3", "WASL", "WIPF1", "CDC42", "CTTN",
  "CFL1", "DIAPH1", "DIAPH3", "DNM2", "FSCN1", "CSRP2", "MMP14",
  "SH3PXD2A", "ENAH", "ITGB1", "ITGB3", "VCL", "PXN", "ZYX", "ILK"
)

mal_cells <- colnames(sce)[sce$malignancy == "Malignant"]
mal_obj <- subset(sce, cells = mal_cells)
counts <- GetAssayData(mal_obj, assay = "RNA", layer = "counts")
frequency <- Matrix::rowMeans(counts > 0)
signature <- intersect(signature23, names(frequency)[frequency >= 0.05 & frequency <= 0.95])

mal_obj <- AddModuleScore(mal_obj, features = list(signature), assay = "RNA",
                          name = "Invadopodia", seed = 42)
cutoff <- quantile(mal_obj$Invadopodia1, 0.75)

sce$inva_score <- NA_real_
sce$inva_score[match(colnames(mal_obj), colnames(sce))] <- mal_obj$Invadopodia1
sce$state <- as.character(sce$celltype)
sce$state[sce$celltype == "Epithelial"] <- "Normal_Epithelial"
sce$state[match(colnames(mal_obj), colnames(sce))] <- ifelse(
  mal_obj$Invadopodia1 >= cutoff, "Hi_inva_Epithelial", "Low_inva_Epithelial"
)

# Sample-level grouping
meta <- sce[[]]
sample_table <- meta %>%
  group_by(sample_id) %>%
  summarise(
    n_epi = sum(state %in% c("Hi_inva_Epithelial", "Low_inva_Epithelial", "Normal_Epithelial")),
    n_hi = sum(state == "Hi_inva_Epithelial"),
    n_caf = sum(state == "Fibroblast"),
    n_non_epi = sum(!state %in% c("Hi_inva_Epithelial", "Low_inva_Epithelial", "Normal_Epithelial")),
    .groups = "drop"
  ) %>%
  mutate(inva_fraction = n_hi / n_epi, caf_fraction = n_caf / n_non_epi)

inva_cut <- median(sample_table$inva_fraction[sample_table$n_epi >= 20], na.rm = TRUE)
caf_cut <- median(sample_table$caf_fraction[sample_table$n_caf >= 3], na.rm = TRUE)
sample_table$inva_group <- ifelse(sample_table$n_epi >= 20 & sample_table$inva_fraction > inva_cut, "High", "Low")
sample_table$caf_group <- ifelse(sample_table$n_caf >= 3 & sample_table$caf_fraction > caf_cut, "High", "Low")

sce$inva_group <- sample_table$inva_group[match(sce$sample_id, sample_table$sample_id)]
sce$caf_group <- sample_table$caf_group[match(sce$sample_id, sample_table$sample_id)]

write.table(sample_table, "results/sample_groups.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(sce, output_file)
