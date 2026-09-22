# Spatial CAFI pseudobulk, DESeq2 and TIDE input
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(DESeq2)})

spatial_list <- readRDS("results/spatial_CAFI_regions.rds")
sample_summary <- read.delim("results/CAFI_sample_summary.tsv")
keep_samples <- sample_summary$sample_id[sample_summary$interface_spots > 10]
spatial_list <- spatial_list[keep_samples]

genes <- Reduce(intersect, lapply(spatial_list, function(x) rownames(x[["Spatial"]])))
pseudobulk <- list(); metadata <- list(); interface_detected <- setNames(numeric(length(genes)), genes)
n_interface <- 0

for (sample_id in names(spatial_list)) {
  obj <- spatial_list[[sample_id]]
  counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")[genes, ]
  groups <- ifelse(obj$is_interface, "Interface", "Non_interface")
  for (group in c("Interface", "Non_interface")) {
    idx <- which(groups == group)
    if (length(idx) >= 10) {
      name <- paste(sample_id, group, sep = "__")
      pseudobulk[[name]] <- Matrix::rowSums(counts[, idx, drop = FALSE])
      metadata[[name]] <- data.frame(profile = name, sample_id, region = group)
    }
  }
  interface_detected <- interface_detected + Matrix::rowSums(counts[, obj$is_interface, drop = FALSE] > 0)
  n_interface <- n_interface + sum(obj$is_interface)
}

count_matrix <- do.call(cbind, pseudobulk)
meta <- do.call(rbind, metadata)
rownames(meta) <- meta$profile
# Retain sections containing both Interface and Non-interface profiles
paired_samples <- names(which(table(meta$sample_id) == 2))
meta <- meta[meta$sample_id %in% paired_samples, , drop = FALSE]
count_matrix <- count_matrix[, meta$profile, drop = FALSE]
meta$sample_id <- factor(meta$sample_id)
meta$region <- factor(meta$region, levels = c("Non_interface", "Interface"))

dds <- DESeqDataSetFromMatrix(round(count_matrix), meta[colnames(count_matrix), ],
                              design = ~ sample_id + region)
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("region", "Interface", "Non_interface")))
res$gene <- rownames(res)
res$detection_fraction <- interface_detected[res$gene] / n_interface
res$signature <- res$padj < 0.05 & res$log2FoldChange > 0.5 & res$detection_fraction > 0.5
write.table(res, "results/CAFI_spatial_DE.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(res$gene[res$signature], "results/CAFI_signature.txt")

# TIDE input: sum raw counts within Interface, Proximal and Distal regions
tide_profiles <- list()
for (sample_id in names(spatial_list)) {
  obj <- spatial_list[[sample_id]]
  counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")[genes, ]
  for (region in c("Interface", "Proximal", "Distal")) {
    idx <- which(obj$region == region)
    if (length(idx) > 0) {
      tide_profiles[[paste(sample_id, region, sep = "__")]] <- Matrix::rowSums(counts[, idx, drop = FALSE])
    }
  }
}
tide_counts <- do.call(cbind, tide_profiles)
log_cpm <- log2(sweep(tide_counts, 2, colSums(tide_counts), "/") * 1e6 + 1)
tide_input <- log_cpm - rowMeans(log_cpm)
write.table(tide_input, "results/spatial_TIDE_input.tsv", sep = "\t", col.names = NA, quote = FALSE)
