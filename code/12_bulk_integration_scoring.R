# Platform-specific bulk integration, ComBat and signature scoring
suppressPackageStartupMessages({library(sva); library(IOBR)})
set.seed(42)

bulk <- readRDS("data/bulk_cohorts.rds")
platform_table <- read.delim("data/bulk_platforms.tsv")
gene_sets <- read.delim("data/gene_sets.tsv")

combat_data <- list()
for (platform in c("RNA-seq", "Microarray")) {
  datasets <- platform_table$dataset[platform_table$platform == platform]
  matrices <- bulk[datasets]
  common_genes <- Reduce(intersect, lapply(matrices, rownames))
  matrices <- lapply(matrices, function(x) x[common_genes, , drop = FALSE])
  merged <- do.call(cbind, matrices)
  batch <- rep(datasets, vapply(matrices, ncol, integer(1)))
  combat_data[[platform]] <- ComBat(merged, batch = batch, mod = NULL, prior.plots = FALSE)
  write.table(combat_data[[platform]], paste0("results/bulk_ComBat_", platform, ".tsv"),
              sep = "\t", col.names = NA, quote = FALSE)
}
saveRDS(combat_data, "results/bulk_ComBat.rds")

# Expression-matched background module scores
score_list <- list()
for (platform in names(combat_data)) {
  expr <- combat_data[[platform]]
  nbin <- min(24, floor(nrow(expr) / 100))
  order_gene <- order(rowMeans(expr), rownames(expr))
  bins <- setNames(floor((seq_len(nrow(expr)) - 1) * nbin / nrow(expr)) + 1, rownames(expr)[order_gene])
  bins <- bins[rownames(expr)]

  platform_scores <- data.frame(sample_id = colnames(expr))
  for (set_name in c("CAFI interface", "Invadopodiahi malignant cell")) {
    genes <- intersect(gene_sets$gene[gene_sets$set_name == set_name], rownames(expr))
    controls <- unique(unlist(lapply(genes, function(gene) {
      sample(names(bins)[bins == bins[gene]], 100, replace = FALSE)
    })))
    score <- colMeans(expr[genes, , drop = FALSE]) - colMeans(expr[controls, , drop = FALSE])
    platform_scores[[make.names(set_name)]] <- score
  }

  # Estimate tumour purity from the same ComBat-corrected matrix
  estimate <- deconvo_tme(expr, method = "estimate", arrays = platform == "Microarray")
  purity <- estimate$TumorPurity_estimate[match(platform_scores$sample_id, estimate$ID)]
  platform_scores$tumor_purity <- purity
  platform_scores$CAFI_purity_weighted <- platform_scores$CAFI.interface * purity
  platform_scores$InvaHi_purity_weighted <- platform_scores$Invadopodiahi.malignant.cell * purity
  score_list[[platform]] <- platform_scores
}

write.table(do.call(rbind, score_list), "results/bulk_signature_scores.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
