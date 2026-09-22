# Four-layer ligand–receptor evidence integration
suppressPackageStartupMessages({library(CellChat); library(dplyr)})

ligand_sc <- read.delim("results/DEG_CAF_High_vs_Low.tsv")
receptor_sc <- read.delim("results/DEG_Hi_vs_other_epithelial.tsv")
spatial <- read.delim("results/CAFI_spatial_DE.tsv")
db <- readRDS("results/CellChat_secreted_database.rds")

candidates <- subset(db$interaction, annotation == "Secreted Signaling")
ligand_genes <- ligand_sc$gene

# Retain pairs whose ligand subunits all have CAF differential-expression estimates
keep <- logical(nrow(candidates))
for (i in seq_len(nrow(candidates))) {
  ligand <- candidates$ligand[i]
  subunits <- if (ligand %in% rownames(db$complex)) {
    x <- unlist(db$complex[ligand, ], use.names = FALSE)
    unique(x[!is.na(x) & x != ""])
  } else ligand
  keep[i] <- all(subunits %in% ligand_genes)
}
candidates <- candidates[keep, ]

get_evidence <- function(molecule, table, fc_col, p_col) {
  subunits <- if (molecule %in% rownames(db$complex)) {
    x <- unlist(db$complex[molecule, ], use.names = FALSE)
    unique(x[!is.na(x) & x != ""])
  } else molecule
  x <- table[match(subunits, table$gene), ]
  if (any(!is.finite(x[[fc_col]])) || any(!is.finite(x[[p_col]]))) {
    return(c(fc = NA_real_, padj = NA_real_))
  }
  c(fc = min(x[[fc_col]]), padj = max(x[[p_col]]))
}

rows <- lapply(seq_len(nrow(candidates)), function(i) {
  ligand <- candidates$ligand[i]; receptor <- candidates$receptor[i]
  c(
    interaction_name = candidates$interaction_name[i],
    pathway_name = candidates$pathway_name[i],
    get_evidence(ligand, ligand_sc, "avg_log2FC", "p_val_adj"),
    get_evidence(receptor, receptor_sc, "avg_log2FC", "p_val_adj"),
    get_evidence(ligand, spatial, "log2FoldChange", "padj"),
    get_evidence(receptor, spatial, "log2FoldChange", "padj")
  )
})
evidence <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(evidence) <- c("interaction", "pathway", "ligand_sc_fc", "ligand_sc_p",
                     "receptor_sc_fc", "receptor_sc_p", "ligand_spatial_fc",
                     "ligand_spatial_p", "receptor_spatial_fc", "receptor_spatial_p")

score_cols <- list(
  c("ligand_sc_fc", "ligand_sc_p"), c("receptor_sc_fc", "receptor_sc_p"),
  c("ligand_spatial_fc", "ligand_spatial_p"), c("receptor_spatial_fc", "receptor_spatial_p")
)
score_matrix <- matrix(NA_real_, nrow(evidence), 4)
complete <- rep(TRUE, nrow(evidence))
for (k in seq_along(score_cols)) {
  fc <- as.numeric(evidence[[score_cols[[k]][1]]])
  p <- as.numeric(evidence[[score_cols[[k]][2]]])
  known <- is.finite(fc) & is.finite(p)
  weight <- rep(0, length(p))
  weight[known & p <= 0.05] <- 1
  weight[known & p > 0.05] <- pmax(0.05, -log10(p[known & p > 0.05]) / -log10(0.05))
  raw <- rep(0.001, length(p))
  raw[known] <- pmax(0.001, pmax(fc[known], 0) * weight[known])
  score_matrix[, k] <- pmin(1, pmax(0.001, raw / quantile(raw, 0.95)))
  complete <- complete & known & fc > 0 & p < 0.05
}
harmonic <- 4 / rowSums(1 / score_matrix)
geometric <- exp(rowMeans(log(score_matrix)))
evidence$balanced_score <- 100 * sqrt(harmonic * geometric)
evidence$complete_evidence <- complete
write.table(evidence, "results/LR_four_layer_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Pathway integration: prioritise complete pairs and use other pairs for the soft score
pathway_rows <- lapply(split(evidence, evidence$pathway), function(x) {
  all_scores <- sort(x$balanced_score, decreasing = TRUE)
  strict_scores <- sort(x$balanced_score[x$complete_evidence], decreasing = TRUE)

  n_soft <- min(3, length(all_scores))
  soft_weights <- c(0.6, 0.3, 0.1)[seq_len(n_soft)]
  soft_score <- sum(all_scores[seq_len(n_soft)] * soft_weights) / sum(soft_weights)

  n_complete <- length(strict_scores)
  strict_score <- 0
  if (n_complete > 0) {
    n_strict <- min(2, n_complete)
    strict_weights <- c(0.7, 0.3)[seq_len(n_strict)]
    strict_score <- sum(strict_scores[seq_len(n_strict)] * strict_weights) / sum(strict_weights)
  }

  integrated <- if (n_complete > 0) {
    min(100, strict_score * (1 + 0.15 * min(n_complete, 2) / 2) + 0.25 * soft_score)
  } else {
    min(30, 0.25 * soft_score)
  }

  data.frame(
    pathway = x$pathway[1],
    n_pairs = nrow(x),
    n_complete = n_complete,
    strict_score = strict_score,
    soft_score = soft_score,
    integrated_score = integrated
  )
})
pathway_scores <- do.call(rbind, pathway_rows)
pathway_scores <- pathway_scores[order(
  -(pathway_scores$n_complete > 0), -pathway_scores$integrated_score,
  -pathway_scores$n_complete, -pathway_scores$n_pairs
), ]
write.table(pathway_scores, "results/LR_pathway_scores.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
