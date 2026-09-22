# CAFI interface and spatial region definitions
suppressPackageStartupMessages({library(Seurat); library(RANN)})

spatial_list <- readRDS("results/spatial_RCTD.rds")
labels <- c(
  "Hi_inva_Epithelial", "Low_inva_Epithelial", "Normal_Epithelial",
  "CD4T", "Treg", "CD8T", "NKT", "NK", "B", "Plasma", "Monocyte",
  "Macrophage", "DC", "Neutrophil", "Mast", "Fibroblast", "Endothelial"
)

summary_list <- list()
for (sample_id in names(spatial_list)) {
  obj <- spatial_list[[sample_id]]
  weights <- t(as.matrix(GetAssayData(obj, assay = "RCTD", layer = "data")))
  colnames(weights) <- gsub("-", "_", colnames(weights))
  weights <- weights[, labels]
  coords <- as.matrix(GetTissueCoordinates(obj)[rownames(weights), c("x", "y")])

  hi <- weights[, "Hi_inva_Epithelial"]
  caf <- weights[, "Fibroblast"]
  low <- weights[, "Low_inva_Epithelial"]
  normal <- weights[, "Normal_Epithelial"]
  dominant <- colnames(weights)[max.col(weights)]
  baseline <- 1 / 17
  joint_score <- 2 * hi * caf / (hi + caf + 1e-6)

  score_candidate <-
    hi > 3 * baseline & caf > 3 * baseline &
    joint_score > quantile(joint_score, 0.75) &
    dominant %in% c("Hi_inva_Epithelial", "Fibroblast") &
    hi > 3 * low & hi > 3 * normal

  nn <- nn2(coords, k = 19)$nn.idx
  spatial_candidate <- logical(nrow(weights))
  for (i in seq_len(nrow(weights))) {
    neighbors <- head(nn[i, nn[i, ] != i], 18)
    spatial_candidate[i] <-
      (dominant[i] == "Hi_inva_Epithelial" & any(dominant[neighbors] == "Fibroblast")) |
      (dominant[i] == "Fibroblast" & any(dominant[neighbors] == "Hi_inva_Epithelial"))
  }

  interface <- score_candidate & spatial_candidate
  region <- rep(NA_character_, nrow(weights))
  distance <- rep(NA_real_, nrow(weights))
  region[interface] <- "Interface"
  distance[interface] <- 0

  if (any(interface)) {
    distance[!interface] <- nn2(coords[interface, ], coords[!interface, ], k = 1)$nn.dists[, 1]
    q25 <- quantile(distance[!interface], 0.25)
    region[!interface] <- ifelse(distance[!interface] <= q25, "Proximal", "Distal")
  }

  obj$joint_score <- joint_score
  obj$is_interface <- interface
  obj$region <- region
  obj$distance_to_interface <- distance
  spatial_list[[sample_id]] <- obj

  summary_list[[sample_id]] <- data.frame(
    sample_id,
    total_spots = ncol(obj),
    score_spots = sum(score_candidate),
    spatial_spots = sum(spatial_candidate),
    interface_spots = sum(interface),
    regional_eligible = sum(interface) > 10
  )
}

write.table(do.call(rbind, summary_list), "results/CAFI_sample_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(spatial_list, "results/spatial_CAFI_regions.rds")
