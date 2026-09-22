# BANKSY domains and the epithelial-side invasive front
suppressPackageStartupMessages({
  library(Seurat); library(Banksy); library(SpatialExperiment); library(Matrix); library(RANN)
})
set.seed(42)

spatial_list <- readRDS("results/spatial_CAFI_regions.rds")
for (sample_id in names(spatial_list)) {
  obj <- spatial_list[[sample_id]]
  counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")
  coords <- as.matrix(GetTissueCoordinates(obj)[colnames(obj), c("x", "y")])
  lib <- Matrix::colSums(counts)
  normcounts <- counts %*% Diagonal(x = median(lib) / lib)

  spe <- SpatialExperiment(assays = list(normcounts = normcounts), spatialCoords = coords)
  spe <- computeBanksy(spe, assay_name = "normcounts", compute_agf = TRUE, k_geom = 15)
  spe <- runBanksyPCA(spe, use_agf = TRUE, lambda = c(0, 0.2), seed = 42)
  spe <- clusterBanksy(spe, use_agf = TRUE, lambda = c(0, 0.2), resolution = 1, seed = 42)

  cluster_col <- grep("lam0.2.*res1", colnames(colData(spe)), value = TRUE)[1]
  cluster <- as.character(colData(spe)[[cluster_col]])

  weights <- t(as.matrix(GetAssayData(obj, assay = "RCTD", layer = "data")))
  colnames(weights) <- gsub("-", "_", colnames(weights))
  category_weights <- cbind(
    Epithelial = rowSums(weights[, c("Hi_inva_Epithelial", "Low_inva_Epithelial", "Normal_Epithelial")]),
    Stromal = rowSums(weights[, c("Fibroblast", "Endothelial")]),
    Immune = rowSums(weights[, setdiff(colnames(weights), c(
      "Hi_inva_Epithelial", "Low_inva_Epithelial", "Normal_Epithelial", "Fibroblast", "Endothelial"
    ))])
  )
  cluster_mean <- aggregate(category_weights, list(cluster = cluster), mean)
  cluster_type <- setNames(colnames(category_weights)[max.col(cluster_mean[, -1])], cluster_mean$cluster)
  spot_type <- unname(cluster_type[cluster])

  nn <- nn2(coords, k = min(8, nrow(coords)))
  radius <- median(nn$nn.dists[, 2]) * 1.75
  invasive_front <- logical(nrow(coords))
  for (i in seq_len(nrow(coords))) {
    neighbors <- nn$nn.idx[i, nn$nn.idx[i, ] != i & nn$nn.dists[i, ] <= radius]
    invasive_front[i] <- spot_type[i] == "Epithelial" & any(spot_type[neighbors] == "Stromal")
  }

  obj$BANKSY_cluster <- cluster
  obj$BANKSY_category <- spot_type
  obj$BANKSY_invasive_front <- invasive_front
  spatial_list[[sample_id]] <- obj
}

saveRDS(spatial_list, "results/spatial_BANKSY.rds")
