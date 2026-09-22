# RCTD spatial deconvolution
suppressPackageStartupMessages({library(Seurat); library(spacexr); library(Matrix)})

reference_sce <- readRDS("results/single_cell_states.rds")
spatial_list <- readRDS("data/spatial_objects.rds")

reference_counts <- GetAssayData(reference_sce, assay = "RNA", layer = "counts")
reference_labels <- factor(reference_sce$state)
names(reference_labels) <- colnames(reference_sce)
reference <- Reference(reference_counts, reference_labels)

dir.create("results/RCTD", recursive = TRUE, showWarnings = FALSE)
for (sample_id in names(spatial_list)) {
  obj <- spatial_list[[sample_id]]
  counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")
  coords <- GetTissueCoordinates(obj)[colnames(counts), c("x", "y")]
  puck <- SpatialRNA(coords, counts, Matrix::colSums(counts))

  fit <- create.RCTD(puck, reference, max_cores = 4)
  fit <- run.RCTD(fit, doublet_mode = "full")
  weights <- normalize_weights(fit@results$weights)

  obj <- obj[, rownames(weights)]
  obj[["RCTD"]] <- CreateAssayObject(data = t(as.matrix(weights)))
  obj$dominant_celltype <- colnames(weights)[max.col(weights)]
  spatial_list[[sample_id]] <- obj
  saveRDS(fit, file.path("results/RCTD", paste0(sample_id, ".rds")))
}

saveRDS(spatial_list, "results/spatial_RCTD.rds")
