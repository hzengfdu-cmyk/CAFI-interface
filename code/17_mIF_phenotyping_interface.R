# mIF marker positivity, cell phenotypes, CAFI interface and density
suppressPackageStartupMessages({library(dplyr); library(RANN); library(reticulate)})

cells <- read.delim("data/mIF_final_cell_table.tsv")
pixel_size_um <- 0.325
radius_um <- 220
min_opposite <- 3

sk_filters <- import("skimage.filters")
np <- import("numpy")
markers <- c("EPCAM", "CD8", "PDGFRB", "TKS5", "CTTN", "GZMB", "PD1", "TGFB1")

# Pool cells across all cores and use the upper threshold to define positivity
for (marker in markers) {
  x <- asinh(pmax(cells[[marker]], 0) / 5)
  thresholds <- as.numeric(sk_filters$threshold_multiotsu(
    np$array(x[is.finite(x)], dtype = "float64"), classes = 3L
  ))
  cells[[paste0(marker, "_asinh5")]] <- x
  cells[[paste0(marker, "_positive")]] <- x > thresholds[2]
}

# Mutually exclusive lineages and imaging-defined invadopodiahi cells
lineage_count <- cells$EPCAM_positive + cells$CD8_positive + cells$PDGFRB_positive
cells$celltype <- "Other"
cells$celltype[lineage_count == 1 & cells$EPCAM_positive] <- "Other_Epithelial"
cells$celltype[lineage_count == 1 & cells$CD8_positive] <- "CD8T"
cells$celltype[lineage_count == 1 & cells$PDGFRB_positive] <- "Fibroblast"
cells$celltype[lineage_count == 1 & cells$EPCAM_positive &
                 cells$TKS5_positive & cells$CTTN_positive] <- "Hi_inva_Epithelial"

cells$is_interface <- FALSE
cells$distance_to_interface_um <- NA_real_
core_metrics <- list()

for (core in unique(cells$core_id)) {
  idx <- which(cells$core_id == core)
  hi <- idx[cells$celltype[idx] == "Hi_inva_Epithelial"]
  caf <- idx[cells$celltype[idx] == "Fibroblast"]
  xy <- as.matrix(cells[idx, c("centroid_x", "centroid_y")])

  if (length(hi) > 0 && length(caf) > 0) {
    hi_xy <- as.matrix(cells[hi, c("centroid_x", "centroid_y")])
    caf_xy <- as.matrix(cells[caf, c("centroid_x", "centroid_y")])
    hi_n <- rowSums(nn2(caf_xy, hi_xy, searchtype = "radius",
                        radius = radius_um / pixel_size_um, k = nrow(caf_xy))$nn.idx > 0)
    caf_n <- rowSums(nn2(hi_xy, caf_xy, searchtype = "radius",
                         radius = radius_um / pixel_size_um, k = nrow(hi_xy))$nn.idx > 0)
    interface <- c(hi[hi_n >= min_opposite], caf[caf_n >= min_opposite])
    cells$is_interface[interface] <- TRUE

    if (length(interface) > 0) {
      nearest <- nn2(as.matrix(cells[interface, c("centroid_x", "centroid_y")]), xy, k = 1)
      cells$distance_to_interface_um[idx] <- nearest$nn.dists[, 1] * pixel_size_um
    }
  }

  xy_um <- xy * pixel_size_um
  hull <- chull(xy_um[, 1], xy_um[, 2])
  hull_area <- abs(sum(
    xy_um[hull, 1] * c(xy_um[hull, 2][-1], xy_um[hull, 2][1]) -
      c(xy_um[hull, 1][-1], xy_um[hull, 1][1]) * xy_um[hull, 2]
  )) / 2

  cd8 <- idx[cells$celltype[idx] == "CD8T"]
  core_metrics[[core]] <- data.frame(
    core_id = core,
    footprint_mm2 = hull_area / 1e6,
    interface_density = sum(cells$is_interface[idx]) / hull_area * 1e6,
    CD8_density = length(cd8) / hull_area * 1e6,
    GZMB_positive_fraction = mean(cells$GZMB_positive[cd8]),
    PD1_positive_fraction = mean(cells$PD1_positive[cd8]),
    TGFB1_positive_fraction = mean(cells$TGFB1_positive[caf], na.rm = TRUE)
  )
}

metrics <- do.call(rbind, core_metrics)
cutoff <- median(metrics$interface_density, na.rm = TRUE)
metrics$CAFI_group <- ifelse(metrics$interface_density >= cutoff, "High", "Low")
write.table(cells, "results/mIF_cell_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(metrics, "results/mIF_core_metrics.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
