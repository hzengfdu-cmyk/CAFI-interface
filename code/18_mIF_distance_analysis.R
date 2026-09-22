# Within-core mIF distance binning and Spearman correlations
suppressPackageStartupMessages(dplyr)

cells <- read.delim("results/mIF_cell_table.tsv")
results <- list(); correlations <- list()

targets <- list(
  GZMB = c("CD8T", "GZMB_asinh5"),
  PD1 = c("CD8T", "PD1_asinh5"),
  TGFB1 = c("Fibroblast", "TGFB1_asinh5")
)

for (marker in names(targets)) {
  target_celltype <- targets[[marker]][1]
  intensity <- targets[[marker]][2]
  data <- cells %>%
    filter(.data$celltype == target_celltype, !is_interface, is.finite(distance_to_interface_um)) %>%
    group_by(core_id) %>%
    mutate(distance_bin = ntile(distance_to_interface_um, 5)) %>%
    group_by(core_id, distance_bin) %>%
    summarise(
      n_cells = n(),
      mean_distance_um = mean(distance_to_interface_um),
      mean_intensity = mean(.data[[intensity]]),
      .groups = "drop"
    ) %>%
    filter(n_cells >= 5)

  test <- cor.test(data$mean_distance_um, data$mean_intensity,
                   method = "spearman", exact = FALSE)
  data$marker <- marker
  results[[marker]] <- data
  correlations[[marker]] <- data.frame(
    marker,
    n_bins = nrow(data),
    n_cores = n_distinct(data$core_id),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}

write.table(do.call(rbind, results), "results/mIF_distance_bins.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, correlations), "results/mIF_distance_correlations.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
