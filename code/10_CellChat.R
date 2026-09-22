# CellChat inference and probability-weighted Ro/e
suppressPackageStartupMessages({library(Seurat); library(CellChat); library(data.table)})
set.seed(42)

sce <- readRDS("results/single_cell_states.rds")
data(CellChatDB.human)
db <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
saveRDS(db, "results/CellChat_secreted_database.rds")

for (group_name in c("All", "High", "Low")) {
  cells <- if (group_name == "All") colnames(sce) else colnames(sce)[sce$inva_group == group_name]
  obj <- subset(sce, cells = cells)
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")
  meta <- data.frame(celltype = obj$state, row.names = colnames(obj))

  chat <- createCellChat(expr, meta = meta, group.by = "celltype")
  chat@DB <- db
  chat <- subsetData(chat)
  chat <- identifyOverExpressedGenes(chat)
  chat <- identifyOverExpressedInteractions(chat)
  chat <- computeCommunProb(
    chat,
    type = "truncatedMean",
    trim = 0.1,
    raw.use = TRUE,
    population.size = FALSE,
    seed.use = 42
  )
  chat <- filterCommunication(chat, min.cells = 10)
  chat <- computeCommunProbPathway(chat)
  chat <- aggregateNet(chat)
  saveRDS(chat, paste0("results/CellChat_", group_name, ".rds"))
  write.table(subsetCommunication(chat), paste0("results/CellChat_", group_name, ".tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Probability-weighted Ro/e from CAFs to the three epithelial states
high <- read.delim("results/CellChat_High.tsv"); high$group <- "High"
low <- read.delim("results/CellChat_Low.tsv"); low$group <- "Low"
records <- rbind(high, low)
records <- subset(records, source == "Fibroblast" &
  target %in% c("Hi_inva_Epithelial", "Low_inva_Epithelial", "Normal_Epithelial") &
  prob > 1e-4)

observed_long <- aggregate(prob ~ interaction_name + group, records, sum)
observed <- dcast(as.data.table(observed_long), interaction_name ~ group,
                  value.var = "prob", fill = 0)
rownames(observed) <- observed$interaction_name
observed <- as.matrix(observed[, -1])
expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
roe <- observed / expected
write.table(roe, "results/CellChat_probability_ROE.tsv", sep = "\t", col.names = NA, quote = FALSE)
