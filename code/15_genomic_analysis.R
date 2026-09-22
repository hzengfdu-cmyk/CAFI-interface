# TCGA and IMvigor210 genomic event processing
suppressPackageStartupMessages(library(dplyr))

mutation <- read.delim("data/somatic_mutations.tsv")
cna <- read.delim("data/copy_number.tsv")
sample_info <- read.delim("data/genomic_sample_info.tsv")

functional_classes <- c(
  "Missense_Mutation", "Splice_Site", "Splice_Region", "Frame_Shift_Del",
  "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Nonsense_Mutation",
  "Nonstop_Mutation", "Translation_Start_Site"
)

mutation <- mutation %>% filter(variant_class %in% functional_classes)
cna <- cna %>% filter(cna %in% c(-2, 2))

# Define CAFI groups using the cohort-specific median
sample_info <- sample_info %>%
  group_by(cohort) %>%
  mutate(CAFI_group = ifelse(CAFI_score >= median(CAFI_score), "High", "Low")) %>%
  ungroup()

mutation_frequency <- mutation %>%
  inner_join(sample_info, by = c("cohort", "sample_id")) %>%
  distinct(cohort, sample_id, gene, CAFI_group) %>%
  count(cohort, gene, CAFI_group, name = "altered")

cna_frequency <- cna %>%
  inner_join(sample_info, by = c("cohort", "sample_id")) %>%
  mutate(event = ifelse(cna == 2, "Amplification", "Deep_deletion")) %>%
  distinct(cohort, sample_id, gene, event, CAFI_group) %>%
  count(cohort, gene, event, CAFI_group, name = "altered")

write.table(mutation_frequency, "results/mutation_frequency.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cna_frequency, "results/CNA_frequency.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
