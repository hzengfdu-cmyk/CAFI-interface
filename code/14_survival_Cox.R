# Kaplan–Meier and final Cox models
suppressPackageStartupMessages({library(survival); library(survminer)})

clinical <- read.delim("data/clinical_survival.tsv")

# Kaplan–Meier groups use the cohort-specific median CAFI score
clinical$KM_group <- ave(clinical$CAFI_score, clinical$cohort,
                         FUN = function(x) ifelse(x >= median(x, na.rm = TRUE), "High", "Low"))

for (cohort_name in unique(clinical$cohort)) {
  d <- subset(clinical, cohort == cohort_name)
  endpoints <- if (cohort_name == "NMIBC") "RFS" else if (cohort_name == "XQ_ICB") c("OS", "PFS") else "OS"
  for (endpoint in endpoints) {
    d$time <- d[[paste0(endpoint, "_months")]]
    d$event <- d[[paste0(endpoint, "_event")]]
    fit <- survfit(Surv(time, event) ~ KM_group, data = d)
    saveRDS(fit, paste0("results/KM_", cohort_name, "_", endpoint, ".rds"))
  }
}

# Cox models use surv_cutpoint for continuous variables
make_group <- function(d, variable) {
  tmp <- data.frame(time = d$time, event = d$event, value = d[[variable]])
  cp <- surv_cutpoint(tmp, time = "time", event = "event", variables = "value")
  cutoff <- cp$cutpoint["value", "cutpoint"]
  factor(ifelse(d[[variable]] > cutoff, "High", "Low"), levels = c("Low", "High"))
}

# Complete-case MIBC model
mibc <- subset(clinical, cohort == "MIBC")
mibc$time <- mibc$OS_months; mibc$event <- mibc$OS_event
mibc$CAFI_group <- make_group(mibc, "CAFI_score")
mibc$Age_group <- make_group(mibc, "Age_num")
mibc_fit <- coxph(
  Surv(time, event) ~ CAFI_group + Age_group + Sex_binary + pT_binary + pN_binary + pM_binary,
  data = mibc
)

# Complete-case NMIBC model
nmibc <- subset(clinical, cohort == "NMIBC")
nmibc$time <- nmibc$RFS_months; nmibc$event <- nmibc$RFS_event
nmibc$CAFI_group <- make_group(nmibc, "CAFI_score")
nmibc$Age_group <- make_group(nmibc, "Age_num")
nmibc_fit <- coxph(
  Surv(time, event) ~ CAFI_group + Age_group + Sex_binary + pT_binary,
  data = nmibc
)

# Complete six-variable IMvigor210 model
imv <- subset(clinical, cohort == "IMvigor210")
imv$time <- imv$OS_months; imv$event <- imv$OS_event
imv$CAFI_group <- make_group(imv, "CAFI_score")
imv$FMOne_group <- make_group(imv, "FMOne_mutation_burden_per_MB")
imv_fit <- coxph(
  Surv(time, event) ~ CAFI_group + FMOne_group + 
    IC_level_binary + TC_level_binary + Immune_inflamed_vs_other,
  data = imv
)

saveRDS(list(MIBC = mibc_fit, NMIBC = nmibc_fit, IMvigor210 = imv_fit),
        "results/final_Cox_models.rds")
