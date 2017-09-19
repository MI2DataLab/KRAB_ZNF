library(RTCGA)          # Functions for querying the data 
library(RTCGA.rnaseq)   # Gene expression data
library(RTCGA.clinical) # Clinical data
library(survival)       # Survival analysis tools
library(survminer)      # Survival analysis tools
library(tidyverse)
library(knitr)
library(pheatmap)
library(plotly)

load("survival_data.Rdata")
load("survival_data_cat_median.Rdata")
load("survival_data_cat_optimal.Rdata")

# Calculate p-values of log-rank tests for each pair [gene, cohort] using optimal cutpoint
genes <- colnames(survival_data_cat_optimal)[3:16]
cohorts <- unique(survival_data_cat_optimal$dataset)
pvals_optimal <- data.frame()

for (i in cohorts) {
  for(j in genes) {
    try(
      {  
        test <- survdiff(as.formula(paste("Surv(times, patient.vital_status) ~ ", j)), 
                         data = survival_data_cat_optimal[survival_data_cat_optimal$dataset == i, ])
        p_value <- 1 - pchisq(test$chisq, 1)
        pvals_optimal[i, j] <- p_value
      }, silent = TRUE)
  }
}
pvals_optimal$cohort <- rownames(pvals_optimal)
save(pvals_optimal, file = "pvals_optimal.Rdata")

# Calculate p-values of log-rank tests for each pair [gene, cohort] using median cutpoint
pvals_median <- data.frame()

for (i in cohorts) {
  for(j in genes) {
    try(
      {  
        test <- survdiff(as.formula(paste("Surv(times, patient.vital_status) ~ ", j)), 
                         data = survival_data_cat_median[survival_data_cat_median$dataset == i, ])
        p_value <- 1 - pchisq(test$chisq, 1)
        pvals_median[i, j] <- p_value
      }, silent = TRUE)
  }
}
pvals_median$cohort <- rownames(pvals_median)
save(pvals_median, file = "pvals_median.Rdata")
