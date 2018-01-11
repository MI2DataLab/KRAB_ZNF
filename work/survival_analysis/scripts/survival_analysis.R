library(RTCGA)          # Functions for querying the data 
library(RTCGA.rnaseq)   # Gene expression data
library(RTCGA.clinical) # Clinical data
library(survival)       # Survival analysis tools
library(survminer)      # Survival analysis tools
library(tidyverse)
library(knitr)
library(pheatmap)
library(plotly)

load("data/survival_data.RData")
load("data/survival_data_cat_median.RData")
load("data/survival_data_cat_optimal.RData")

# Calculate p-values of log-rank tests for each pair [gene, cohort] using optimal cutpoint
genes <- colnames(survival_data_cat_median)[5:ncol(survival_data_cat_median)]
cohorts <- unique(survival_data_cat_optimal$dataset)
pvals_optimal <- data.frame()
test_results_optimal_list <- list()
k <- 1
for (i in cohorts) {
  for(j in genes) {
    try(
      {  
        test <- survdiff(as.formula(paste("Surv(times, patient.vital_status) ~ ", j)), 
                         data = survival_data_cat_optimal[survival_data_cat_optimal$dataset == i, ])
        p_value <- NULL
        p_value <- 1 - pchisq(test$chisq, 1)
        pvals_optimal[i, j] <- p_value
        n_high <- survival_data_cat_optimal %>% filter(dataset == i) %>% select(j) %>% filter(.[[1]] == "high") %>% nrow()
        n_low <- survival_data_cat_optimal %>% filter(dataset == i) %>% select(j) %>% filter(.[[1]] == "low") %>% nrow()
        hazard_ratio <- (test$obs[1]/test$exp[1])/(test$obs[2]/test$exp[2])
        test_results_optimal_list[[k]] <- c(cohort = i, gene = j, n_low = n_low, n_high = n_high, p_value = p_value, hazard_ratio = hazard_ratio) %>% t() %>% data.frame()
        k <- k+1
      }, silent = TRUE)
  }
}
pvals_optimal$cohort <- rownames(pvals_optimal)
save(pvals_optimal, file = "pvals_optimal.Rdata")
test_results_optimal <- bind_rows(test_results_optimal_list)
save(test_results_optimal, file = "test_results_optimal.Rdata")
# Calculate p-values of log-rank tests for each pair [gene, cohort] using median cutpoint
pvals_median <- data.frame()
test_results_median_list <- list()
k <- 1
for (i in cohorts) {
  for(j in genes) {
    try(
      {  
        test <- survdiff(as.formula(paste("Surv(times, patient.vital_status) ~ ", j)), 
                         data = survival_data_cat_median[survival_data_cat_median$dataset == i, ])
        p_value <- 1 - pchisq(test$chisq, 1)
        pvals_median[i, j] <- p_value
        n_high <- survival_data_cat_median %>% filter(dataset == i) %>% select(j) %>% filter(.[[1]] == "high") %>% nrow()
        n_low <- survival_data_cat_median %>% filter(dataset == i) %>% select(j) %>% filter(.[[1]] == "low") %>% nrow()
        hazard_ratio <- (test$obs[1]/test$exp[1])/(test$obs[2]/test$exp[2])
        test_results_median_list[[k]] <- c(cohort = i, gene = j, n_low = n_low, n_high = n_high, p_value = p_value, hazard_ratio = hazard_ratio) %>% t() %>% data.frame()
        k <- k+1
      }, silent = TRUE)
  }
}
#pvals_median$cohort <- rownames(pvals_median)
save(pvals_median, file = "pvals_median.RData")
test_results_median <- bind_rows(test_results_median_list)
save(test_results_median, file = "test_results_median.RData")
