library(RTCGA)          # Functions for querying the data 
library(RTCGA.rnaseq)   # Gene expression data
library(RTCGA.clinical) # Clinical data
library(survival)       # Survival analysis tools
library(survminer)      # Survival analysis tools
library(tidyverse)

# Genenes selected for investigation
selected_genes <- c("ZNF695",
                    "ZNF643",
                    "ZNF789",
                    "ZNF320",
                    "ZNF273",
                    "ZNF707",
                    "ZNF205",
                    "ZNF468",
                    "ZNF714",
                    "ZNF485",
                    "ZNF525",
                    "ZNF267",
                    "ZNF282",
                    "ZNF114")

# extractiong column names for selected genes
selected_genes_colnames <- sapply(selected_genes, grep, colnames(BRCA.rnaseq), value = TRUE)

# downloading expressions 
expressions <- expressionsTCGA(ACC.rnaseq,
                               BLCA.rnaseq,
                               BRCA.rnaseq,
                               CESC.rnaseq,
                               CHOL.rnaseq,
                               COAD.rnaseq,
                               COADREAD.rnaseq,
                               DLBC.rnaseq,
                               ESCA.rnaseq,
                               #FPPP.rnaseq,
                               GBM.rnaseq,
                               GBMLGG.rnaseq,
                               HNSC.rnaseq,
                               KICH.rnaseq,
                               KIPAN.rnaseq,
                               KIRC.rnaseq,
                               KIRP.rnaseq,
                               LAML.rnaseq,
                               LGG.rnaseq,
                               LIHC.rnaseq,
                               LUAD.rnaseq,
                               LUSC.rnaseq,
                               #MESO.rnaseq,
                               OV.rnaseq,
                               PAAD.rnaseq,
                               PCPG.rnaseq,
                               PRAD.rnaseq,
                               READ.rnaseq,
                               SARC.rnaseq,
                               SKCM.rnaseq,
                               STAD.rnaseq,
                               STES.rnaseq,
                               TGCT.rnaseq,
                               THCA.rnaseq,
                               THYM.rnaseq,
                               UCEC.rnaseq,
                               UCS.rnaseq,
                               UVM.rnaseq, 
                               extract.cols = selected_genes_colnames) %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # expressions for type "01"
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) # cleaning patient barcode

# downloading survival times
survival_times <- survivalTCGA(ACC.clinical,
                               BLCA.clinical,
                               BRCA.clinical,
                               CESC.clinical,
                               CHOL.clinical,
                               COAD.clinical,
                               COADREAD.clinical,
                               DLBC.clinical,
                               ESCA.clinical,
                               FPPP.clinical,
                               GBM.clinical,
                               GBMLGG.clinical,
                               HNSC.clinical,
                               KICH.clinical,
                               KIPAN.clinical,
                               KIRC.clinical,
                               KIRP.clinical,
                               LAML.clinical,
                               LGG.clinical,
                               LIHC.clinical,
                               LUAD.clinical,
                               LUSC.clinical,
                               MESO.clinical,
                               OV.clinical,
                               PAAD.clinical,
                               PCPG.clinical,
                               PRAD.clinical,
                               READ.clinical,
                               SARC.clinical,
                               SKCM.clinical,
                               STAD.clinical,
                               STES.clinical,
                               TGCT.clinical,
                               THCA.clinical,
                               THYM.clinical,
                               UCEC.clinical,
                               UCS.clinical,
                               UVM.clinical) %>%
                               distinct()  #There are some duplicates

# Joining survival times with gene expressions
survival_data <- survival_times %>% 
                 left_join(expressions, by = "bcr_patient_barcode") %>%
                 filter(times >= 0) %>% # omit strange cases of negative times
                 filter(!is.na(dataset))

# Check for join multiplicates
survival_data %>% group_by(bcr_patient_barcode, dataset) %>% 
                  summarise(count = n()) %>% 
                  filter(count > 1)

# We remove SKCM.rnaseq, because there is only one observation.
survival_data <- survival_data %>% filter(dataset != "SKCM.rnaseq")

# And shorten gene names
colnames(survival_data)[5:18] <- substr(colnames(survival_data)[5:18], 1, 6)

survival_data$dataset <- gsub(".rnaseq", "", survival_data$dataset)

# Data overview
survival_data %>% group_by(dataset) %>% summarise(n_of_observations = n())



# Categorisation of gene expressions
# In order to use log-rank test for comparing two survival curves we need to split 
# gene expressions variables into categorical variables (high/low expression). 

# Cathegorization of gene expression using the maximally selected rank statistics, 
# for each cohort separately.
categorized_list <- list()
surv_cutpoint_list <- list()
cutpoints_list <- list()

for(i in unique(survival_data$dataset)) {
  surv_cutpoint_list[[i]] <-  surv_cutpoint(
    survival_data[survival_data$dataset == i,],
    time = "times",
    event = "patient.vital_status",
    variables = c(colnames(survival_data)[5:18], "dataset")
  )
  cutpoints_list[[i]] <- surv_cutpoint_list[[i]]$cutpoint[, 1]
  categorized_list[[i]] <- surv_categorize(surv_cutpoint_list[[i]])
}

cutpoints_optimal_matrix <- bind_cols(cutpoints_list)
rownames(cutpoints_matrix) <- colnames(survival_data)[5:18]
cutpoints__optimal_matrix <- t(cutpoints_optimal_matrix)
survival_data_cat_optimal <- bind_rows(categorized_list) # Merging lists
survival_data_cat_optimal$dataset <- gsub(".rnaseq", "", survival_data_cat_optimal$dataset)

# We can also use median as a cutpoint for gene expressions
expression_medians <- survival_data[, 4:18] %>% group_by(dataset) %>% summarise_each(funs(median))
survival_data_cat_median <- survival_data

survival_data_split <- split(survival_data, survival_data$dataset)


median_cut <- function(x) {
  for(i in colnames(x)[5:18]) {
    x[, i] <- cut(x[, i], 
                  c(min(x[, i]), median(x[, i]), max(x[, i])), 
                  include.lowest = TRUE,
                  labels = c("low", "high")
                  )
  }
  return(x)
}

survival_data_cat_median <- bind_rows(lapply(survival_data_split, median_cut))
survival_data_cat_median <- survival_data_cat_median[, -2]
survival_data_cat_median$dataset <- gsub(".rnaseq", "", survival_data_cat_median$dataset)

save(survival_data_cat_optimal, file = "survival_data_cat_optimal.Rdata")
save(survival_data_cat_median, file = "survival_data_cat_median.Rdata")
save(survival_data, file = "survival_data.Rdata")
save(cutpoints_matrix, file = "cutpoints_matrix.Rdata")

