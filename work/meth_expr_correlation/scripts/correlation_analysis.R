library(RTCGA.clinical)
library(RTCGA.rnaseq)
library(RTCGA.methylation)
library(RTCGA.rnaseq)
library(dplyr)
library(ggplot2)
library(feather)


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

# loading expressions previously downloaded from TCGA via RTCGA package
load("data/expressions.Rdata")

# cleaning cohort name and gene name
expressions$dataset <- substr(expressions$dataset, 1, nchar(expressions$dataset) - 7)
colnames(expressions)[3:16] <- substr(colnames(expressions)[3:16], 1, 6)

# splitting expressions 
expressions_split <- list()
for(i in unique(expressions$dataset))
  expressions_split[[i]] <- expressions %>% filter(dataset == paste(i)) %>% as.data.frame()

# calculating ntiles
k <- 10

for(i in unique(expressions$dataset))
  for(j in selected_genes)
    expressions_split[[i]][, j] <- ntile(expressions_split[[i]][, j], k)

# downloading methylation
methylation <- list(
  BRCA = BRCA.methylation,
  COAD = COAD.methylation,
  COADREAD = COADREAD.methylation,
  GBM = GBM.methylation,
  GBMLGG = GBMLGG.methylation,
  KIPAN = KIPAN.methylation,
  KIRC = KIRC.methylation,
  KIRP = KIRP.methylation,
  LAML = LAML.methylation,
  LUAD = LUAD.methylation,
  LUSC = LUSC.methylation,
  OV = rbind(OV.methylation1, OV.methylation2),
  READ = READ.methylation,
  STAD = STAD.methylation,
  STES = STES.methylation,
  UCEC = UCEC.methylation
)


# cleaning methylation data and subseting only cancer cells
for(i in names(methylation)) {
  methylation[[i]] <- methylation[[i]] %>% 
    filter(substr(bcr_patient_barcode, 14, 15) == "01") %>%
    mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))
}

# saving data
# for(i in names(methylation)) {
#   if(nrow(methylation[[i]]) > 10) {
#   selected_methylation <- methylation[[i]]#[, !is.na(methylation[[i]][1, ])]
#   splitted_methylation <- split(selected_methylation, rep(1:ceiling(nrow(selected_methylation)/10), each = 10, length.out = nrow(selected_methylation)))
#   #write_feather(selected_methylation, path = paste("data/methylation/", i, "methylation.feather", sep = ""))
#   for(j in 1:length(splitted_methylation)){
#     methylation_part <- splitted_methylation[[j]]
#     save(methylation_part, file = paste("data/methylation/", i, "methylation_pt_", j, ".Rdata", sep = ""))
#   }
#   }
# }
files <- list.files(path = "data/methylation", pattern = "BRCA", full.names = TRUE)
selected_methylation <- data.frame()
for(file in files) {
  load(file)
  selected_methylation <- rbind(selected_methylation, methylation_part)
}

# Joining methylation with expression
methylation_expression <- list()

for(i in names(methylation)) {
  try({
    methylation_expression[[i]] <- expressions_split[[i]] %>% left_join(methylation[[i]])
  })
}
save(methylation_expression, file = "data/meth_expr.Rdata")

# load if necessary previously saved data
# load("data/meth_expr.Rdata")
# load("data/methylation.Rdata")

# extracting names of cpg islands
cpg_islands <- colnames(methylation_expression[[1]])[17:length(colnames(methylation_expression[[1]]))]
n_of_combinations <- length(cpg_islands)*length(selected_genes)*length(methylation)


# For each [cohort, gene, cpg island] we calculate (if possible) a t-test comparing two
# groups of observations (lowest and highest 10% with respect to gene expression level)
# Each cohort saved in a separate file
for(m in 1:length(names(methylation_expression))) {
  row_n <- 1
  pvals_list <- list()
  
  
  for(i in names(methylation_expression)[m]) {
    for(j in selected_genes) {
      for(k in cpg_islands) {
        low_expr <- methylation_expression[[i]][methylation_expression[[i]][, j] == 1 , k]
        high_expr <- methylation_expression[[i]][methylation_expression[[i]][, j] == 10 , k]
        pvalue <- NA
        try({
          test <- t.test(low_expr, high_expr)
          pvalue <- test$p.value
        }, silent = TRUE)
        pvals_list[[row_n]] <- data.frame("cohort" = i, 
                                          "gene" = j, 
                                          "cpg_island" = k, 
                                          "pvalue" = pvalue, 
                                          "low_expr_meth_mean" = mean(low_expr, na.rm = TRUE), 
                                          "high_expr_meth_mean" = mean(high_expr, na.rm = TRUE))
        row_n <- row_n + 1
      }
    }
  }
  save(pvals_list, file = paste("correlation_data/pvals_list_", m, ".Rdata"))
}


