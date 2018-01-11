# converting data to csv
library(dplyr)

load("illumina_humanmethylation_27_data.rda")
files <- list.files(pattern = "pvals_list")

# We load each cohort file separately, adjust p-values via False Discovery Rate correction,
# join additional columns and ve in csv both whole file and filtered out by pvalue_adjusted < 0.05 
# and difference in means > 0.25 

for(i in files) {
load(i)
pvals <- bind_rows(pvals_list)
pvals$pvalue <- as.numeric(pvals$pvalue)
pvals$pvalue_adjusted <- p.adjust(pvals$pvalue, method = "fdr")

pvals_extended <- pvals %>% arrange(pvalue_adjusted) %>%
  mutate(log_odds_ratio =  log2((1-high_expr_meth_mean)*low_expr_meth_mean/((1-low_expr_meth_mean)*high_expr_meth_mean))) %>%
  merge(illumina_humanmethylation_27_data[, c("Name", 
                                              "Gene_ID", 
                                              "Symbol", 
                                              "TSS_Coordinate", 
                                              "Distance_to_TSS",
                                              "CPG_ISLAND",
                                              "CHR",	
                                              "MAPINFO", 
                                              "UCSC_REFGENE_NAME",	
                                              "UCSC_REFGENE_GROUP",
                                              "RELATION_TO_UCSC_CPG_ISLAND",
                                              "Gene_Strand",
                                              "Accession",
                                              "Product")], 
        by.x = "cpg_island", by.y = "Name") 

pvals_extended_signif <- pvals_extended[pvals_extended$pvalue_adjusted < 0.05 & !is.na(pvals_extended$pvalue_adjusted), ] %>%
  filter(abs(high_expr_meth_mean - low_expr_meth_mean) > 0.25) %>%
  arrange(pvalue_adjusted)

write.csv(pvals_extended, file = paste("csv/", pvals_extended$cohort[1], "all.csv"))
write.csv(pvals_extended_signif, file = paste("csv/", pvals_extended_signif$cohort[1], "significant.csv"))
remove(pvals_list)
}