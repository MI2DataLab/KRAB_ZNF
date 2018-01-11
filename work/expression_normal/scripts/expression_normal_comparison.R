library(tidyverse)
library(xlsx)
library(RTCGA.rnaseq)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(plotly)
library(d3heatmap)

krab_corsinotti <- list.files("data/source")
krab_corsinotti_xlsx <- read.xlsx(paste("data/source/", krab_corsinotti, sep = ""), 1)
krab_corsinotti_names <- krab_corsinotti_xlsx$Gene.symbol
#krab_corsinotti_names <- krab_corsinotti_xlsx %>% filter(substr(Gene.symbol, 1, 3) == "ZNF") 
selected_genes <- droplevels(krab_corsinotti_names)

selected_genes_original <- c("ZNF695",
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

all_tcga_genes <- sapply(strsplit(colnames(BRCA.rnaseq), split = "\\|"), function(x) x[[1]]) %>% sort()
selected_colnames_logical <- sapply(strsplit(colnames(ACC.rnaseq), split = "\\|"), function(x) x[[1]]) %in% selected_genes
selected_colnames <- colnames(BRCA.rnaseq)[selected_colnames_logical]
found_genes <- sapply(strsplit(selected_colnames, split = "\\|"), function(x) x[[1]])
missing_krabs_tcga <- setdiff(selected_genes, found_genes)

selected_genes_all <- found_genes
save(selected_genes_all, file = "data/selected_genes_all.RData")
save(missing_krabs_tcga, file = "data/missing_krabs_tcga.RData")

expressions_normal <- expressionsTCGA(BRCA.rnaseq, extract.cols = selected_colnames) %>% filter(substr(bcr_patient_barcode, 14, 15) == "11")
# downloading expressions
expressions_normal <- expressionsTCGA(ACC.rnaseq,
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
                               extract.cols = selected_colnames) %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "11")  # expressions for type "11", normal
  #mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) # cleaning patient barcode
#save(expressions_normal, file = "data/expressions_normal_raw.RData")

colnames(expressions_normal)[3:ncol(expressions_normal)] <-
  sapply(strsplit(colnames(expressions_normal)[3:ncol(expressions_normal)], split = "\\|"), function(x)
    x[[1]])

expressions_normal$dataset <- gsub(pattern = ".rnaseq", replacement = "", x = expressions_normal$dataset)
save(expressions_normal, file = "data/expressions_normal.RData")


load("data/expressions_normal.RData")
expressions_long <- melt(expressions_normal)
colnames(expressions_long) <- c("patient_code", "cohort", "gene", "expression")
expressions_long$cohort <- expressions_long$cohort %>% tolower() %>%
  gsub(pattern = "acc", replacement = "Adenoid cystic carcinoma (ACC)") %>%
  gsub(pattern = "blca", replacement = "Bladder (BLCA)") %>%
  gsub(pattern = "brca", replacement = "Breast (BRCA)") %>%
  gsub(pattern = "cesc", replacement = "Cervical squamous cell (CESC)") %>%
  gsub(pattern = "chol", replacement = "Cholangiocarcinoma (CHOL)") %>%
  gsub(pattern = "coadread", replacement = "Colon (COAD, COADREAD)") %>%
  gsub(pattern = "coad", replacement = "Colon (COAD, COADREAD)") %>%
  gsub(pattern = "dlbc", replacement = "B-Cells (DLBC)") %>%
  gsub(pattern = "esca", replacement = "Esophagus (ESCA)") %>%
  gsub(pattern = "gbm", replacement = "Brain (GBM, LGG)") %>%
  gsub(pattern = "hnsc", replacement = "Head and neck (HNSC)") %>%
  gsub(pattern = "kich", replacement = "Kidney (KICH, KIPAN, KIRC, KIRP)") %>%
  gsub(pattern = "kipan", replacement = "Kidney (KICH, KIPAN, KIRC, KIRP)") %>%
  gsub(pattern = "kirc", replacement = "Kidney (KICH, KIPAN, KIRC, KIRP)") %>%
  gsub(pattern = "kirp", replacement = "Kidney (KICH, KIPAN, KIRC, KIRP)") %>%
  gsub(pattern = "lgg", replacement = "Brain (GBM, LGG)") %>%
  gsub(pattern = "lihc", replacement = "Liver (LIHC)") %>%
  gsub(pattern = "luad", replacement = "Lung (LUAD, LUSC)") %>%
  gsub(pattern = "lusc", replacement = "Lung (LUAD, LUSC)") %>%
  gsub(pattern = "ov", replacement = "Ovary (OV)") %>%
  gsub(pattern = "paad", replacement = "Pancreas (PAAD)") %>%
  gsub(pattern = "pcpg", replacement = "Pheochromocytoma and Paraganglioma (PCPG)") %>%
  gsub(pattern = "prad", replacement = "Prostate (PRAD)") %>%
  gsub(pattern = "read", replacement = "Rectum (READ)") %>%
  gsub(pattern = "sarc", replacement = "Pheochromocytoma and Paraganglioma (PCPG)") %>%
  gsub(pattern = "sarc", replacement = "Bone (SARC)") %>%
  gsub(pattern = "skcm", replacement = "Skin (SKCM)") %>%
  gsub(pattern = "stad", replacement = "Stomach (STAD, STES)") %>%
  gsub(pattern = "stes", replacement = "Stomach (STAD, STES)") %>%
  gsub(pattern = "tgct", replacement = "Testicles (TGCT)") %>%
  gsub(pattern = "thca", replacement = "Thyroid (THCA)") %>%
  gsub(pattern = "thym", replacement = "Thymus (THYM)") %>%
  gsub(pattern = "ucec", replacement = "Uterus (UCEC, UCS)") %>%
  gsub(pattern = "ucs", replacement = "Uterus (UCEC, UCS)") %>%
  gsub(pattern = "uvm", replacement = "Uvea (UVM)")
#expressions_long$cohort <- substr(expressions_long$cohort, 1, regexpr("\\.", expressions_long$cohort) - 1) %>% as.factor()
cohorts_summary <- expressions_long %>% group_by(cohort) %>% summarise(number_of_obs = n()/length(unique(expressions_long$gene)))
expressions_means <- expressions_long %>% group_by(cohort, gene) %>% summarise(mean = mean(expression))
expressions_means <- spread(expressions_means, gene, mean) %>% data.frame()
rownames(expressions_means) <- expressions_means$cohort
expressions_means <- expressions_means[, -1]
save(expressions_means, file = "data/expression_means.RData")

#heatmap
plot_ly(z = as.matrix(expressions_means[, selected_genes_original]), 
        x = colnames(expressions_means[, selected_genes_original]), 
        y = rownames(expressions_means[, selected_genes_original]), type = "heatmap", colors = "Greys")

pheatmap(t(log(expressions_means + 0.001)), cellwidth = 10, cellheight = 10)

expressions_long %>% filter(cohort %in% cohort) %>% ggplot(aes(x = gene, y = expression)) + geom_boxplot()

t(expressions_means[, selected_genes_original]) %>% as.matrix() + 1 %>% log()

RColorBrewer::brewer.pal.info


split(methylation[[1]], rep(1:ceiling(nrow(methylation[[1]])/10), each = 10, length.out = nrow(methylation[[1]])))

