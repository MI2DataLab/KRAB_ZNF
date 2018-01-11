# Methylation-expression correlation analysis
Important note: all files for this part of analysis are available at:
https://drive.google.com/drive/folders/0B9C7gUUOHbNpaXl0R2tCR013dEU
due to Github file size limitations

## Data
Methylation and expression data downloaded from TCGA via RTCGA package.

Following methylation datasets were available in RTCGA package:
BRCA.methylation,
COAD.methylation,
COADREAD.methylation,
GBM.methylation,
GBMLGG.methylation,
KIPAN.methylation,
KIRC.methylation,
KIRP.methylation,
LAML.methylation,
LUAD.methylation,
LUSC.methylation,
OV.methylation1, 
OV.methylation2
READ.methylation,
STAD.methylation,
STES.methylation,
UCEC.methylation

The downloaded copy of the data can be found in a list object in data/mathylation.Rdata. 


While for expressions there were:
ACC.rnaseq,
BLCA.rnaseq,
BRCA.rnaseq,
CESC.rnaseq,
CHOL.rnaseq,
COAD.rnaseq,
COADREAD.rnaseq,
DLBC.rnaseq,
ESCA.rnaseq,
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
UVM.rnaseq

The downloaded copy of the expression data can be found in data/expression.Rdata

Additonal data about cpg island we take from the file illumina_humanmethylation_27

##Analysis
Expression and methylation data is later merged by patient barcode.
Then, for each [cohort, gene, cpg island] we calculate (if possible, because some in some places in the methylation data there are many NA values causing not sufficient amount of observations for t-test) a t-test comparing two
groups of observations (lowest and highest 10% with respect to gene expression level).

Results of t-test are saved in .Rdata files (list objects), each cohort in a separate file.

Script correlation_data/data_to_csv.R joins these results of test with information about genes like "Distance_to_TSS" etc. taken form illumina_humanmethylation_27_data.rda and transforms them into csv files, stored in folder csv/.

Note that each cohort has 2 csv files:
1. All results of tests stored in for example "BRCA all.csv"
2. Filtered significant results (p value < 0.05 and difference in methylation means > 0.25) stored for example in " BRCA significant.csv"

  