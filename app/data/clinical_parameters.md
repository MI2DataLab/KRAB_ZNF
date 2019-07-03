---
title: "Clinical Parameters Description"
output: 
  html_document:
    toc: true
---

# Description for parametrers used in [KRAB ZNF Explorer](http://mi2.mini.pw.edu.pl:8080/KRAB_ZNF/).

## UCEC

### Clinical stage 

Combination of TNM results for each patient, based on patient history, physical
examination, and any imaging done before initiation of treatment.

Stage 0 - cancer in situ,

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

### Histology_type 

Histologic subtype of Uterine Corpus Endometrial Carcinoma submitted for TCGA

msi_status_7-marker_call - MSI (microsatellite instability) in cancer genome:

MSI-H – MSI-high,

MSI-L – MSI-low,

MSS – MS-stable.

### Tumor_grade 

Classification of the microscopic cell appearance abnormality and deviations in their rate of growth with the goal of predicting developments at tissue level.

GX: Grade cannot be assessed (undetermined grade).

G1: Well differentiated (low grade).

G2: Moderately differentiated (intermediate grade).

G3: Poorly differentiated (high grade).

G4: Undifferentiated (high grade).

###  X_PANCAN_DNAMethyl_UCEC

Unsupervised clustering of DNA methylation data generated from Illumina InfiniumDNAmethylation
arrays revealed four unique subtypes (based on PMID: 23636398).

Cluster 1 - heavily methylated subtype reminiscent of the CpG island methylator phenotype(CIMP)
described in colon cancers and glioblastomas, associated with the MSI subtype and attributable to
promoter hypermethylation of MLH1.

Cluster 3 – a serous-like cluster with minimal DNA methylation changes, composed primarily of
serous tumours and some endometrioid tumours .

## THCA

### histological_type

Histologic subtype of Thyroid Cancer submitted for TCGA.

### meth_Cluster 

Unsupervised clustering of DNA methylation data generated from Illumina InfiniumDNAmethylation arrays revealed four unique subtypes (based onPMID: 25417114 ):

CpG island methylated – hypermethylation of a large number of CpG sites in islands and shores,

Follicular – few methylation changes compared to normal thyroid.

### mRNA_Cluster_number – Different THCA subtypes based on mRNA expression profiling (based on
PMID: 25417114).

### Pathologic_M – Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### Pathologic_T – characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### Pathologic_N –characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### Pathologic_stage 

Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.Stage IV - advanced or metastatic cancer. 

## PRAD

### histological_type

Histologic subtypes of Prostate Adenocarcinoma submitted for TCGA.

### clinical_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### methylation_cluster 

Unsupervised hierarchical clustering of the most variably hypermethylated

CpGs identified four epigenetically distinct groups of prostate cancers (based on PMID: 26544944).

### mRNA_cluster 

Different PRAD subtypes based on mRNA expression profiling (based on PMID:
26544944).

### mRNA_subtype 

Molecular subtypes of prostate cancer based on known and novel genetic drivers
of the disease; four are characterized by gene fusions ,three are defined by gene mutations (based on
PMID: 26544944).

### Pathologic_N 
Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 

Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

## LUSC
### expression_subtypes_LUSC 

Whole-transcriptome expression profiles generated by RNA sequencing
and by microarrays (based on PMID: 22960745).

### histological_type

Histologic subtypes of tumors submitted for TCGA.PANCAN_Cluster_Cluster_PANCAN– Unsupervised clustering of DNA methylation data generated
from Illumina InfiniumDNAmethylation arrays.

### pathologic_N 
Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 
Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### pathologic_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage 
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

### Smoking_history

1 –Lifelong Non-smoker

2 –Current smoker

3 –Current reformed smoker for > 15 years

4 –Current reformed smoker for < or = 15 years

5 –Current Reformed Smoker, Duration Not Specified

## LUAD

### expression_subtypes_LUAD 
Lung adenocarcinoma subtypes based on mRNA expression (based on
PMID: 25079552).

### histological_type
Histologic subtypes of tumors submitted for TCGA.

### pathologic_N 
Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 
Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### pathologic_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage 
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

### Smoking_history

1 –Lifelong Non-smoker

2 –Current smoker

3 –Current reformed smoker for > 15 years

4 –Current reformed smoker for < or = 15 years

5 –Current Reformed Smoker, Duration Not Specified

## LIHC

### HBV-consensus 
Clinical or molecular evidence of HBV infection.

### HCV_consensus 
Serological and/or molecular markers of HCV infection.

### Histological_type 
Histologic subtypes of tumors submitted for TCGA.

### Hypermethylation.Cluster.Laird.group 
Unsupervised clustering of HCC using CpG sites showing
cancer-specific DNA hypermethylation (based on PMID: 28622513).

### Hypomethylation.Cluster.Laird.group 

Unsupervised clustering of HCC using CpG sites showing
cancer-specific DNA hypomethylation (based on PMID: 28622513).

### mRNA.clusters.5.group.NMF.Hoadley.group
Liver hepatocellular carcinoma subtypes based on mRNA expression.
### pathologic_N 
Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 
Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.
### pathologic_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage 
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

## KIRP

### DNA_methylation_subtype 

Molecular subtyping by means of a DNA methylation platform revealed
three subtypes of papillary renal-cell carcinoma (PRCC), one of which showed widespread DNA
hypermethylation patterns characteristic of CIMP-associated tumors (the other subtypes are
identified as cluster 1 and cluster 2) (based on PMID: 26536169).

### histological_subtype 
Histologic subtypes of tumors submitted for TCGA.

### mRNA_subtype
Kidney Renal Papillary Cell Carcinoma subtypes based on mRNA expression (based
on PMID: 26536169).

### pathologic_N 
Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 
Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### pathologic_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage 
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

### Smoking_history

1 –Lifelong Non-smoker

2 –Current smoker

3 –Current reformed smoker for > 15 years

4 –Current reformed smoker for < or = 15 years

5 –Current Reformed Smoker, Duration Not Specified

## KIRC
### mRNA_cluster 
Unsupervised clustering methods identified four stable subsets in mRNA expression
data sets (based on PMID: 23792563).
### pathologic_N 
Characterization of the regional lymph nodes.NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 
Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.
### pathologic_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage 
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

### Smoking_history

1 –Lifelong Non-smoker

2 –Current smoker

3 –Current reformed smoker for > 15 years

4 –Current reformed smoker for < or = 15 years

5 –Current Reformed Smoker, Duration Not Specified

### tumor_grade 
Classification of the microscopic cell appearance abnormality and deviations in their
rate of growth with the goal of predicting developments at tissue level.GX: Grade cannot be assessed (undetermined grade)

G1: Well differentiated (low grade)

G2: Moderately differentiated (intermediate grade)

G3: Poorly differentiated (high grade)

G4: Undifferentiated (high grade)

## KICH

### histological_type_eosinophilic.1_classic.0 
Histologic subtypes of chromophobe renal cell
carcinoma samples submitted for TCGA. 1 – eosinophilic, 0 – classic.
### pathologic_N 
Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 

Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### pathologic_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage 
combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.Stage IV - advanced or metastatic cancer.

### Smoking_history

1 –Lifelong Non-smoker

2 –Current smoker

3 –Current reformed smoker for > 15 years

4 –Current reformed smoker for < or = 15 years

5 –Current Reformed Smoker, Duration Not Specified

## HNSC

### clinical_N 
Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### clinical_T 
Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### clinical_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### clinical_stage 
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.Stage IV - advanced or metastatic cancer.

### histlogical_type 
Histologic subtypes of Head and Neck Squamous Cell Carcinoma samples
submitted for TCGA.

### Hpv.status.ish 
HPV testing based on HPV16 in situ hybridization (ISH).

### Hpv.status.p16 
HPV testing based on p16 immunohistochemistry.
### histological_type
Histologic subtypes of Head and Neck Squamous Cell Carcinoma samples
submitted for TCGA.
### Methylation_subype 
Unsupervised clustering of DNA methylation data generated from Illumina
InfiniumDNAmethylation arrays revealed four unique subtypes (based on PMID: 25631445).
### RNA_subtype 
Head and Neck subtypes based on mRNA expression (based on PMID: 25631445).
### pathologic_N 
characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 
Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### pathologic_M

Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage 
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

### Smoking_history

1 –Lifelong Non-smoker

2 –Current smoker

3 –Current reformed smoker for > 15 years

4 –Current reformed smoker for < or = 15 years

5 –Current Reformed Smoker, Duration Not Specified

## ESCA

### columnar_metaplasia_present 
Replacement of the normal stratified squamous epithelium lining of
the esophagus by simple columnar epithelium with goblet cells.

### columnar_mucosa_dysplasia 
Apre-malignant lesion the esophagus associated with Barrett's
esophagus; considered the precursor of esophageal adenocarcinoma.

### histological_type 
Histologic subtypes of Esophageal Carcinoma samples submitted for TCGA.

### histologic_grade 
Classification of the microscopic cell appearance abnormality and deviations in
their rate of growth with the goal of predicting developments at tissue level.

GX: Grade cannot be assessed (undetermined grade)

G1: Well differentiated (low grade)

G2: Moderately differentiated (intermediate grade)

G3: Poorly differentiated (high grade)

G4: Undifferentiated (high grade)

### H.PYLORI-Infection 
status of patient’s 
Helicobacter pylori infection.

### pathologic_N 

Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 

Characterizationof the size and extent of the main tumor.TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### pathologic_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

### Smoking_history

1 –Lifelong Non-smoker

2 –-Current smoker

3 –Current reformed smoker for > 15 years

4 –Current reformed smoker for < or = 15 years

5 –Current Reformed Smoker, Duration Not Specified

## BRCA

### ER.status 

The presence of estrogen receptors on the surface of the cancer cell

### HER2.Final.Status 
The presence of a growth-promoting protein (HER2) on the surface of the cancer
cell

### histological_type 
Histologic subtypes of Breast Cancer samples submitted for TCGA.

### Integrated.Clusters.with.PAM50 
Breast Cancer subtypes based on mRNA expression integrated
with PAM50 tumor profiling test.Methylation.Clusters – genome-wide DNA methylation pattern within different Breast Cancer
samples submitted for TCGA (based on PMID: 23000897).

### Cluster 3 
A hyper-methylated phenotype significantly enriched for Luminal B mRNA-subtype and
under-represented for PIK3CA and MAP3K1/MAP2K4 mutations.

### Cluster 5 
The lowest levels of DNA methylation, overlapped with the Basal-like mRNA-subtype, and
a high frequency of TP53 mutations.

### pathologic_N 
Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 
Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### pathologic_M 
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage 
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

### PR.status 
The presence of progesterone receptors on the surface of the cancer cell.

### BLCAhistological_subtype
Histologic subtypes of Bladder Cancer samples submitted for TCGA.

### neoplasm_histologic_grade 
Measure of anaplasia in a sampled tumors.

### pathologic_N 
Characterization of the regional lymph nodes.

NX: Cancer in nearby lymph nodes cannot be measured.

N0: There is no cancer in nearby lymph nodes.

N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the
number after the N, the more lymph nodes that contain cancer.

### pathologic_T 
Characterizationof the size and extent of the main tumor.

TX: Main tumor cannot be measured.

T0: Main tumor cannot be found.

T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T,
the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to
provide more detail, such as T3a and T3b.

### pathologic_M
Characterization of the distant metastasis.

MX: Metastasis cannot be measured.

M0: Cancer has not spread to other parts of the body.

M1: Cancer has spread to other parts of the body.

### pathologic_stage 
Combination of TNM results for each patient.

Stage 0 - cancer in situ.

Stage I -early-stage cancer; a small tumor without spreading to the lymph nodes or other parts of the
body.

Stage II and III - larger cancers or tumors with possible spreading to lymph nodes but not to other
parts of the body.

Stage IV - advanced or metastatic cancer.

### Smoking_history

1 –Lifelong Non-smoker

2 –Current smoker

3 –Current reformed smoker for > 15 years

4 –Current reformed smoker for < or = 15 years5 –Current Reformed Smoker, Duration Not Specified

### X_PANCAN_Cluster_Cluster_PANCAN
genome-wide DNA methylation pattern within different

### BLCA samples 
submitted for TCGA.
