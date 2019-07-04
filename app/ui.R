library("plotly")
library("shiny")
library("shinythemes")

# precalculated data
# Gene expression data
load("data/selected_genes_all.RData")
load("data/selected_genes_all_restricted.RData")
# Methylation data
source('data/krabmen/R/read_data.R')

# inicialize dictionaries
tumor.names <- load.tumor.names()
gene.list <- load.gene.list()
gene.names <- gene.list

selected_group <- c(
  "ZNF695",
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
  "ZNF114"
)
all_genes <- selected_genes_all
cohorts <- c(
  "ACC",
  "BLCA",
  "BRCA",
  "CESC",
  "CHOL",
  "COAD",
  "COADREAD",
  "DLBC",
  "ESCA",
  "GBM",
  "GBMLGG",
  "HNSC",
  "KICH",
  "KIPAN",
  "KIRC",
  "KIRP",
  "LGG",
  "LIHC",
  "LUAD",
  "LUSC",
  "OV",
  "PAAD",
  "PCPG",
  "PRAD",
  "READ",
  "SARC",
  "STAD",
  "STES",
  "TGCT",
  "THCA",
  "THYM",
  "UCEC",
  "UCS",
  "UVM"
)

cohorts_methylation <- c(
  "BRCA",
  "COAD",
  "COADREAD",
  "GBM",
  "GBMLGG",
  "KIPAN",
  "KIRC",
  "KIRP",
  "LAML",
  "LUAD",
  "LUSC",
  "OV",
  "READ",
  "STAD",
  "STES",
  "UCEC"
)

# navbar part of the interface
navbarPage(
  title = "KRAB ZNF explorer",
  theme = shinytheme("yeti"),
  tabPanel("About",
           wellPanel(helpText(
             HTML("<b>ABOUT</b>")
           ),
           helpText(
             HTML(
               '<b> KRAB-ZNF explorer - expression of KRAB-ZNF- transcription factors in The Cancer Genome Atlas study</b>
               <br>
               KRAB-ZNF explorer is a web-based application for analyzing the expression of KRAB-ZNF transcription factors in the data from The Cancer Genome Atlas study.
               <br>
               Key functionalities cover:
               <br>
               1) comparative analysis of KRAB-ZNF genes between normal and tumor samples, <br>
               2) correlation of KRAB-ZNF expression with various clinico-pathological parameters, <br>
               3) analysis and visualization of the link between patient survival and expression of KRAB-ZNF genes, <br>
               4) analyses of the association between KRAB-ZNF expression and CpG methylation status,  <br>
               5) analyses of isoform expressions for KRAB-ZNF genes in normal and tumor samples, <br>
               4) comparative analysis of KRAB-ZNF expression in normal tissues <br>
               <br>
               KRAB-ZNF explorer is developed by Rafal Cylwa and collaborators: Urszula Oleksiewicz, Marta Machnik, Kornel Kielczewski and Przemyslaw Biecek. <br>
               Exploration is based on a shiny library [1], plots are avaliable through plotly [2], pheatmap [3] and ggplot2 [4] libraries. Survival analysis is done with survminer [5] package.
               Analyses are performed based on <a href="https://github.com/RTCGA/RTCGA">The Cancer Genome Atlas (TCGA) data</a> accessed with RTCGA package [6] . <br>
               R sources of the application are in this <a href="https://github.com/MI2DataLab/KRAB_ZNF">Github repository</a>. Example application is presented in [7].<br>
               <br>
               References:<br>
                [1] Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2019). <b>shiny: Web Application Framework for R.</b> R package version 1.3.2. https://CRAN.R-project.org/package=shiny<br>
                [2] Carson Sievert (2018) <b>plotly for R.</b> https://plotly-r.com<br>
                [3] Raivo Kolde (2019). <b>pheatmap: Pretty Heatmaps.</b> R package version 1.0.12. https://CRAN.R-project.org/package=pheatmap<br>
                [4] H. Wickham. <b>ggplot2: Elegant Graphics for Data Analysis.</b> Springer-Verlag New York, 2016<br>
                [5] Alboukadel Kassambara and Marcin Kosinski (2019). <b>survminer: Drawing Survival Curves using `ggplot2`.</b> R package version 0.4.4. https://CRAN.R-project.org/package=survminer<br>
                [6] Marcin Kosinski and Przemyslaw Biecek (2019). <b>RTCGA: The Cancer Genome Atlas Data Integration.</b> R package version 1.14.0. https://rtcga.github.io/RTCGA<br>
                [7] Marta Machnik, Rafał Cylwa, Kornel Kiełczewski, Przemysław Biecek, Triantafillos Liloglou, Andrzej Mackiewicz, and Urszula Oleksiewicz (2019). <b>The expression signature of cancer-associated KRAB-ZNF factors identified in TCGA pan-cancer transcriptomic data.</b>. Mol Oncol. 2019 Apr;13(4):701-724. doi: 10.1002/1878-0261.12407.
               '
             )
           )),
           br(),
           br(),
           tags$a(img(src = "MI2_logo.jpg", width = "100px", height = "100px", align = "left"), href = "http://mi2.mini.pw.edu.pl/index.php/mi2-datalab/"),
           tags$a(img(src = "tcga_logo.svg", width = "400px", height = "100px", align = "left"), href = "https://cancergenome.nih.gov/")
  ),
  tabPanel("Expression in Tumor vs Normal",
           sidebarLayout(
             sidebarPanel(
               wellPanel(helpText(
                 HTML(
                   "This panel helps to compare genes expression between tumor and healthy sample for selected cohort and selected gene. Parameters in the left panel help to select suitable cohort/gene and tune the graphical elements in the plot.<br><br>
                   Boxplots compare expression levels of a selected gene between normal and cancer tissue. An bottom table shows descriptive statistics and t-test results of the same comparisons."
                 )
               )),
               selectInput('select.tumor', label = 'Tumor',
                           choices = tumor.names,
                           selected = tumor.names[1], multiple = T),

               hr(),

               selectInput('select.gene', label = 'Gene',
                           choices = gene.list,
                           selected = gene.list[1], multiple = T),

               hr(),

               selectInput('select.scale', label = 'Scale',
                           choices = c("linear", "log"),
                           selected = 1),

               numericInput('scale.min', 'Scale minimum', 1, min = 0, max = NA, step = 0.1),

               uiOutput('breaks'),

               hr(),

               numericInput('input.font.size', label = 'Font size',
                            value = 24, min = 1, max = 1024),

               selectInput('select.file.type', label = 'Plot file type',
                           choices = c('tiff', 'pdf', 'eps'),
                           selected = 1),

               checkboxInput('greyscale', label = 'Greyscale', value = F),

               downloadButton('download.plot'),
               hr(),
               selectInput('select.table.file.type', label = 'Table file type',
                           choices = c('csv', 'txt'),
                           selected = 1),
               downloadButton('pvalues.download')
             ),
             mainPanel(
               plotOutput('distPlot'),
               dataTableOutput(outputId = 'pvalues'),
               HTML("Legend: <br>
                    se - standard error, <br>
                    25% - 1st quartile, 25 centile,<br>
                    50% - median, 2nd quartile, 50 centile,<br>
                    75% - 3rd quartile, 75 centile.")
             )
           )),
  tabPanel("Clinical Parameters",
           sidebarLayout(
             sidebarPanel(
               wellPanel(helpText(
                 HTML(
                   "This panel helps to understand differences in gene expression between different histological types of a selected tumor. If some types are linked with significantly low or high gene expression, then corresponding boxplots will be shifted.<br><br>
                   Each boxplot shows median, and quartile expression (boundies od the box), min and max expression and outliers. This plot compares of expression levels of a selected gene across selected subtypes. An bottom table shows descriptive statistics and t-test results of the same comparisons."
                 )
               )),
               HTML("Clinical parameters description in "), a(href="https://raw.githack.com/MI2DataLab/KRAB_ZNF/master/app/data/clinical_parameters.html", "html"),
               downloadButton(outputId = "clinical.parameters.description", label = "pdf"),
               selectInput('tumor.name.subtypes', label = 'Tumor',
                           choices = setdiff(tumor.names, c("STES", "KIPAN", "CHOL")),
                           selected = tumor.names[1],
                           multiple = F),

               selectInput('gene.names.subtypes', label = 'Gene',
                           choices = gene.names,
                           selected = gene.names[1], multiple = T),

               uiOutput('subtypes'),

               uiOutput('subtype.values'),

               hr(),

               selectInput('select.scale.subtypes', label = 'Scale',
                           choices = c('linear', 'log'),
                           selected = 1),

               numericInput('scale.min.subtypes', 'Scale minimum', 1, min = 0, max = NA, step = 0.1),

               uiOutput("breaks.subtypes"),

               hr(),

               numericInput('input.font.size.subtypes', label = 'Font size',
                            value = 24, min = 1, max = 1024),

               checkboxInput('rotate.x.axis.subtypes', label = 'Rotate x axis text',
                             value = F),

               hr(),

               selectInput('select.file.type.subtypes', label = 'File type',
                           choices = c('tiff', 'pdf', 'eps'),
                           selected = 2),

               downloadButton(outputId = "plot.subtypes.download", label = "Download plot"),
               hr(),
               selectInput('select.table.subtypes.file.type', label = 'Table file type',
                           choices = c('csv', 'txt'),
                           selected = 1),

               downloadButton(outputId = "table.subtypes.download", label = "Download table")

             ),

             mainPanel(plotOutput('distPlot.subtypes'),
                       dataTableOutput('pvalues.subtypes'))
           )),
  navbarMenu(
    "Survival Analysis",
    tabPanel("Survival plots",
             sidebarLayout(
               sidebarPanel(
                 wellPanel(helpText(
                   HTML(
                     "This panel helps to understand if a expression of a given gene is linked with patients survival.<br><br>
                      Kaplan-Meier estimates of survival curves for a chosen gene and cohort along with the at-risk table and p-values. This plot is generated with the survminer package."
                   )
                 )),
                 selectInput("gene", "Choose a gene:",
                             choices = selected_genes_all_restricted),
                 selectInput(
                   "cohort",
                   "Choose a cohort",
                   choices = cohorts,
                   selected = "KIPAN"
                 ),
                 selectInput(
                   "curve_cutpoint",
                   "Choose a method for cutpoint:",
                   c("median", "max-rank statistic")
                 ),
                 checkboxInput(
                   "show_conf_int",
                   "Show confidence interval",
                   value = TRUE,
                   width = NULL
                 ),
                 selectInput(
                   "conf_int_style",
                   "Choose confidence interval style",
                   c("ribbon", "step"),
                   selected = "ribbon"
                 ),
                 sliderInput(
                   "survival_font_main",
                   "Font main:",
                   min = 5,
                   max = 30,
                   value = 15
                 ),
                 sliderInput(
                   "survival_font_tickslab",
                   "Font tickslab:",
                   min = 5,
                   max = 30,
                   value = 12
                 ),
                 sliderInput(
                   "survival_font_xlab",
                   "Font xlab:",
                   min = 5,
                   max = 30,
                   value = 15
                 ),
                 sliderInput(
                   "survival_font_ylab",
                   "Font ylab:",
                   min = 5,
                   max = 30,
                   value = 15
                 ),
                 sliderInput(
                   "survival_font_legend",
                   "Font legend:",
                   min = 5,
                   max = 30,
                   value = 15
                 ),
                 sliderInput(
                   "survival_height_px",
                   "Downloaded plot height in px:",
                   min = 200,
                   max = 2000,
                   value = 600,
                   step = 10
                 ),
                 sliderInput(
                   "survival_width_px",
                   "Downloaded plot width in px:",
                   min = 200,
                   max = 2000,
                   value = 800,
                   step = 10
                 ),
                 downloadButton(outputId = "survival_plot_png", label = "Download PNG plot"),
                 downloadButton(outputId = "survival_plot_pdf", label = "Download PDF plot"),
                 downloadButton(outputId = "survival_plot_eps", label = "Download EPS plot"),
                 downloadButton(outputId = "survival_plot_tiff", label = "Download TIFF plot"),
                 width = 3
               ),
               mainPanel(plotOutput("survival_curve", width = "90%", height = 500))
             )),
    tabPanel(
      "Cutpoints for expressions",
      sidebarLayout(
        sidebarPanel(
          wellPanel(helpText(
            HTML(
              "This panel helps to understand how one should discretize expression of a selected gene.<br><br>
              Visualisation of cutpoints for the expression of a selected gene using maximally selected rank statistics is generated with the survminer package."
            )
          )),
          selectInput("gene2", "Choose a gene:",
                      choices = selected_genes_all_restricted),
          selectInput(
            "cohort2",
            "Choose a cohort",
            choices = cohorts,
            selected = "KIPAN"
          ),
          sliderInput(
            "cutpoint_font_main",
            "Font main:",
            min = 5,
            max = 30,
            value = 15
          ),
          sliderInput(
            "cutpoint_font_tickslab",
            "Font tickslab:",
            min = 5,
            max = 30,
            value = 12
          ),
          sliderInput(
            "cutpoint_font_xlab",
            "Font xlab:",
            min = 5,
            max = 30,
            value = 15
          ),
          sliderInput(
            "cutpoint_font_ylab",
            "Font ylab:",
            min = 5,
            max = 30,
            value = 15
          ),
          sliderInput(
            "cutpoint_font_legend",
            "Font legend:",
            min = 5,
            max = 30,
            value = 15
          ),
          sliderInput(
            "cutpoint_height_px",
            "Downloaded plot height in px:",
            min = 200,
            max = 2000,
            value = 600,
            step = 10
          ),
          sliderInput(
            "cutpoint_width_px",
            "Downloaded plot width in px:",
            min = 200,
            max = 2000,
            value = 800,
            step = 10
          ),
          downloadButton(outputId = "cutpoint_plot_png", label = "Download PNG plot"),
          downloadButton(outputId = "cutpoint_plot_pdf", label = "Download PDF plot"),
          downloadButton(outputId = "cutpoint_plot_eps", label = "Download EPS plot"),
          downloadButton(outputId = "cutpoint_plot_tiff", label = "Download TIFF plot"),
          width = 3
        ),

        # Show a summary of the dataset and an HTML table with the
        # requested number of observations
        mainPanel(plotOutput("cutpoint_plot"))
      )
    ),
    tabPanel(
      "Heatmaps of p-values",
      sidebarLayout(
        sidebarPanel(
          wellPanel(helpText(
            HTML(
              "This panel helps to understand which genes are linked with survival in particular tumors.<br><br>
              Heatmap representing p-values for the log-rank test for each pair of gene and cohort. A low p-value indicates a significant difference in survival. To customize the gene list delete all the genes from the box (Ctrl+A+Del) and select preferred ones via dropdown menu.
              p-values are after the FDR correction. The darkest color corresponds to the p-value < 0.1. Additional dendrograms help to understand which genes and which tumors are more similar to each other."
            )
          )),
          selectInput(
            "method",
            "Choose cutpoint method",
            choices = c("median", "max-rank statistic")
          ),
          selectInput(
            "scale",
            "Choose scale for p-values",
            choices = c("normal", "logarithmic")
          ),
          selectInput(
            "heatmap_log_rank_scale",
            "Select color scale",
            choices = rownames(RColorBrewer::brewer.pal.info),
            selected = "Purples"
          ),
          checkboxInput("all_genes_logrank", "Select all genes", value = TRUE),
          selectInput(
            "selected_genes_manually_log_rank",
            "Select genes manually:",
            choices = selected_genes_all_restricted,
            selected = selected_genes_all_restricted,#selected_group,
            multiple = TRUE
          ),
          sliderInput(
            "log_rank_heatmap_fontsize",
            "Font size:",
            min = 5,
            max = 30,
            value = 15
          ),
          sliderInput(
            "log_rank_heatmap_fontsize_row",
            "Row font size:",
            min = 5,
            max = 30,
            value = 10
          ),
          sliderInput(
            "log_rank_heatmap_fontsize_col",
            "Column font size:",
            min = 5,
            max = 30,
            value = 10
          ),
          selectInput('select_file_type_log_rank_heatmap', label = 'File type',
                      choices = c('tiff', 'pdf', 'eps'),
                      selected = 1),
          downloadButton('download_heatmap_log_rank'),
          width = 3
        ),
        mainPanel(plotOutput("heatmap_log_rank_plot", width = "100%"))
      )
    ),
    tabPanel(
      "Log-rank test tables",
      sidebarLayout(
        sidebarPanel(
          wellPanel(helpText(
            HTML(
              "This panel helps to understand which genes are linked with survival in particular tumors. It supplements the heatmap presented in previous tab.<br><br>
              A table containing results of a log-rank test (null hypothesis = the two groups have identical survival and hazard functions) comparing survival of patients from two groups (a low and high group of expression of a given gene for each cohort separately chosen according to the median or optimal cutpoint). For each cohort and gene, we provide p-value of a log-rank test, numbers of observations in both groups and hazard ratio."
            )
          )),
          selectInput(
            "log_rank_table_method",
            "Choose cutpoint method",
            choices = c("median", "max-rank statistic")
          ),
          downloadButton(outputId = "log_rank_table_csv", label = "Download CSV table"),
          downloadButton(outputId = "log_rank_table_txt", label = "Download TXT table")
        ),
        mainPanel(dataTableOutput("log_rank_table"))
      )
    )
  ),
  navbarMenu(
    "Methylation and expression",
    tabPanel(
      "High/low expressions tests",
      sidebarLayout(
        sidebarPanel(
          wellPanel(helpText(
            HTML(
              "This panel helps to understand which genes and which cpg islands are linked with differences in gene expression.<br><br>
              The table shows results of t-test (p-values are FDR adjusted) for methylation level difference between a group with high expression (top 10%) and low expression (bottom 10%) of a given KRAB-ZNF factor conducted for each cohort and gene separately. <br/> Hit 'Generate output' to run calculations. Note that they may be time-consuming as we need to process all cpg island."
            )
          )),
          checkboxInput(
            "only_significant_results",
            "Show only significant results (FDR adjusted p-value < 0.05, fold-change > 4)",
            value = FALSE
          ),
          selectInput("gene6", "Choose a gene", choices = all_genes),
          selectInput("cohort6", "Choose a cohort", choices = cohorts_methylation),
          actionButton("perform_methylation", "Generate output"),
          hr(),
          downloadButton(outputId = "methylation_expression_correlation_csv", label = "Download CSV"),
          width = 3
        ),
        mainPanel(
          dataTableOutput("methylation_expression_correlation_table")
        )
      )
    )
  ),
  navbarMenu(
    "Isoforms expression",
    tabPanel(
      "Isoforms expression normal vs cancer test",
      sidebarLayout(
        sidebarPanel(
          wellPanel(helpText(
            HTML(
              "This panel helps to understand which isoforms are expressed differently. It supplements boxplots presented in the next tab.<br><br>
              A table shows results of a t-test of gene isoform expression level between a group of normal tissue and cancer tissue conducted for each gene and cohort separately with a joined summary statistics for each group. Click Generate Output to view and download the table. Empty table means no data in TCGA available for selected variables. P-value equal to NA indicates an insufficient amount of data to conduct the test."
            )
          )),
          selectInput("cohort_isoform_table", "Choose cohort", choices = cohorts, selected = "BRCA"),
          selectInput("gene_isoform_table", "Choose gene", choices = all_genes),
          actionButton("perform_isoform_table", "Generate output"),
          hr(),
          downloadButton(outputId = "isoform_test_all_csv", label = "Download CSV table"),
          hr(),
          downloadButton(outputId = "isoform_test_all_txt", label = "Download TXT table")
        ),
        mainPanel(dataTableOutput("normal_vs_cancer_isoforms_table"),
                  HTML("Legend: <br>
                    se - standard error, <br>
                    25% - 1st quartile, 25 centile,<br>
                    50% - median, 2nd quartile, 50 centile,<br>
                    75% - 3rd quartile, 75 centile."))
      )
    ),
    tabPanel(
      "Isoforms expression plots",
      sidebarLayout(
        sidebarPanel(
          wellPanel(helpText(
            HTML(
              "This panel helps to understand which isoforms are expressed differently. It supplements the table presented in the previouse tab.<br><br>
                Boxplots present the expression of different isoforms of a selected gene in normal and cancer tissues. "
            )
          )),
          selectInput("gene3", "Choose a gene:",
                      choices = all_genes, selected = "ZNF695"),
          selectInput(
            "cohort3",
            "Choose a cohort",
            choices = cohorts,
            selected = "BRCA"
          ),
          checkboxInput("boxplot_log_scale", label = "Log scale for the plot", value = TRUE),
          checkboxInput("boxplot_grey_scale", label = "Grey scale for the plot", value = TRUE),
          uiOutput('breaks.isoforms'),
          selectInput("select_isoform_boxplot_xlab_font", label = "Select font size", choices = 1:25, selected = 12),
          selectInput('select.file.type.normal.cancer', label = 'File type',
                      choices = c('tiff', 'pdf', 'eps'),
                      selected = 1),
          downloadButton('download.plot.normal.cancer'),
          width = 3
        ),
        mainPanel(plotOutput("normal_vs_cancer_isoforms_plot"))
      )
    ),
    tabPanel(
      "Cancer tissue isoform expressions",
      sidebarLayout(
        sidebarPanel(
          wellPanel(helpText(
            HTML(
              "This panel helps to understand which isoforms are more common in particular tumor. It supplements the table and boxplots presented in the previouse tabs.<br><br>
                Each bar in the plot represents a percentage of total gene expression for a particular isoform of a selected gene."
            )
          )),
          selectInput("gene4", "Choose a gene:",
                      choices = all_genes),
          selectInput(
            "cohort4",
            "Choose a cohort",
            choices = cohorts,
            selected = "BRCA"
          ),
          selectInput("select_isoform_expression_font", label = "Select font size", choices = 1:25, selected = 12),
          selectInput('select.file.type.isoform.expression', label = 'File type',
                      choices = c('tiff', 'pdf', 'eps'),
                      selected = 1),
          downloadButton('download_isoform_expression_plot'),
          width = 3
        ),
        mainPanel(
          plotOutput("isoform_expressions_plot"),
          dataTableOutput("isoforms_expression_table")
        )
      )
    ),
    tabPanel(
      "Isoform expressions across cohorts",
      sidebarLayout(
        sidebarPanel(
          wellPanel(helpText(
            HTML(
              "This panel helps to understand which isoforms are more common in particular tumor. It supplements plots presented in the previous tabs but shows all tumors in a single plot. Isoformes are color coded. This way it is easier to notice differences between tumors in isoformes composition.<br><br>
                Plots present, for each gene, percentages of total expression of each isoform across all cancer cohorts."
            )
          )),
          selectInput("gene5", "Choose a gene:",
                      choices = all_genes),
          selectInput("select_all_cohorts_isoform_expression_font", label = "Select font size", choices = 1:25, selected = 12),
          selectInput('select.file.type.all.cohorts.isoform.expression', label = 'File type',
                      choices = c('tiff', 'pdf', 'eps'),
                      selected = 1),
          downloadButton('download_all_cohorts_isoform_expression_plot'),
          width = 3
        ),
        mainPanel(
          plotOutput("all_cohorts_isoform_expressions_plot"),
          dataTableOutput("all_cohorts_isoforms_expression_table")
        )
      )
    )
  ),
  navbarMenu(
    "Expression in normal tissue",
    tabPanel("Heatmap",
             sidebarLayout(
               sidebarPanel(
                 wellPanel(helpText(
                   HTML(
                     "This panel helps to understand which genes (rows) are differentially expressed in different tumors (columns).  Different colors in rows identifies genes with different expressions across tumors. This plots supplements the table in the next panel.<br><br>
                    Heatmap presents mean (or log of mean) expression of selected genes. To customize the gene list delete all the genes from the box (Ctrl+A+Del) and select preferred ones via dropdown menu."
                   )
                 )),
                 #checkboxInput("restrict", "Restrict to selected genes", value = TRUE),
                 checkboxInput("log_scale", "Log scale", value = TRUE),
                 selectInput(
                   "heatmap_color_scale1",
                   "Select color scale",
                   choices = rownames(RColorBrewer::brewer.pal.info),
                   selected = "Greys"
                 ),
                 checkboxInput("all_genes_normal", "Select all genes", value = TRUE),
                 selectInput(
                   "selected_genes_manually",
                   "Or select genes manually:",
                   choices = all_genes,
                   selected = all_genes,
                   multiple = TRUE
                 ),
                 sliderInput(
                   "normal_heatmap_fontsize",
                   "Font size:",
                   min = 5,
                   max = 30,
                   value = 15
                 ),
                 sliderInput(
                   "normal_heatmap_fontsize_row",
                   "Row font size:",
                   min = 5,
                   max = 30,
                   value = 10
                 ),
                 sliderInput(
                   "normal_heatmap_fontsize_col",
                   "Column font size:",
                   min = 5,
                   max = 30,
                   value = 10
                 ),
                 selectInput('select_file_type_expression_heatmap', label = 'File type',
                             choices = c('tiff', 'pdf', 'eps'),
                             selected = 1),
                 downloadButton('download_expression_heatmap1'),
                 width = 3
               ),
               mainPanel(plotOutput("expression_heatmap1"))
             )),
    tabPanel("Expressions boxplots and table",
             sidebarLayout(
               sidebarPanel(
                 wellPanel(helpText(
                   HTML(
                     "This panel helps to understand which genes (rows) are differentially expressed in different tumors (columns). This plots supplements the heatmap presented in the previous panel.<br><br>
                    Boxplot presents a natural log of expression across all normal tissues along with a summary table with the data from normal tissues."
                   )
                 )),
                 selectInput("gene_normal_expression", "Select gene", choices = all_genes, selected = all_genes[1]),
                 selectInput("select_normal_boxplot_font", label = "Select font size", choices = 1:25, selected = 12),
                 selectInput('select_file_type_expression_boxplot', label = 'File type',
                             choices = c('tiff', 'pdf', 'eps'),
                             selected = 1),
                 downloadButton("normal_expression_boxplot_download", label = "Download plot"),
                 hr(),
                 selectInput('select.normal.table.file.type', label = 'File type',
                             choices = c('csv', 'txt'),
                             selected = 1),
                 downloadButton("normal_expression_pvalues_table_download", label = "Download table")
               ),
               mainPanel(plotOutput("normal_expression_boxplot", height = 600),
                         dataTableOutput("normal_expression_pvalues_table_render"),
                         HTML("Legend: <br>
                    se - standard error, <br>
                    25% - 1st quartile, 25 centile,<br>
                    50% - median, 2nd quartile, 50 centile,<br>
                    75% - 3rd quartile, 75 centile."))
             ))
  )
)
