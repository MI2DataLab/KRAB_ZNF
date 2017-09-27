library(plotly)
library(shinythemes)

load("data/expression_means.Rdata")
selected_group <- c("ZNF695",
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
all_genes <- colnames(expressions_means)
cohorts <- c("ACC",
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
             "LAML",
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
             "SKCM",
             "STAD",
             "STES",
             "TGCT",
             "THCA",
             "THYM",
             "UCEC",
             "UCS",
             "UVM")
cohorts_methylation <- c("BRCA",
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
                         "UCEC")

navbarPage("KRAB ZNF", theme = shinytheme("yeti"),
           tabPanel("About", 
             wellPanel(
             helpText(HTML("<b>ABOUT</b>")),
             helpText(HTML('The aim of this application is to analise and visualize survival data of cancer patients with respect to expressions of some seleted genes. 
                           <br><br>
                           It allows user to divide gene expressions into high/low expressions using different cutpoints and plot and download Kaplan-Meier curves.
                           <br><br>                           
                           The data is obtained from The Cancer Genome Atlas (TCGA).
                           <br><br>
                           App code as well as R scripts used for data preparation and analysis are accesible via <a href="https://github.com/rcylwa/survival_analysis_App">Github repository</a>.'))
             )),
           navbarMenu("Survival Analysis", 
           tabPanel("Survival plots",
                    sidebarLayout(
                      sidebarPanel(
                        selectInput("gene", "Choose a gene:", 
                                    choices = selected_group),
                        selectInput("cohort", "Choose a cohort",
                                    choices = cohorts, selected = "KIPAN"),
                        selectInput("curve_cutpoint", "Choose a method for cutpoint:",
                                    c("median", "max-rank statistic")),
                        checkboxInput("show_conf_int", "Show confidence interval", value = FALSE, width = NULL),
                        selectInput("conf_int_style", "Choose confidence interval style",
                                    c("ribbon", "step"), selected = "step"),
                      downloadButton(outputId = "survival_plot_png", label = "Download PNG plot"),
                      downloadButton(outputId = "survival_plot_pdf", label = "Download PDF plot"),
                      downloadButton(outputId = "survival_plot_eps", label = "Download EPS plot"),
                      width = 3),
                      
                      mainPanel(
                        plotOutput("survival_curve", width = "90%", height = 500)
                        )
                    )),
           tabPanel("Cutpoints for expressions",
                    sidebarLayout(
                      sidebarPanel(
                        wellPanel(
                          helpText(HTML("Visualisation of cutpoints for gene expression selected using maximally selected rank statistics."))
                          ),
                        selectInput("gene2", "Choose a gene:", 
                                    choices = selected_group),
                        selectInput("cohort2", "Choose a cohort",
                                    choices = cohorts, selected = "KIPAN"),
                        downloadButton(outputId = "cutpoint_plot_png", label = "Download PNG plot"),
                        downloadButton(outputId = "cutpoint_plot_pdf", label = "Download PDF plot"),
                        downloadButton(outputId = "cutpoint_plot_eps", label = "Download EPS plot"),
                        width = 3
                      ),
                      
                      # Show a summary of the dataset and an HTML table with the 
                      # requested number of observations
                      mainPanel(
                        plotOutput("cutpoint_plot")
                      ))
                    ),
           tabPanel("Heatmaps of p-values",
                    sidebarLayout(
                      sidebarPanel(
                        wellPanel(
                          helpText(HTML("Heatmap representing p-values for log-rank test 
                                        for each pair of gene and cohort. Dark blue indicates 
                                        significance in survial difference."))
                        ),
                        selectInput("method", "Choose cutpoint method",
                                    choices = c("median", "max-rank statistic")),
                        selectInput("scale", "Choose scale for p-values",
                                    choices = c("normal", "logarithmic")),
                        width = 3
                      ),
                      mainPanel(plotlyOutput("heathmap"))
                    ))
           ),
           navbarMenu("Methylation and expression",
                      tabPanel("High/low expressions tests", 
                               sidebarLayout(
                                 sidebarPanel(
                                   wellPanel(
                                     helpText(HTML("Table shows results of t-test of methylaion level between
                                                   a group of high expression and low expression 
                                                   conducted for each cohort and gene separately."))
                                     ),
                                   checkboxInput("only_significant_results", "Show only significant results", value = FALSE),
                                   selectInput("gene6", "Choose a gene", choices = all_genes),
                                   selectInput("cohort6", "Choose a cohort", choices = cohorts_methylation),
                                   downloadButton(outputId = "methylation_expression_correlation_csv", label = "Download CSV"), width = 3),
                                 mainPanel(dataTableOutput("methylation_expression_correlation_table"))
                                 ))
                      ),
           navbarMenu("Isoforms expression", 
           tabPanel("Isoforms expression normal vs cancer test", 
                    sidebarLayout(
                      sidebarPanel(
                        wellPanel(
                          helpText(HTML("Table shows results of t-test of gene isoform expression level between
                                        a group of normal tissue (group 11) and cancer tissue (group 01) 
                                        conducted for each cohort separately. Note that in the downloaded CSV \
                                        some p_values are \"NA\". That indicates that there were no observations 
                                        of one two groups (01 or 11)"))
                          ),
                        downloadButton(outputId = "isoform_test_all_csv", label = "Download CSV table")),
                      mainPanel(dataTableOutput("normal_vs_cancer_isoforms_table"))
                        )),
           tabPanel("Isoforms expression plots", 
                    sidebarLayout(
                      sidebarPanel(
                        wellPanel(
                          helpText(HTML("Plots present comparison between normal (11) and cancer (01) tissue 
                                        of expression levels of different isoforms 
                                        for selected gene and cancer cohort."))
                        ),
                        selectInput("gene3", "Choose a gene:", 
                                    choices = selected_group),
                        selectInput("cohort3", "Choose a cohort",
                                    choices = cohorts, selected = "BRCA"),
                        checkboxInput("boxplot_log_scale", label = "Log scale for the plot", value = TRUE),
                        width = 3
                      ),
                      mainPanel(plotOutput("normal_vs_cancer_isoforms_plot"))
                    )),
           tabPanel("Cancer tissue isoform expressions",
                    sidebarLayout(
                      sidebarPanel(
                        wellPanel(
                          helpText(HTML("Plots present expression levels of different isoforms 
                                        for selected gene and cancer cohort. Table below the plot 
                                        shows exact values of expression means for selected gene 
                                        and cohort."))
                          ),
                        selectInput("gene4", "Choose a gene:", 
                                    choices = selected_group),
                        selectInput("cohort4", "Choose a cohort",
                                    choices = cohorts, selected = "BRCA"),
                        width = 3
                          ),
                      mainPanel(plotOutput("isoform_expressions_plot"),
                                dataTableOutput("isoforms_expression_table"))
                        )),
           tabPanel("Gene isoform expressions across cohorts",
                    sidebarLayout(
                      sidebarPanel(
                        wellPanel(
                          helpText(HTML("Plots present, for each gene, percentages of total expression 
                                        of each isoform across all cancer cohorts."))
                          ),
                        selectInput("gene5", "Choose a gene:", 
                                    choices = selected_group),
                        width = 3
                          ),
                      mainPanel(plotOutput("all_cohorts_soform_expressions_plot"),
                                dataTableOutput("all_cohorts_isoforms_expression_table"))
                        ))
           ),
           navbarMenu("Expression in normal tissue",
                      tabPanel("Heatmap1",
                               sidebarLayout(
                                 sidebarPanel(
                                   wellPanel(
                                     helpText(HTML("You can view heatmap of mean expression either restricted to selected genes or select them manually."))
                                   ),
                                   selectInput("selected_genes_manually", "Or select genes manually:", choices = all_genes, selected = all_genes, multiple = TRUE),
                                   width = 3
                                 ),
                                 mainPanel(plotOutput("expression_heatmap"))
                               )),
                      tabPanel("Heatmap2",
                               sidebarLayout(
                                 sidebarPanel(
                                   wellPanel(
                                     helpText(HTML("You can view heatmap of mean expression either restricted to selected genes or select them manually."))
                                   ),
                                   checkboxInput("restrict", "Restrict to selected genes", value = TRUE),
                                   checkboxInput("log_scale", "Log scale", value = TRUE),
                                   selectInput("heatmap_color_scale", "Select color scale", 
                                               choices = rownames(RColorBrewer::brewer.pal.info), selected = "Greys"),
                                   selectInput("selected_genes_manually", "Or select genes manually:", choices = all_genes, selected = all_genes, multiple = TRUE),
                                   width = 3
                                 ),
                                 mainPanel(plotlyOutput("expression_heatmap2"))
                               )))
)
