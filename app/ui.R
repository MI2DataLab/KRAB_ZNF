library(plotly)

navbarPage("Survival analysis",
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
           
           tabPanel("Survival plots",
                    sidebarLayout(
                      sidebarPanel(
                        selectInput("gene", "Choose a gene:", 
                                    choices = c("ZNF695",
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
                                                "ZNF114")),
                        selectInput("cohort", "Choose a cohort",
                                    choices = c("ACC",
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
                                                "UVM"), selected = "KIPAN"),
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
                                    choices = c("ZNF695",
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
                                                "ZNF114")),
                        selectInput("cohort2", "Choose a cohort",
                                    choices = c("ACC",
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
                                                "UVM"), selected = "KIPAN"),
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
                    )),
           tabPanel("Data",
                    dataTableOutput("survival_data_table")),
           tabPanel("Links")
)