library(shiny)
library(RTCGA)
library(plotly)
library(pheatmap)
library(ggplot2)
library(ggthemes)
library(survminer)
library(survival)
library(plyr)

load("data/survival_data.Rdata")
load("data/survival_data_cat_median.Rdata")
load("data/survival_data_cat_optimal.Rdata")
load("data/pvals_optimal.Rdata")
load("data/pvals_median.Rdata")

shinyServer(function(input, output) {
  
  survival_data_selected <- function(){
    if(input$curve_cutpoint == "max-rank statistic") {
      selected_data <- survival_data_cat_optimal[survival_data_cat_optimal$dataset == input$cohort, 
                                                 c("times", "patient.vital_status", paste(input$gene))]
    }
    else {
      selected_data <- survival_data_cat_median[survival_data_cat_median$dataset == input$cohort,
                                                c("times", "patient.vital_status", paste(input$gene))]
    }
    colnames(selected_data)[3] <- "expression"
    return(selected_data)
  }
  
  survival_curve_plot <- reactive({
    selected_data <- survival_data_selected()
    #plot(selected_data[, "x2"])
    fit <- survfit(Surv(times, patient.vital_status) ~ expression, data = selected_data)

    selected_gene <- input$gene
    limit_time <- min(max(selected_data[selected_data[, "expression"] == "low", "times"]),
                      max(selected_data[selected_data[, "expression"] == "high", "times"]))
    n_low <- length(selected_data[selected_data[, "expression"] == "low", "times"])
    n_high <- length(selected_data[selected_data[, "expression"] == "high", "times"])
    par(oma=c(0,0,2,0))
    ggsurvplot(
      fit,                     # survfit object with calculated statistics.
      data = selected_data,
      risk.table = TRUE,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      conf.int = input$show_conf_int,         # show confidence intervals for
      conf.int.style = input$conf_int_style,
      title = paste("Cohort = ", input$cohort, ", gene = ", input$gene,
                    ", n_low = ", n_low, ", n_high = ", n_high),
      break.y.by = 0.1,
      xlim = c(0, limit_time),        # present narrower X axis, but not affect
      xlab = "Time (days)",
      # survival estimates.
      break.time.by = 1000,    # break X axis in time intervals by 500.
      palette = c("#9b9696", "#000000"),
      ggtheme = theme_bw(), # customize plot and risk table with a theme.
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = FALSE, # show bars instead of names in text annotations
      # in legend of risk table
      legend = "top",
      legend.title = NULL,
      surv.scale = "percent",
      break.x.by = 365,
      #risk.table.fontsize = 15,
      font.main = c(15, "bold"),
      font.tickslab = 12,
      font.x = 15,
      font.y = 15,
      font.legend = 15
    )
  }
  )
  cutpoint_plot <- function(){
    max_stat_cutpoint <- surv_cutpoint(
      survival_data[survival_data$dataset == input$cohort2, ],
      time = "times",
      event = "patient.vital_status",
      variables = c(input$gene2)
    )
    
    plot(max_stat_cutpoint, xlab = paste(input$gene2, " mRNA expression level"), 
         title = paste("Cohort: ", input$cohort2))
  }
  
  output$survival_curve <- renderPlot({
    survival_curve_plot()
  })
  
  output$cutpoint_plot <- renderPlot({
    cutpoint_plot()
   
  })
  
  output$survival_plot_png <- downloadHandler(
    filename = function() { paste(input$gene, input$cohort, '.png', sep='') },
    content = function(file) {
      png(file, width=800, height=600)
      par(oma=c(0,0,4,0))
      print(survival_curve_plot())
      dev.off()
    }
  )
  
  output$survival_plot_pdf <- downloadHandler(
    filename = function() { paste(input$gene, input$cohort, '.pdf', sep='') },
    content = function(file) {
      cairo_pdf(file)
      print(survival_curve_plot())
      dev.off()
    }
  )
  
  output$survival_plot_eps <- downloadHandler(
    filename = function() { paste(input$gene, input$cohort, '.eps', sep='') },
    content = function(file) {
      #setEPS()
      cairo_ps(file, width=8, height=6, fallback_resolution = 600)
      par(oma=c(0,0,4,0))
      print(survival_curve_plot())
      dev.off()
    }
  )
  
  output$cutpoint_plot_png <- downloadHandler(
    filename = function() { paste(input$gene, input$cohort, '.png', sep='') },
    content = function(file) {
      png(file)
      par(oma=c(2,2,2,2))
      print(cutpoint_plot())
      dev.off()
    }
  )
  
  output$cutpoint_plot_pdf <- downloadHandler(
    filename = function() { paste(input$gene, input$cohort, '.pdf', sep='') },
    content = function(file) {
      cairo_pdf(file)
      print(cutpoint_plot())
      dev.off()
    }
  )
  
  output$cutpoint_plot_eps <- downloadHandler(
    filename = function() { paste(input$gene, input$cohort, '.eps', sep='') },
    content = function(file) {
      cairo_ps(file, width=8, height=6, fallback_resolution = 600)
      #par(oma=c(0,0,4,0))
      print(cutpoint_plot())
      dev.off()
    }
  )
  
  # Generate a summary of the dataset
  output$summary <- renderPrint({
    print(input$gene)
  })
  
  output$heathmap <- renderPlotly({
    if(input$method == "median")
      m <- as.matrix(pvals_median[, 1:14])
    else
      m <- as.matrix(pvals_optimal[, 1:14])
    if(input$scale == "logarithmic")
      m <- log(m + 0.0001)
    plot_ly(
      x = colnames(m), y = rownames(m),
      z = m, type = "heatmap"
    ) %>%
      layout(margin = list(l = 130, b = 70), autosize = F, height = 600) %>%
      config(displayModeBar = T)
  })
  
  output$survival_data_table <- renderDataTable(
    survival_data
  )
  
})