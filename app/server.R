library(shiny)
library(RTCGA)
library(plotly)
library(pheatmap)
library(ggplot2)
library(ggthemes)
library(survminer)
library(survival)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(httpuv)
#library(RTCGA.methylation)
library(heatmaply)



shinyServer(function(input, output) {
  Sys.setlocale("LC_ALL","English")
  load("data/survival_data.Rdata")
  load("data/survival_data_cat_median.Rdata")
  load("data/survival_data_cat_optimal.Rdata")
  load("data/pvals_optimal.Rdata")
  load("data/pvals_median.Rdata")
  load("data/all_isoforms.rda")
  load("data/test_cancer_vs_normal.rda")
  load("data/expression_means.Rdata")
  load("data/expressions_all.Rdata")
  load("data/illumina_humanmethylation_27_data.rda")
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
  
  output$methylation_expression_correlation_table <- renderDataTable({
    methylation_expression_correlation()
  })
  
  output$methylation_expression_correlation_csv <- downloadHandler(
    filename = function() { paste('file.csv', sep='') },
    content = function(file) {
      content = methylation_expression_correlation()
      write.csv(content, file = file)
    }
  )
  
  
  
  methylation_expression_correlation <- reactive({
    gene <- input$gene6
    cohort <- input$cohort6
    k <- 10
    expressions_selected <- expressions_all[, c("bcr_patient_barcode", "dataset", gene)] %>% filter(dataset == cohort) %>% as.data.frame()
    expressions_selected[, gene] <- ntile(expressions_selected[, gene], k)
    
    # load files from prats
    files <- list.files(path = "data/methylation", pattern = cohort, full.names = TRUE)
    selected_methylation <- data.frame()
    for(file in files) {
      load(file)
      selected_methylation <- rbind(selected_methylation, methylation_part)
    }
    
    selected_methylation <- selected_methylation[, !is.na(selected_methylation[1, ])]
    methylation_expression <- expressions_selected %>% filter(.[[3]] %in% c(1, 10)) %>%inner_join(selected_methylation)
    cpg_islands <- colnames(selected_methylation)[-1]
    
    row_n <- 1
    pvals_list <- list()
    for(k in cpg_islands) {
      low_expr <- methylation_expression[methylation_expression[, 3] == 1 , k]
      high_expr <- methylation_expression[methylation_expression[, 3] == 10 , k]
      pvalue <- NA
      try({
        test <- t.test(low_expr, high_expr)
        pvalue <- test$p.value
      }, silent = TRUE)
      pvals_list[[row_n]] <- data.frame("cohort" = cohort,
                                        "gene" = gene,
                                        "cpg_island" = k,
                                        "pvalue" = pvalue,
                                        "low_expr_meth_mean" = mean(low_expr, na.rm = TRUE),
                                        "high_expr_meth_mean" = mean(high_expr, na.rm = TRUE))
      row_n <- row_n + 1
    }
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

    if(input$only_significant_results)
      pvals_extended_signif
    else
      pvals_extended
    })
  
  output$normal_vs_cancer_isoforms_plot <- renderPlot({
    plot_data <- all_isoforms[all_isoforms$cohort == input$cohort3 & 
                           all_isoforms$gene == input$gene3 & 
                           all_isoforms$group %in% c("01", "11"), ] %>% 
                  mutate(group = sub(group, pattern = "01", replacement = "Cancer")) %>%
                  mutate(group = sub(group, pattern = "11", replacement = "Normal"))
    ylim1 = boxplot.stats(plot_data$value)$stats[c(1, 5)]
    plot_data$group <- as.factor(plot_data$group)
    group_colors <- c("Cancer" = "#333BFF", "Normal" = "#CC6600")
    if(input$boxplot_log_scale) {
    ggplot(plot_data, aes(isoform, log(value))) +
      geom_boxplot(aes(colour = group)) +
      #coord_cartesian(ylim = ylim1*4) +
      ylab("log(expression level)") + 
      scale_fill_manual(values = group_colors)
    }
    else {
      ggplot(plot_data, aes(isoform, value)) +
        geom_boxplot(aes(colour = group)) +
        #coord_cartesian(ylim = ylim1*4) +
        ylab("expression level") + 
        scale_fill_manual(values = group_colors)
    }
    
  })
  


  
  normal_vs_cancer_isoforms_table_reactive <- reactive({
    isoform_name_mapping <- read.table(file = "data/knownToRefSeq.txt")
    colnames(isoform_name_mapping) <- c("isoform", "GenBank_name")
    isoform_name_mapping$isoform_cut <- substr(isoform_name_mapping$isoform, 1, 8)
    isoform_summary_table_tmp <- all_isoforms %>% filter(group %in% c("01", "11")) %>%
      mutate(group = sub(x = group, "01", "cancer tissue")) %>%
      mutate(group = sub(x = group, "11", "normal tissue")) %>%
      group_by(cohort, gene, isoform, group) %>%
      dplyr::summarize(mean_expression = mean(value), 
                       `25%`= quantile(value, probs=0.25, na.rm = TRUE),
                       `50%`= quantile(value, probs=0.5),
                       `75%`= quantile(value, probs=0.75),
                       n_onservations = n()) %>%
      mutate(isoform_cut = substr(isoform, 1, 8))
    test_cancer_vs_normal %>% filter(p_value != "NA") %>% 
    left_join(isoform_summary_table_tmp) %>%
    left_join(isoform_name_mapping[, -1])
    })
  
  output$normal_vs_cancer_isoforms_table <- renderDataTable({
    normal_vs_cancer_isoforms_table_reactive()
  }, options = list(pageLength = 10))
  
  output$isoform_test_all_csv <- downloadHandler(
    filename = function() { paste('normal_vs_caner_tests.csv', sep='') },
    content = function(file) {
      write.csv(normal_vs_cancer_isoforms_table_reactive(), file = file)
    }
  )
  
  output$isoform_expressions_plot <- renderPlot({
    all_isoforms %>% filter(group == "01") %>%
                     group_by(cohort, gene, isoform) %>%
                     summarise(sum_of_expression = sum(value)) -> a
    
    all_isoforms %>% filter(group == "01") %>%
                     group_by(cohort, gene) %>%
                     summarise(total_expression = sum(value)) -> b
    
    a %>% left_join(b) %>%
          mutate(expr_percentage = 100*sum_of_expression/total_expression) %>%
          filter(cohort == input$cohort4, gene == input$gene4) %>%
          ggplot(aes(reorder(isoform, expr_percentage), expr_percentage)) + 
          geom_bar(stat = "identity") + 
          coord_flip() + 
          xlab("isoform") + 
          ggtitle(paste(input$cohort4, ", ", input$gene4, sep = "")) + 
          scale_fill_manual(values="Greys")
  })
  
  output$all_cohorts_soform_expressions_plot <- renderPlot({
    all_isoforms %>% filter(group == "01") %>% 
      group_by(cohort, gene, isoform) %>% 
      summarise(sum_of_expression = sum(value)) -> a
    
    all_isoforms %>% filter(group == "01") %>%
      group_by(cohort, gene) %>%
      summarise(total_expression = sum(value)) -> b
    
    a %>% left_join(b) %>%
      mutate(expr_percentage = 100*sum_of_expression/total_expression) %>%
      filter(gene == input$gene5) %>%
      ggplot(aes(cohort, expr_percentage)) + 
      geom_bar(aes(fill = isoform), stat = "identity") + 
      theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
      ggtitle(paste(input$gene5))
  })
  
  
  output$expression_heatmap <- renderPlot({
    if(input$restrict)
      pheatmap(log(t(expressions_means[, selected_group])), cellwidth = 25, cellheight = 25)
    else
      pheatmap(log(t(expressions_means[, input$selected_genes_manually]) + 1), cellwidth = 20, cellheight = 20)
    }, height = function() {
      if(!input$restrict) length(input$selected_genes_manually)*22 + 10
      else "auto"})

  output$expression_heatmap2 <- renderPlotly({
    if(input$restrict) {
      if(input$log_scale)
        z <- log(t(expressions_means[, selected_group]) + 1) %>% as.matrix()
      else
        z <- t(expressions_means[, selected_group]) %>% as.matrix()
      height <- 400
    }
    else {
      if(input$log_scale)
        z <- log(t(expressions_means[, input$selected_genes_manually] + 1)) %>% as.matrix()
      else
        z <- t(expressions_means[, input$selected_genes_manually]) %>% as.matrix()
      height <- length(input$selected_genes_manually)*22 + 10
      }
      # plot_ly(z = z,
      #         x = colnames(z),
      #         y = rownames(z),
      #         type = "heatmap",
      #         height = height,
      #         colors = input$heatmap_color_scale
      #         )
    color_palette <- colorRampPalette(brewer.pal(nrow(z), input$heatmap_color_scale))
    heatmaply(x = z, 
              xlab = colnames(z),
              ylab = rownames(z),
              colors = color_palette
              )

  })
  
}) 