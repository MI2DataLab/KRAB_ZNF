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
library(heatmaply)
library(tidyr)
#library(RTCGA.methylation)


shinyServer(function(input, output) {
  # Sys.setlocale("LC_ALL","English")
  load("data/test_results_median.RData")
  load("data/test_results_optimal.RData")
  load("data/survival_data.RData")
  load("data/survival_data_cat_median.RData")
  load("data/survival_data_cat_optimal.RData")
  load("data/pvals_optimal.RData")
  load("data/pvals_median.RData")
  load("data/all_isoforms.rda")
  #load("data/test_cancer_vs_normal.rda")
  load("data/expression_means.RData")
  load("data/expressions_all.RData")
  load("data/expressions_normal.RData")
  load("data/illumina_humanmethylation_27_data.rda")
  load("data/selected_genes_all.RData")
  load("data/selected_genes_all_restricted.RData")
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
  
  survival_data_selected <- function() {
    if (input$curve_cutpoint == "max-rank statistic") {
      selected_data <-
        survival_data_cat_optimal[survival_data_cat_optimal$dataset == input$cohort,
                                  c("times", "patient.vital_status", paste(input$gene))]
    }
    else {
      selected_data <-
        survival_data_cat_median[survival_data_cat_median$dataset == input$cohort,
                                 c("times", "patient.vital_status", paste(input$gene))]
    }
    colnames(selected_data)[3] <- "expression"
    selected_data$times <- selected_data$times / 365
    return(selected_data)
  }
  
  survival_curve_plot <- reactive({
    selected_data <- survival_data_selected()
    #plot(selected_data[, "x2"])
    fit <-
      survfit(Surv(times, patient.vital_status) ~ expression, data = selected_data)
    
    selected_gene <- input$gene
    limit_time <-
      max(max(selected_data[selected_data[, "expression"] == "low", "times"]),
          max(selected_data[selected_data[, "expression"] == "high", "times"]))
    n_low <-
      length(selected_data[selected_data[, "expression"] == "low", "times"])
    n_high <-
      length(selected_data[selected_data[, "expression"] == "high", "times"])
    par(oma = c(0, 0, 2, 0))
    ggsurvplot(
      fit,
      # survfit object with calculated statistics.
      data = selected_data,
      risk.table = TRUE,
      # show risk table.
      pval = TRUE,
      # show p-value of log-rank test.
      conf.int = input$show_conf_int,
      # show confidence intervals for
      conf.int.style = input$conf_int_style,
      title = paste(
        "Cohort = ",
        input$cohort,
        ", gene = ",
        input$gene,
        ", n_low = ",
        n_low,
        ", n_high = ",
        n_high
      ),
      break.y.by = 0.1,
      xlim = c(0, limit_time),
      # present narrower X axis, but not affect
      xlab = "Time (years)",
      # survival estimates.
      #break.time.by = 1000,    # break X axis in time intervals by 500.
      palette = c("#9b9696", "#000000"),
      ggtheme = theme_bw(),
      # customize plot and risk table with a theme.
      risk.table.y.text.col = T,
      # colour risk table text annotations.
      risk.table.y.text = FALSE,
      # show bars instead of names in text annotations
      # in legend of risk table
      legend = "top",
      legend.title = NULL,
      surv.scale = "percent",
      break.x.by = 1,
      #risk.table.fontsize = 15,
      font.main = c(input$survival_font_main, "bold"),
      font.tickslab = input$survival_font_tickslab,
      font.x = input$survival_font_xlab,
      font.y = input$survival_font_ylab,
      font.legend = input$survival_font_legend
    )
  })
  cutpoint_plot <- function() {
    max_stat_cutpoint <- surv_cutpoint(
      survival_data[survival_data$dataset == input$cohort2,],
      time = "times",
      event = "patient.vital_status",
      variables = c(input$gene2)
    )
    
    plot(
      max_stat_cutpoint,
      xlab = paste(input$gene2, " mRNA expression level"),
      title = paste("Cohort: ", input$cohort2),
      palette = c("#f00909", "#000000"),
      ggtheme = theme_bw(),
      font.main = c(input$survival_font_main, "bold"),
      font.title = 20,
      font.tickslab = input$cutpoint_font_tickslab,
      font.x = input$cutpoint_font_xlab,
      font.y = input$cutpoint_font_ylab,
      font.legend = input$cutpoint_font_legend
    )
    
  }
  
  output$survival_curve <- renderPlot({
    survival_curve_plot()
  })
  
  output$cutpoint_plot <- renderPlot({
    cutpoint_plot()
  })
  
  
  
  output$survival_plot_png <- downloadHandler(
    filename = function() {
      paste(input$gene, input$cohort, '.png', sep = '')
    },
    content = function(file) {
      png(
        file,
        width = input$survival_width_px,
        height = input$survival_height_px
      )
      par(oma = c(0, 0, 4, 0))
      print(survival_curve_plot())
      dev.off()
    }
  )
  
  output$survival_plot_pdf <- downloadHandler(
    filename = function() {
      paste(input$gene, input$cohort, '.pdf', sep = '')
    },
    content = function(file) {
      cairo_pdf(
        file,
        width = input$survival_width_px / 100,
        height = input$survival_height_px / 100
      )
      print(survival_curve_plot())
      dev.off()
    }
  )
  
  output$survival_plot_eps <- downloadHandler(
    filename = function() {
      paste(input$gene, input$cohort, '.eps', sep = '')
    },
    content = function(file) {
      #setEPS()
      cairo_ps(
        file,
        width = input$survival_width_px / 100,
        height = input$survival_height_px / 100,
        fallback_resolution = 600
      )
      par(oma = c(0, 0, 4, 0))
      print(survival_curve_plot())
      dev.off()
    }
  )
  
  output$survival_plot_tiff <- downloadHandler(
    filename = function() {
      paste(input$gene, input$cohort, '.tiff', sep = '')
    },
    content = function(file) {
      tiff(
        file,
        width = input$survival_width_px,
        height = input$survival_height_px
      )
      par(oma = c(0, 0, 4, 0))
      print(survival_curve_plot())
      dev.off()
    }
  )
  
  output$cutpoint_plot_png <- downloadHandler(
    filename = function() {
      paste(input$gene, input$cohort, '.png', sep = '')
    },
    content = function(file) {
      png(
        file,
        width = input$cutpoint_width_px,
        height = input$cutpoint_height_px
      )
      par(oma = c(2, 2, 2, 2))
      print(cutpoint_plot())
      dev.off()
    }
  )
  
  output$cutpoint_plot_tiff <- downloadHandler(
    filename = function() {
      paste(input$gene, input$cohort, '.tiff', sep = '')
    },
    content = function(file) {
      tiff(
        file,
        width = input$cutpoint_width_px,
        height = input$cutpoint_height_px
      )
      par(oma = c(2, 2, 2, 2))
      print(cutpoint_plot())
      dev.off()
    }
  )
  
  output$cutpoint_plot_pdf <- downloadHandler(
    filename = function() {
      paste(input$gene, input$cohort, '.pdf', sep = '')
    },
    content = function(file) {
      cairo_pdf(
        file,
        width = input$cutpoint_width_px / 100,
        height = input$cutpoint_height_px / 100
      )
      print(cutpoint_plot())
      dev.off()
    }
  )
  
  output$cutpoint_plot_eps <- downloadHandler(
    filename = function() {
      paste(input$gene, input$cohort, '.eps', sep = '')
    },
    content = function(file) {
      cairo_ps(
        file,
        width = input$cutpoint_width_px / 100,
        height = input$cutpoint_height_px / 100,
        fallback_resolution = 600
      )
      #par(oma=c(0,0,4,0))
      print(cutpoint_plot())
      dev.off()
    }
  )
  
  # Generate a summary of the dataset
  output$summary <- renderPrint({
    print(input$gene)
  })
  
  heatmap_log_rank_reactive <- function() {
    genes_log_rank <- input$selected_genes_manually_log_rank
    if(input$all_genes_logrank)
      genes_log_rank <- selected_genes_all_restricted
    if (input$method == "median")
      m <- as.matrix(pvals_median[, genes_log_rank])
    else
      m <- as.matrix(pvals_optimal[, genes_log_rank])
    if (input$scale == "logarithmic")
      m <- log(m + 1)
    hm <- pheatmap(
      t(m),
      silent = TRUE,
      cellwidth = 20,
      cellheight = 20,
      fontsize = input$log_rank_heatmap_fontsize,
      fontsize_row = input$log_rank_heatmap_fontsize_row,
      fontsize_col =  input$log_rank_heatmap_fontsize_col,
      treeheight_row = 15,
      treeheight_col = 15,
      color = colorRampPalette(rev(
        brewer.pal(n = 7, name = input$heatmap_log_rank_scale)
      ))(10)
    )
    #grid::grid.newpage()
    grid::grid.draw(hm$gtable)
  }
  
  
  
  
  output$heatmap_log_rank_plot <- renderPlot({
    par(mar = c(5, 6, 4, 1) + .1)
    heatmap_log_rank_reactive()
  }, height = function() {
    length(input$selected_genes_manually_log_rank) * 20 + 200
  }, execOnResize = TRUE)
  
  output$download_heatmap_log_rank <- downloadHandler(
    filename = function() { 
      sprintf('%s.%s', "heatmap", input$select_file_type_log_rank_heatmap)
      # paste("heatmap.", input$select_file_type_expression_heatmap, sep = "")
    },
    content = function(file) {
      if (input$select_file_type_log_rank_heatmap == 'tiff') {
        device <- function(..., width, height) {
          grDevices::tiff(..., width = width, height = height, res = 600, units = "in")
        }
      } else if (input$select_file_type_log_rank_heatmap == 'pdf'){
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height, onefile = FALSE)
        }
      } else {
        device <- function(..., width, height) {
          postscript(..., width = width, height = height)
        }
      }
      ggsave(file, 
             plot = heatmap_log_rank_reactive(), 
             device = device(height = length(input$selected_genes_manually_log_rank)/2 + 2, width = 12))
    })
  
  output$heatmap_log_rank_plot_png <- downloadHandler(
    filename = function() {
      paste("heatmap", '.png', sep = '')
    },
    content = function(file) {
      heatmap_log_rank_file(file)
    }
  )
  
  output$heatmap_log_rank_plot_pdf <- downloadHandler(
    filename = function() {
      paste("heatmap", '.pdf', sep = '')
    },
    content = function(file) {
      heatmap_log_rank_file(file)
    }
  )
  
  output$heatmap_log_rank_plot_eps <- downloadHandler(
    filename = function() {
      paste("heatmap", '.eps', sep = '')
    },
    content = function(file) {
      cairo_ps(file, fallback_resolution = 600)
      par(oma = c(0, 0, 4, 0))
      print(heatmap_log_rank_reactive())
      dev.off()
    }
  )
  
  output$heatmap_log_rank_plot_tiff <- downloadHandler(
    filename = function() {
      paste("heatmap", '.tiff', sep = '')
    },
    content = function(file) {
      heatmap_log_rank_file(file)
    }
  )
  
  output$log_rank_table <- renderDataTable({
    if (input$log_rank_table_method == "median")
      test_results_median
    else
      test_results_optimal
  })
  
  output$log_rank_table_csv <- downloadHandler(
    filename = function() {
      paste(input$log_rank_table_method, 'log_rank_test.csv', sep = '')
    },
    content = function(file) {
      content = test_results_optimal
      if (input$log_rank_table_method == "median")
        content = test_results_median
      write.csv(content, file = file)
    }
  )
  
  output$log_rank_table_txt <- downloadHandler(
    filename = function() {
      paste(input$log_rank_table_method, 'log_rank_test.txt', sep = '')
    },
    content = function(file) {
      content = test_results_optimal
      if (input$log_rank_table_method == "median")
        content = test_results_median
      write.table(content, file = file)
    }
  )
  
  output$survival_data_table <- renderDataTable(survival_data)
  
  output$methylation_expression_correlation_table <-
    renderDataTable({
      methylation_expression_correlation()
    })
  
  output$methylation_expression_correlation_csv <- downloadHandler(
    filename = function() {
      paste('file.csv', sep = '')
    },
    content = function(file) {
      content = methylation_expression_correlation()
      write.csv(content, file = file)
    }
  )
  
  
  
  methylation_expression_correlation <- eventReactive(input$perform_methylation, {
    withProgress(message = "Calculating output. Please Wait...", value = 0,{
    gene <- input$gene6
    cohort <- input$cohort6
    k <- 10
    expressions_selected <-
      expressions_all[, c("bcr_patient_barcode", "dataset", gene)] %>% filter(dataset == cohort) %>% as.data.frame()
    expressions_selected[, gene] <-
      ntile(expressions_selected[, gene], k)
    
    # load files from prats
    files <-
      list.files(path = "data/methylation",
                 pattern = cohort,
                 full.names = TRUE)
    selected_methylation <- data.frame()
    for (file in files) {
      load(file)
      selected_methylation <-
        rbind(selected_methylation, methylation_part)
    }
    
    selected_methylation <-
      selected_methylation[,!is.na(selected_methylation[1,])]
    methylation_expression <-
      expressions_selected %>% filter(.[[3]] %in% c(1, 10)) %>% mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% inner_join(selected_methylation)
    cpg_islands <- colnames(selected_methylation)[-1]
    
    row_n <- 1
    pvals_list <- list()
    for (k in cpg_islands) {
      low_expr <-
        methylation_expression[methylation_expression[, 3] == 1 , k]
      high_expr <-
        methylation_expression[methylation_expression[, 3] == 10 , k]
      pvalue <- NA
      try({
        test <- t.test(low_expr, high_expr)
        pvalue <- test$p.value
      }, silent = TRUE)
      pvals_list[[row_n]] <- data.frame(
        "cohort" = cohort,
        "gene" = gene,
        "cpg_island" = k,
        "pvalue" = pvalue,
        "low_expr_meth_mean" = mean(low_expr, na.rm = TRUE),
        "high_expr_meth_mean" = mean(high_expr, na.rm = TRUE)
      )
      row_n <- row_n + 1
      incProgress(1/length(cpg_islands), detail = paste("Processing cpg island:", k))
    }
    suppressWarnings(pvals <- bind_rows(pvals_list))
    
    pvals$pvalue <- as.numeric(pvals$pvalue)
    pvals$pvalue_adjusted <- p.adjust(pvals$pvalue, method = "fdr")
    
    pvals_extended <- pvals %>% arrange(pvalue_adjusted) %>%
      mutate(log_odds_ratio =  log2((1 - high_expr_meth_mean) * low_expr_meth_mean /
                                      ((1 - low_expr_meth_mean) * high_expr_meth_mean)
      )) %>%
      merge(illumina_humanmethylation_27_data[, c(
        "Name",
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
        "Product"
      )],
      by.x = "cpg_island",
      by.y = "Name")
    
    pvals_extended_signif <-
      pvals_extended[pvals_extended$pvalue_adjusted < 0.05 &
                       !is.na(pvals_extended$pvalue_adjusted),] %>%
      filter(abs(high_expr_meth_mean - low_expr_meth_mean) > 0.25) %>%
      arrange(pvalue_adjusted)
    
    if (input$only_significant_results)
      pvals_extended_signif
    else
      pvals_extended
  })
    })
  
  normal_vs_cancer_isoforms_plot_reactive <- reactive({
    plot_data <- all_isoforms[all_isoforms$cohort == input$cohort3 &
                                all_isoforms$gene == input$gene3 &
                                all_isoforms$group %in% c("01", "11"), ] %>%
      mutate(group = sub(group, pattern = "01", replacement = "Cancer")) %>%
      mutate(group = sub(group, pattern = "11", replacement = "Normal"))
    #ylim1 = boxplot.stats(plot_data$value)$stats[c(1, 5)]
    plot_data$group <- as.factor(plot_data$group)
    levels(plot_data$group) = c("Cancer", "Normal")
    # if(input$boxplot_log_scale)
    #   plot_data$value <- log(plot_data$value + 1)
    group_colors <-
      c("Cancer" = "#eeeeee", "Normal" = "#777777")  #c('#eeeeee', '#777777')
    if(!input$boxplot_grey_scale)
    group_colors <- c("Cancer" = "#f00505", "Normal" = "#000000")
    
    p0 <- ggplot(plot_data, aes(isoform, value, fill = group)) +
      geom_boxplot() +
      #coord_cartesian(ylim = ylim1*4) +
      ylab(ifelse(input$boxplot_log_scale, "mRNA log(x+1) scale",  "mRNA")) +
      scale_fill_manual(values = group_colors) +
      theme(
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA),
        text = element_text(size = as.integer(input$select_isoform_boxplot_xlab_font))
      )
    
    p0 + scale.isoforms()
    
  })
  
  output$breaks.isoforms <- renderUI({
    data <- all_isoforms[all_isoforms$cohort == input$cohort3 &
                                all_isoforms$gene == input$gene3 &
                                all_isoforms$group %in% c("01", "11"), ] %>%
      mutate(group = sub(group, pattern = "01", replacement = "Cancer")) %>%
      mutate(group = sub(group, pattern = "11", replacement = "Normal"))
    b <- max(data$value)
    if(!input$boxplot_log_scale)
      breaks <- round(seq(input$scale.min, b, length.out = 5), digits = 0)
    else
      breaks <- round(10^(seq(0, log10(b), length.out = 5)), digits = 0)
    textInput('breaks.isoforms', 'Breaks', value = toString(breaks))
  })
  
  
  
  scale.isoforms <- reactive({
    
    if (is.null(input$breaks.isoforms) || '' == input$breaks.isoforms) {
      b <- waiver()
    } else {
      b <- as.numeric(as.list(strsplit(input$breaks.isoforms, ',')[[1]]))
    }
    
    
    
    if (!input$boxplot_log_scale) {
      m <- max(input$scale.min, 0)
      scale <- scale_y_continuous(limits = c(m, NA), breaks = b, labels = function(x) format(x, big.mark = "", scientific = F))
    } else {
      m <- 0 #max(input$scale.min, 0.01)
      b <- b[which(b >= m)]
      scale <- scale_y_continuous(limits = c(m, NA), breaks = b, labels = function(x) format(x, big.mark = "", scientific = F), trans = "log1p")
      # scale <- scale_y_log10(limits = c(NA, NA), breaks = b, labels = function(x) {
      #   c(format(x[which(x < 1)], nsmall = 2), format(x[which(x >= 1)], nsmall = 0))
      #})
    }
  })
  
  output$normal_vs_cancer_isoforms_plot <- renderPlot({
    normal_vs_cancer_isoforms_plot_reactive()
  })
  
  
  output$download.plot.normal.cancer <- downloadHandler(
    filename = function() { 
      tumors <- paste(input$cohort3, collapse = ' ')
      sprintf('%s-%s.%s', tumors, input$gene3, input$select.file.type.normal.cancer)
    },
    content = function(file) {
      if (input$select.file.type.normal.cancer == 'tiff') {
        device <- function(..., width, height) {
          grDevices::tiff(..., width = width, height = height, res = 300, units = "in")
        }
      } else if (input$select.file.type.normal.cancer == 'pdf'){
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height)
        }
      } else {
        device <- function(..., width, height) {
          postscript(..., width = width, height = height)
        }
      }
      ggsave(file, plot = normal_vs_cancer_isoforms_plot_reactive(), device = device)
    })
  
  normal_vs_cancer_isoforms_table_reactive <-
    eventReactive(
      input$perform_isoform_table, {
      withProgress(message = "Calculating output. Please Wait...", value = 0, {
        isoform_name_mapping <- read.table(file = "data/knownToRefSeq.txt")
        colnames(isoform_name_mapping) <- c("isoform", "GenBank_name")
        isoform_name_mapping$isoform_cut <-
          substr(isoform_name_mapping$isoform, 1, 8)
        isoform_summary_table_tmp <-
          all_isoforms %>% filter(cohort == input$cohort_isoform_table & gene == input$gene_isoform_table) %>%
          filter(group %in% c("01", "11")) %>%
          mutate(group = sub(x = group, "01", "cancer tissue")) %>%
          mutate(group = sub(x = group, "11", "normal tissue")) %>%
          group_by(cohort, gene, isoform, group) %>%
          dplyr::summarize(
            mean_expression = mean(value),
            se = sd(value) / sqrt(n()),
            median_expression = median(value),
            `25%` = quantile(value, probs = 0.25, na.rm = TRUE),
            `50%` = quantile(value, probs = 0.5),
            `75%` = quantile(value, probs = 0.75),
            n_onservations = n()
          ) %>%
          mutate(isoform_cut = substr(isoform, 1, 8))
        counter <- 1
        test_cancer_vs_normal_list <- list()
          j <- input$gene_isoform_table
          i <- input$cohort_isoform_table
          genes <- unique(all_isoforms$gene)
          
          isoforms <-
            unique(all_isoforms[all_isoforms$cohort == i & all_isoforms$gene == j, "isoform"])
          normal_group <-
            all_isoforms %>% filter(cohort == i,
                                    gene == j,
                                    substr(bcr_patient_barcode, 14, 15) == "11")
          cancer_group <-
            all_isoforms %>% filter(cohort == i,
                                    gene == j,
                                    substr(bcr_patient_barcode, 14, 15) == "01")
          
          for (k in isoforms) {
            p_value <- NA
            test_cancer_vs_normal_list[[counter]] <-
              data.frame(i, j, k, p_value)
            try({
              normal_group_isoform <- normal_group %>% filter(isoform == k) 
              cancer_group_isoform <- cancer_group %>% filter(isoform == k)   
              p_value <- t.test(normal_group_isoform$value, cancer_group_isoform$value)$p.value
              test_cancer_vs_normal_list[[counter]][1, 4] <- p_value
            }, silent = TRUE)
            counter <- counter + 1
            
          }
          
        
        
        test_cancer_vs_normal <- bind_rows(test_cancer_vs_normal_list)
        colnames(test_cancer_vs_normal) <- c("cohort", "gene", "isoform", "pvalue")
        test_cancer_vs_normal %>% filter(p_value != "NA") %>%
          left_join(isoform_summary_table_tmp) %>%
          left_join(isoform_name_mapping[, -1])
      })
    })
  
  output$normal_vs_cancer_isoforms_table <- renderDataTable({
    normal_vs_cancer_isoforms_table_reactive()
  }, options = list(pageLength = 10))
  
  output$isoform_test_all_csv <- downloadHandler(
    filename = function() {
      paste('normal_vs_caner_tests.csv', sep = '')
    },
    content = function(file) {
      write.csv(normal_vs_cancer_isoforms_table_reactive(), file = file)
    }
  )
  
  output$isoform_test_all_txt <- downloadHandler(
    filename = function() {
      paste('normal_vs_caner_tests.txt', sep = '')
    },
    content = function(file) {
      write.table(normal_vs_cancer_isoforms_table_reactive(), file = file)
    }
  )
  
  isoform_expressions_plot_reactive <- reactive({
    all_isoforms %>% filter(group == "01" ) %>%
      filter(cohort == input$cohort4, gene == input$gene4) %>%
      #filter(cohort == "BRCA", gene == "ZNF695") %>%
      group_by(cohort, gene, isoform) %>%
      dplyr::summarise(sum_of_expression = sum(value)) -> a
    
    all_isoforms %>% filter(group == "01") %>%
      filter(cohort == input$cohort4, gene == input$gene4) %>%
      group_by(cohort, gene) %>%
      dplyr::summarise(total_expression = sum(value)) -> b
    
    a %>% left_join(b) %>%
      mutate(expr_percentage = 100 * sum_of_expression / total_expression) %>%
      #filter(cohort == input$cohort4, gene == input$gene4) %>%
      ggplot(aes(reorder(isoform, expr_percentage), expr_percentage)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      xlab("Isoforms") +
      ylab("Expression percentage") +
      ggtitle(paste(input$gene4, " in ", input$cohort4, sep = "")) + 
      theme(panel.background = element_rect(fill = 'white'),
            panel.border = element_rect(colour = 'black', fill = NA),
            text = element_text(size = as.integer(input$select_isoform_expression_font)))
      
  })
  
  output$isoform_expressions_plot <- renderPlot({
    isoform_expressions_plot_reactive()
  })
  
  output$download_isoform_expression_plot <- downloadHandler(
    filename = function() { 
      tumors <- paste(input$cohort4, collapse = ' ')
      sprintf('%s-%s.%s', tumors, input$gene4, input$select.file.type.isoform.expression)
    },
    content = function(file) {
      if (input$select.file.type.isoform.expression == 'tiff') {
        device <- function(..., width, height) {
          grDevices::tiff(..., width = width, height = height, res = 300, units = "in")
        }
      } else if (input$select.file.type.isoform.expression == 'pdf'){
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height)
        }
      } else {
        device <- function(..., width, height) {
          postscript(..., width = width, height = height)
        }
      }
      ggsave(file, plot = isoform_expressions_plot_reactive(), device = device)
    })
  
  
  
  all_cohorts_isoform_expressions_plot_reactive <- reactive({
    all_isoforms %>% filter(group == "01") %>%
      group_by(cohort, gene, isoform) %>%
      dplyr::summarise(sum_of_expression = sum(value)) -> a
    
    all_isoforms %>% filter(group == "01") %>%
      group_by(cohort, gene) %>%
      dplyr::summarise(total_expression = sum(value)) -> b
    
    a %>% left_join(b) %>%
      mutate(expr_percentage = 100 * sum_of_expression / total_expression) %>%
      filter(gene == input$gene5) %>%
      ggplot(aes(cohort, expr_percentage)) +
      geom_bar(aes(fill = isoform), stat = "identity") +
      theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
      ggtitle(paste(input$gene5)) + 
      ylab("Expression percentage") +
      theme(panel.background = element_rect(fill = 'white'),
            panel.border = element_rect(colour = 'black', fill = NA),
            text = element_text(size = as.integer(input$select_all_cohorts_isoform_expression_font)))
  })
  
  output$all_cohorts_isoform_expressions_plot <- renderPlot({
    all_cohorts_isoform_expressions_plot_reactive()
  })
  
  output$download_all_cohorts_isoform_expression_plot <- downloadHandler(
    filename = function() { 
      sprintf('%s.%s', input$gene5, input$select.file.type.all.cohorts.isoform.expression)
    },
    content = function(file) {
      if (input$select.file.type.all.cohorts.isoform.expression == 'tiff') {
        device <- function(..., width, height) {
          grDevices::tiff(..., width = width, height = height, res = 300, units = "in")
        }
      } else if (input$select.file.type.all.cohorts.isoform.expression == 'pdf'){
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height)
        }
      } else {
        device <- function(..., width, height) {
          postscript(..., width = width, height = height)
        }
      }
      ggsave(file, plot = all_cohorts_isoform_expressions_plot_reactive(), device = device)
    })
  
  output$all_cohorts_isoform_expressions_plot_png <-
    downloadHandler(
      filename = function() {
        paste("all_cohorts_isoform_expressions_plot", '.png', sep = '')
      },
      content = function(file) {
        png(file)
        par(oma = c(2, 2, 2, 2))
        print(all_cohorts_isoform_expressions_plot_reactive())
        dev.off()
      }
    )
  
  output$all_cohorts_isoform_expressions_plot_pdf <-
    downloadHandler(
      filename = function() {
        paste("all_cohorts_isoform_expressions_plot_plot",
              '.pdf',
              sep = '')
      },
      content = function(file) {
        cairo_pdf(file)
        print(all_cohorts_isoform_expressions_plot_reactive())
        dev.off()
      }
    )
  
  
  expression_heatmap1_reactive <- function() {
    m <- list(
      l = 65,
      r = 40,
      b = 200,
      t = 200,
      pad = 4
    )
    genes <- input$selected_genes_manually
    if(input$all_genes_normal)
      genes <- selected_genes_all
    m <- as.matrix(t(expressions_means))
    if (input$log_scale)
      m <- log(m + 1)
    m <- m[genes,]
    
    p <- pheatmap(
      m,
      cellwidth = 20,
      cellheight = 20,
      fontsize = input$normal_heatmap_fontsize,
      fontsize_row = input$normal_heatmap_fontsize_row,
      fontsize_col =  input$normal_heatmap_fontsize_col,
      margin = m,
      silent = TRUE,
      color = colorRampPalette(rev(
        brewer.pal(n = 7, name = input$heatmap_color_scale1)
      ))(10)
    )
    #grid::grid.newpage()
    grid::grid.draw(p$gtable)
    
  }
  
  output$expression_heatmap1 <- renderPlot({
    expression_heatmap1_reactive()
  }, height = function() {
    if(input$all_genes_normal) 
      length(selected_genes_all) * 20 + 300
    else
      length(input$selected_genes_manually) * 20 + 300})
  
  output$download_expression_heatmap1 <- downloadHandler(
    filename = function() { 
     sprintf('%s.%s', "heatmap", input$select_file_type_expression_heatmap)
     # paste("heatmap.", input$select_file_type_expression_heatmap, sep = "")
    },
    content = function(file) {
      if (input$select_file_type_expression_heatmap == 'tiff') {
        device <- function(..., width, height) {
          grDevices::tiff(..., width = width, height = height, res = 600, units = "in")
        }
      } else if (input$select_file_type_expression_heatmap == 'pdf'){
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height, onefile = FALSE)
        }
      } else {
        device <- function(..., width, height) {
          postscript(..., width = width, height = height)
        }
      }
      ggsave(file, 
             plot = expression_heatmap1_reactive(), 
             device = device(height = length(input$selected_genes_manually)/2 + 2, width = length(input$selected_genes_manually)/2 + 2))
    })
  
  output$expression_heatmap2 <- renderPlotly({
    if (input$restrict2) {
      if (input$log_scale2)
        z <-
          log(t(expressions_means[, selected_group]) + 1) %>% as.matrix()
      else
        z <- t(expressions_means[, selected_group]) %>% as.matrix()
      height <- 400
    }
    else {
      if (input$log_scale2)
        z <-
          log(t(expressions_means[, input$selected_genes_manually2] + 1)) %>% as.matrix()
      else
        z <-
          t(expressions_means[, input$selected_genes_manually2]) %>% as.matrix()
      height <- length(input$selected_genes_manually2) * 22 + 10
    }
    m <- list(
      l = 65,
      r = 40,
      b = 100,
      t = 50,
      pad = 4
    )
    plot_ly(
      z = z,
      x = colnames(z),
      y = rownames(z),
      type = "heatmap",
      height = height,
      colors = input$heatmap_color_scale2
    ) %>% layout(margin = m)
    # color_palette <- colorRampPalette(brewer.pal(nrow(z), input$heatmap_color_scale))
    # heatmaply(x = z,
    #           xlab = colnames(z),
    #           ylab = rownames(z),
    #           colors = color_palette,
    #           margins = c(60,100,40,20),
    #           subplot_heights = c(height, height)
    #           )
    
  })
  
  normal_expression_pvalues_table <- reactive({
    expressions_normal_long <- expressions_normal %>% gather(gene, expression, HKR1:ZNF99, factor_key=TRUE)
    expressions_normal_long$dataset <- expressions_normal_long$dataset %>%
      tolower() %>%
      gsub(pattern = "acc", replacement = "Adenoid cystic carcinoma (ACC)") %>%
      gsub(pattern = "blca", replacement = "Bladder (BLCA)") %>%
      gsub(pattern = "brca", replacement = "Breast (BRCA)") %>%
      gsub(pattern = "cesc", replacement = "Cervix (CESC)") %>%
      gsub(pattern = "chol", replacement = "Bile duct (CHOL)") %>%
      #gsub(pattern = "coadread", replacement = "remove") %>%
      gsub(pattern = "coad", replacement = "Colon (COAD)") %>%
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
      gsub(pattern = "pcpg", replacement = "Neuroendocrine glands (PCPG)") %>%
      gsub(pattern = "prad", replacement = "Prostate (PRAD)") %>%
      gsub(pattern = "read", replacement = "Rectum (READ)") %>%
      ## gsub(pattern = "sarc", replacement = "Pheochromocytoma and Paraganglioma (PCPG)") %>%
      gsub(pattern = "sarc", replacement = "Bone (SARC)") %>%
      gsub(pattern = "skcm", replacement = "Skin (SKCM)") %>%
      gsub(pattern = "stad", replacement = "Stomach (STAD)") %>%
      ## gsub(pattern = "stes", replacement = "Stomach (STAD)") %>%
      gsub(pattern = "tgct", replacement = "Testicles (TGCT)") %>%
      gsub(pattern = "thca", replacement = "Thyroid (THCA)") %>%
      gsub(pattern = "thym", replacement = "Thymus (THYM)") %>%
      gsub(pattern = "ucec", replacement = "Uterus (UCEC, UCS)") %>%
      gsub(pattern = "ucs", replacement = "Uterus (UCEC, UCS)") %>%
      gsub(pattern = "uvm", replacement = "Uvea (UVM)")
      
    
    expressions_normal_summary <- expressions_normal_long %>%
      filter(dataset != "STES") %>%
      group_by(dataset, gene) %>% 
      dplyr::summarize(
        mean_expression = mean(expression),
        se = sd(expression) / sqrt(n()),
        median_expression = median(expression),
        `25%` = quantile(expression, probs = 0.25, na.rm = TRUE),
        `50%` = quantile(expression, probs = 0.5),
        `75%` = quantile(expression, probs = 0.75),
        n_onservations = n()
      ) 
    
  })
  
  output$normal_expression_pvalues_table_render <- renderDataTable({
    normal_expression_pvalues_table()
  })
  
  output$normal_expression_pvalues_table_download <- downloadHandler(
    filename = function() { 
      sprintf('%s.%s', "pvalues", input$select.normal.table.file.type)
    },
    content = function(file) {
      if(input$select.table.file.type == "csv")
        write.csv(x = normal_expression_pvalues_table(), file = file, sep = ",")
      if(input$select.table.file.type == "txt")
        write.table(x = normal_expression_pvalues_table(), file = file)
    })
  
  normal_expression_boxplot_reactive <- reactive({
    selected_gene_expressions <- expressions_normal %>% 
      select("dataset", paste(input$gene_normal_expression)) %>%
      filter(dataset != "COADREAD") %>%
      filter(dataset != "STES")
    selected_gene_expressions$dataset <- selected_gene_expressions$dataset %>%
      tolower() %>%
      gsub(pattern = "acc", replacement = "Adenoid cystic carcinoma (ACC)") %>%
      gsub(pattern = "blca", replacement = "Bladder (BLCA)") %>%
      gsub(pattern = "brca", replacement = "Breast (BRCA)") %>%
      gsub(pattern = "cesc", replacement = "Cervix (CESC)") %>%
      gsub(pattern = "chol", replacement = "Bile duct (CHOL)") %>%
      #gsub(pattern = "coadread", replacement = "remove") %>%
      gsub(pattern = "coad", replacement = "Colon (COAD)") %>%
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
      gsub(pattern = "pcpg", replacement = "Neuroendocrine glands (PCPG)") %>%
      gsub(pattern = "prad", replacement = "Prostate (PRAD)") %>%
      gsub(pattern = "read", replacement = "Rectum (READ)") %>%
      ## gsub(pattern = "sarc", replacement = "Pheochromocytoma and Paraganglioma (PCPG)") %>%
      gsub(pattern = "sarc", replacement = "Bone (SARC)") %>%
      gsub(pattern = "skcm", replacement = "Skin (SKCM)") %>%
      gsub(pattern = "stad", replacement = "Stomach (STAD)") %>%
      ## gsub(pattern = "stes", replacement = "Stomach (STAD)") %>%
      gsub(pattern = "tgct", replacement = "Testicles (TGCT)") %>%
      gsub(pattern = "thca", replacement = "Thyroid (THCA)") %>%
      gsub(pattern = "thym", replacement = "Thymus (THYM)") %>%
      gsub(pattern = "ucec", replacement = "Uterus (UCEC, UCS)") %>%
      gsub(pattern = "ucs", replacement = "Uterus (UCEC, UCS)") %>%
      gsub(pattern = "uvm", replacement = "Uvea (UVM)")
    
    
    #ylim1 = boxplot.stats(selected_gene_expressions[, 2]) $stats[c(1, 5)]
    ggplot(selected_gene_expressions, aes(dataset, log(selected_gene_expressions[, 2] + 1))) +
      geom_boxplot() +
      #coord_cartesian(ylim = ylim1*4) +
      ylab("log(expression + 1)") +
      scale_fill_manual(values = c("#000000", "#000000")) +
      theme(panel.background = element_rect(fill = 'white'),
            panel.border = element_rect(colour = 'black', fill = NA)) + 
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = input$select_normal_boxplot_font),
            axis.text.y = element_text(size = input$select_normal_boxplot_font))
      
  })
  
  output$normal_expression_boxplot <- renderPlot({
    normal_expression_boxplot_reactive()
  })
  
  output$normal_expression_boxplot_download = downloadHandler(
    filename = function() { 
      sprintf('%s.%s', input$gene_normal_expression, input$select_file_type_expression_boxplot)
    },
    content = function(file) {
      if (input$select_file_type_expression_boxplot == 'tiff') {
        device <- function(..., width, height) {
          grDevices::tiff(..., width = width, height = height, res = 300, units = "in")
        }
      } else if (input$select_file_type_expression_boxplot == 'pdf'){
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height)
        }
      } else {
        device <- function(..., width, height) {
          postscript(..., width = width, height = height)
        }
      }
      ggsave(file, plot = normal_expression_boxplot_reactive(), device = device)
    })
  
  ## Kornel's part #####
  source('data/krabmen/R/read_data.R')
  
  ## Boxplots ##########
  load.datasets()
  
  scale <- reactive({
    
    if (is.null(input$breaks) || '' == input$breaks) {
      b <- waiver()
    } else {
      b <- as.numeric(as.list(strsplit(input$breaks, ',')[[1]]))
    }
    
    
    
    if (input$select.scale == 'linear') {
      m <- max(input$scale.min, 0)
      scale <- scale_y_continuous(limits = c(m, NA), breaks = b, labels = function(x) format(x, big.mark = "", scientific = F))
    } else {
      m <- max(input$scale.min, 0.01)
      b <- b[which(b >= m)]
      scale <- scale_y_log10(limits = c(m, NA), breaks = b, labels = function(x) {
        c(format(x[which(x < 1)], nsmall = 2), format(x[which(x >= 1)], nsmall = 0))
      })
    }
  })
  
  output$breaks <- renderUI({
    data <- filtered.data()
    b <- max(data$expression)
    if(input$select.scale == "linear")
      breaks <- round(seq(input$scale.min, b, length.out = 5), digits = 0)
    else
      breaks <- round(10^(seq(0, log10(b), length.out = 5)), digits = 0)
    textInput('breaks', 'Breaks', value = toString(breaks))
  })
  
  output$value <- renderPrint({ input$select.tumor })
  
  pvalues.reactive <- function(){
    vals <- pvalues2[pvalues2$gene %in% input$select.gene & pvalues2$tumor %in% input$select.tumor,c('padj', 'tumor', 'gene')]
    vals$padj <- sprintf("%.4f", vals$padj)
    
    vals %>% left_join(summary.all.tumors.long())
  }
  
  output$pvalues <- renderDataTable({ 
   pvalues.reactive()
  })
  
  summary.all.tumors.long <- reactive({
    all.tumors.long2 %>% 
      filter(tumor %in% input$select.tumor) %>%
      filter(gene %in% input$select.gene) %>%
      group_by(tumor, type, gene) %>% 
      dplyr::summarize(
      mean_expression = mean(expression),
      se = sd(expression) / sqrt(n()),
      median_expression = median(expression),
      `25%` = quantile(expression, probs = 0.25, na.rm = TRUE),
      `50%` = quantile(expression, probs = 0.5),
      `75%` = quantile(expression, probs = 0.75),
      n_onservations = n()
    ) 
  })
  
  output$pvalues.download <- downloadHandler(
    filename = function() { 
            sprintf('%s.%s', "pvalues", input$select.table.file.type)
    },
    content = function(file) {
      if(input$select.table.file.type == "csv")
        write.csv(x = pvalues.reactive(), file = file, sep = ",")
      if(input$select.table.file.type == "txt")
        write.table(x = pvalues.reactive(), file = file)
    })
  
  
  filtered.data <- reactive({
    filtered <- all.tumors.long2[all.tumors.long2$tumor %in% input$select.tumor, ]
    filtered[filtered$gene %in% input$select.gene, ]
  })
  
  
  plotGeneExpression <- function() {
    
    if (input$greyscale) {
      colors <- scale_fill_manual(values = c('#eeeeee', '#777777'))
    } else {
      colors <- scale_fill_manual(values = c('#aaaaaa', '#ff0000'))
    }
    
    ggplot(filtered.data(), aes(y = expression,
                                x = tumor, 
                                fill = factor(type, labels = c('Healthy', 'Tumor')))) + 
      geom_boxplot(alpha = 1, position = position_dodge(0.8)) +
      colors +
      ggtitle('Boxplots for selected tumors') +
      ylab(paste0(input$select.gene, ' expression')) +
      theme(legend.title=element_blank()) + 
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0))) +
      theme(axis.text = element_text(colour = 'black')) + 
      theme(plot.title = element_text(margin = margin(0, 0, 20, 0))) + 
      theme(text = element_text(size = input$input.font.size)) +
      theme(panel.background = element_rect(fill = 'white'),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = 'black', fill = NA)) + 
      scale() + facet_wrap(~gene)
  }
  
  output$distPlot <- renderPlot({
    plotGeneExpression()
  })
  
  output$download.plot = downloadHandler(
    filename = function() { 
      tumors <- paste(input$select.tumor, collapse = ' ')
      sprintf('%s-%s.%s', tumors, input$select.gene, input$select.file.type)
    },
    content = function(file) {
      if (input$select.file.type == 'tiff') {
        device <- function(..., width, height) {
          grDevices::tiff(..., width = width, height = height, res = 300, units = "in")
        }
      } else if (input$select.file.type == 'pdf'){
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height)
        }
      } else {
        device <- function(..., width, height) {
          postscript(..., width = width, height = height)
        }
      }
      ggsave(file, plot = plotGeneExpression(), device = device)
    })
  
  ## Boxplots subtypes##
  tumor.data <- .tumor.data()
  
  # output$breaks.subtypes <- renderUI({
  #   data <- filtered.data.subtypes()
  #   b <- max(data$expression)
  #   if(input$select.scale.subtypes == "linear")
  #     breaks <- round(seq(input$scale.min.subtypes, b, length.out = 5), digits = 0)
  #   else
  #     breaks <- round(10^(seq(0, log10(b), length.out = 5)), digits = 0)
  #   textInput('breaks.subtypes', 'Breaks', value = toString(breaks))
  # })
  
  breaks.subtypes <- function(){
    data2 <- filtered.data.subtypes()

    b <- max(na.omit(data2$expression))
    if(input$select.scale.subtypes == "linear")
      breaks <- round(seq(input$scale.min.subtypes, b, length.out = 5), digits = 0)
    else
      breaks <- round(10^(seq(0, log10(b), length.out = 5)), digits = 0)
  }
  
  scale.subtypes <- reactive({
    
    # if (is.null(input$breaks.subtypes) || '' == input$breaks.subtypes) {
    #   b <- waiver()
    # } else {
    #   b <- as.numeric(as.list(strsplit(input$breaks.subtypes, ',')[[1]]))
    # }
    b <- breaks.subtypes()
    
    
    
    if (input$select.scale.subtypes == 'linear') {
      m <- max(input$scale.min.subtypes, 0)
      scale <- scale_y_continuous(limits = c(m, NA), breaks = b, labels = function(x) format(x, big.mark = "", scientific = F))
    } else {
      m <- max(input$scale.min.subtypes, 0.01)
      b <- b[which(b >= m)]
      scale <- scale_y_log10(limits = c(m, NA), breaks = b, labels = function(x) {
        c(format(x[which(x < 1)], nsmall = 2), format(x[which(x >= 1)], nsmall = 0))
      })
    }
  })
  
  generatePlot.subtypes <- function(){
    if (is.null(input$breaks.subtypes) || '' == input$breaks.subtypes) {
      b <- waiver()
    } else {
      b <- as.numeric(as.list(strsplit(input$breaks.subtypes, ',')[[1]]))
    }
    
    if (input$select.scale.subtypes == 'linear') {
      m <- max(input$scale.min.subtypes, 0)
      scale <- scale_y_continuous(limits = c(m, NA), breaks = b, labels = function(x) format(x, big.mark = "", scientific = F))
    } else {
      # m <- max(input$scale.min.subtypes, 0.01)
      # b <- b[which(b >= m)]
      #
      # if (max(filtered.data.subtypes()$expression) < max(b)) {
      #   limit <- max(b)
      # } else {
      #   limit <- NA
      # }
      #
      # scale <- scale_y_log10(limits = c(m, limit), breaks = b, labels = function(x) {
      #   c(format(x[which(x < 1)], nsmall = 2), format(x[which(x >= 1)], nsmall = 0))
      # })
      scale <- scale_y_log10()
    }
    scale <- scale.subtypes()
    if (input$rotate.x.axis.subtypes == T) {
      rotate <- theme(axis.text.x = element_text(angle = 90, hjust = 1))  
    } else {
      rotate <- theme()
    }
    
    data <- filtered.data.subtypes()
    data[, c(input$subtypes)] <- gsub(x = data[, c(input$subtypes)], " ", " \n")
    
    ggplot(data, aes_string(y = 'expression', x = input$subtypes)) + 
      geom_boxplot(alpha = 1, position = position_dodge(0.8)) + 
      ggtitle(input$subtypes) +
      ylab(paste0(input$gene.names.subtypes, ' expression')) +
      theme(legend.title = element_blank()) + 
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0))) +
      theme(axis.text = element_text(colour = 'black')) + 
      theme(plot.title = element_text(margin = margin(0, 0, 20, 0))) + 
      theme(text = element_text(size = input$input.font.size.subtypes)) + 
      rotate +
      scale + 
      theme(panel.background = element_rect(fill = 'white'),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = 'black', fill = NA)) + 
      facet_wrap(~gene) 
    #theme(text = element_text(size=20))
  }
  
  renderPvaluesTable.subtypes <- reactive({
    path <- file.path('data/krabmen/boxplots-subtypes/resources', input$tumor.name.subtypes, paste0(input$subtypes, '.RData'))
    
    df <- local(get(load(path)))
    path_all <- file.path('data/krabmen/boxplots-subtypes/resources', paste0(input$tumor.name.subtypes, '.RData'))
    summarised <- local(get(load(path_all))) %>% 
      select(tumor, gene, expression, input$subtypes) %>%
      mutate(subtype = .[[4]]) %>%
      group_by(tumor, gene, subtype) %>%
      dplyr::summarize(
        mean_expression = mean(expression),
        se = sd(expression) / sqrt(n()),
        median_expression = median(expression),
        `25%` = quantile(expression, probs = 0.25, na.rm = TRUE),
        `50%` = quantile(expression, probs = 0.5),
        `75%` = quantile(expression, probs = 0.75),
        n_onservations = n()
      ) 
      
    
    df[df$gene %in% input$gene.names.subtypes,] %>% cbind(input$tumor.name.subtypes) %>% left_join(summarised)
  })
  
  output$subtypes <- renderUI({
    
    df <- tumor.data[[input$tumor.name.subtypes]]
    
    choices <- setdiff(colnames(df), c('gene', 'expression', 'tumor', 'overall_survival', 'disease_status'))
    if(input$tumor.name.subtypes== "ESCA")
      choices <- setdiff(choices, c("clinical_M"))
    if(input$tumor.name.subtypes== "BLCA")
      choices <- setdiff(choices, c("clinical_T"))
    if(input$tumor.name.subtypes== "KIRP")
      choices <- setdiff(choices, c( "pathologic_N2", "clinical_M"))
    if(input$tumor.name.subtypes== "LIHC")
      choices <- setdiff(choices, c( "Hypermethylation.cluster.laird.groups.", 
                                     "Hypermethylation.cluster.laird.groups.1", 
                                     "Hypomethylation.cluster.laird.groups.", 
                                     "Hypomethylation.cluster.laird.groups.1"))
    if(input$tumor.name.subtypes== "THCA")
      choices <- setdiff(choices, c("meth_Cluster_number", "mRNA_Cluster_number_Feb2014", "meth_Cluster_number_Feb2014"))
    if(input$tumor.name.subtypes== "UCEC")
      choices <- setdiff(choices, c("IntegrativeCluster", "mrna_expression_cluster"))
    
    selectInput('subtypes', label = 'Subtypes', choices = levels(factor(choices)), selected = choices[1])
  })
  
  
  
  output$subtype.values <- renderUI({
    if (length(input$subtypes) > 0) {
      df <- tumor.data[[input$tumor.name.subtypes]]
      
      # TODO: fires on change of input$tumor.name before input$subtypes is refreshed
      if (input$subtypes %in% colnames(df)) {
        choices <- setdiff(unique(df[,c(input$subtypes)]), 
                           c("NA", "Unknown",  "Not evaluated", "Performed but not available", "not performer", "indeterminate, equivocal", "N/A", "Notassigned",
                             "Not Available"))
        
        
        selectInput('subtype.values', label = 'Subtype values', choices = levels(factor(choices)), selected = choices, multiple = T)
      }
    }
    
  })
  
  filtered.data.subtypes <- reactive({
    df <- tumor.data[[input$tumor.name.subtypes]]
    
    # Selected subtype 
    df <- df[, c('expression', 'gene', input$subtypes)]
    
    # Filter genes
    df <- df[df$gene %in% input$gene.names.subtypes, ]
    
    # Filter subtype values
    df <- df[df[, c(input$subtypes)] %in% input$subtype.values, ]
    
    df[, c(input$subtypes)] <- ordered(df[, c(input$subtypes)], levels = input$subtype.values)
    
    df
  })
  
  output$apply.subtypes <- renderUI({
    actionButton('apply.subtypes', 'Apply')
  })
  
  output$pvalues.subtypes <- renderDataTable({
    renderPvaluesTable.subtypes()
  })
  
  output$distPlot.subtypes <- renderPlot({
    generatePlot.subtypes()
  })
  
  output$plot.subtypes.download = downloadHandler(
    filename = function() {
      sprintf(
        'plot.%s',
        # input$tumor.name.subtypes,
        # input$subtypes,
        input$select.file.type.subtypes
      )
    },
    content = function(file) {
      if (input$select.file.type.subtypes == 'tiff') {
        device <- function(..., width, height) {
          grDevices::tiff(
            ...,
            width = width,
            height = height,
            res = 300,
            units = "in"
          )
        }
      } else if (input$select.file.type.subtypes == 'pdf') {
        device <- function(..., width, height) {
          grDevices::pdf(..., width, height)
        }
      } else {
        device <- function(..., width, height) {
          postscript(..., width = width, height = height)
        }
      }
      
      ggsave(
        file,
        plot = generatePlot.subtypes(),
        # width = input$plot.width.subtypes,
        # height = input$plot.height.subtypes,
        device = device,
        units = 'in'
      )
    }
  )
  
  output$table.subtypes.download <- downloadHandler(
    filename = function() { 
      sprintf('%s.%s', "pvalues", input$select.table.subtypes.file.type)
    },
    content = function(file) {
      if(input$select.table.subtypes.file.type == "csv")
        write.csv(x =renderPvaluesTable.subtypes(), file = file, sep = ",")
      if(input$select.table.subtypes.file.type == "txt")
        write.table(x = renderPvaluesTable.subtypes(), file = file)
    })
  
  
  output$clinical.parameters.description <- downloadHandler(
    filename = "clinical_parameters.pdf",
    content = function(file) {
      file.copy("data/clinical_parameters.pdf", file)
    }
  )
  
  ######################
  
}) 