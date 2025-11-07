#!/usr/bin/env Rscript

# Copyright (C) 2025 Rishvanth Prabakar
#
# Authors: Rish Prabakar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

suppressMessages(library("viridis"))
suppressMessages(library("optparse"))
suppressMessages(library("pheatmap"))
suppressMessages(library("tidyverse"))

main <- function() {

  parser <- OptionParser()
  parser <- add_option(parser, c("-t", "--tssNormFile"),
              help = "single cell flagstat QC file")
  parser <- add_option(parser, c("-c", "--clusterIdFile"),
              help="single cell cluster ID TSV file")
  parser <- add_option(parser, c("-o", "--outPrefix"),
              help = "Outfile prefix")
  opt <- parse_args(parser)

  if (is.null(opt$tssNormFile) | is.null(opt$outPrefix)) {
    print_help(parser)
    quit(status = 1)
  }


  # read data
  tss <- read_tsv(opt$tssNormFile)
  
  # remove unwanted columns
  tss.plot <- tss %>%
    select(!barcode)
  tss.plot <- as.matrix(tss.plot)

  ## TSS enrichment score is the max of the normalized TSS
  tss.ess <- tibble(barcode = tss$barcode, 
    ess = apply(tss[, 2:ncol(tss)], 1, max))

  # heatmap for all cells
  col.breaks <- seq(quantile(tss.plot, 0.01),
    quantile(tss.plot, 0.99), length.out=101)

  pdf(sprintf("%s_bc_tss_enrichment.pdf", opt$outPrefix), 
    height = 6, width = 8)
  pheatmap(tss.plot[order(tss.ess$ess, decreasing = TRUE),], 
    cluster_cols = FALSE, cluster_rows = FALSE, color = viridis(101),
    show_rownames = FALSE, show_colnames = FALSE, breaks = col.breaks)
  dev.off()


  ## TSS enrichment score is the max of the normalized TSS
  
  ggplot(data = tss.ess) +
    geom_histogram(mapping = aes(ess), binwidth = 0.01) +
    labs(x = "TSS enrichment score", y = "Number of cells",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()
  ggsave(sprintf("%s_bc_ess.pdf", opt$outPrefix), 
    height = 5, width = 5)



  # heatmap with cluster ids
  if (!is.null(opt$clusterIdFile)) {
   
    # read the cluster IDs
    id <- read_tsv(opt$clusterIdFile, col_names = FALSE)
    names(id) <- c("barcode", "cluster")
    id$cluster <- as.factor(id$cluster)  
  
    # keep only the cells that have a cluster ID
    tss <- tss[tss$barcode %in% id$barcode, ]
  
    tss.ess <- tss.ess[tss.ess$barcode %in% id$barcode, ]
    tss.ess <- tss.ess %>% 
      left_join(id)
  
 
    # format for plotting
    tss.plot <- column_to_rownames(tss, var = "barcode")
    tss.plot <- as.matrix(tss.plot)
  
 
    anno <- data.frame(cluster = id$cluster)
    rownames(anno) <- id$barcode

    # plot with cluster annotation
    col.breaks <- seq(quantile(tss.plot, 0.01),
      quantile(tss.plot, 0.99), length.out=101)

    pdf(sprintf("%s_cluster_tss_enrichment.pdf", opt$outPrefix), 
      height = 6, width = 8)
    pheatmap(tss.plot[order(tss.ess$ess, decreasing = TRUE),], 
      cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = anno, 
      show_rownames = FALSE, show_colnames = FALSE,
      color = viridis(101), breaks = col.breaks)
    dev.off()

    # plot the TSS enrichment score 
    ggplot(data = tss.ess) +
      geom_boxplot(mapping = aes(x = cluster, y = ess)) +
      labs(x = "Cluster", y = "TSS enrichment score",
           title = sprintf("%s", opt$outPrefix)) +
      theme_bw()
    ggsave(sprintf("%s_cluster_ess.pdf", opt$outPrefix), 
      height = 4, width = 6)
 
  }
 
}

main()
