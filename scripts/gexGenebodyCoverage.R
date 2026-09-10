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
suppressMessages(library("cowplot"))
suppressMessages(library("optparse"))
suppressMessages(library("pheatmap"))
suppressMessages(library("tidyverse"))
suppressMessages(library("RColorBrewer"))

main <- function() {

  parser <- OptionParser()
  parser <- add_option(parser, c("-g", "--genebodyCovFile"),
              help = "single cell flagstat QC file")
  parser <- add_option(parser, c("-c", "--clusterIdFile"),
              help="single cell cluster ID TSV file")
  parser <- add_option(parser, c("-o", "--outPrefix"),
              help = "Outfile prefix")
  opt <- parse_args(parser)

  if (is.null(opt$genebodyCovFile) | is.null(opt$outPrefix)) {
    print_help(parser)
    quit(status = 1)
  }


  # read data
  gb <- read_tsv(opt$genebodyCovFile)
  
  # remove unwanted columns
  gb.plot <- gb %>%
    select(!barcode)
 
  # normalize
  gb.plot <- (gb.plot / rowSums(gb.plot)) * 100
 
  # heatmap for all cells
  pdf(sprintf("%s_bc_genebody_cov.pdf", opt$outPrefix), 
    height = 6, width = 8)
  pheatmap(gb.plot, cluster_cols = FALSE, color = viridis(100),
    show_rownames = FALSE, show_colnames = FALSE)
  dev.off()


  # heatmap with cluster ids
  if (!is.null(opt$clusterIdFile)) {
   
    # read the cluster IDs
    id <- read_tsv(opt$clusterIdFile, col_names = F)
    names(id) <- c("barcode", "cluster")
    id$cluster <- as.factor(id$cluster)  
  
    # keep only the cells that have a cluster ID
    gb <- gb[gb$barcode %in% id$barcode, ]
  
    # plot summary over all cells in a cluster
    gb.cluster <- gb %>% 
      left_join(id) %>% 
      select(barcode, cluster, everything())

    gb.summary <- gb.cluster %>% 
      group_by(cluster) %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      ungroup() 
  
    gb.summary[, 2:ncol(gb.summary)] <- 
      (gb.summary[, 2:ncol(gb.summary)] / rowSums(gb.summary[, 2:ncol(gb.summary)])) * 100

    gb.summary <- gb.summary %>% 
      pivot_longer(!cluster, names_to="gene.body", values_to = "count")

    cluster.col <- brewer.pal(length(unique(id$cluster)), "Set1")
    names(cluster.col) <- levels(id$cluster)


    ## plot cluster summary
    ggplot(data = gb.summary) + 
      geom_line(mapping = aes(x = as.numeric(gene.body), y = count, 
        group = cluster, color = cluster)) +
      labs(x = "Gene body (5' - 3')", y = "Fraction of reads",
        title = sprintf("%s", opt$outPrefix)) +
      scale_colour_manual(values = cluster.col) +
      theme_cowplot(10)
    ggsave(sprintf("%s_cluster_genebody_cov_summary.pdf", opt$outPrefix), 
      height = 5, width = 5)
 
    # format for plotting
    gb <- column_to_rownames(gb, var = "barcode")
    gb <- (gb / rowSums(gb)) * 100
 
    anno <- data.frame(cluster = id$cluster)
    rownames(anno) <- id$barcode

    # plot with cluster annotation
    cluster.col <- list(cluster = cluster.col)
    pdf(sprintf("%s_cluster_genebody_cov.pdf", opt$outPrefix), 
      height = 6, width = 8)
    pheatmap(gb, cluster_cols = FALSE, annotation_row = anno,
      annotation_color = cluster.col,  border_color = NA, 
      show_rownames = FALSE, show_colnames = FALSE,
      color = viridis(100))
    dev.off() 
  }
 
}

main()
