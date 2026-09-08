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
suppressMessages(library("RColorBrewer"))

main <- function() {

  parser <- OptionParser()
  parser <- add_option(parser, c("-f", "--fraglenFile"),
              help = "single cell flagstat QC file")
  parser <- add_option(parser, c("-s", "--sampleFraglenFile"),
              help = "single cell flagstat QC file")
  parser <- add_option(parser, c("-c", "--clusterIdFile"),
              help="single cell cluster ID TSV file")
  parser <- add_option(parser, c("-m", "--minFraglen"), default = 25,
              help = "min fragment length to plot [default: %default]")
  parser <- add_option(parser, c("-M", "--maxFraglen"), default = 300,
              help = "max fragment length to plot [default: %default]")
  parser <- add_option(parser, c("-o", "--outPrefix"),
              help = "Outfile prefix")
  opt <- parse_args(parser)

  if (is.null(opt$fraglenFile) | is.null(opt$sampleFraglenFile) | 
      is.null(opt$outPrefix)) {
    print_help(parser)
    quit(status = 1)
  }

  ## sample fragment length distribution
  # read the data
  sample.fl <- read_tsv(opt$sampleFraglenFile)
  
  # plot
  ggplot(data = sample.fl) +
    geom_line(mapping = aes(x = frag_len, y = norm_count)) +
    labs(x = "Fragment length", y = "Normalized fragment count",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()
  ggsave(sprintf("%s_sample_frag_len.pdf", opt$outPrefix), 
    height = 5, width = 5)
    

  ## single cell fragment length distribution
  # read data
  fl <- read_tsv(opt$fraglenFile)
  
  # remove unwanted columns
  fl.len <- fl %>%
    select(!barcode) %>%
    select(seq(opt$minFraglen, opt$maxFraglen))
  
  # heatmap for all cells
  pdf(sprintf("%s_bc_frag_len.pdf", opt$outPrefix), 
    height = 6, width = 8)
  pheatmap(fl.len, cluster_cols = FALSE, color = viridis(100), 
    show_rownames = FALSE, show_colnames = FALSE)
  dev.off()


  # heatmap with cluster ids
  if (!is.null(opt$clusterIdFile)) {
   
    # read the cluster IDs
    id <- read_tsv(opt$clusterIdFile, col_names = F)
    names(id) <- c("barcode", "cluster")
    id$cluster <- as.factor(id$cluster)  
  
    # keep only the cells that have a cluster ID
    fl <- fl[fl$barcode %in% id$barcode, ]
  
    # keep only the needed fragment lenghts 
    fl.len <- fl %>%
      select(barcode, seq(opt$minFraglen, opt$maxFraglen))


    ## summary plot
    fl.cluster <- fl.len %>% 
      left_join(id) %>% 
      select(barcode, cluster, everything())
    print(fl.cluster) 

    fl.summary <- fl.cluster %>% 
      group_by(cluster) %>% 
      summarise(across(where(is.numeric), mean)) %>% 
      ungroup() 
  
    fl.summary <- fl.summary %>% 
      pivot_longer(!cluster, names_to="frag.len", values_to = "count")

    cluster.col <- brewer.pal(length(unique(id$cluster)), "Set1")
    names(cluster.col) <- levels(id$cluster)

    ## plot cluster summary
    ggplot(data = fl.summary) + 
      geom_line(mapping = aes(x = as.numeric(frag.len), y = count, 
        group = cluster, color = cluster)) +
      labs(x = "Fragment length", y = "Normalized fragment count",
        title = sprintf("%s", opt$outPrefix)) +
      scale_colour_manual(values = cluster.col) +
      theme_bw()
    ggsave(sprintf("%s_cluster_frag_len_summary.pdf", opt$outPrefix), 
      height = 5, width = 5)

 
    # format for plotting
    fl.len <- column_to_rownames(fl.len, var = "barcode")
 
    anno <- data.frame(cluster = id$cluster)
    rownames(anno) <- id$barcode

    # plot with cluster annotation
    cluster.col <- list(cluster = cluster.col)
    pdf(sprintf("%s_cluster_frag_len.pdf", opt$outPrefix), 
      height = 6, width = 8)
    pheatmap(fl.len, cluster_cols = FALSE, annotation_row = anno, 
      show_rownames = FALSE, show_colnames = FALSE,
      annotation_color = cluster.col, color = viridis(100))
    dev.off() 
  }
 
}

main()
