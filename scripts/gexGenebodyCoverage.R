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

suppressMessages(library("optparse"))
suppressMessages(library("pheatmap"))
suppressMessages(library("tidyverse"))

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
  gb <- gb %>%
    select(!barcode)
 
  # normalize
  gb <- (gb / rowSums(gb)) * 100
 
  # heatmap for all cells
  pdf(sprintf("%s_bc_genebody_cov.pdf", opt$outPrefix), 
    height = 6, width = 8)
  pheatmap(gb, cluster_cols = FALSE, 
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
  
 
    # format for plotting
    gb <- column_to_rownames(gb, var = "barcode")
 
    anno <- data.frame(cluster = id$cluster)
    rownames(anno) <- id$barcode

    # plot with cluster annotation
    pdf(sprintf("%s_cluster_genebody_cov.pdf", opt$outPrefix), 
      height = 6, width = 8)
    pheatmap(gb, cluster_cols = FALSE, annotation_row = anno, 
      show_rownames = FALSE, show_colnames = FALSE)
    dev.off() 
  }
 
}

main()
