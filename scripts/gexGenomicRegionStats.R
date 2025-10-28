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
suppressMessages(library("tidyverse"))

main <- function() {

  parser <- OptionParser()
  parser <- add_option(parser, c("-f", "--featureCounts"),
              help="Cell barcoded feature count QC file")
  parser <- add_option(parser, c("-c", "--clusterIdFile"),
              help="single cell cluster ID TSV file")
  parser <- add_option(parser, c("-o", "--outPrefix"),
              help="Outfile prefix")
  opt <- parse_args(parser)

  if (is.null(opt$featureCounts) | is.null(opt$outPrefix)) {
    print_help(parser)
    quit(status = 1)
  }
 
  #### counts for each feature
  ## Read the feature qc files
  counts <- read_tsv(opt$featureCounts)
 
  ## as fraction of aligned bases
  counts[, 2:ncol(counts)] <- (counts[, 2:ncol(counts)] / counts$bases) * 100
  print(counts) 
 
  ## format for ggplot
  row.order <- c("pc_cds", "pc_utr", "pc_intron", 
    "intergenic", "chrM_protein_coding", "lncRNA")

  counts <- counts %>% 
    select(c("barcode", all_of(row.order))) %>% 
    pivot_longer(!barcode, names_to = "region", values_to = "value")
  counts$region <- factor(counts$region, 
    labels = row.order, levels = row.order)
  print(counts)  

  ## plot the data
  ggplot(data = counts) + 
    geom_boxplot(mapping = aes(x = region, y = value), outlier.size = 0.5) +
    labs(x = "Genomic Region", y = "Percent of bases",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()
  ggsave(sprintf("%s_gex_feature_counts.pdf", opt$outPrefix), 
    height = 4, width = 6)

  ## plot by cluster
  if (!is.null(opt$clusterIdFile)) {
    # read the cluster IDs
    id <- read_tsv(opt$clusterIdFile, col_names = F)
    names(id) <- c("barcode", "cluster")
    id$cluster <- as.factor(id$cluster)  

    # keep only the cells that have a cluster ID   
    counts <- counts[counts$barcode %in% id$barcode, ]
 
    # join the cluster IDs
    counts <- counts %>%
      left_join(id)
 
    # plot 
    ggplot(data = counts) + 
      geom_boxplot(mapping = aes(x = cluster, y = value), outlier.size = 0.5) +
      facet_wrap(vars(region), scales = "free_y") +
      labs(x = "Cluster", y = "Percent of bases",
           title = sprintf("%s", opt$outPrefix)) +
      theme_bw()
    ggsave(sprintf("%s_cluster_gex_feature_counts.pdf", opt$outPrefix), 
      height = 4, width = 6)
  }

}

main()
