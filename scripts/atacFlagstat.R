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
  parser <- add_option(parser, c("-f", "--flagstatFile"),
              help="single cell flagstat QC file")
  parser <- add_option(parser, c("-c", "--clusterIdFile"),
              help="single cell cluster ID TSV file")
  parser <- add_option(parser, c("-o", "--outPrefix"),
              help="Outfile prefix")
  opt <- parse_args(parser)

  if (is.null(opt$flagstatFile) | is.null(opt$outPrefix)) {
    print_help(parser)
    quit(status = 1)
  }


  ## Read the input
  fstat <- read_tsv(opt$flagstatFile)
 
  ## Normalize
  fstat.norm <- tibble(barcode = fstat$barcode, 
    primary = fstat$primary)

  fstat.norm$primary_duplicates <- fstat$primary_duplicates / fstat$primary
  fstat.norm$primary_mapped <- fstat$primary_mapped / fstat$primary
  fstat.norm$properly_paired <- fstat$properly_paired / fstat$primary



  # format data
  fstat.norm <- fstat.norm %>%
    # select(!primary) %>% 
    pivot_longer(!barcode, names_to = "type", values_to = "value")
  
  # plot sample level stats  
  ggplot(data = fstat.norm) +
    geom_histogram(mapping = aes(value), bins = 100) +
    facet_wrap(vars(type), scales = "free") +
    labs(x = "Number/Fraction of primary alignments", 
         y = "Number of cells",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()
  ggsave(sprintf("%s_sample_flagstat.pdf", opt$outPrefix), 
    height = 4, width = 6)


  ## if the cluster id file is provided, make cluster-wise plots
  if (!is.null(opt$clusterIdFile)) {
    # read the cluster IDs
    id <- read_tsv(opt$clusterIdFile, col_names = F)
    names(id) <- c("barcode", "cluster")
    id$cluster <- as.factor(id$cluster)  
  
    # keep only the cells that have a cluster ID
    fstat.norm <- fstat.norm[fstat.norm$barcode %in% id$barcode, ]
    
    # join the cluster IDs
    fstat.norm <- fstat.norm %>%
      left_join(id)

    # plot cluster level stats
    ggplot(data = fstat.norm) +
      geom_boxplot(mapping = aes(x = cluster, y = value), outliers = FALSE) +
      facet_wrap(vars(type), scales = "free_y") +
      labs(x = "Cluster", y = "Number/Fraction of primary alignments",
           title = sprintf("%s", opt$outPrefix)) +
      theme_bw()
    ggsave(sprintf("%s_cluster_flagstat.pdf", opt$outPrefix), 
      height = 4, width = 6)
  

  }
}


main()
