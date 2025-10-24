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
  parser <- add_option(parser, c("-c", "--featureCounts"),
              help="Cell barcoded feature count QC file")
  parser <- add_option(parser, c("-a", "--featureAlignLen"),
              help="Cell barcoded aligned fragment QC file")
  parser <- add_option(parser, c("-f", "--featureFragLen"),
              help="Cell barcoded fragment length QC file")
  parser <- add_option(parser, c("-o", "--outPrefix"),
              help="Outfile prefix")
  opt <- parse_args(parser)

  if (is.null(opt$featureCounts) | is.null(opt$featureAlignLen) |
      is.null(opt$featureFragLen) | is.null(opt$outPrefix)) {
    print_help(parser)
    quit(status = 1)
  }
 
  ## Read the feature qc files
  counts <- read_tsv(opt$featureCounts)
 
  ## as fraction of aligned bases
  counts[, 2:ncol(counts)] <- (counts[, 2:ncol(counts)] / counts$count) * 100
  
  ## format for ggplot
  row.order <- colnames(counts[, 2:ncol(counts)])
  
  counts <- counts %>% 
    select(!c(count)) %>% 
    pivot_longer(!barcode, names_to = "region", values_to = "value")
  counts$region <- factor(counts$region, 
    labels = row.order, levels = row.order)
  
  ## plot the data
  ggplot(data = counts) + 
    geom_boxplot(mapping = aes(x = region, y = value), outlier.size = 0.5) +
    labs(x = "Genomic Region", y = "Percent of reads",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()

  ggsave(sprintf("%s_genome_feature_counts.pdf", opt$outPrefix), 
    height = 4, width = 6)




  #### Aligned region of entire fragment 
  ## Read the feature qc files
  align <- read_tsv(opt$featureAlignLen)

  ## as fraction of aligned bases
  align[, 2:ncol(align)] <- (align[, 2:ncol(align)] / align$bases) * 100

  ## format for ggplot
  align <- align %>% 
    select(!c(bases)) %>% 
    pivot_longer(!barcode, names_to = "region", values_to = "value")
  align$region <- factor(align$region, 
    labels = row.order, levels = row.order)

  ## plot the data
  ggplot(data = align) + 
    geom_boxplot(mapping = aes(x = region, y = value), outlier.size = 0.5) +
    labs(x = "Genomic Region", y = "Percent of aligned bases",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()

  ggsave(sprintf("%s_genome_feature_align_len.pdf", opt$outPrefix), 
    height = 4, width = 6)



  #### mean fragment lenth for each feature
  ## read QC file
  frag.len <- read_tsv(opt$featureFragLen)


  ## format for ggplot
  frag.len <- frag.len %>% 
    pivot_longer(!barcode, names_to = "region", values_to = "value")
  frag.len$region <- factor(frag.len$region, 
    labels = row.order, levels = row.order)


  ## plot the data
  ggplot(data = frag.len) + 
    geom_boxplot(mapping = aes(x = region, y = value), outlier.size = 0.5) +
    labs(x = "Genomic Region", y = "Mean fragment length",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()

  ggsave(sprintf("%s_genome_feature_frag_len.pdf", opt$outPrefix), 
    height = 4, width = 6)

}

main()
