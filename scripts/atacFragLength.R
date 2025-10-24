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
  parser <- add_option(parser, c("-f", "--fraglenFile"),
              help = "single cell flagstat QC file")
  parser <- add_option(parser, c("-s", "--sampleFraglenFile"),
              help = "single cell flagstat QC file")
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
  fl <- fl %>%
    select(!barcode) %>%
    select(seq(opt$minFraglen, opt$maxFraglen))

 
  print(fl)
 
  # plot
  pdf(sprintf("%s_bc_frag_len.pdf", opt$outPrefix), 
    height = 6, width = 8)
  pheatmap(fl, cluster_cols = FALSE, 
    show_rownames = FALSE, show_colnames = FALSE)
  dev.off()
 
}

main()
