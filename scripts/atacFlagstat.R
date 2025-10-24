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
  print(fstat.norm)

  ## Plots
  ggplot(data = fstat.norm) +
    geom_histogram(mapping = aes(primary), binwidth = 0.01) +
    scale_x_log10() +
    labs(x = "Number of primary alignments", y = "Number of cells",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()
  ggsave(sprintf("%s_primary_alignments.pdf", opt$outPrefix), 
    height = 5, width = 5)

  ggplot(data = fstat.norm) +
    geom_histogram(mapping = aes(primary_duplicates), binwidth = 0.01) +
    labs(x = "Fraction of PCR duplicates", y = "Number of cells",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()
  ggsave(sprintf("%s_pcr_duplicates.pdf", opt$outPrefix), 
    height = 5, width = 5)

  ggplot(data = fstat.norm) +
    geom_histogram(mapping = aes(primary_mapped), binwidth = 0.01) +
    labs(x = "Fraction of primary mapped alignments", y = "Number of cells",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()
  ggsave(sprintf("%s_mapped.pdf", opt$outPrefix), 
    height = 5, width = 5)

  ggplot(data = fstat.norm) +
    geom_histogram(mapping = aes(properly_paired), binwidth = 0.01) +
    labs(x = "Fraction of properly paired alignments", y = "Number of cells",
         title = sprintf("%s", opt$outPrefix)) +
    theme_bw()
  ggsave(sprintf("%s_proper_pair.pdf", opt$outPrefix), 
    height = 5, width = 5)

}


main()
