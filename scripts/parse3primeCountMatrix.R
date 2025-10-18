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

suppressMessages(library("tidyverse"))
suppressMessages(library("optparse"))
suppressMessages(library("Seurat"))


main <- function() {

  parser <- OptionParser()
  parser <- add_option(parser, c("-r", "--rawCountDir"),
    help = "Directory containing STAR count files")
  parser <- add_option(parser, c("-m", "--matrixFile"), default = "matrix.mtx",
    help = "Count matrix file name [default: %default]")
  parser <- add_option(parser, c("-f", "--filteredBCs"),
    help = "Filtered barcode list")
  parser <- add_option(parser, c("-n", "--minCount"), default = 500,
    help = "Min. count per barcode to keep (used if -f is NULL) 
                [default: %default]")
  parser <- add_option(parser, c("-a", "--barcode1"),
    help = "Parse barcode 1 CSV file")
  parser <- add_option(parser, c("-b", "--barcode2"),
    help = "Parse barcode 2 and 3 CSV file")
  parser <- add_option(parser, c("-l", "--subLibrary"),
    help = "Sub-library ID")
  parser <- add_option(parser, c("-s", "--sampleIDs"),
    help = "CSV file with barcode 1 sample IDs")
  parser <- add_option(parser, c("-o", "--outPrefix"),
    help = "Output prefix")
  opt <- parse_args(parser)

  if (is.null(opt$rawCountDir) |
      is.null(opt$barcode1) | is.null(opt$barcode2) |
      is.null(opt$subLibrary) | is.null(opt$sampleIDs) | 
      is.null(opt$outPrefix)) {
    print_help(parser)
    quit(status = 1)
  }


  # read raw count matrix
  counts <- ReadMtx(
    mtx = paste0(opt$rawCountDir, "/", opt$matrixFile),
    cells = paste0(opt$rawCountDir, "/barcodes.tsv"),
    features = paste0(opt$rawCountDir, "/features.tsv")) 
  print(dim(counts))

 

  if (!is.null(opt$filteredBCs)) {
    # read filtered barcode list
    filtered.cells <- read_tsv(opt$filteredBCs, col_names = FALSE)

    # subset the raw matrix to have only the filtered barcodes
    counts <- counts[, filtered.cells$X1]
    print(dim(counts))
  }
  else {
    # keep the barcodes with min count
    counts <- counts[, colSums(counts) >= opt$minCount]
    print(dim(counts))
      
  }

  # read parse barcodes. Parse has 2 barcode files: one for the
  # first barcode that can either be a 3' or random hexmer for
  # each cell. And one file for the second and third barcodes.
  bc1 <- read_csv(opt$barcode1)
  bc2 <- read_csv(opt$barcode2) 

  # Separate out the count matrix based on 3' barcodes and
  # random hexamer barcodes
  bc1.t <- bc1 %>% filter(stype == "T")
  counts.t <- counts[, substr(colnames(counts), 19, 19+7) %in% bc1.t$sequence]
  print(bc1.t)
  print(dim(counts.t))


  # rename barcodes to match the parse splitpile style by
  # renamint the barcode with the well number
  bc.names <- tibble(bc = colnames(counts.t))
  bc.names$bc3 <- substr(bc.names$bc, 1, 8)
  bc.names$bc2 <- substr(bc.names$bc, 10, 17)
  bc.names$bc1 <- substr(bc.names$bc, 19, 26)

  bc1.t <- bc1.t %>% select(c(bci, sequence, well))
  bc2 <- bc2 %>% select(c(bci, sequence))

  bc.names <- bc.names %>% 
    left_join(bc2, join_by(bc3 == sequence)) %>% rename(bci3 = bci)
  bc.names <- bc.names %>% 
    left_join(bc2, join_by(bc2 == sequence)) %>% rename(bci2 = bci)
  bc.names <- bc.names %>% 
    left_join(bc1.t, join_by(bc1 == sequence)) %>% rename(bci1 = bci)

  print(bc.names)

  # demux the samples based on the first barcode
  sample.id <- read_csv(opt$sampleIDs, col_names = FALSE)
  names(sample.id) <- c("well", "sample")
  
  bc.names <- bc.names %>% left_join(sample.id) 
  print(table(bc.names$sample))
 
  if (sum(bc.names$bc == colnames(counts.t)) == ncol(counts.t)) {
    colnames(counts.t) <- paste(
      str_pad(bc.names$bci1, 2, pad = "0"), 
      str_pad(bc.names$bci2, 2, pad = "0"), 
      str_pad(bc.names$bci3, 2, pad = "0"), sep = "_")
    print(head(colnames(counts.t)))
  
    # add the sub-library barcode
    colnames(counts.t) <- paste0(colnames(counts.t), "__", opt$subLibrary)

    # create seurat objet
    counts.t <- CreateSeuratObject(counts = counts.t, names.field = 0)
    # add sample metadata
    counts.t$sample <- bc.names$sample
     
    # write the 3' count matrix
    SaveSeuratRds(counts.t, sprintf("%s_3prime.Rds", opt$outPrefix))
    print(counts.t) 
  }
  else {
    # sanity check
    print("Colname mismatch between counts.t and bc.names")
  }
  
  
}

main()
