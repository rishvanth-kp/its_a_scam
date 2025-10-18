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

  bc1.r <- bc1 %>% filter(stype == "R")
  counts.r <- counts[, substr(colnames(counts), 19, 19+7) %in% bc1.r$sequence]
  print(bc1.r)
  print(dim(counts.r))

  bc1.t$rsequence <- bc1.r$sequence
  print(bc1.t)

  # random hexamer collapsing
  counts.c <- as.matrix(counts.t)
  print(counts.c[1:5, 1:5]) 
  
  n.match <- 0
  n.gt <- 0 

  for (i in 1:ncol(counts.c)) { 
    t.id <- colnames(counts.c)[i]
  
    # replace the 3' barcode with the random barcode
    r.id <- paste0(substr(t.id, 1, 18), 
      bc1.t[bc1.t$sequence == substr(t.id, 19, 19+7), ]$rsequence)

    # make sure it exists only once
    if (sum(colnames(counts.r) %in% r.id) == 1) {
      # add the random count to the 3' count
      # print(head(counts.c[, t.id]))
      # print(head(counts.r[, r.id]))
      counts.c[, t.id] <- counts.c[, t.id] + counts.r[, r.id]
      n.match <- n.match + 1
      # print(head(counts.c[, t.id]))
      # print(dim(counts.c))
    }
    else if (sum(colnames(counts.r) %in% r.id) > 1) {
      n.gt <- n.gt + 1
    }


  }
  print(n.match)
  print(n.gt)
  

  # rename barcodes to match the parse splitpile style by
  # renamint the barcode with the well number
  bc.names <- tibble(bc = colnames(counts.c))
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

    counts.t <- CreateSeuratObject(counts = counts.t, names.field = 0)
    counts.t$sample <- bc.names$sample
  }
  else {
    # sanity check
    print("Colname mismatch between counts.t and bc.names")
  }
  
  if (sum(bc.names$bc == colnames(counts.c)) == ncol(counts.c)) {
    colnames(counts.c) <- paste(
      str_pad(bc.names$bci1, 2, pad = "0"), 
      str_pad(bc.names$bci2, 2, pad = "0"), 
      str_pad(bc.names$bci3, 2, pad = "0"), sep = "_")
    print(head(colnames(counts.c)))
    
    # add the sub-library barcode
    colnames(counts.c) <- paste0(colnames(counts.c), "__", opt$subLibrary)
    
    counts.c <- CreateSeuratObject(counts = counts.c, names.field = 0)
    counts.c$sample <- bc.names$sample
  }
  else {
    # sanity check
    print("Colname mismatch between counts.c and bc.names")
  }
  
   
  # write the 3' count matrix
  SaveSeuratRds(counts.t, sprintf("%s_3prime.Rds", opt$outPrefix))
  print(counts.t) 

  # write the random hexamer collapsed count matrix
  SaveSeuratRds(counts.c, sprintf("%s_collapse.Rds", opt$outPrefix))
  print(counts.c) 
    

}

main()
