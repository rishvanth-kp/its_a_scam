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

 
suppressMessages(library("Seurat"))
suppressMessages(library("Signac"))
suppressMessages(library("optparse"))
suppressMessages(library("GenomicRanges"))

main <- function() {

  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--inFragFile"),
    help = "Fragmets file with index")
  parser <- add_option(parser, c("-o", "--outfileName"),
    help = "Name of output Rds file")
  parser <- add_option(parser, c("-m", "--minFrags"), default = 1000,
    help = "Number number of fragmets per cell [default: %default]")
  parser <- add_option(parser, c("-c", "--minCells"), default = 0,
    help = "Min. cells to keep a feature [default: %default]")
  parser <- add_option(parser, c("-f", "--minFeatures"), default = 0,
    help = "Min. features to keep a cell [default: %default]")
  parser <- add_option(parser, c("-p", "--macsPath"), default = "macs3",
    help = "Min. features to keep a cell [default: %default]")
  opt <- parse_args(parser)

  if (is.null(opt$inFragFile) || is.null(opt$outfileName)) {
    print_help(parser)
    quit(status = 1)
  }

  # compute the number of fragmets per barcode
  frag.count <- CountFragments(opt$inFragFile)
  print(dim(frag.count)) 
  print(head(frag.count))

  # get a list of barcodes that meet the min frags
  keep.bc <- frag.count[frag.count$frequency_count >= opt$minFrags, ]$CB
  print(length(keep.bc))

  # create a fragment object
  frags <- CreateFragmentObject(path = opt$inFragFile, cells = keep.bc)   

  # call peaks
  peaks <- CallPeaks(frags, macs2.path = opt$macsPath)
  print(dim(peaks))
  print(head(peaks))

  # create cell-peak matrix
  count.matrix <- FeatureMatrix(fragments = frags, 
    features = peaks, cells = keep.bc, sep = c(":", "-"))
  print(dim(count.matrix))

  # create chromatic assay
  chrom.assay <- CreateChromatinAssay(counts = count.matrix, 
    sep = c(":", "-"), fragments = opt$inFragFile, 
    min.cells = opt$minCells, min.features = opt$minFeatures)

  # create seruat assay
  sample <- CreateSeuratObject(counts = chrom.assay, assay = "peaks")
  print(sample)

  # save seurat object
  SaveSeuratRds(sample, file = opt$outfileName)

}   

main()
