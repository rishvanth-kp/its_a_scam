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
suppressMessages(library("cowplot"))

main <- function() {
  
  ## parse command line arguments
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--inputMetagene"),
              help = "Metagene input CSV file")
  parser <- add_option(parser, c("-g", "--geneList"),
              help = "List of genes to include in the analysis")
  parser <- add_option(parser, c("-c", "--groupList"),
              help = "List of cell groups to include")
  parser <- add_option(parser, c("-o", "--outPrefix"),
              help = "Outfile prefix") 
  opt <- parse_args(parser)

  if (is.null(opt$inputMetagene) || is.null(opt$outPrefix)) {
    print_help(parser)
    quit(status = 1)
  }

  
  ## read the input 
  metagene <- read_tsv(opt$inputMetagene)

  ## Optionally, keep only the needed groups
  if (!is.null(opt$groupList)) {
    group.list <- read_csv(opt$groupList, col_names = FALSE)
    print(group.list)
    metagene <- metagene %>%
      filter(metagene$group %in% group.list$X1)
    print(metagene)
  }


  ## Nomalize the data
  metagene <-  metagene %>% group_by(group) %>% 
    mutate(across(where(is.numeric)), 
    (across(where(is.numeric)) * 10^6) / sum(across(where(is.numeric))) ) %>%
    ungroup() 

  ## Optionally, keep only the needed genes
  if (!is.null(opt$geneList)) {
    gene.list <- read_tsv(opt$geneList, col_names = FALSE)
    metagene <- metagene %>% 
      filter(metagene$feature %in% gene.list$X1)
    print(metagene)
  }

  ## get the mean of the merged metagene
  metagene.mean <- metagene %>% group_by(group) %>% 
    summarise(across(where(is.numeric), mean))  %>% 
    ungroup() %>%
    pivot_longer(cols = !c("group"), names_to = "pct", values_to = "count" )
  
  metagene.mean$group <- as.factor(metagene.mean$group)
  ## plot the merged metagene 
  ggplot(metagene.mean, aes(x = as.numeric(pct), y = as.numeric(count), 
      color = group)) + 
    geom_line(linewidth=0.5) + 
    labs(x = "Genomic position", y = "Bases per million") + 
    theme_bw()
  ggsave(sprintf("%s_metagene.pdf", opt$outPrefix), height = 4, width = 5)


}

main()
