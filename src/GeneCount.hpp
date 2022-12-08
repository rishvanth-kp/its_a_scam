/*
* GeneCount: Determine gene counts of aligned reads
* Copyright (C) 2022 Rishvanth Prabakar
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef GENE_COUNT_HPP
#define GENE_COUNT_HPP

#include <string>
#include <unordered_map>

#include "SamEntry.hpp"
#include "GtfReader.hpp"
#include "FeatureVector.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStepVector.hpp"

class GeneCount {
public:
  GeneCount();
  ~GeneCount();

  void preprocess_gff(const std::string &gff_file);

  void add(const SamEntry &e);
  void add(const SamEntry &e1, const SamEntry &e2);

  void get_gene_counts(vector<pair<string, size_t>> &counts);


private:
  GenomicStepVector<FeatureVector<string>> gene_locations;
  std::unordered_map<string, size_t> gene_count;
};

#endif
