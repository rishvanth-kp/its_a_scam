/*
* AlignedGenomicFeature: Determine feature of aligned reads
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

#ifndef ALIGNED_GENOMIC_FEATURE_HPP
#define ALIGNED_GENOMIC_FEATURE_HPP

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "SamEntry.hpp"
#include "GtfReader.hpp"
#include "FeatureVector.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStepVector.hpp"

class AlignedGenomicFeature {
public:
  AlignedGenomicFeature();
  ~AlignedGenomicFeature();

  void preprocess_gff(const std::string& gff_file);

  void process_barcodes(const std::string& bc_file);

  void add(const SamEntry &e);

  void add(const SamEntry &e, const std::string &bc);

  void add(const SamEntry &e, const std::string &bc, 
           const std::string &umi);

  void clear_counts();

  void get_feature_counts(std::vector<pair<std::string, size_t>> &counts,
                          size_t &n_bases);
  void feature_count_to_file(const std::string& file) const;


private:
  GenomicStepVector<std::string> genomic_features;
  
  // for bulk samples that do not have barcodes
  // This is slightly ineffiecnet since these structures
  // are created even the sample is barcoded. 
  std::unordered_map<std::string, size_t> feature_count;
  size_t match_bases;

  // for barcoded samples
  std::unordered_map<std::string, size_t> feature_index;
  std::unordered_map<std::string, size_t> bc_index;

  std::vector<std::vector<size_t>> bc_feature_count;
  std::vector<size_t> bc_counted_bases;

  // for barcoded and UMI samples
  std::unordered_set<std::string> seen_cb_umi;
};


#endif
