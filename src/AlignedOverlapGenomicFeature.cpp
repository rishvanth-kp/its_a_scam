/*
* AlignedOverlapGenomicFeature: Determine feature of aligned reads
*   allowing a read to align to multiple overlapping features
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

#include <sstream>
#include <fstream>
#include <unordered_set>

#include "GtfReader.hpp"
#include "BedReader.hpp"
#include "AlignedOverlapGenomicFeature.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::unordered_set;

AlignedOverlapGenomicFeature::AlignedOverlapGenomicFeature() {
  bc_counter = 0;
  feature_counter = 0;
}

AlignedOverlapGenomicFeature::~AlignedOverlapGenomicFeature() {

}


void
AlignedOverlapGenomicFeature::add_gtf_features(const string& gtf_file) {

  GtfReader gtf_reader(gtf_file);
  GencodeGtfEntry entry;

  GenomicStepVector<FeatureVector<string>> gtf_features;
  unordered_set<string> chroms;

  while (gtf_reader.read_gencode_gtf_line(entry)) {
    // keep track of all the chromosomes
    chroms.insert(entry.name);

    // keep track of genes, UTR, CDS. Need to keep genes to figure out
    // introns later.
    if (entry.feature == "gene") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{entry.gene_name});
    }
    else if (entry.feature == "UTR") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{"UTR"});
    }
    else if (entry.feature == "CDS") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{"CDS"});
    }
  }

  for (auto chrom_it = chroms.begin(); chrom_it != chroms.end(); ++chrom_it) {

    vector<pair<GenomicRegion, FeatureVector<string>>> out;
    gtf_features.at(*chrom_it, out);
    for (auto it = out.begin(); it != out.end(); ++it) {
      bool cds = false;
      bool utr = false;
      bool gene = false;

      for (size_t j = 0; j < it->second.size(); ++j) {
        if (it->second.at(j) == "CDS") {
          cds = true;
        }
        else if (it->second.at(j) == "UTR") {
          utr = true;
        }
        else {
          gene = true;
        }
      }

      // build reference
      // A region can be associated with multiple features.
      // If a region covers an UTR or CDS, it is assined as that.
      // If a region covers a gene, but does not have an UTR or CDS,
      // it is assigned as intronic.

      if (cds) {
        genomic_features.add(it->first, FeatureVector<string>{"CDS"});
      }
      if (utr) {
        genomic_features.add(it->first, FeatureVector<string>{"UTR"});
      }
      if (!cds & !utr & gene) {
        genomic_features.add(it->first, FeatureVector<string>{"intron"});
      }

    }
  }

  // keep track of the feature types
  feature_index["CDS"] = feature_counter++;
  feature_index["UTR"] = feature_counter++;
  feature_index["intron"] = feature_counter++;
  feature_index["intergene"] = feature_counter++;
}

void
AlignedOverlapGenomicFeature::add_bed_features(const std::string& bed_file,
                                               const std::string feature_name) {

  BedReader bed_reader(bed_file);
  GenomicRegion bed_region;
  while (bed_reader.read_bed3_line(bed_region)) {
    genomic_features.add(bed_region, FeatureVector<string>{feature_name});
  }

  // added to index
  feature_index[feature_name] = feature_counter++;

}

void
AlignedOverlapGenomicFeature::process_barcodes(const std::string& bc_file) {

  std::ifstream bc_in(bc_file);
  if (!bc_in) {
    throw std::runtime_error("Cannot open " + bc_file);
  }
  string line;

  // keep crate an index for all barcodes
  while (getline(bc_in, line)) {
    string barcode;
    std::istringstream iss(line);
    iss >> barcode;

    bc_index[barcode] = bc_counter++;
  }
  bc_in.close();

  // initialize count matrix
  feature_counts.resize(bc_index.size());
  for (size_t i = 0; i < feature_counts.size(); ++i) {
    feature_counts[i].resize(feature_index.size());
  }

  // initialize counted_bases vector
  counted_bases.resize(bc_index.size());
}

void
AlignedOverlapGenomicFeature::add(const SamEntry &e1,
                                  const SamEntry &e2,
                                  const string &bc) {

  unordered_map<string, size_t>::iterator bc_it;
  // search for valid barcode
  bc_it = bc_index.find(bc);
  if (bc_it != bc_index.end()) {
    size_t entry_bc_index = bc_it->second;

    // find the start and end location of the reads
    size_t e1_start = e1.pos - 1;
    size_t e1_end = SamCigar::reference_end_pos(e1);
    size_t e2_start = e2.pos - 1;
    size_t e2_end = SamCigar::reference_end_pos(e2);

    // find the stand and end location of the fragments
    // TODO: verify that these match the TLEN
    size_t frag_start = e1_start;
    if (e2_start < e1_start)
      frag_start = e2_start;

    size_t frag_end = e1_end;
    if (e2_end < e1_end)
      frag_end = e2_end;

    // find overlapping regions
    vector<pair<GenomicRegion, FeatureVector<string>>> out;
    genomic_features.at(GenomicRegion(e1.rname, frag_start, frag_end),
                        out, true);

    for (auto it = out.begin(); it != out.end(); ++it) {
      const size_t match_len = it->first.end - it->first.start;
      // add match length for the sample.
      counted_bases[entry_bc_index] += match_len;

      // intergenic, if not aligned to any other feature
      if (it->second.size() == 0) {
        size_t index = feature_index["intergene"];
        feature_counts[entry_bc_index][index] += match_len;
      }

      for (size_t j = 0; j < it->second.size(); ++j) {
        size_t index = feature_index[it->second.at(j)];
        feature_counts[entry_bc_index][index] += match_len;
      }

    }
  }

}

void AlignedOverlapGenomicFeature::feature_counts_to_file(
      const std::string &file_name) const {

  std::ofstream out(file_name);
  if (!out)
    throw std::runtime_error("Cannot open " + file_name);

  // write header
  out << "barcode" << "\tbases";
  vector<pair<string, size_t>> features;
  for (auto it = feature_index.begin(); it != feature_index.end(); ++it) {
    features.push_back(std::make_pair(it->first, it->second));
    out << "\t" << it->first;
  }
  out << endl;


  for (auto it = bc_index.begin(); it != bc_index.end(); ++it) {
    out << it->first;
    out << "\t" << counted_bases[it->second];
    for (size_t j = 0; j < features.size(); ++j) {
      out << "\t" << feature_counts[it->second][features[j].second];
    }
    out << endl;
  }

  out.close();
}
