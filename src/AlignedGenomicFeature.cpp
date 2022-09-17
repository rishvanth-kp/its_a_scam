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

#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <vector>

#include "AlignedGenomicFeature.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::unordered_set;

AlignedGenomicFeature::AlignedGenomicFeature() {
  match_bases = 0;
}

AlignedGenomicFeature::~AlignedGenomicFeature() {

}


void
AlignedGenomicFeature::preprocess_gff(const string& gff_file) {

  GtfReader gff_reader(gff_file);
  GencodeGtfEntry entry;

  GenomicStepVector<FeatureVector<string>> gtf_features;
  unordered_set<string> chroms;
  unordered_set<string> features;

  while (gff_reader.read_gencode_gtf_line(entry)) {
    chroms.insert(entry.name);

    if (entry.feature == "gene") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{entry.gene_type});
    }
    else if (entry.feature == "UTR") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{"UTR"});
    }
    else if (entry.feature == "CDS") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{"CDS"});
    }
    else if (entry.feature == "exon") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{"exon"});
    }
  }

  for (auto chrom_it = chroms.begin(); chrom_it != chroms.end(); ++chrom_it) {
    vector<pair<GenomicRegion, FeatureVector<string>>> out;
    gtf_features.at(*chrom_it, out);
    for (auto it = out.begin(); it != out.end(); ++it) {
      size_t gene_count{0};
      bool exon{false};
      bool cds{false};
      bool utr{false};
      string gene_type;

      /*
      cout << it->first << "\t";
      for (size_t j = 0; j < it->second.size(); ++j) {
        cout << it->second.at(j) << "\t";
      }
      cout << endl;
      */

      for (size_t j = 0; j < it->second.size(); ++j) {
        if (it->second.at(j) == "exon") {
          exon = true;
        }
        else if (it->second.at(j) ==  "CDS") {
          cds = true;
        }
        else if (it->second.at(j) == "UTR") {
          utr = true;
        }
        else {
          ++gene_count;
          gene_type = it->second.at(j);
        }
      }

      // build reference
      if (it->first.name == "chrM") {
        genomic_features.add(it->first, "chrM_" + gene_type);
        features.insert("chrM_" + gene_type);
      }
      else {
        if (gene_count > 1) {
          genomic_features.add(it->first, "ambiguous");
          features.insert("ambiguous");
        }
        else if (gene_type == "protein_coding") {
          if (!exon) {
            genomic_features.add(it->first, "pc_intron");
            features.insert("pc_intron");
          }
          else if (cds & !utr) {
            genomic_features.add(it->first, "pc_cds");
            features.insert("pc_cds");
          }
          else if (!cds & utr) {
            genomic_features.add(it->first, "pc_utr");
            features.insert("pc_utr");
          }
          else if (cds & utr) {
            genomic_features.add(it->first, "pc_both");
            features.insert("pc_utr");
          }
          else {
            genomic_features.add(it->first, "pc_none");
            features.insert("pc_utr");
          }

        }
        else {
          genomic_features.add(it->first, gene_type);
          features.insert(gene_type);
        }
      }
    }

  }

  features.insert("intergenic");
  for (auto it = features.begin(); it != features.end(); ++it) {
    feature_count.insert(make_pair(*it, 0));
  }
/*
  for (auto chrom_it = chroms.begin(); chrom_it != chroms.end(); ++chrom_it) {
    vector<pair<GenomicRegion, string>> out;
    genomic_features.at(*chrom_it, out);
    for (auto it = out.begin(); it != out.end(); ++it) {
      cout << it->first << "\t" << it->second << endl;
    }
  }
*/

}


void
AlignedGenomicFeature::add(const SamEntry &e) {
  SamCigar::CigarRegions ref_regions;
  SamCigar::cigar_to_reference_regions(e, ref_regions);

  vector<pair<GenomicRegion, string>> out;
  for (auto it = ref_regions.begin(); it != ref_regions.end(); ++it) {
    if (it->first == SamCigar::Cigar::aln_match ||
        it->first == SamCigar::Cigar::seq_match ||
        it->first == SamCigar::Cigar::seq_mismatch) {

      genomic_features.at(it->second, out, true);
      for (auto jt = out.begin(); jt != out.end(); ++jt) {
        const size_t match_len = jt->first.end - jt->first.start;
        match_bases += match_len;
        if (jt->second == string{}) {
          feature_count["intergenic"] += match_len;
        }
        else {
          feature_count[jt->second] += match_len;
        }
      }
    }
  }
}


void
AlignedGenomicFeature::clear_counts() {
  for (auto it = feature_count.begin(); it != feature_count.end(); ++it) {
    it->second = 0;
  }
  match_bases = 0;
}


static bool
sort_features(pair<string, size_t> a, pair<string, size_t> b) {
  return (a.first < b.first);
}

void
AlignedGenomicFeature::get_feature_counts(vector<pair<string, size_t>> &counts,
                                          size_t &n_bases) {
  counts.clear();
  for (auto it = feature_count.begin(); it != feature_count.end(); ++it) {
    counts.push_back(std::make_pair(it->first, it->second));
  }
  std::sort(counts.begin(), counts.end(), sort_features);
  n_bases = match_bases;
}

void
AlignedGenomicFeature::feature_count_to_file(const string& file_name) const{
  std::ofstream out_file(file_name);

  vector<pair<string, size_t>> counts;
  for (auto it = feature_count.begin(); it != feature_count.end(); ++it) {
    counts.push_back(std::make_pair(it->first, it->second));
  }
  std::sort(counts.begin(), counts.end(), sort_features);

  for (auto it = counts.begin(); it != counts.end(); ++it) {
    const float feature_percent =
      (static_cast<float>(it->second) / static_cast<float>(match_bases)) * 100;
    out_file << it->first << "\t"
             << it->second << "\t"
             << feature_percent << endl;
  }

  out_file.close();
}
