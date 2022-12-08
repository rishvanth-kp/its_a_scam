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

#include "GeneCount.hpp"

#include <unordered_set>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::unordered_map;
using std::unordered_set;

GeneCount::GeneCount() {

}

GeneCount::~GeneCount() {

}


void
GeneCount::preprocess_gff(const string &gff_file) {

  GtfReader gff_reader(gff_file);
  GencodeGtfEntry entry;

  while (gff_reader.read_gencode_gtf_line(entry)) {
    if (entry.feature == "exon") {
      gene_locations.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>(entry.gene_id));
      gene_count.insert(make_pair(entry.gene_id, 0));
    }
    // keep track of metadata for each gene
    else if (entry.feature == "gene") {

    }
  }

}


void
GeneCount::add(const SamEntry &e) {
  SamCigar::CigarRegions ref_regions;
  SamCigar::cigar_to_reference_regions(e, ref_regions);

  // cout << "here" << endl;
  // cout << e.qname << "\t" << e.pos << endl;

  unordered_set<string> aligned_genes;

  vector<pair<GenomicRegion, FeatureVector<string>>> out;
  for (auto it = ref_regions.begin(); it != ref_regions.end(); ++it) {
    if (it->first == SamCigar::Cigar::aln_match ||
        it->first == SamCigar::Cigar::seq_match ||
        it->first == SamCigar::Cigar::seq_mismatch) {
      // cout << it->second << endl;
      gene_locations.at(it->second, out);
      for (auto jt = out.begin(); jt != out.end(); ++jt) {
        // cout << jt->first << "\t";
        for (size_t k = 0; k < jt->second.size(); ++k) {
          // cout << jt->second.at(k) << "\t";
          aligned_genes.insert(jt->second.at(k));
        }
      }
      // cout << endl;
    }
  }
/*
  cout << "gene count: " << aligned_genes.size() << endl;
  for (auto it = aligned_genes.begin(); it != aligned_genes.end(); ++it) {
    cout << *it << endl;
  }
  cout << endl;
*/
  if (aligned_genes.size() == 1) {
    ++gene_count[*aligned_genes.cbegin()];
  }
}


void
GeneCount::add(const SamEntry &e1, const SamEntry &e2) {

  unordered_set<string> aligned_genes;
  vector <pair<GenomicRegion, FeatureVector<string>>> out;
  
  // process first read
  SamCigar::CigarRegions ref_regions1;
  SamCigar::cigar_to_reference_regions(e1, ref_regions1);
  for (auto it = ref_regions1.begin(); it != ref_regions1.end(); ++it) {
    if (it->first == SamCigar::Cigar::aln_match || 
        it->first == SamCigar::Cigar::seq_match ||
        it->first == SamCigar::Cigar::seq_mismatch) {
      
      gene_locations.at(it->second, out);
      for (auto jt = out.begin(); jt != out.end(); ++jt) {
        for (size_t k = 0; k < jt->second.size(); ++k) {
          aligned_genes.insert(jt->second.at(k));
        }
      }
    }
  } 
 
  // process second read
  SamCigar::CigarRegions ref_regions2;
  SamCigar::cigar_to_reference_regions(e2, ref_regions2);
  for (auto it = ref_regions2.begin(); it != ref_regions2.end(); ++it) {
    if (it->first == SamCigar::Cigar::aln_match || 
        it->first == SamCigar::Cigar::seq_match ||
        it->first == SamCigar::Cigar::seq_mismatch) {
      
      gene_locations.at(it->second, out);
      for (auto jt = out.begin(); jt != out.end(); ++jt) {
        for (size_t k = 0; k < jt->second.size(); ++k) {
          aligned_genes.insert(jt->second.at(k));
        }
      }
    }
  } 


  // add gene count
  if (aligned_genes.size() == 1) {
    ++gene_count[*aligned_genes.cbegin()];
  } 
}

void
GeneCount::get_gene_counts(vector<pair<string, size_t>> &counts) {
  counts.clear();
  for (auto it = gene_count.begin(); it != gene_count.end(); ++it) {
    cout << it->first << "\t" << it->second << endl;
  }
}
