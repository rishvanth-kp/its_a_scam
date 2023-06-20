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

}

AlignedOverlapGenomicFeature::~AlignedOverlapGenomicFeature() {

}
  

void
AlignedOverlapGenomicFeature::add_gtf_features(const string& gtf_file) {
  
  GtfReader gtf_reader(gtf_file);
  GencodeGtfEntry entry;

  GenomicStepVector<FeatureVector<string>> gtf_features;
  unordered_set<string> chroms;

  cout << "prerocessing" << endl;
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
    cout << *chrom_it << endl;

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
  features.insert("CDS");
  features.insert("UTR");
  features.insert("intron");
}
  
void 
AlignedOverlapGenomicFeature::add_bed_features(const std::string& bed_file, 
                                               const std::string feature_name) {

  BedReader bed_reader(bed_file);
  GenomicRegion bed_region;
  while (bed_reader.read_bed3_line(bed_region)) {
    genomic_features.add(bed_region, FeatureVector<string>{feature_name});
  }

  features.insert(feature_name);
}
