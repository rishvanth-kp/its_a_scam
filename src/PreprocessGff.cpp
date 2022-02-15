/*
* PreprocessGff: parse features from GFF file 
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

#include "PreprocessGff.hpp"

using std::cout;
using std::cerr;
using std::endl;

PreprocessGff::PreprocessGff(const string &chrom_size_file, 
                             bool verbose) {
  VERBOSE = verbose;
  read_chrom_sizes(chrom_size_file);
}

PreprocessGff::~PreprocessGff() {
};

void
PreprocessGff::read_chrom_sizes(const string &chrom_size_file) {

  if (VERBOSE)
    cerr << "READING CHROM SIZES FROM " << chrom_size_file << endl;
  
  std::ifstream in(chrom_size_file);
  if (!in) 
    throw std::runtime_error("cannot open " + chrom_size_file);     

  string line;
  while (getline(in, line)) {
    size_t split = line.find('\t');
    string chrom = line.substr(0, split);
    size_t sz = std::stoi(line.substr(split));
    chrom_sizes[chrom] = sz; 
  }

  in.close();
  
  if (VERBOSE)
    cerr << chrom_sizes.size() << " CHROM SIZES INSERTED" << endl;
}

struct GenomeFeatures {
  uint8_t gene_count = 0;
  string gene_type;
  bool exon;
  bool cds;
  bool utr;
};

static void
add_feature (const GencodeGtfEntry &entry, 
             vector<GenomeFeatures> &chr_features) {

  if (entry.feature == "gene") {
    for (size_t i = entry.start; i <= entry.end; ++i) {
      ++chr_features[i].gene_count;
      chr_features[i].gene_type = entry.gene_type;
    }
  }

  else if (entry.feature == "exon") {
    for (size_t i = entry.start; i <= entry.end; ++i) {
      chr_features[i].exon = true;
    }
  }

  else if (entry.feature == "CDS") {
    for (size_t i = entry.start; i <= entry.end; ++i) {
      chr_features[i].cds = true;
    }
  }

  else if (entry.feature == "UTR") {
    for (size_t i = entry.start; i <= entry.end; ++i) {
      chr_features[i].utr = true;
    }
  }

}


static void
process_chr_features (const string &chrom,
                      const vector<GenomeFeatures> &chr_features) {
  
  size_t i = 0;
  size_t start = 0;
  string old_feature;
  if (chr_features[0].gene_count == 0) { 
    old_feature = "intergenic";
  }
  else {
    old_feature = chr_features[0].gene_type; 
  }
  string new_feature = old_feature;

  while (i < chr_features.size()) {
    if (chr_features[i].gene_count == 0) {
      new_feature = "intergenic";
    }

    else if (chr_features[i].gene_count > 1) {
      new_feature = "ambiguous";
    }
  
    else if (chr_features[i].gene_type == "protein_coding") {
      if (!chr_features[i].exon) {
        new_feature = "pc_intron";
      }
      else if (chr_features[i].cds & !chr_features[i].utr) {
        new_feature = "pc_cds";
      }
      else if (!chr_features[i].cds & chr_features[i].utr) {
        new_feature = "pc_utr";
      }
      else if (chr_features[i].cds & chr_features[i].utr) {
        new_feature = "pc_both";
      }
      else if (!chr_features[i].cds & !chr_features[i].utr) {
        new_feature = "pc_none";
      }
    }
    
    else {
      new_feature = chr_features[i].gene_type;
    }

    if (old_feature != new_feature) {
      cout << chrom << "\t"
           << start << "\t"
           << i << "\t"
           << old_feature << endl;

      start = i;
      old_feature = new_feature;
    }

  ++i;
  }
  cout << chrom << "\t"
       << start << "\t"
       << i << "\t"
       << old_feature << endl;

}


void 
PreprocessGff::parse_genome_features(const string &gff_file) {

  if (VERBOSE)
    cerr << "READING GFF FILE " << gff_file << endl;
  
  GtfReader gff_reader(gff_file); 
  GencodeGtfEntry entry;
  string curr_chr = "";
  vector<GenomeFeatures> chr_features;
  while (gff_reader.read_gencode_gtf_line(entry)) {
    if (entry.name != curr_chr) {
      if (VERBOSE) {
        cerr << "PROCESSING CHROM: " << entry.name << endl;
      }

      if (!curr_chr.empty())
        process_chr_features(curr_chr, chr_features);

      curr_chr = entry.name;
      chr_features.clear();
      chr_features.resize(chrom_sizes[curr_chr]);
      // cout << chr_features.size() << endl;  
    }
    else {
      add_feature(entry, chr_features);  
    }
  } 
  process_chr_features(curr_chr, chr_features);

}
