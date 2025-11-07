/*
* Copyright (C) 2023 Rishvanth Prabakar
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

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>

#include "gcatlib/SamEntry.hpp"
#include "gcatlib/SamReader.hpp"
#include "gcatlib/BedReader.hpp"
#include "gcatlib/GenomicRegion.hpp"
#include "gcatlib/FeatureVector.hpp"
#include "gcatlib/GenomicStepVector.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::unordered_map;
using std::unordered_set;


struct TssMetadata {
  string chrom;
  size_t start;
  size_t end;
  string tss_id;
  bool strand;
};

static void
split_string (const string &in, vector<string> &tokens,
              const char delim = ':') {

  tokens.clear();
  size_t start = 0;
  size_t end = in.find(delim);
  while (end != string::npos) {
    tokens.push_back(in.substr(start, end - start));
    start = ++end;
    end = in.find(delim, start);
  }
  tokens.push_back(in.substr(start));
}

static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-r TSS region bed file [required]" << endl
      << "\t-s side distance to add to the ends of the TSS [default: 1000]" 
          << endl
      << "\t-b barcode list (if empty, treated as bulk sample)" << endl
      << "\t-o outfile prefix [required]" << endl
      << "\t-n min fragment length [default: 10]" << endl
      << "\t-x max fragment length [default: 1000]" << endl
      << "\t-d name split delimeter " 
          << "[default: \":\"; ignored if -t is provided]" << endl
      << "\t-c barcode field in name " 
          << "[default: 7 (0 based); ignored if -t is provided]" << endl
      << "\t-t barcode tag in SAM/BAM file [default \"\"]" << endl 
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all the flags are present [default: 3]" << endl
      << "\t-F only include if none of the flags are present [default: 3340]" 
          << endl
      << "\t-v verbose [default: false]" << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {

    string aln_file;
    string tss_file;
    string bc_file; 
    string out_prefix;

    size_t side_dist = 1000;

    size_t min_frag_len = 1;
    size_t max_frag_len = 1000;

    char bc_delim = ':';
    uint8_t bc_col = 7;
    string bc_tag;

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;    

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:r:s:b:o:n:x:d:c:t:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'r')
        tss_file = optarg;
      else if (opt == 's')
        side_dist = std::stoi(optarg);
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'n')
        min_frag_len = std::stoi(optarg);
      else if (opt == 'x')
        max_frag_len = std::stoi(optarg);
      else if (opt == 'd')
        bc_delim = optarg[0];
      else if (opt == 'c')
        bc_col = std::stoi(optarg);
      else if (opt == 't')
        bc_tag = optarg;
      else if (opt == 'q')
        min_mapq = std::stoi(optarg);
      else if (opt == 'f')
        include_all = std::stoi(optarg);
      else if (opt == 'F')
        include_none = std::stoi(optarg);
      else if (opt == 'v')
        VERBOSE = true;
      else 
        throw std::runtime_error(print_usage(argv[0]));
    } 

    if (aln_file.empty() || tss_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }


    // process barcodes
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;
    
    // to index the TSS enhancement martix
    unordered_map<string, size_t> bc_index;
    size_t bc_counter = 0;
    vector<string> bc_metadata;    
    bool bulk_sample = false;  

    if (bc_file.empty()) {
      // if a barocde file is not proivded, treat it as a bulk sample
      // all the reads in the sam file will be considerd as one sample
      bulk_sample = true;
      bc_index[out_prefix] = bc_counter++;
      bc_metadata.push_back(out_prefix);
    }
    else {
      // if a barcode file is provided, each barcode is a separate sample
      std::ifstream bc_in(bc_file);
      string line;
      while (getline(bc_in, line)) {
        vector<string> tokens;
        split_string(line, tokens, '\t');
        bc_index[tokens[0]] = bc_counter++;
        bc_metadata.push_back(tokens[0]);
      } 
      bc_in.close();
    }

    // create matrix to store the coverage around TSS for each sample
    vector<vector<size_t>> tss_coverage(bc_counter,
                                        vector<size_t>(2*side_dist + 1, 0));


    // check if the cell barcode is in the name or a tag
    // if a tag parameter is not provided, then the barcode is in the name
    bool bc_in_tag = false;
    if (!bc_tag.empty()) {
      bc_in_tag = true;
    }

    if (VERBOSE) {
      if (bulk_sample) {
        cerr << "\tBulk sample: " << out_prefix << endl;
      }
      else {
        cerr << "\tNumber of barocdes: " << bc_counter << endl;
      }
    }
      

    // process TSS file
    if (VERBOSE)
      cerr << "[PROCESSING TSS REGIONS]" << endl;

    GenomicStepVector<FeatureVector<string>> tss;
    unordered_map<string, TssMetadata> tss_metadata;

    size_t tss_counter = 0;
  
    BedReader tss_reader(tss_file);
    GenomicRegion bed_region;
    vector<string> bed_fields;
    while (tss_reader.read_bed_line(bed_region, bed_fields)) {

      // store the tss regrion
      // storing just the TSS name. The distance from the tss can be 
      // determined from the TSS start stored in the metadata.
      tss.add(bed_region.name, bed_region.start - side_dist, 
              bed_region.end + side_dist + 1, 
              FeatureVector<string>(bed_fields[0])); 
    
 
      // store TSS metadata
      TssMetadata metadata;
      metadata.chrom = bed_region.name;
      metadata.start = bed_region.start;
      metadata.end = bed_region.end;
      metadata.tss_id = bed_fields[0];
      if (bed_fields[2][0] == '-')
        metadata.strand = false;
      else
        metadata.strand = true;
      tss_metadata[bed_fields[0]] = metadata;

      // keep track of number of tss
      ++tss_counter;
    }   

    if (VERBOSE) 
      cerr << "\tNumber of TSS: " << tss_counter << endl;


    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;
    
    SamReader sam_reader(aln_file);
    SamEntry entry1, entry2;
    
    size_t aln_count = 0;
    size_t aln_pass_count = 0;
    size_t in_bc_count = 0;
  
    while (sam_reader.read_pe_sam(entry1, entry2)) {
      ++aln_count;
      if (VERBOSE) {
        if (!(aln_count % 1000000)) {
          cerr << "\tprocessed " << aln_count << " pairs" << endl;
        }
      }

      if (entry1.mapq >= min_mapq && entry2.mapq >= min_mapq &&
          SamFlags::is_all_set(entry1.flag, include_all) &&
          SamFlags::is_all_set(entry2.flag, include_all) &&
          !SamFlags::is_any_set(entry1.flag, include_none) &&
          !SamFlags::is_any_set(entry2.flag, include_none)) {

        ++aln_pass_count;

        // get the cell barcode
        string cell_bc;
        if (bc_in_tag) {
          // read bc from the appropriate tag
          SamTags::get_tag(entry1.tags, bc_tag, cell_bc);
        }
        else {
          // parser the sam name to get the bc
          vector<string> tokens;
          split_string(entry1.qname, tokens, bc_delim);
          if (tokens.size() > bc_col) {
            cell_bc = tokens[bc_col];
          }
        }
      
        // check if the barcode is in list and get the barcode index
        bool bc_found = bulk_sample;
        size_t tss_row_index = 0;
        if (!bulk_sample) {
          unordered_map<string, size_t>::iterator bc_it;
          bc_it = bc_index.find(cell_bc);
          if (bc_it != bc_index.end()) {
            bc_found = true;
            tss_row_index = bc_it->second;
          }
        }
 
        // proceed if barcode is found
        if (bc_found) {
          ++in_bc_count;

          // find the start and end postion of the reads
          size_t e1_start = entry1.pos - 1;
          size_t e1_end = SamCigar::reference_end_pos(entry1);
          size_t e2_start = entry2.pos - 1;
          size_t e2_end = SamCigar::reference_end_pos(entry2);

          // find the start and postion of the fragment
          size_t frag_start = e1_start;
          if (e2_start < e1_start)
            frag_start = e2_start;

          size_t frag_end = e1_end;
          if (e2_end > e1_end) 
            frag_end = e2_end;

          size_t frag_len = frag_end - frag_start;
          if ((frag_len >= min_frag_len) && (frag_len <= max_frag_len)) {

            // check if the fragment is in a TSS region
            unordered_set<string> aligned_tss;
            vector<pair<GenomicRegion, FeatureVector<string>>> out;
            tss.at(GenomicRegion(entry1.rname, frag_start, frag_end), out); 
            for (auto it = out.begin(); it != out.end(); ++it) {
              for (size_t j = 0; j < it->second.size(); ++j) {
                aligned_tss.insert(it->second.at(j));
              }
            }

            for (auto it = aligned_tss.begin(); it != aligned_tss.end(); ++it) {
              TssMetadata metadata = tss_metadata[*it];
              for (size_t j = frag_start; j < frag_end; ++j) {
                int tss_offset = j - metadata.start + side_dist;
                if (!metadata.strand)
                  tss_offset = (2*side_dist) - tss_offset;
                if (tss_offset >= 0 && tss_offset <= 2*side_dist) {
                  ++tss_coverage[tss_row_index][tss_offset];
                } 
              }
            }

          }


        }
 
      }
    }  

    // normalize tss enrichment matrix
    vector<size_t> norm_factor(bc_metadata.size(), 0);
    const size_t norm_len = 100;
    for (size_t i = 0; i < bc_metadata.size(); ++i) {
      for (size_t j = 0; j < norm_len; ++j) {
        norm_factor[i] += tss_coverage[i][j];
        norm_factor[i] += tss_coverage[i][(2*side_dist) - j];
      }
    }

    // write output
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;

    std::ofstream out_en(out_prefix + "_tss_enrichmet.txt");
    std::ofstream out_counts(out_prefix + "_tss_raw_counts.txt");
    std::ofstream out_norm(out_prefix + "_tss_norm_factor.txt");
    for (size_t i = 0; i < bc_metadata.size(); ++i) {
      out_en << bc_metadata[i];
      out_counts << bc_metadata[i];
      out_norm << bc_metadata[i] << "\t" << norm_factor[i];
      if (norm_factor[i] > 0) {
        for (size_t j = 0; j < tss_coverage[i].size(); ++j) {
          out_en << "\t" << static_cast<float>(tss_coverage[i][j]) /
                            static_cast<float>(norm_factor[i]);
          out_counts << "\t" << tss_coverage[i][j];
        }
      }
      else {
        for (size_t j = 0; j < tss_coverage[i].size(); ++j) {
          out_en << "\t" <<  0;
          out_counts << "\t" << tss_coverage[i][j];
        }
      }
      out_en << endl;
      out_counts << endl;
      out_norm << endl;
    }
    out_en.close();
    out_counts.close();
    out_norm.close();

  } 
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl; 
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
