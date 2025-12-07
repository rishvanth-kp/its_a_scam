/*
* Copyright (C) 2024 Rishvanth Prabakar
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

#include "gcatlib/SamEntry.hpp"
#include "gcatlib/Metagene.hpp"
#include "gcatlib/SamReader.hpp"

using std::cout;
using std::cerr;
using std::string;
using std::vector;
using std::unordered_map;

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
      << "\t-b barcode list file [required]" << endl
      << "\t-r regions bed file [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-n number of divisions for a region [default: 100]" << endl
      << "\t-m min. fragment length [default: 0]" << endl
      << "\t-M max. fragment length [default: 1024]" << endl
      << "\t-d name split delimeter"
          << "[default: \":\"]; ignored if -t is provided" << endl
      << "\t-c barcode field in name"
          << "[default: 7 (0 based); ignored if -t is provided]" << endl
      << "\t-t barcode tag in SAM file [default \"\"]" << endl
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

    // parse args
    string aln_file;
    string regions_file;
    string bc_file;

    string out_prefix;

    size_t n_divisions = 100;

    char bc_delim = ':';
    size_t bc_col = 7;
    string bc_tag;

    size_t min_frag_len = 0;
    size_t max_frag_len = 1024;

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:r:o:n:d:c:t:m:M:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'r')
        regions_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'n')
        n_divisions = std::stoi(optarg);
      else if (opt == 'd')
        bc_delim = optarg[0];
      else if (opt == 'c')
        bc_col = std::stoi(optarg);
      else if (opt == 't')
        bc_tag = optarg;
      else if (opt == 'm')
        min_frag_len = std::stoi(optarg);
      else if (opt == 'M')
        max_frag_len = std::stoi(optarg);
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

    if (aln_file.empty() || bc_file.empty() || regions_file.empty() ||
        out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }



    // process the barcodes and pseudobulk info
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;

    unordered_map<string, string> bc_group;
    unordered_map<string, size_t> group_index;
    vector<string> group_names;
    size_t n_groups = 0;
    size_t n_bcs = 0;

    std::ifstream bc_in(bc_file);
    string line;
    while (getline(bc_in, line)) {
      vector<string> tokens;
      split_string(line, tokens, '\t');

      // keep track of the group for each barcode
      bc_group[tokens[0]] = tokens[1];
      ++n_bcs;

      // if a new group is encountered, create a new index
      auto it = group_index.find(tokens[1]);
      if (it == group_index.end()) {
        group_index[tokens[1]] = n_groups++;
        group_names.push_back(tokens[1]);
      }
    }
    bc_in.close();




    // process the regions
    if (VERBOSE)
      cerr << "[CREATING METAGENE OF REGIONS]" << endl;
    Metagene metagene(regions_file, n_divisions);
    // get the names of the regions and compute a region index
    vector<string> feature_names;
    metagene.get_feature_names(feature_names);
    unordered_map<string, size_t> feature_index;
    for (size_t i = 0; i < feature_names.size(); ++i) {
      feature_index[feature_names[i]] = i;
    }


    // create count matrix
    if (VERBOSE)
      cerr << "[INITIALIZING COUNT MATRIX]" << endl;
    size_t n_features = metagene.get_n_features();
    vector<vector<size_t>> metagene_matrix(n_groups * n_features,
                                            vector<size_t>(n_divisions, 0));

    if (feature_names.size() != n_features) {
      throw std::runtime_error("feature counts do not match");
    }

    // create vectors to store normalization factors
    // For tracking all the reads 
    vector<size_t> group_frag_counts(n_groups, 0);
    vector<size_t> group_base_counts(n_groups, 0); 
    // For tracking reads that pass the fragment length criteria
    vector<size_t> passed_frag_counts(n_groups, 0);
    vector<size_t> passed_base_counts(n_groups, 0);
    // For tacking reads that fall in the desired region
    vector<size_t> region_frag_counts(n_groups, 0);
    vector<size_t> region_base_counts(n_groups, 0);
    

    if (VERBOSE) {
      cerr << "\tNumber of barcodes: " << n_bcs << endl;
      cerr << "\tNumber of groups: " << n_groups << endl;
      cerr << "\tNumber of features: " << n_features << endl;
    }


    // initalize sam reader
    if (VERBOSE)
      cerr << "[INITIALIZING SAM/BAM READER]" << endl;

    SamReader reader(aln_file);
    SamEntry entry1, entry2;


    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;

    size_t aln_count = 0;
    while (reader.read_pe_sam(entry1, entry2)) {
      ++aln_count;
      if (VERBOSE) {
        if (!(aln_count % 1000000)) {
          cerr << "\tprocessed " << aln_count << " fragments" << endl;
        }
      }


      // make sure the fragment passed qc
      if (entry1.mapq >= min_mapq && entry2.mapq >= min_mapq &&
          SamFlags::is_all_set(entry1.flag, include_all) &&
          SamFlags::is_all_set(entry2.flag, include_all) &&
          !SamFlags::is_any_set(entry1.flag, include_none) &&
          !SamFlags::is_any_set(entry2.flag, include_none)) {

        // get the cell barcode, either from a tag or the name
        string cell_bc;
        if (!bc_tag.empty()) {
          // read bc from the appropriate tag
          SamTags::get_tag(entry1.tags, bc_tag, cell_bc);
        }
        else {
          // parse the name to get the bc
          vector<string> tokens;
          split_string(entry1.qname, tokens, bc_delim);
          cell_bc = tokens[bc_col];
        }

        // check if the barcode needs to be processed
        unordered_map<string, string>::iterator bc_it;
        bc_it = bc_group.find(cell_bc);
        if (bc_it != bc_group.end()) {

          // find the start and end location of each read
          size_t e1_start = entry1.pos - 1;
          size_t e1_end = SamCigar::reference_end_pos(entry1);
          size_t e2_start = entry2.pos - 1;
          size_t e2_end = SamCigar::reference_end_pos(entry2);

          // fing the start and end location of the fragment
          size_t frag_start = e1_start;
          if (e2_start < e1_start)
            frag_start = e2_start;

          size_t frag_end = e1_end;
          if (e2_end >= e1_end)
            frag_end = e2_end;

          size_t frag_len = frag_end - frag_start;


          string cell_group = bc_it->second;
          size_t cell_group_index = group_index[cell_group];
         
          // for group normalization
          group_frag_counts[cell_group_index] += 1;
          group_base_counts[cell_group_index] += frag_len; 

          // query the metagene, if it is within the desirefd frag len
          if (frag_len >= min_frag_len && frag_len <= max_frag_len) {
          
            // for normalization of reads the pass the fragment length
            passed_frag_counts[cell_group_index] += 1;
            passed_base_counts[cell_group_index] += frag_len; 

            GenomicRegion entry_in;

            vector<pair<GenomicRegion,
                      FeatureVector<pair<string, size_t>>>> out;

            entry_in.name = entry1.rname;
            entry_in.start = frag_start;
            entry_in.end = frag_end;

            bool in_region = false;
            size_t region_size = 0;

            metagene.at(entry_in, out);
            for(size_t i = 0; i < out.size(); ++i) {
              // for normalization
              in_region = true;
              region_size += (out[i].first.end - out[i].first.start);

              for (size_t j = 0; j < out[i].second.size(); ++j) {

                size_t cell_feature_index =
                  feature_index[out[i].second.at(j).first];
                size_t row_index =
                  (cell_group_index * n_features) + cell_feature_index;

                  metagene_matrix[row_index][out[i].second.at(j).second] +=
                    (out[i].first.end - out[i].first.start);

              }
            }
            
            // for normalization of reads that are in the region
            if (in_region) {
              region_frag_counts[cell_group_index] += 1;
              region_base_counts[cell_group_index] += region_size; 
            }
            /*
            for (size_t i = 0; i < regions.size(); ++i) {
              size_t cell_feature_index = feature_index[regions[i]];
              size_t row_index = (cell_group_index * n_features) + cell_feature_index;
              for (size_t j = first[i]; j <= last[i]; ++j) {
                ++metagene_matrix[row_index][j];
              }
            }
            */

          }


        }
      }
    }


    // write output
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;

    std::ofstream counts_file(out_prefix + "_metagene.txt");

    // write the header
    counts_file << "group\tfeature";
    for (size_t i = 0; i < n_divisions; ++i) {
      counts_file << "\t" << i;
    }
    counts_file << endl;

    for (size_t i = 0; i < n_groups; ++i) {
      for (size_t j = 0; j < n_features; ++j) {
        counts_file << group_names[i] << "\t" << feature_names[j];
        for (size_t k = 0; k < n_divisions; ++k) {
          counts_file << "\t"
                      << metagene_matrix[(i * n_features) + j][k];
        }
        counts_file << endl;
      }
    }

    counts_file.close();


    // write the normalization
    std::ofstream norm_file(out_prefix + "_metagene_normalization.txt");
    
    // write the header
    norm_file << "group"
              << "\tgroup_frag_count\tgroup_base_count"
              << "\tpassed_frag_count\tpassed_base_count"
              << "\tregion_frag_count\tregion_base_count" << endl; 

    // write the group counts
    for (size_t i = 0; i < n_groups; ++i) {
      norm_file << group_names[i] 
                << "\t" << group_frag_counts[i]
                << "\t" << group_base_counts[i]
                << "\t" << passed_frag_counts[i]
                << "\t" << passed_base_counts[i]
                << "\t" << region_frag_counts[i]
                << "\t" << region_base_counts[i]
                << endl;
    }

    norm_file.close();

  }
  catch (const std::exception  &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
