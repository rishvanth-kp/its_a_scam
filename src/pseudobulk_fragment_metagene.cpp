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

    char bc_delim = ':';
    uint8_t bc_col = 7;
    string bc_tag;

    size_t min_frag_len = 0;
    size_t max_frag_len = 1024;

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:r:o:d:c:t:m:M:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b') 
        bc_file = optarg;
      else if (opt == 'r')
        regions_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
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
    size_t group_counter = 0;
    size_t bc_counter = 0;

    std::ifstream bc_in(bc_file);
    string line;
    while (getline(bc_in, line)) {
      vector<string> tokens;
      split_string(line, tokens, '\t');

      // keep track of the group for each barcode
      bc_group[tokens[0]] = tokens[1];
      ++bc_counter;

      // if a new group is encountered, create a new index
      auto it = group_index.find(tokens[1]);
      if (it == group_index.end()) {
        group_index[tokens[1]] = group_counter++;
      }
    }
    bc_in.close();
    
    if (VERBOSE) {
      cerr << "\tNumber of barcodes: " << bc_counter << endl;
      cerr << "\tNumber of groups: " << group_counter << endl;
    }
     


    // process the regions
    if (VERBOSE)
      cerr << "[CREATING METAGENE OF REGIONS]" << endl;
    Metagene metagene(regions_file);
   
    // create count matrix
    if (VERBOSE) 
      cerr << "[INITIALIZING COUNT MATRIX]" << endl;
    size_t n_features = 19;
    vector<vector<size_t>> metagene_matrix(group_counter * n_features,
                                            vector<size_t>(101, 0)); 
   


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
          string cell_group = bc_it->second;

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
 

          // query the metagene
          vector<string> regions;
          vector<size_t> first, last;
          GenomicRegion entry_in;
      
          entry_in.name = entry1.rname;
          entry_in.start = frag_start;
          entry_in.end = frag_end;

          metagene.at(entry_in, regions, first, last);
          for (size_t i = 0; i < regions.size(); ++i) {
            for (size_t j = first[i]; j <= last[i]; ++j) {
              
            }
          }
 
        }
      }
    }



    // write output  



  }
  catch (const std::exception  &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
