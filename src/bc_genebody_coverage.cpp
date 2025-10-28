/*
* bc_genebody_coverage.cpp: get the coverage stats over the genebody
*   region
*
* Copyright (C) 2025 Rishvanth Prabakar
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
#include "unistd.h"
#include <unordered_map>
#include <unordered_set>

#include "gcatlib/Metagene.hpp"
#include "gcatlib/SamEntry.hpp"
#include "gcatlib/SamReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
using std::unordered_set;


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
      << "\t-n number of divisions to for a region [default: 100]" << endl
      << "\t-c cell barcode tag in SAM/BAM file [default: CB]" << endl
      << "\t-u UMI tag in SAM/BAM file [default: \"\"]" << endl
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all the flags are present [default: 0]" << endl
      << "\t-F only include if none of the flags are present [default: 2308]"
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

    string bc_tag = "CB";
    string umi_tag = "";
    
    size_t min_mapq = 0;
    size_t include_all = 0;
    size_t include_none = 0x0904;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:r:n:o:c:u:q:f:F:v")) != -1) {
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
      else if (opt == 'c')
        bc_tag = optarg;
      else if (opt == 'u')
        umi_tag = optarg;
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
 

    // process the barcodes
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;

    // index into the count matrix
    unordered_map<string, size_t> bc_index;
    size_t bc_counter = 0;
    vector<string> bc_metadata;
    
    std::ifstream bc_in(bc_file);
    string line;
    while (getline(bc_in, line)) {
      vector<string> tokens;
      split_string(line, tokens, '\t');
      bc_index[tokens[0]] = bc_counter++;
      bc_metadata.push_back(tokens[0]);
    }
    bc_in.close();

    // check if umi needs to be prcoressed
    bool process_umi = false;
    if (!umi_tag.empty()) {
      process_umi = true;
    }
    // set to keep track of bc+umi pairs
    // will be used only if process_umi is true
    unordered_set<string> seen_cb_umi;
    


    // process the regions
    if (VERBOSE)
      cerr << "[CREATING METAGENE OF REGIONS]" << endl;
    Metagene genebody_metagene(regions_file, n_divisions);
    
    // create the count matrix
    if (VERBOSE) {
      cerr << "[INITIALIZING COUNT MATRIX]" << endl;
      cerr << "\tNumber of barcodes: " << bc_counter << endl;
      cerr << "\tNumber of genes:" << genebody_metagene.get_n_features() 
              << endl;
    }

    vector<vector<size_t>> genebody_mat(bc_counter, 
                                        vector<size_t>(n_divisions, 0));

    // initialize sam/bam reader
    if (VERBOSE)
      cerr << "[INITIALIZING SAM/BAM READER]" << endl;
    
    SamReader sam_reader(aln_file);
    SamEntry sam_entry;


    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;

    size_t aln_count = 0;
    while (sam_reader.read_sam_line(sam_entry)) {
      ++aln_count;
      if (VERBOSE) {
        if (!(aln_count % 1000000)) {
          cerr << "\tprocessed " << aln_count << " fragments" << endl;
        }
      }

      // make sure the read passes QC
      if (sam_entry.mapq >= min_mapq &&
          SamFlags::is_all_set(sam_entry.flag, include_all) &&
          !SamFlags::is_any_set(sam_entry.flag, include_none)) {


        // get cell barcode
        bool process_cell = false;
        string cell_bc;
        SamTags::get_tag(sam_entry.tags, bc_tag, cell_bc);          
        unordered_map<string, size_t>::iterator bc_it;
        bc_it = bc_index.find(cell_bc);
        if (bc_it != bc_index.end()) {
          process_cell = true;
        }

        // if needed, get cell umi
        string cell_umi;
        if (process_cell && process_umi) {
          SamTags::get_tag(sam_entry.tags, umi_tag, cell_umi);
          // check if the BC+UMI combo has already been processed
          if (seen_cb_umi.find(cell_bc + cell_umi) != seen_cb_umi.end()) {
            // has already been seen
            process_cell = false;
          }
          else {
            // has not been seen, so add to list
            seen_cb_umi.insert(cell_bc + cell_umi);
          }
        } 
    
        // get the gene body regions
        if (process_cell) {

          // query the metagene index
          GenomicRegion query;
          query.name = sam_entry.rname;
          query.start = sam_entry.pos - 1;
          query.end = sam_entry.pos;

          vector<string> feature;
          vector<size_t> first;
          vector<size_t> last;

          genebody_metagene.at_ends(query, feature, first, last);
          for (size_t i = 0; i < feature.size(); ++i) {
            ++genebody_mat[bc_it->second][first[i]];
          }

        } 

      }
    }


    // write output
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;
   
    std::ofstream cov_file(out_prefix + "_genebody_coverage.txt"); 
    // write header
    cov_file << "barcode";
    for (size_t i = 0; i < n_divisions; ++i) {
      cov_file << "\t" << i;
    }
    cov_file << endl;

    // write the content
    for (size_t i = 0; i < bc_metadata.size(); ++i) {
      cov_file << bc_metadata[i];
      // sanity check
      if (bc_index.find(bc_metadata[i])->second != i) {
        throw std::runtime_error("bc metadata does not match index");
      }
      for (size_t j = 0; j < n_divisions; ++j) {
        cov_file << "\t" << genebody_mat[i][j];
      }
      cov_file << endl;
    }
    cov_file.close();

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
