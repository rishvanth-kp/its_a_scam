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
#include "gcatlib/SamReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
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
      << "\t-b barcode list [required]" << endl
      << "\t-o outfile prefix [required]" << endl
      << "\t-d name split delimeter " 
          << "[default: \":\"; ignored if -t is provided]" << endl
      << "\t-c barcode field in name " 
          << "[default: 7 (0 based); ignored if -t is provided]" << endl
      << "\t-t barcode tag in SAM file [default: \"\"]" << endl
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-v verbose [default: false]" << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {
  
    string aln_file;
    string bc_file;
    string out_prefix;

    char bc_delim = ':';
    uint8_t bc_col = 7;
    string bc_tag;

    size_t min_mapq = 0;
    
    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:o:d:c:t:q:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'd')
        bc_delim = optarg[0];
      else if (opt == 'c')
        bc_col = std::stoi(optarg);
      else if (opt == 't')
        bc_tag = optarg;
      else if (opt == 'q')
        min_mapq = std::stoi(optarg);
      else if (opt == 'v')
        VERBOSE = true;
      else 
        throw std::runtime_error(print_usage(argv[0]));
    }
 
    if (aln_file.empty() || bc_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    
    // process barcodes
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;  
    
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

    if (VERBOSE) {
      cerr << "\tNumber of barcodes: " << bc_counter << endl;
    }

        
    // initialize sam/bam reader
    if (VERBOSE)
      cerr << "[INITIALIZING SAM/BAM READER]" << endl;

    SamReader reader(aln_file);
    SamEntry entry;

    // initialize array for storing flagstat
    // 1 row for each barcode
    // column info: 
    //  0: primary
    //  1: secondary
    //  2: supplementary
    //  3: duplicates
    //  4: primary duplicates
    //  5: mapped
    //  6: primary mapped
    //  7: paired in sequencing
    //  8: read 1
    //  9: read 2
    //  10: properly paired
    //  11: with itself and mate mapped
    //  12: singletons

    vector<vector<size_t>> bc_flagstat(bc_counter, vector<size_t>(13, 0));

    
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;
    
    size_t aln_count = 0;
  
    while (reader.read_sam_line(entry)) {
      ++ aln_count;
      if (VERBOSE) {
        if (!(aln_count % 1000000)) {
          cerr << "\tprocessed " << aln_count << " alignments" << endl;
        }
      }

      if (entry.mapq >= min_mapq) {
        
        // get cell barcode, eiter from the tag or the name
        string cell_bc;
        if (!bc_tag.empty()) {
          // read bc from the appropriate tag
          SamTags::get_tag(entry.tags, bc_tag, cell_bc);
        }
        else {
          vector<string> tokens;
          split_string(entry.qname, tokens, bc_delim);
          cell_bc = tokens[bc_col];
        }

        // check if barcode needs to be processed
        unordered_map<string, size_t>::iterator bc_it;
        bc_it = bc_index.find(cell_bc);
        if (bc_it != bc_index.end()){
          size_t entry_bc_index = bc_it->second;
          //  0: primary
          if (!is_set(entry.flag, SamFlags::Flag::not_primary_aln) && 
              !is_set(entry.flag, SamFlags::Flag::supplementary_aln)) {
            ++bc_flagstat[entry_bc_index][0];
          }
          //  1: secondary
          if (is_set(entry.flag, SamFlags::Flag::not_primary_aln)) {
            ++bc_flagstat[entry_bc_index][1];
          }
          //  2: supplementary
          if (is_set(entry.flag, SamFlags::Flag::supplementary_aln)) {
            ++bc_flagstat[entry_bc_index][2];
          }
          //  3: duplicates
          if (is_set(entry.flag, SamFlags::Flag::pcr_duplicate)) {
            ++bc_flagstat[entry_bc_index][3];
          }
          //  4: primary duplicates
          if (is_set(entry.flag, SamFlags::Flag::pcr_duplicate) &&
              !is_set(entry.flag, SamFlags::Flag::not_primary_aln) && 
              !is_set(entry.flag, SamFlags::Flag::supplementary_aln)) {
            ++bc_flagstat[entry_bc_index][4];
          }
          //  5: mapped
          if (!is_set(entry.flag, SamFlags::Flag::read_unmapped)) {
            ++bc_flagstat[entry_bc_index][5];
          }
          //  6: primary mapped
          if (!is_set(entry.flag, SamFlags::Flag::read_unmapped) &&
              !is_set(entry.flag, SamFlags::Flag::not_primary_aln) &&
              !is_set(entry.flag, SamFlags::Flag::supplementary_aln)) {
            ++bc_flagstat[entry_bc_index][6];
          }
          //  7: paired in sequencing
          if (is_set(entry.flag, SamFlags::Flag::read_paired)) {
            ++bc_flagstat[entry_bc_index][7];
          }
          //  8: read 1
          if (is_set(entry.flag, SamFlags::Flag::read_paired) &&
              is_set(entry.flag, SamFlags::Flag::first_in_pair)) {
            ++bc_flagstat[entry_bc_index][8];
          }
          //  9: read 2
          if (is_set(entry.flag, SamFlags::Flag::read_paired) &&
              is_set(entry.flag, SamFlags::Flag::first_in_pair)) {
            ++bc_flagstat[entry_bc_index][9];
          }
          //  10: properly paired
          if (is_set(entry.flag, SamFlags::Flag::read_paired) &&
              is_set(entry.flag, SamFlags::Flag::proper_pair) &&
              !is_set(entry.flag, SamFlags::Flag::read_unmapped)) { 
            ++bc_flagstat[entry_bc_index][10];
          }
          //  11: with itself and mate mapped
          if (is_set(entry.flag, SamFlags::Flag::read_paired) &&
              !is_set(entry.flag, SamFlags::Flag::read_unmapped) &&
              !is_set(entry.flag, SamFlags::Flag::mate_unmapped)) {
            ++bc_flagstat[entry_bc_index][11];
          }
          //  12: singletons
          if (is_set(entry.flag, SamFlags::Flag::read_paired) &&
              !is_set(entry.flag, SamFlags::Flag::read_unmapped) &&
              is_set(entry.flag, SamFlags::Flag::mate_unmapped)) {
            ++bc_flagstat[entry_bc_index][12];
          }


        }
      }
    }

    // write output
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;
    
    std::ofstream out(out_prefix + "_bc_flagstat.txt");

    // write header
    out << "barcode" 
        << "\tprimary"
        << "\tsecondary"
        << "\tsupplementary"
        << "\tduplicates"
        << "\tprimary_duplicates"
        << "\tmapped"
        << "\tprimary_mapped"
        << "\tpaired_in_sequencing"
        << "\tread_1"
        << "\tread_2"
        << "\tproperly_paired"
        << "\twith_itself_and_mate_mapped"
        << "\tsingletons" << endl;

    // write flagstat for each barcode
    for (size_t i = 0; i < bc_flagstat.size(); ++i) {
      out << bc_metadata[i];
      for (size_t j = 0; j < bc_flagstat[i].size(); ++j) {
        out << "\t" << bc_flagstat[i][j];
      }
      out << endl;
    }

    out.close();

  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
