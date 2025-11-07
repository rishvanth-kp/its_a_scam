/*
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
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <unordered_set>

#include "gcatlib/SamEntry.hpp"
#include "gcatlib/SamReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
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
      << "\t-b barcode list file [default = \"\"]" << endl
      << "\t-o out file prefix [default = \"\"]" << endl
      << "\t-l Number of bases to shift left [default: 4]" << endl
      << "\t-r Number of bases to shift right [default: -5]" << endl
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
    string bc_file;
    string out_prefix;
    
    int shift_left = 4;
    int shift_right = -5;

    char bc_delim = ':';
    uint8_t bc_col = 7; 
    string bc_tag;

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:o:l:r:d:c:t:q:f:F:v")) != -1) {
      if (opt == 'a') 
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'l')
        shift_left = std::stof(optarg);
      else if (opt == 'r')
        shift_right = std::stof(optarg);
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
    

    if (aln_file.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    // if a barcode file is provided, keep track of them
    unordered_set<string> barcodes;
    if (!bc_file.empty()) {
      if (VERBOSE) 
        cerr << "[PROCESSING BARCODES]" << endl;

      std::ifstream bc_in(bc_file);
      string line;
      while (getline(bc_in, line)) {
        vector<string> tokens;
        split_string(line, tokens, '\t');
        barcodes.insert(tokens[0]);
      }
      bc_in.close();
    }

    // initialize sam reader
    if (VERBOSE)
      cerr << "[INITIALIZING SAM/BAM READER]" << endl;

    SamReader sam_reader(aln_file);
    SamEntry entry1, entry2;



    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;
    
    size_t aln_count = 0;

    string outfile_name = "fragments.tsv";
    if (!out_prefix.empty())
      outfile_name = out_prefix + "_fragments.tsv";

    std::ofstream frag_file(outfile_name);

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

        // get cell barcode, either from a tag or the name
        string cell_bc;
        if (!bc_tag.empty()) {
          // read barcode from the tag
          SamTags::get_tag(entry1.tags, bc_tag, cell_bc);
        }
        else {
          // parse the name to get the bc
          vector<string> tokens;
          split_string(entry1.qname, tokens, bc_delim);
          cell_bc = tokens[bc_col];
        }

        // check if the barcode needs to be added
        unordered_set<string>::const_iterator it;
        it = barcodes.find(cell_bc);
        if ((barcodes.empty() && !cell_bc.empty()) || (it != barcodes.end())) {

          // find the start and end location of each read
          size_t e1_start = entry1.pos - 1;
          size_t e1_end = SamCigar::reference_end_pos(entry1);
          size_t e2_start = entry2.pos - 1;
          size_t e2_end = SamCigar::reference_end_pos(entry2);

          // find the start and end location of the fragment
          size_t frag_start = e1_start;
          if (e2_start < e1_start)
            frag_start = e2_start;

          size_t frag_end = e1_end;
          if (e2_end >= e1_end) 
            frag_end = e2_end;

          // write to fragments file
          frag_file << entry1.rname 
                    << "\t" << frag_start + shift_left
                    << "\t" << frag_end + shift_right
                    << "\t" << cell_bc
                    << "\t1" << endl; 

        }

      }

    }

    frag_file.close();

  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
