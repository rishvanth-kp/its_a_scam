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

#include "gcatlib/SamEntry.hpp"
#include "gcatlib/SamReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
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
      << "\t-o out file prefix [required]" << endl
      << "\t-m max fragment length to track [default: 600]" << endl
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
    string bc_file;
    string out_prefix;

    char bc_delim = ':';
    char bc_col = 7;
    string bc_tag;

    size_t max_frag_len = 600;

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:o:m:d:c:t:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'm')
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

    if (aln_file.empty() || bc_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }


    // process barcodes
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;

    // to index into the count matrix
    unordered_map<string, size_t> bc_index;
    size_t bc_counter = 0;
    vector<string> bc_metadata;

    std::ifstream bc_in(bc_file);
    string line;
    while(getline(bc_in, line)) {
      vector<string> tokens;
      split_string(line, tokens, '\t');
      bc_index[tokens[0]] = bc_counter++;
      bc_metadata.push_back(tokens[0]);
    }
    bc_in.close();

    if (VERBOSE) {
      cerr << "\tNumber of barcodes: " << bc_counter << endl;
    }

    // create frag dist matrix for individual sample
    vector<vector<size_t>> bc_frag_dist(bc_counter,
                                      vector<size_t>(max_frag_len, 0));
    // create distribution for all samples
    vector<size_t> frag_dist(max_frag_len, 0);


    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;

    SamReader sam_reader(aln_file);
    SamEntry entry1, entry2;

    size_t aln_count = 0;

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


        
        // get the cell barcode, either from a tag or the name
        string cell_bc;
        if (!bc_tag.empty()) {
          // read barcode form the tag
          SamTags::get_tag(entry1.tags, bc_tag, cell_bc);
        }
        else {
          // parse the name to get the bc
          vector<string> tokens;
          split_string(entry1.qname, tokens, bc_delim);
          cell_bc = tokens[bc_col];
        }

        unordered_map<string, size_t>::iterator bc_it;
        bc_it = bc_index.find(cell_bc);
        if (bc_it != bc_index.end()) {

          size_t entry_bc_index = bc_it->second;

          // find the start and end location of the reads
          size_t e1_start = entry1.pos - 1;
          size_t e1_end = SamCigar::reference_end_pos(entry1);
          size_t e2_start = entry2.pos - 1;
          size_t e2_end = SamCigar::reference_end_pos(entry2);

          // find the stand and end location of the fragments
          size_t frag_start = e1_start;
          if (e2_start < e1_start)
            frag_start = e2_start;

          size_t frag_end = e1_end;
          if (e2_end > e1_end)
            frag_end = e2_end;

          size_t frag_len = frag_end - frag_start;

          // keep track of stuff
          if (frag_len <= max_frag_len) {
            ++bc_frag_dist[entry_bc_index][frag_len - 1];
            ++frag_dist[frag_len - 1];
          }

        }
      }
    }


    // write output
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;

    // write the per barcode fragment length distribution
    std::ofstream bc_out(out_prefix + "_bc_frag_len.txt");
    // write header line
    bc_out << "barcode";
    for (size_t i = 1; i <= max_frag_len; ++i) {
      bc_out << "\t" << i;
    }
    bc_out << endl;

    // write frag len matrix
    for (auto it = bc_index.begin(); it != bc_index.end(); ++it) {
      bc_out << it->first;
      size_t bc_frag_count = 0;
      // get the total number of fragemnts for the barcode for normalization
      for (size_t j = 0; j < max_frag_len; ++j) {
        bc_frag_count += bc_frag_dist[it->second][j];
      }
      // write the normalized frag length distribution
      for (size_t j = 0; j < max_frag_len; ++j) {
        bc_out << "\t" << static_cast<float>(bc_frag_dist[it->second][j]) /
                          static_cast<float>(bc_frag_count);
      }
      bc_out << endl;
    }
    bc_out.close();

    // write the sample fragment length distribution
    std::ofstream sample_out(out_prefix + "_sample_frag_len.txt");
    // write header line
    sample_out << "frag_len\tcount\tnorm_count" << endl;
    size_t sample_frag_count = 0;
    for (size_t i = 0; i < max_frag_len; ++i) {
      sample_frag_count += frag_dist[i];
    }
    for (size_t i = 0; i < max_frag_len; ++i) {
      sample_out << i+1 << "\t"
                 << frag_dist[i] << "\t"
                 << static_cast<float>(frag_dist[i]) /
                    static_cast<float>(sample_frag_count) << endl;
    }
    sample_out.close();


  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
