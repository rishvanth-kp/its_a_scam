/*
* barcode_count: count barcodes in a sam file
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
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <algorithm>
#include <unordered_map>

#include "SamReader.hpp"
#include "SamEntry.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::unordered_map;

void
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


static bool
sort_bc_count(pair<string, size_t> a, pair<string, size_t> b) {
  return (a.second > b.second);
}


static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-m minimun barcode count to output [default: 1000]" << endl
      << "\t-d name split delimeter" 
          << "[default: \":\"; ignored if -t is provided]" << endl
      << "\t-c barcode field in name" 
          << "[default: 7 (0 based); ignored if -t is provided]" << endl
      << "\t-t barcode tag in SAM file [default: \"\"]" << endl
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
    string out_prefix;

    char bc_delim = ':';
    size_t bc_col = 7;
    size_t min_bc_count = 1000;
    string bc_tag;

    size_t min_mapq = 0;
    uint16_t include_all = 0x0003;
    uint16_t include_none = 0x0D0C;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:o:m:d:c:t:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'm')
        min_bc_count = std::stoi(optarg);
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

    if (aln_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    SamReader sam_reader(aln_file);
    unordered_map<string, size_t> bc_counter;
    SamEntry entry1, entry2;

    size_t aln_count = 0;
    size_t aln_pass_count = 0;

    bool bc_in_tag = false;
    if (!bc_tag.empty()) {
      bc_in_tag = true;
    }

    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;

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


        string cell_bc;
        if (bc_in_tag) {
          SamTags::get_tag(entry1.tags, bc_tag, cell_bc);
        }
        else { 
          vector<string> tokens;
          split_string(entry1.qname, tokens, bc_delim);
          if (tokens.size() > bc_col)
              cell_bc = tokens[bc_col];
        }
 
        if (!cell_bc.empty()) {
          unordered_map<string, size_t>::iterator it;
          it = bc_counter.find(cell_bc);
          if (it == bc_counter.end()) {
            bc_counter[cell_bc] = 1;
          }
          else {
            ++bc_counter[cell_bc];
          }
        }

      }
    }

    // keep only barcodes meeting the min count
    if (VERBOSE)
      cerr << "[FILTERING LOW BARCODE COUNTS]" << endl;
    vector<pair<string, size_t>> bc_output;
    for (auto it = bc_counter.begin(); it != bc_counter.end(); ++it) {
      if (it->second >= min_bc_count) {
        bc_output.push_back(std::make_pair(it->first, it->second));
      }
    }

    // sort so that the highest counts appear first
    if (VERBOSE)
      cerr << "[SORTING BARCODES]" << endl;
    std::sort(bc_output.begin(), bc_output.end(), sort_bc_count);

    // write output
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;
    std::ofstream out_counts (out_prefix + "_bc_counts.txt");
    for (auto it = bc_output.begin(); it != bc_output.end(); ++it) {
      out_counts << it->first << "\t" << it->second <<endl;
    }
    out_counts.close();

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
