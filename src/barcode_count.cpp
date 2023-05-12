/*
* barcode_count: count barcodes in a sam file
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

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <unistd.h>
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

static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all the flags are present [default: 3]" << endl
      << "\t-F only include if none of the flags are present [default: 2316]" 
        << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {
    
    string aln_file;
    string out_prefix;
    
    size_t min_mapq = 0;
    uint16_t include_all = 0x0003;
    uint16_t include_none = 0x0804;

    int opt;
    while ((opt = getopt(argc, argv, "a:o:q:f:F:")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'q')
        min_mapq = std::stoi(optarg);
      else if (opt == 'f')
        include_all = std::stoi(optarg);
      else if (opt == 'F')
        include_none = std::stoi(optarg);
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (aln_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    SamReader sam_reader(aln_file);
    unordered_map<string, size_t> bc_counter;
    SamEntry entry1, entry2;
    while (sam_reader.read_pe_sam(entry1, entry2)) {

      if (entry1.mapq >= min_mapq && entry2.mapq >= min_mapq &&
          SamFlags::is_all_set(entry1.flag, include_all) &&
          SamFlags::is_all_set(entry2.flag, include_all) &&
          !SamFlags::is_any_set(entry1.flag, include_none) &&
          !SamFlags::is_any_set(entry2.flag, include_none)) {

        vector<string> tokens;
        split_string(entry1.qname, tokens);
         
        unordered_map<string, size_t>::iterator it;
        it = bc_counter.find(tokens[7]);
        if (it == bc_counter.end()) {
          bc_counter[tokens[7]] = 1;
        }      
        else {
          bc_counter[tokens[7]] += 1;
        }
      }
    }
    for (auto it = bc_counter.begin(); it != bc_counter.end(); ++it) {
      cout << it->first << "\t" << it->second << endl;
    } 

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
