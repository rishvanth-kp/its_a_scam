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
#include <sstream>
#include <fstream>
#include <unistd.h>

#include "gcatlib/SamEntry.hpp"
#include "gcatlib/Metagene.hpp"
#include "gcatlib/SamReader.hpp"

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


    // process the regions
    Metagene metagene(regions_file);
  }
  catch (const std::exception  &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
