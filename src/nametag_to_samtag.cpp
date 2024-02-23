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
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-o out SAM file [required]" << endl
      << "\t-d name split delimeter [default: \":\"]" << endl
      << "\t-c barcode field in name [default: 7]" << endl
      << "\t-t SAM tag to add barcode to [default: CB]" << endl
      << "\t-z SAM tag type [default: z]" << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {
    
    string aln_file;
    string out_file;
    
    char bc_delim = ':';
    size_t bc_col = 7;

    string bc_tag = "CB";
    char bc_tag_type = 'z';
  
    int opt;
    while ((opt = getopt(argc, argv, "a:o:d:c:t:z")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'o')
        out_file = optarg;
      else if (opt == 'd')
        bc_delim = optarg[0];
      else if (opt == 'c')
        bc_col = std::stoi(optarg);
      else if (opt == 't')
        bc_tag = optarg;
      else if (opt == 'z')
        bc_tag_type = optarg[0];
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (aln_file.empty() || out_file.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

  }
  
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
