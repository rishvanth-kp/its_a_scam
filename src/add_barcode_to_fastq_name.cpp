/*
* add_barcode_to_fastq_name: read barcode fastq file and applend it to
*   fastq name
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
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;


static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-1 Fastq file 1 [required]" << endl
      << "\t-2 Fastq file 2 [required for paired-end data]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-c space delimated column to append barcode to [default: 0]" 
        << endl
      << "\t-s barcode read split position [default: 0]" << endl
      << "\t-d barcode read split delimeter [default: '+']"  << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {

    string in_file_1;
    string in_file_2;
    string out_prefix;

    size_t bc_col = 0;
    size_t bc_split_pos = 0;
    char bc_split_delim = '+';
    
    int opt;
    while ((opt = getopt(argc, argv, "1:2:o:f:")) != -1) {
      if (opt == '1') 
        in_file_1 = optarg;
      else if (opt == '2')
        in_file_2 = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'c')
        bc_col = std::stoi(optarg);
      else if (opt == 's')
        bc_split_pos = std::stoi(optarg);
      else if (opt == 'd')
        bc_split_delim = optarg[0];
      else 
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (in_file_1.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

  }

  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
