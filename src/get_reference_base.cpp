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

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;


static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-g reference genome file [required]" << endl
      << "\t-b bed file with positions to fetch [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-p number of previous base to feteh [default: 0]" << endl
      << "\t-n number of next bases to fetch [default: 0]" << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {
  
    string genome_file;
    string bed_file;
    string out_file;    

    size_t prev_len = 0;
    size_t next_len = 0;

    int opt;
    while((opt = getopt(argc, argv, "g:b:o:p:n")) != -1) {
      if (opt ==  'g')
        genome_file = optarg;
      else if (opt == 'b')
        bed_file = optarg;
      else if (opt == 'o')
        out_file = optarg;
      else if (opt == 'p')
        prev_len = std::stoi(optarg);
      else if (opt == 'n')
        next_len = std::stoi(optarg);
      else 
        throw std::runtime_error(print_usage(argv[0]));
    } 

    if (genome_file.empty() || bed_file.empty() || out_file.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} 
