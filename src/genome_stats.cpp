/*
* genome_stats: calculate basic genome statistics
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
#include <unistd.h>

#include "gcatlib/GenomeReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-g genome file name [required]" << endl
      << "\t-v verbose [default: false]" << endl;
  
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {
   
    string genome_file;
    
    bool VERBOSE = false;
    
    int opt;
    while ((opt = getopt(argc, argv, "g:v")) != -1) {
      if (opt == 'g')
        genome_file = optarg;
      else if (opt == 'v')
        VERBOSE = true;
      else
        throw std::runtime_error(print_usage(argv[0]));
    } 
 
    if (genome_file.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    GenomeReader genome_reader(genome_file);
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
