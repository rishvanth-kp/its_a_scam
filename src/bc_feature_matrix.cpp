/*
* bc_feature_matrix: count barcodes aligning to genomic features
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
#include <unistd.h>

#include "AlignedOverlapGenomicFeature.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;

static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-b barcode list file [required]" << endl
      << "\t-g gtf file [optional]" << endl
      << "\t-t tsr bed file [optional]" << endl
      << "\t-e enhancers bed file [optional]" << endl
      << "\t-v verbose [default: false]" << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {
   
    string aln_file;
    string bc_file;
    
    string gtf_file;
    string tsr_file;
    string enhancer_file;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:g:t:e:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'g')
        gtf_file = optarg;
      else if (opt == 't')
        tsr_file = optarg;
      else if (opt == 'e')
        enhancer_file = optarg;
      else if (opt == 'v')
        VERBOSE = true;
      else 
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (VERBOSE)
      cerr << "[PROCESSING REGIONS]" << endl;

    AlignedOverlapGenomicFeature aligned_feature;
    aligned_feature.add_gtf_features(gtf_file);
 
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
