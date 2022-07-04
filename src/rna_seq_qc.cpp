/*
* rna_seq_qc: determines rna-seq QC metrics
* Copyright (C) 2021 Rishvanth Prabakar
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
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <unistd.h>

#include "GtfReader.hpp"
#include "SamReader.hpp"
#include "SamEntry.hpp"
#include "PreprocessGff.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStepVector.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::map;
using std::vector;
using std::string;

static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-g gtf file (required)" << endl
      << "\t-c chrom size file (required)" << endl
      << "\t-v verbose (default: false)" << endl;
  return oss.str();
}

int
main(int argc, char *argv[]) {
  try {

    string gtf_file;
    string chrom_size_file;
    bool VERBOSE{false};

    int opt;
    while ((opt = getopt(argc, argv, "g:c:v")) != -1) {
      if (opt == 'g')
        gtf_file = optarg;
      else if (opt == 'c')
        chrom_size_file = optarg;
      else if (opt == 'v')
        VERBOSE = true;
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (gtf_file.empty() || chrom_size_file.empty())
      throw std::runtime_error(print_usage(argv[0]));

    // PreprocessGff gff_processor(chrom_size_file, VERBOSE);
    // gff_processor.parse_genome_features(gtf_file);

    GenomicStepVector<size_t> step_vec;
    step_vec.add("chr1", 15, 20, 1);
    step_vec.add("chr1", 3, 5, 2);
    step_vec.add("chr1", 20, 21, 3);
    step_vec.add("chr1", 21, 22, 3);
    step_vec.add("chr1", 22, 23, 3);
    step_vec.add("chr1", 23, 24, 3);
    step_vec.add("chr1", 24, 31, 3);
    vector<std::pair<GenomicRegion, size_t>> out;
    step_vec.at(GenomicRegion("chr1", 0, 35), out, true);
    for (size_t i = 0; i < out.size(); ++i) {
      cout << out[i].first << "\t"
           << out[i].second << endl;
    }
  

  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
