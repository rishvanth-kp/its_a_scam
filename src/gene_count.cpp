/*
* gene_count: count bases/reads aligning to each gene
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

#include "gcatlib/GeneCount.hpp"
#include "gcatlib/SamReader.hpp"
#include "gcatlib/SamEntry.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-g GTF feature file [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-s single end reads? [default: false]" << endl
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all the flags are present [default: 0]" << endl
      << "\t-F only incldue if none of the flags are present [default: 2316]"
        << endl;
  return oss.str();
}

int
main(int argc, char* argv[]) {
  try {

    string aln_file;
    string out_prefix;
    string gtf_file;

    bool se = false;
    size_t min_mapq = 0;
    uint16_t include_all = 0;
    uint16_t include_none = 0x0804;

    int opt;
    while ((opt = getopt(argc, argv, "a:g:q:f:F:o:s")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'g')
        gtf_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'q')
        min_mapq = std::stoi(optarg);
      else if (opt == 'f')
        include_all = std::stoi(optarg);
      else if (opt == 'F')
        include_none = std::stoi(optarg);
      else if (opt == 's')
        se = true;
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (aln_file.empty() || gtf_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    GeneCount gene_counter;
    gene_counter.preprocess_gff(gtf_file);

    SamReader sam_reader(aln_file);

    if (se) {
      SamEntry entry;
      while (sam_reader.read_sam_line(entry)) {
        if (entry.mapq >= min_mapq &&
            SamFlags::is_all_set(entry.flag, include_all) &&
            !SamFlags::is_any_set(entry.flag, include_none)) {

          gene_counter.add(entry);
        }
      }
    }
    else {
      SamEntry entry1, entry2;
      while (sam_reader.read_pe_sam(entry1, entry2)) {
        if (entry1.mapq >= min_mapq && entry2.mapq >= min_mapq &&
            SamFlags::is_all_set(entry1.flag, include_all) &&
            SamFlags::is_all_set(entry2.flag, include_all) &&
            !SamFlags::is_any_set(entry1.flag, include_none) &&
            !SamFlags::is_any_set(entry2.flag, include_none)) {

          gene_counter.add(entry1, entry2);
        }
      }
    }

    // vector<pair<string, size_t>> gene_counts;
    // gene_counter.get_gene_counts(gene_counts);
    gene_counter.gene_counts_to_file(out_prefix + ".txt");

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
