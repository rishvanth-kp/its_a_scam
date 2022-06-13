/*
* alignment_gc_distribution: compute histogram of GC content
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
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <unistd.h>

#include "SamReader.hpp"
#include "SamEntry.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

static float
gc_content(const string &seq, const size_t start, const size_t end) {

  float gc = 0;
  size_t i = start;
  while (i < end) {
    if (seq[i] == 'G' || seq[i] == 'g' ||
        seq[i] == 'C' || seq[i] == 'c') {
      ++gc;
    }
    ++i;
  }

  return ((gc / static_cast<float>(end - start)) * 100);
}

static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all of the flags are present [default: 0]"
          << endl
      << "\t-F only include if none of the flags are present [default: 2052]"
          << endl;
  return oss.str();
}

int
main(int argc, char* argv[]) {
  try {

    string aln_file;
    string out_prefix;
    size_t min_mapq = 0;

    uint16_t include_all = 0;
    uint16_t include_none = 0x0804;

    int opt;
    while ((opt = getopt(argc, argv, "a:q:o:")) != -1) {
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

    SamReader reader(aln_file);
    SamEntry entry;

    vector<size_t> gc_dist(101, 0);

    while (reader.read_sam_line(entry)) {
      if (entry.mapq >= min_mapq &&
          SamFlags::is_all_set(entry.flag, include_all) &&
          !SamFlags::is_any_set(entry.flag, include_none)) {

        size_t start_pos = 0;
        size_t end_pos = entry.seq.length();

        SamCigar::CigarTuples tuples;
        SamCigar::string_to_tuple(entry, tuples);

        SamCigar::CigarTuples::iterator it = tuples.begin();
        if (it->first == 'H')
          ++it;
        if (it->first == 'S')
          start_pos += it->second;

        SamCigar::CigarTuples::reverse_iterator rit = tuples.rbegin();
        if (rit->first == 'H')
          ++rit;
        if (rit->first == 'S')
          end_pos -= rit->second;

        ++gc_dist[std::round(gc_content(entry.seq, start_pos, end_pos))];
      }
    }

    std::ofstream hist_file(out_prefix + "_gc_dist.txt");
    hist_file << "gc" << "\t" << "reads" << endl;
    for (size_t i = 0; i < gc_dist.size(); ++i) {
      hist_file << i << "\t"
                << gc_dist[i] << endl;
    }
    hist_file.close();

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
