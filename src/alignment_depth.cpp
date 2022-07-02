/*
* alignment_depth: compute depth of aligned regions
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
#include <sstream>
#include <fstream>
#include <unistd.h>

#include "SamEntry.hpp"
#include "SamReader.hpp"
#include "GenomicStepVector.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;


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
    uint16_t include_none = 0x804;

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

    SamReader reader(aln_file);
    SamEntry entry;
    GenomicStepVector<size_t> coverage;

    // parse genomic regions from sam file
    vector<GenomicRegion> ref_chroms;
    string header;
    reader.read_sam_header(header);
    get_seq_lengths(header, ref_chroms);

    // determine coverage
    while (reader.read_sam_line(entry)) {
      if (entry.mapq >= min_mapq &&
          SamFlags::is_all_set(entry.flag, include_all) &&
          !SamFlags::is_any_set(entry.flag, include_none)) {

        // sam positons are 1-based
        size_t aln_pos = entry.pos - 1;
        SamCigar::CigarTuples tuples;
        SamCigar::string_to_tuple(entry, tuples);
        for (auto it = tuples.begin(); it != tuples.end(); ++it) {
          // matches and mismatches consumes query and referene, and
          // contributes to coverage
          if (it->first == SamCigar::Cigar::aln_match ||
              it->first == SamCigar::Cigar::seq_match ||
              it->first == SamCigar::Cigar::seq_mismatch) {

            coverage.add(entry.rname, aln_pos, aln_pos + it->second, 1);
            aln_pos += it->second;
          }
          // deletion and skipped regions consumes only reference, and
          // does not contribute to coverage
          else if (it->first == SamCigar::Cigar::ref_del ||
                   it->first == SamCigar::Cigar::ref_skip) {
            aln_pos += it->second;
          }
          // insertions, clipped, and padded bases does not consume referece
          // nothing to be done
        }
      }
    }

    // write output
    std::ofstream depth_file(out_prefix + "_depth.txt");
    vector<pair<GenomicRegion, size_t>> out;
    for (size_t i = 0; i < ref_chroms.size(); ++i) {
      coverage.at(ref_chroms[i], out);
      for (size_t j = 0; j < out.size(); ++j) {
        depth_file << out[j].first << "\t"
                   << out[j].second << endl;
      }
    }
    depth_file.close();
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
