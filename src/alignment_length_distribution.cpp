/*
* alignment_length_distribution: compute histogram of alignment lengths
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

#include "gcatlib/SamReader.hpp"
#include "gcatlib/SamEntry.hpp"

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
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all of the flags are present [default: 0]"
          << endl
      << "\t-F only include if none of the flags are present [default: 2052]"
          << endl
      << "\t max read length [default: 200]"
      << "\t-o out file prefix [required]" << endl;
  return oss.str();
}

int
main(int argc, char* argv[]) {
  try {

    string aln_file;
    string out_prefix;
    size_t min_mapq = 0;
    size_t max_len = 200;

    uint16_t include_all = 0;
    uint16_t include_none = 0x0804;


    int opt;
    while ((opt = getopt(argc, argv, "a:q:f:F:m:o:")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'q')
        min_mapq = std::stoi(optarg);
      else if (opt == 'f')
        include_all = std::stoi(optarg);
      else if (opt == 'F')
        include_none = std::stoi(optarg);
      else if (opt == 'm')
        max_len = std::stoi(optarg);
      else if (opt == 'o')
        out_prefix = optarg;
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (aln_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    SamReader reader(aln_file);
    SamEntry entry;

    vector<size_t> read_len(max_len, 0);
    vector<size_t> aln_len(max_len, 0);

    while (reader.read_sam_line(entry)) {
      if (entry.mapq >= min_mapq &&
          SamFlags::is_all_set(entry.flag, include_all) &&
          !SamFlags::is_any_set(entry.flag, include_none)) {

        const size_t seq_len = entry.seq.length();
        size_t clip_len = 0;

        SamCigar::CigarTuples tuples;
        SamCigar::string_to_tuple(entry, tuples);

        // 'H' can only be present as the first or last operation.
        // 'H' bases are not present in SEQ, skip any 'H'.
        // 'S' can only have 'H' between them and the end.
        // 'S' bases are present in SEQ, so remove them from aligned len.
        SamCigar::CigarTuples::iterator it = tuples.begin();
        if (it->first == SamCigar::Cigar::hard_clip)
          ++it;
        if (it->first == SamCigar::Cigar::soft_clip)
          clip_len += it->second;

        SamCigar::CigarTuples::reverse_iterator rit = tuples.rbegin();
        if (rit->first == SamCigar::Cigar::hard_clip)
          ++rit;
        if (rit->first == SamCigar::Cigar::soft_clip)
          clip_len += rit->second;

        ++read_len[seq_len];
        ++aln_len[seq_len - clip_len];
      }
    }

    std::ofstream hist_file(out_prefix + "_aln_len_dist.txt");
    hist_file << "length" << "\t"
              << "seq_bases" << "\t"
              << "aln_bases" << endl;
    for (size_t i = 0; i < read_len.size(); ++i) {
      hist_file << i << "\t"
                << read_len[i] << "\t"
                << aln_len[i] << endl;
    }
    hist_file.close();

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
