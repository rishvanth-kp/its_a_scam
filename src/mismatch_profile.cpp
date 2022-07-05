/*
* mismatch_profile: report mismatch location ans statistics
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
#include <unordered_map>

#include "SamReader.hpp"
#include "SamEntry.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;


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


    unordered_map<string, size_t> mismatch_profile;

    vector<string> nucs {"A", "T", "G", "C"};
    for (size_t i = 0; i < nucs.size(); ++i) {
      for (size_t j = 0; j < nucs.size(); ++j) {
        if (i != j) {
          mismatch_profile.insert(make_pair(nucs[i] + nucs[j], 0));
        }
      }
    }



    size_t total_aln_length = 0;
    while (reader.read_sam_line(entry)) {
      if (entry.mapq >= min_mapq &&
          SamFlags::is_all_set(entry.flag, include_all) &&
          !SamFlags::is_any_set(entry.flag, include_none)) {

        SamCigar::CigarTuples cigar_tuple;
        SamCigar::string_to_tuple(entry, cigar_tuple);

        // calculte alignment length
        size_t clip_len = 0;
        SamCigar::CigarTuples::iterator it = cigar_tuple.begin();
        if (it->first == SamCigar::Cigar::hard_clip)
          ++it;
        if (it->first == SamCigar::Cigar::soft_clip)
          clip_len += it->second;

        SamCigar::CigarTuples::reverse_iterator rit = cigar_tuple.rbegin();
        if (rit->first == SamCigar::Cigar::hard_clip)
          ++rit;
        if (rit->first == SamCigar::Cigar::soft_clip)
          clip_len += rit->second;

        const size_t seq_len = entry.seq.length();
        total_aln_length += (seq_len - clip_len);

        string md_tag;
        vector<pair<size_t, string>> md_tuple;
        SamTags::get_tag(entry.tags, "MD", md_tag);
        SamTags::md_to_tuple(md_tag, md_tuple);

        size_t ref_offset = 0;
        for (auto it = md_tuple.begin(); it != md_tuple.end(); ++it) {
          ref_offset += it->first;
          if (it->second[0] == '^')
            ref_offset += (it->second.length() - 1);
          else if (it->second != "") {
            ++ref_offset;
            size_t ref_pos = 0;
            size_t query_pos = 0;
            SamCigar::move_in_reference(cigar_tuple, ref_offset,
              ref_pos, query_pos);
            mismatch_profile[it->second + entry.seq[query_pos]]++;
          }
        }

      }
    }

    // write output
    std::ofstream mm_stats_file(out_prefix + "_mismatch_stats.txt");
    for (size_t i = 0; i < nucs.size(); ++i) {
      for (size_t j = 0; j < nucs.size(); ++j) {
        if (i != j) {
          size_t mm_count = mismatch_profile[nucs[i] + nucs[j]];
          float mm_frac = static_cast<float>(mm_count) /
                          static_cast<float>(total_aln_length);
          mm_stats_file << nucs[i] + nucs[j] << "\t"
                        << mm_count << "\t"
                        << mm_frac <<  endl;
        }
      }
    }
    mm_stats_file.close();

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
