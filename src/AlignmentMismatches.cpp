/*
* AlignmentMismatches: class of keeping track of mismatches
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

#include "AlignmentMismatches.hpp"

using std::string;

AlignmentMismatch::AlignmentMismatch() {

}

AlignmentMismatch::~AlignmentMismatch() {

}

void
AlignmentMismatch::add(const SamEntry &sam_entry) {

  // parse cigar string
  SamCigar::CigarTuples cigar_tuple;
  SamCigar::string_to_tuple(sam_entry, cigar_tuple);

  // keep track of coverage for all the bases (not just for
  // mismatched bases). We cannot know if a base is going to be
  // mismatched in a read that is not yet seen.
  size_t aln_pos = sam_entry.pos - 1;
  for (auto it = cigar_tuple.begin(); it != cigar_tuple.end(); ++it) {
    // matches and mismatches consumes query and referene, and
    // contributes to coverage
    if (it->first == SamCigar::Cigar::aln_match ||
        it->first == SamCigar::Cigar::seq_match ||
        it->first == SamCigar::Cigar::seq_mismatch) {

      coverage.add(sam_entry.rname, aln_pos, aln_pos + it->second, 1);
      aln_pos += it->second;
    }
    // deletion and skipped regions consumes only reference, and
    // does not contribute to coverage
    else if (it->first == SamCigar::Cigar::ref_del ||
             it->first == SamCigar::Cigar::ref_skip) {
      aln_pos += it->second;
    }
  }

  // parse MD tag
  string md_tag;
  vector<pair<size_t, string>> md_tuple;
  SamTags::get_tag(sam_entry.tags, "MD", md_tag);
  SamTags::md_to_tuple(md_tag, md_tuple);

  // figure out the mismatch locations and keep track of them
  size_t ref_offset = 0;
  for (auto it = md_tuple.begin(); it != md_tuple.end(); ++it) {
    ref_offset += it->first;
    if (it->second[0] == '^')
      ref_offset += (it->second.length() - 1);
    else if (it->second != "") {
      ++ref_offset;
      size_t ref_pos = sam_entry.pos - 1;
      size_t query_pos = 0;
      SamCigar::move_in_reference(cigar_tuple, ref_offset,
        ref_pos, query_pos);
      mismatches.add(sam_entry.rname, ref_pos, ref_pos + 1,
        string{sam_entry.seq[query_pos]});
    }
  }
}


void
AlignmentMismatch::at(const GenomicRegion &in_region,
          vector<GenomicRegion> &out_region,
          vector<size_t> &out_aln_count,
          vector<string> &out_mismatches) {

  out_region.clear();
  out_aln_count.clear();
  out_mismatches.clear();

  vector<pair<GenomicRegion, string>> mm_out;
  vector<pair<GenomicRegion, size_t>> cov_out;
  mismatches.at(in_region, mm_out);
  for (auto it = mm_out.begin(); it != mm_out.end(); ++it) {
    for (size_t j = it->first.start; j < it->first.end; ++j) {
      out_region.push_back(GenomicRegion{it->first.name,
                            j, j + 1});
      out_mismatches.push_back(it->second);

      coverage.at(out_region.back(), cov_out);
      out_aln_count.push_back(cov_out[0].second);
    }
  }
}
