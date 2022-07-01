/*
* SamEntry: class to store SAM entry
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

#ifndef SAM_ENTRY_HPP
#define SAM_ENTRY_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

using std::string;
using std::string;
using std::vector;
using std::pair;

class SamEntry {
public:

  SamEntry() {};
  SamEntry(const string &line);
  ~SamEntry();

  void parse_entry(const string &line);

  string qname;
  uint16_t flag;
  string rname;
  uint32_t pos;
  uint16_t mapq;
  string cigar;
  string rnext;
  uint32_t pnext;
  int tlen;
  string seq;
  string qual;

  string tag_string;
  vector<string> tags;
};


namespace SamCigar {

  enum class Cigar : char {
    aln_match = 'M',
    ref_insert = 'I',
    ref_del = 'D',
    ref_skip = 'N',
    soft_clip = 'S',
    hard_clip = 'H',
    padding = 'P',
    seq_match = '=',
    seq_mismatch = 'X'
  };

  using CigarTuples = vector<pair<Cigar, size_t>>;

  void
  string_to_tuple(const SamEntry &e, CigarTuples &tuples);
}

namespace SamFlags {

  enum class Flag : uint16_t {
    read_paired = 0x0001,
    proper_pair = 0x0002,
    read_unmapped = 0x0004,
    mate_unmapped = 0x0008,
    read_reverse = 0x0010,
    mate_reverse = 0x0020,
    first_in_pair = 0x0040,
    second_in_pair = 0x0080,
    not_primary_aln = 0x0100,
    fail_qc = 0x0200,
    pcr_duplicate = 0x0400,
    supplementary_aln = 0x0800
  };

  constexpr bool is_any_set(const uint16_t flag, const uint16_t check) {
    return (flag & check);
  }

  constexpr bool is_all_set(const uint16_t flag, const uint16_t check) {
    return ((flag & check) == check);
  }

  constexpr bool is_set(const uint16_t flag, Flag f) {
    return (static_cast<uint16_t>(f) & flag);
  }

  constexpr uint16_t set(const uint16_t flag, Flag f) {
    return (static_cast<uint16_t>(f) | flag);
  }

}

#endif
