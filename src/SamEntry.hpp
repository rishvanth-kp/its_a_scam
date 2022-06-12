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
  using CigarTuples = vector<pair<char, size_t>>;

  void
  string_to_tuple(const SamEntry &e, CigarTuples &tuples);
}

#endif
