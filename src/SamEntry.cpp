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

#include "SamEntry.hpp"

using std::cout;
using std::cerr;
using std::endl;

SamEntry::SamEntry(const string &line) {
  parse_entry(line);
}

SamEntry::~SamEntry() {

}

void
SamEntry::parse_entry(const string &line) {

  std::istringstream iss(line);
  iss >> qname >> flag >> rname >> pos >> mapq >> cigar
      >> rnext >> pnext >> tlen >>  seq >> qual;

  tags.clear();
  string tag;
  while (iss >> tag) {
    tags.push_back(tag);
  }
}



void
SamCigar::string_to_tuple(const SamEntry &e, CigarTuples &tuples) {

  tuples.clear();

  size_t start_pos = 0;
  size_t end_pos = e.cigar.find_first_of("MIDNSHP=X");

  while (end_pos != string::npos) {

    tuples.push_back(std::make_pair(
      static_cast<SamCigar::Cigar>(e.cigar[end_pos]),
      std::stoi(e.cigar.substr(start_pos, end_pos - start_pos))));

    start_pos = ++end_pos;
    end_pos = e.cigar.find_first_of("MIDNSHP=X", end_pos);
  }
}

bool
SamTags::get_tag(const vector<string> &tags, const string &tag, string &value) {
  value.clear();
  bool found = false;
  size_t i = 0;
  // Sam tags have are in TAG:TYPE:VALUE format.
  // TAG is always 2 characters and TYPE is always a single character
  while (!found && (i < tags.size())) {
    if (tags[i].substr(0, 2) == tag) {
      value = tags[i].substr(5);
      found = true;
    }
    ++i;
  }
  return found;
}

void
SamTags::md_to_tuple(const string &md_tag,
                     vector<pair<size_t, string>> &tuples) {

  tuples.clear();

  size_t len;
  string op;

  size_t start_pos = 0;
  size_t end_pos = md_tag.find_first_of("ATGCatgc^");
  len = std::stoi(md_tag.substr(start_pos, end_pos - start_pos));

  while (end_pos != string::npos) {

    start_pos = end_pos;
    end_pos = md_tag.find_first_of("0123456789", end_pos);
    op = md_tag.substr(start_pos, end_pos - start_pos);

    tuples.push_back(std::make_pair(len, op));

    start_pos = end_pos;
    end_pos = md_tag.find_first_of("ATGCatgc^", end_pos);
    len = std::stoi(md_tag.substr(start_pos, end_pos - start_pos));
  }
  op = "";
  tuples.push_back(std::make_pair(len, op));
}
