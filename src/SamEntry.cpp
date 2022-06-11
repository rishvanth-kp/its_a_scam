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

  /*
  cout << qname << endl;
  cout << flag << endl;
  cout << rname << endl;
  cout << pos << endl;
  cout << mapq << endl;
  cout << cigar << endl;
  cout << rnext << endl;
  cout << pnext << endl;
  cout << tlen << endl;
  cout << seq << endl;
  cout << qual << endl;
  for (size_t i = 0; i < tags.size(); ++i)
    cout << tags[i] << endl;
  */
}
