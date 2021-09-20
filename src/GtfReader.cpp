/*
* GtfReader: class to read GTF files 
* Copyright (C) 2021 Rishvanth Prabakar
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

#include "GtfReader.hpp"

using std::cout;
using std::cerr;
using std::endl;

GtfReader::GtfReader(const string &in_file) {
  in.open(in_file);
  if (!in)
    throw std::runtime_error("Cannot open: " + in_file); 
  
  // Gobble comment lines
  string line;
  while (in.peek() == '#') {
    getline(in, line);
    cout << "comm: " << line << endl; 
  }
}

GtfReader::~GtfReader() {
  in.close();
}

bool
GtfReader::read_gtf_line(GtfEntry &g) {
  string line;
  if (getline(in, line)) {
    parse_gtf_line(line, g); 
    return true;
  }
  return false;
}

void
GtfReader::parse_gtf_line(const string &in, GtfEntry &g) {
  cout << "parsing: " << in << endl;
  
}
