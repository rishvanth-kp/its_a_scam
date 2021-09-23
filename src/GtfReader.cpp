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

#include <limits>

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
GtfReader::read_gtf_file(vector<GtfEntry> &g) {
  string line;
  GtfEntry a; 
  while (getline(in, line)) {
    parse_gtf_line(line, a);
    g.push_back(a); 
  }
}

void
GtfReader::parse_gtf_line(const string &in, GtfEntry &g) {
  cout << "parsing: " << in << endl;
  // parse the first 8 required fields
  vector<string> required;
  size_t start = 0, end;
  for (size_t i = 0; i < 8; ++i) {
    end = in.find('\t', start);
    required.push_back(in.substr(start, end - start));
    cout << required.back() << endl;
    start = ++end;
  }

  g.name = required[0];
  g.source = required[1];
  g.feature = required[2];
  g.start = atoi(required[3].c_str());
  g.end = atoi(required[4].c_str());

  if (required[5][0] != '.')
    g.score = atof(required[5].c_str());
  else 
    g.score = std::numeric_limits<float>::min();

  g.strand = required[6][0];

  if (required[7][0] == '.')
    g.frame = std::numeric_limits<size_t>::max();


  // parse attributes, if any
  // attrbutes are ';' separated
  // key and value of an attribute are ' ' separated
  // there is also a space after a ';'
  while (end != in.length()) {
    pair<string, string> attr;

    end = in.find(' ', start);
    attr.first = in.substr(start, end - start);
    start = ++end;
    cout << "key: " << attr.first << endl;

    end = in.find(';', start);
    attr.second = in.substr(start, end - start);
    ++end;
    start = end + 1;
    cout << "value: " << attr.second << endl;

    g.attribute.push_back(attr);
  }
}
