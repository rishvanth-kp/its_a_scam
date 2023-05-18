/*
* BedReader: class to read bed files
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

#include <iostream>
#include <sstream>

#include "BedReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;

BedReader::BedReader (const string &in_file) {
  in.open(in_file);
  if (!in)
    throw std::runtime_error("Cannot open " + in_file);
}

BedReader::~BedReader () {
  in.close();
}


static void
split_string (const string &in, vector<string> &tokens,
              const char delim = '\t') {

  tokens.clear();
  size_t start = 0;
  size_t end = in.find(delim);
  while (end != string::npos) {
    tokens.push_back(in.substr(start, end - start));
    start = ++end;
    end = in.find(delim, start);
  }
  tokens.push_back(in.substr(start));
}

bool
BedReader::read_bed3_line (GenomicRegion &g) {

  string line;
  if (getline(in, line)) {
    vector<string> tokens;
    split_string(line, tokens);
    g.name = tokens[0];
    g.start = atoi(tokens[1].c_str());
    g.end = atoi(tokens[2].c_str());
    return true;
  }

  return false;
}


bool
BedReader::read_bed_line(GenomicRegion &g, vector<string> &fields) {
  string line;
  if (getline(in, line)) {
    std::istringstream iss(line);
    iss >> g.name >> g.start >> g.end;

    fields.clear();
    string field;
    while (iss >> field) {
      fields.push_back(field);
    }

    return true;
  }
  return false;
}


void
BedReader::read_bed3_file (vector<GenomicRegion> &g) {

  string line;
  while (getline(in, line)) {
    vector<string> tokens;
    split_string(line, tokens);
    GenomicRegion a(tokens[0],
                    atoi(tokens[1].c_str()),
                    atoi(tokens[2].c_str()));
    g.push_back(a);
  }
}
