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


#include "BedReader.hpp"

using std::cout;
using std::cerr;
using std::endl;

BedReader::BedReader (const string &in_file) {
  in.open(in_file);
  cout << "Opening file " << in_file << endl;   
  if (!in) 
    throw std::runtime_error("Cannot open " + in_file);
}

BedReader::~BedReader () {
  cout << "closing bed" << endl; 
  in.close();
}


void
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
    g.set_name(tokens[0]);
    g.set_start(atoi(tokens[1].c_str()));
    g.set_end(atoi(tokens[2].c_str()));
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
