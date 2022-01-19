/*
* PreprocessGff: parse features from GFF file 
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

#include "PreprocessGff.hpp"

using std::cout;
using std::cerr;
using std::endl;

PreprocessGff::PreprocessGff(const string &chrom_size_file, 
                             bool verbose) {
  VERBOSE = verbose;
  read_chrom_sizes(chrom_size_file);
}

PreprocessGff::~PreprocessGff() {
};

void
PreprocessGff::read_chrom_sizes(const string &chrom_size_file) {

  if (VERBOSE)
    cerr << "READING CHROM SIZES FROM " << chrom_size_file << endl;
  
  std::ifstream in(chrom_size_file);
  if (!in) 
    throw std::runtime_error("cannot open " + chrom_size_file);     

  string line;
  while (getline(in, line)) {
    size_t split = line.find('\t');
    string chrom = line.substr(0, split);
    size_t sz = std::stoi(line.substr(split));
    chrom_sizes[chrom] = sz; 
  }

  in.close();
  
  if (VERBOSE)
    cerr << chrom_sizes.size() << " CHROM SIZES INSERTED" << endl;
}
