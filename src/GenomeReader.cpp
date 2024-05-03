/*
* GenomeReader: class to read reference genomes
* Copyright (C) 2023 Rishvanth Prabakar
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

#include <algorithm>
#include <iostream>
#include <fstream>

#include "GenomeReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::unordered_map;

GenomeReader::GenomeReader(const string &in_file) {
  read_genome(in_file);
}


void
GenomeReader::read_genome (const string &in_file) {

  std::ifstream in(in_file);
  if (!in)
    throw std::runtime_error("Cannot open " + in_file);

  string line;
  while (getline(in, line)) {
    cout << line << endl;
    if (line[0] == '>') {
      chr_name.push_back(line.substr(1));
      chr_seq.push_back(string{});
      chr_index[line.substr(1)] = ++n_chr;
    }
    else {
      chr_seq.back() += line;
    }
  }

  chr_abs_pos.push_back(0);
  for (size_t i = 0; i < n_chr; ++i) {
    std::transform(chr_seq[i].begin(), chr_seq[i].end(),
                   chr_seq[i].begin(), toupper);

    genome_size += chr_seq[i].length();
    chr_abs_pos.push_back(genome_size);
  }

  cout << "genome size: " << genome_size << endl;
  for (size_t i = 0; i < n_chr; ++i) {
    cerr << "\t" << chr_name[i]
         << "\t" << chr_seq[i].length()
         << "\t" << chr_abs_pos[i]
         << "\t" << chr_seq[i] << endl;
  }

  in.close();
}


size_t
GenomeReader::chrom_count() const {
  return n_chr;
}


size_t
GenomeReader::chrom_len(const string chrom) {
  unordered_map<string, size_t>::iterator it = chr_index.find(chrom);
  if (it != chr_index.end())
    return chr_seq[it->second].length();
  else
    return 0;
}



std::string
GenomeReader::chrom_substr(const string chrom, const size_t start,
                         const size_t end) {

  unordered_map<string, size_t>::iterator it = chr_index.find(chrom);
  if (it != chr_index.end())
    return chr_seq[it->second].substr(start, end - start);
  else
    return string();
}
