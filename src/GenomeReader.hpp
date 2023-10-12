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

#ifndef GENOME_READER_HPP
#define GENOME_READER_HPP

#include <vector>
#include <string>


class GenomeReader {
public:
  GenomeReader(const std::string &in_file);

  size_t chr_count() const;
  size_t chr_len(const size_t i) const;
  std::string chr_tag(const size_t i) const;
  std::string chr_substr(const size_t chr, const size_t pos,
                         const size_t len) const;
 
private:

  std::vector<std::string> chr_name;
  std::vector<std::string> chr_seq;
  std::vector<size_t> chr_abs_pos;
  size_t n_chr;
  size_t genome_size;
  bool VERBOSE = false;

  void read_genome (const std::string &in_file);
};

#endif
