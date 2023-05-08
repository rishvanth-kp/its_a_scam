/*
* FastqReader: class to read fastq files
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

#ifndef FASTQ_READER_HPP
#define FASTQ_READER_HPP

#include <iostream>
#include <fstream>
#include <string>

struct FastqEntry {
  std::string name;
  std::string seq;
  std::string separator;
  std::string quality;
};

class FastqReader {
public:
  FastqReader(const std::string &in_file);
  FastqReader(const std::string &in_file_1, const std::string &in_file_2);

  bool read_se_entry(FastqEntry &e);
  bool read_pe_entry(FastqEntry &e1, FastqEntry &e2);

  ~FastqReader();

private:
  std::ifstream in_1;
  std::ifstream in_2;

  bool is_pe;
};

std::ostream&
operator<<(std::ostream &out, const FastqEntry &e);

#endif
