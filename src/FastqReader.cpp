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

#include "FastqReader.hpp"

using std::string;
using std::cout;
using std::cerr;
using std::endl;


FastqReader::FastqReader (const std::string &in_file) {

  is_pe = false;

  in_1.open(in_file);
  if (!in_1)
    throw std::runtime_error("Cannot open " + in_file);
}


FastqReader::FastqReader (const std::string &in_file_1,
                          const std::string &in_file_2) {

  is_pe = true;

  in_1.open(in_file_1);
  if (!in_1)
    throw std::runtime_error("Cannot open " + in_file_1);

  in_2.open(in_file_2);
  if (!in_2)
    throw std::runtime_error("Cannot open " + in_file_2);

}

FastqReader::~FastqReader () {
  in_1.close();
  if (is_pe)
    in_2.close();
}


bool
FastqReader::read_se_entry (FastqEntry &e) {
  if (getline(in_1, e.name) &&
      getline(in_1, e.seq) &&
      getline(in_1, e.separator) &&
      getline(in_1, e.quality)) {

    return true;
  }

  return false;

}


bool
FastqReader::read_pe_entry(FastqEntry &e1, FastqEntry &e2) {
  if (getline(in_1, e1.name) &&
      getline(in_1, e1.seq) &&
      getline(in_1, e1.separator) &&
      getline(in_1, e1.quality) &&
      getline(in_2, e2.name) &&
      getline(in_2, e2.seq) &&
      getline(in_2, e2.separator) &&
      getline(in_2, e2.quality)) {

    return true;
  }

  return false;
}


std::ostream&
operator<<(std::ostream &out, const FastqEntry &e) {
  out << e.name << endl
      << e.seq << endl
      << e.separator << endl
      << e.quality << endl;
  return out;
}
