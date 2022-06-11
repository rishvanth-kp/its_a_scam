/*
* SamReader: class to read SAM files
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

#include "SamReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

SamReader::SamReader(const string &in_file) {
  cout << "Opening sam file: " << in_file << endl;

  if (!(hts = hts_open(in_file.c_str(), "r")))
    throw std::runtime_error("Cannot open SAM/BAM file: " + in_file);

  if (hts_get_format(hts)->category != sequence_data)
    throw std::runtime_error(in_file + " is not a sequence file");

  if (!(header = sam_hdr_read(hts)))
    throw std::runtime_error("No header in " + in_file);

  if (!(b = bam_init1()))
    throw std::runtime_error("Failed to initialize bam");

  ks_initialize(&sam_entry);
  eof = false;
}

SamReader::~SamReader() {
  hts_close(hts);
  sam_hdr_destroy(header);
  bam_destroy1(b);
  ks_free(&sam_entry);
}

bool
SamReader::read_sam_line(SamEntry &entry) {
  int read_ret = 0;
  if ((read_ret = sam_read1(hts, header, b)) >= 0) {
    int fmt_ret = 0;
    if ((fmt_ret = sam_format1(header, b, &sam_entry)) <= 0) {
      throw std::runtime_error("error reading sam entry");
    }
    entry.parse_entry(sam_entry.s);
  }
  else if (read_ret == -1)
    eof = true;
  else
    throw std::runtime_error("error reading sam entry");


  return !eof;
}
