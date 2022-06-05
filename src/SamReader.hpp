/*
* BedReader: class to read SAM files
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

#ifndef SAM_READER_HPP
#define SAM_READER_HPP

#include <iostream>
#include <string>

#include "htslib/hts.h"
#include "htslib/sam.h"

using std::string;

class SamReader {
public:
  SamReader(const string &in_file);
  ~SamReader();

  bool read_sam_line();

private:

  htsFile *hts;
  sam_hdr_t *header;
  bam1_t *b;
  kstring_t sam_entry;
};

#endif
