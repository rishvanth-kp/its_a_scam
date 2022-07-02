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

#include <sstream>

#include "SamReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

SamReader::SamReader(const string &in_file) {

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


void
SamReader::read_sam_header(string &hdr) {
  hdr = sam_hdr_str(header);
}


void
get_seq_lengths(const string &hdr, vector<GenomicRegion> &out) {
  out.clear();
  std::istringstream iss(hdr);
  string line;
  // @SQ header lines contain the reference sequence dictionaly
  // SN tag contains the reference sequence name
  // LN tag contains the reference sequence length. 1-based.
  while(getline(iss, line)) {
    if (line.substr(0, 3) == "@SQ") {
      std::istringstream token(line);
      string tag;
      string rname;
      size_t rlen;
      while(token >> tag) {
        if (tag.substr(0, 2) == "SN")
          rname = tag.substr(3);
        else if (tag.substr(0, 2) == "LN")
          rlen = std::stoi(tag.substr(3));
      }
      out.push_back(GenomicRegion(rname, 0, rlen));
    }
  }
}
