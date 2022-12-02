/*
* GeneCount: Determine gene counts of aligned reads
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

#include "GeneCount.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::unordered_map;

GeneCount::GeneCount() {

}

GeneCount::~GeneCount() {

}


void
GeneCount::preprocess_gff(const string &gff_file) {
  
  GtfReader gff_reader(gff_file);
  GencodeGtfEntry entry;

  while (gff_reader.read_gencode_gtf_line(entry)) {
    if (entry.feature == "exon") {
      cout << entry.gene_name << endl;
    }
  }
  
}
