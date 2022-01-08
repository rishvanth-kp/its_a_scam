/*
* GenomcRegion: class to store genomic intervals
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

#include "GenomicArray.hpp"

using std::cout;
using std::cerr;
using std::endl;

GenomicArray::GenomicArray() {
  n_chrom = 0;
  n_entry = 0;
}


bool
GenomicArray::add_entry(string chr, size_t start, size_t end, 
  string feature, string name, string id) {
/*
  cout << "Adding: " 
       << start << "\t"
       << end << "\t"
       << feature << "\t"
       << name << "\t"
       << id << endl;
*/

  ArrayEntry in;
  in.end = end;
  in.feature = feature;
  in.name = name;
  in.id = name;

  GenomeArrayType::iterator chr_map;
  chr_map = genomic_array.find(chr);
  if (chr_map == genomic_array.end()) {
    cout << "adding new chromosome " << chr << endl;
    ++n_chrom;
    ++n_entry;
    genomic_array[chr].insert(pair<size_t, ArrayEntry>(start, in));
  }
  else {
    ++n_entry;
    chr_map->second.insert(pair<size_t, ArrayEntry>(start, in));
  }

 
  return true;
}
