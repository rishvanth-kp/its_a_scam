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

#ifndef GENOMIC_ARRAY_HPP
#define GENOMIC_ARRAY_HPP
  
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>

using std::pair;
using std::vector;
using std::string;
using std::multimap;
using std::unordered_map;

struct ArrayEntry {
  size_t end;
  string feature;
  string name;
  string id;
};

class GenomicArray {
public:
  GenomicArray();

  bool add_entry(string chr, size_t start, size_t end, 
                 string feature, string name, string id); 

  size_t chrom_count() const { return n_chrom; }
  size_t entry_count() const { return n_entry; }


private:
  
  typedef multimap<size_t, ArrayEntry> ChromArrayType;
  typedef unordered_map<string, ChromArrayType> GenomeArrayType;

  GenomeArrayType genomic_array;
  bool VERBOSE = false;
  size_t n_chrom;
  size_t n_entry;
};

#endif
