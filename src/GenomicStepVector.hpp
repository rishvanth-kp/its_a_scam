/*
* GenomcRegion: class to store genomic intervals
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

#ifndef GENOMIC_STEP_VECTOR_HPP
#define GENOMIC_STEP_VECTOR_HPP

#include <iostream>
#include <string>
#include <unordered_map>

#include "StepVector.hpp"

using std::pair;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::unordered_map;

template<typename T>
class GenomicStepVector {
public:
  GenomicStepVector();  

  void add(const string chr, const size_t start, const size_t end, const T val);
   

  size_t chrom_count() const { return n_chrom; }
  size_t entry_count() const { return n_entry; }  

private: 
  size_t n_chrom{};
  size_t n_entry{};
 
  unordered_map<string, StepVector<T>> genomic_vector;
};


template<typename T>
GenomicStepVector<T>::GenomicStepVector() {
  cout << "Initializing genomic step vector" << endl;
}

template<typename T>
void
GenomicStepVector<T>::add(const string chr, const size_t start,
                          const size_t end, const T val) {
  cout << "adding entry" << endl;
  cout << chr << "\t" << start << "\t" << end << "\t" << val << endl;

  typename unordered_map<string, StepVector<T>>::iterator chr_map;
  chr_map = genomic_vector.find(chr);
  if (chr_map == genomic_vector.end()) {
    cout << "adding new chrom" << endl;
    genomic_vector[chr].add(start, end, val); 
    ++n_chrom;
    ++n_entry;
    genomic_vector[chr].print_elements();
  } 
  else {
    chr_map->second.add(start, end, val); 
    ++n_entry;
    chr_map->second.print_elements();
  }
}

#endif
