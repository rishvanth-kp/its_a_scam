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
#include "GenomicRegion.hpp"

using std::pair;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::make_pair;
using std::unordered_map;

template<typename T>
class GenomicStepVector {
public:
  GenomicStepVector();

  void add(const string chr, const size_t start, const size_t end, const T val);
  void add(const GenomicRegion &g, const T val);

  void at(const GenomicRegion &g,
          vector<pair<GenomicRegion, T>> &out,
          bool keep_0 = false) const;

  size_t chrom_count() const { return n_chrom; }
  size_t entry_count() const { return n_entry; }

private:
  size_t n_chrom{};
  size_t n_entry{};

  typedef unordered_map<string, StepVector<T>> GenomicStepVectorType;
  GenomicStepVectorType genomic_vector;
};


template<typename T>
GenomicStepVector<T>::GenomicStepVector() {
}

template<typename T>
void
GenomicStepVector<T>::add(const string chr, const size_t start,
                          const size_t end, const T val) {

  typename GenomicStepVectorType::iterator chr_map;
  chr_map = genomic_vector.find(chr);
  if (chr_map == genomic_vector.end()) {
    genomic_vector[chr].add(start, end, val);
    ++n_chrom;
    ++n_entry;
    // genomic_vector[chr].print_elements();
  }
  else {
    chr_map->second.add(start, end, val);
    ++n_entry;
    // chr_map->second.print_elements();
  }
}

template<typename T>
void
GenomicStepVector<T>::add(const GenomicRegion &g, const T val) {
  add(g.name, g.start, g.end, val);
}


template<typename T>
void
GenomicStepVector<T>::at(const GenomicRegion &g,
                         vector<pair<GenomicRegion, T>> &out,
                         bool keep_0) const {

  out.clear();

  typename GenomicStepVectorType::const_iterator chr_map;
  chr_map = genomic_vector.find(g.name);
  if (chr_map != genomic_vector.end()) {
    vector<pair<size_t, T>> step_out;
    chr_map->second.at_range(g.start, g.end, step_out);

    GenomicRegion region;
    for (size_t i = 0; i < step_out.size() - 1; ++i) {
      if (keep_0 || (!keep_0 && step_out[i].second != T{})) {
        region.name = g.name;
        region.start = step_out[i].first;
        region.end = step_out[i + 1].first;
        out.push_back(make_pair(region, step_out[i].second));
      }
    }
  }
}

#endif
