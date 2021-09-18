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

#ifndef GENOMIC_REGION_HPP
#define GENOMIC_REGION_HPP

#include <iostream>
#include <string>
#include <limits>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

class GenomicRegion {
public:
  GenomicRegion(string n="", size_t s=0, 
                size_t e=std::numeric_limits<size_t>::max(), 
                char st='.') {
    name = n;
    start = s;
    end = e;
    strand = st;
    if (!is_valid())
      throw std::runtime_error("Invalid genomic interval"); 
  }
  
  string get_name() const { return name; }
  size_t get_start() const { return start; }
  size_t get_end() const { return end; }
  char get_strand() const { return strand; }
  
  void set_name(const string n) { name = n; }
  void set_start(const size_t s) { start = s; }
  void set_end(const size_t e) { end = e; }
  void set_stand(const size_t st) { strand = st; }

  bool overlaps(GenomicRegion &a) const;

private:
  // Following the BED convntion of half-open intervals
  string name;
  size_t start;
  size_t end;
  char strand;

  bool is_valid() const;
};



#endif
