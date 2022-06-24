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
  explicit GenomicRegion(string n="", size_t s=0,
                size_t e=std::numeric_limits<size_t>::max(),
                char st='.') {
    name = n;
    start = s;
    end = e;
    strand = st;
  }

  // Following the BED convntion of half-open intervals
  string name;
  size_t start;
  size_t end;
  char strand;

};

std::ostream&
operator<<(std::ostream &out, const GenomicRegion &g);

#endif
