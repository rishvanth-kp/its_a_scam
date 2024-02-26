/*
* AlignmentMismatches: class of keeping track of mismatches
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

#ifndef ALIGNMENT_MISMATCHES_HPP
#define ALIGNMENT_MISMATCHES_HPP

#include <string>
#include <vector>

#include "SamEntry.hpp"
#include "GenomicStepVector.hpp"

/**
* \brief Keep track of mismatched bases and thier location.
*
* Parses the MD tag and CIGAR string from a SAM entry and
* tracks the location and count of mismatched bases.
*/
class AlignmentMismatch {
public:
  /**
  * Default constructor
  */
  AlignmentMismatch();
  /**
  * Default destructor
  */
  ~AlignmentMismatch();

  /**
  * Adds the bases covered and mismates in a sam entry.
  * The cigar string is parsed and all the matching and mismatching
  * bases are added to coverage. Then the MD tag is parsed to find the
  * mismatches, and the mismatched bases are recoreded. The reference base
  * is not kept track of.
  *
  * @param [in] sam_entry sam entry to be added.
  */
  void add(const SamEntry &sam_entry);

  /**
  * After all the sam entries are added, this method is used to retrive
  * the mismatched bases and the coverage at that base.
  *
  * @param [in] in_region Region to get the mismatches from.
  * @param [out] out_region One base wide region for each mismatch in the
  *   input region. 'aln_count' and 'mismatches' are parallel arrays
  *   corresponding the location in 'out_region'.
  * @param [out] aln_count Coverage at the position.
  * @param [out] mismatches A string with all the mismatched bases at the
  *   position. The length of this string is equal to the number of
  *   mismatched bases.
  */
  void at(const GenomicRegion &in_region,
          vector<GenomicRegion> &out_region,
          vector<size_t> &aln_count,
          vector<string> &mismatches) const;

private:
  GenomicStepVector<string> mismatches;
  GenomicStepVector<size_t> coverage;
};

#endif
