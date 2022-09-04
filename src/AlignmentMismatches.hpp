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

  

private:
  
}; 

#endif
