/*
* FeatureVector: Vector compatible with GenomicStepVector
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

#ifndef FEATURE_VECTOR_HPP
#define FEATURE_VECTOR_HPP

#include <iostream>
#include <vector>

using std::vector;

template<typename T>
class FeatureVector {
public:
  FeatureVector() = default;
  FeatureVector(const T &in);
  FeatureVector(const FeatureVector<T> &in);


private:
  vector<T> features;
};


#endif
