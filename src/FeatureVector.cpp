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

#include "FeatureVector.hpp"

template<typename T>
FeatureVector<T>::FeatureVector(const T &in) {
  features.push_back(in);
}


template<typename T>
FeatureVector<T>::FeatureVector(const FeatureVector<T> &in) {
  features.clear();
  for (size_t i = 0; i < in.size(); ++i) {
    features.push_back(in[i]);
  }  
}
