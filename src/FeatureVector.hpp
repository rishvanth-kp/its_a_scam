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
  explicit FeatureVector(const T &in);
  FeatureVector(const FeatureVector<T> &in);

  FeatureVector<T> operator+(const FeatureVector<T> &in);
  void operator+=(const FeatureVector<T> in);
  bool operator==(const FeatureVector &in) const;
  bool operator!=(const FeatureVector &in) const;

  void push_back(const T in);
  T at(const size_t i) const;
  size_t size() const;

private:
  vector<T> features;
};



template<typename T>
FeatureVector<T>::FeatureVector(const T &in) {
  features.push_back(in);
}


template<typename T>
FeatureVector<T>::FeatureVector(const FeatureVector<T> &in) {
  features.clear();
  for (size_t i = 0; i < in.size(); ++i) {
    features.push_back(in.at(i));
  }
}


template<typename T>
FeatureVector<T>
FeatureVector<T>::operator+(const FeatureVector<T> &in) {
  FeatureVector<T> out;
  for (size_t i = 0; i < features.size(); ++i) {
    out.push_back(features[i]);
  }
  for (size_t i = 0; i < in.size(); ++i) {
    out.push_back(in.at(i));
  }
  return out;
}


template<typename T>
void
FeatureVector<T>::operator+=(const FeatureVector<T> in) {
  for (size_t i = 0; i < in.size(); ++i) {
    features.push_back(in.at(i));
  }
}


template<typename T>
bool
FeatureVector<T>::operator==(const FeatureVector &in) const {
  bool match = false;
  if (in.size() == features.size()) {
    size_t i = 0;
    while ((i < features.size()) && (in.at(i) == features[i]))
      ++i;
    if (i == features.size())
      match = true;
  }
  return match;
}


template<typename T>
bool
FeatureVector<T>::operator!=(const FeatureVector<T> &in) const {
  return !(*this == in);
}


template<typename T>
void
FeatureVector<T>::push_back(const T in) {
  features.push_back(in);
}

template<typename T>
T
FeatureVector<T>::at(const size_t i) const {
  return features[i];
}

template<typename T>
size_t
FeatureVector<T>::size() const {
  return features.size();
}

#endif
