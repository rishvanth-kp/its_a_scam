/*
* StepVector: 
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

#ifndef STEP_VECTOR_HPP
#define STEP_VECTOR_HPP

#include <iostream>
#include <vector>
#include <map>

using std::map;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::vector;

template<typename T>
class StepVector {
public:
  StepVector();
  
  // add element
  void add(const size_t start, const size_t end, const T val);
  // access elements in range
  void at_range(const size_t start, const size_t end, 
                vector<pair<size_t, T>>& out) const;
  // access value at a location
  T at(const size_t pos) const;

  // print elements
  void print_elements();
  

private:

  map<size_t, T> step_vec;
  bool VERBOSE;
};

template<typename T>
StepVector<T>::StepVector() {
  cout << "Step vector init" << endl;
}
  




template<typename T>
void 
StepVector<T>::add(const size_t start, 
                const size_t end, 
                const T val) {

  if (start < end) {

    // Inserting at the end
    // if the end postion is an an existing element, there is nothing to 
    // be done
    typename map<size_t, T>::iterator end_it = step_vec.lower_bound(end);
    // ending after the end of all existing entries, add default value.
    if (end_it == step_vec.end()) {
      // cout << "End next element at end" << endl;
      T default_val{};
      end_it = step_vec.insert(pair<size_t, T>(end, default_val)).first;
    }
    // ending before the end of exisiting entries, 
    else if (!(end_it->first == end)) {
      T prev_val = (--end_it)->second;
      // cout << "insering end prev val: " << prev_val << endl;
      end_it = step_vec.insert(pair<size_t, T>(end, prev_val)).first;
    }
    // else {
    //   cout << "ending at an exitsting element" << endl;
    // }
    
    // Inserting at the start
    typename map<size_t, T>::iterator start_it = step_vec.lower_bound(start);
    if (start_it->first == start) {
      // If the element exists, update it's value. 
      // cout << "start exists" << endl;
      start_it->second += val;

      // Iterator points to the start element. Increment to point to 
      // next element
      ++start_it; 
    }
    else {
      // If an element does not exist, insert a new element with
      // the current value added to the previous element's value. 
      // cout << "start does not exist" << endl;
      if (start_it == step_vec.begin()) {
        start_it = 
          step_vec.insert(pair<size_t, T>(start, val)).first;
      }
      else {
        T prev_val = (--start_it)->second;
        start_it = 
          step_vec.insert(pair<size_t, T>(start, prev_val + val)).first;
      }

      // Iterator points to the inserted element 
      // Increment to point to next element
     ++start_it;
    }

    // Modifying values in [start + 1. end)
    // cout << "inserting in the middle" << endl;
    // cout << start_it->first << "\t" << end_it->first << endl;
    for (auto it = start_it; it != end_it; ++it) {
      // T prev_val = it->second;
      // cout << "prev val: " << prev_val << endl;
      it->second += val;
    }

  }
}


template<typename T>
T
StepVector<T>::at(const size_t pos) const {
  typename map<size_t, T>::const_iterator it = step_vec.upper_bound(pos);
  if (it == step_vec.begin())
    return T{};
  else 
    return (--it)->second;
  
}

template<typename T>
void
StepVector<T>::at_range(const size_t start, const size_t end,
                        vector<pair<size_t, T>>& out) const {

  out.clear();

  if (start < end) {

    typename map<size_t, T>::const_iterator start_it, end_it;
    start_it = step_vec.upper_bound(start);
    if (start_it == step_vec.begin()) {
      out.push_back(std::make_pair(start, T{}));
    }
    else {
      out.push_back(std::make_pair(start, (--start_it)->second));
      ++start_it;
    }

    end_it = step_vec.upper_bound(end);
    for (auto it = start_it; it != end_it; ++it) {
      out.push_back(std::make_pair(it->first, it->second));
    }
   
    if (end_it == step_vec.begin()) {
      out.push_back(std::make_pair(end, T{}));
    }
    else {
      --end_it;
      if (end_it->first != end) 
        out.push_back(std::make_pair(end, end_it->second));
    }

  } 
}


template<typename T>
void
StepVector<T>::print_elements() {
  
  for (auto it = step_vec.begin(); it != step_vec.end(); ++it) {
    cout << it->first << "\t" << it->second << endl;
  }
  cout << endl;
}


#endif
