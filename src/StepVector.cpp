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

#include "StepVector.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::pair;

StepVector::StepVector() {
  cout << "Step vector init" << endl;
}
  

void 
StepVector::add(const size_t start, 
                const size_t end, 
                const size_t val) {

  cout << endl;
  cout << "Adding: " << start << " " << end << " " << val << endl;

  // Inserting at the end
  map<size_t, size_t>::iterator end_it = step_vec.lower_bound(end);
  if (end_it == step_vec.end()) {
    cout << "End next element at end" << endl;
    end_it = step_vec.insert(pair<size_t, size_t>(end, 0)).first;
  }
  else if (!(end_it->first == end)) {
    size_t prev_val = (--end_it)->second;
    cout << "insering end prev val: " << prev_val << endl;
    end_it = step_vec.insert(pair<size_t, size_t>(end, prev_val)).first;
  }
  else {
    cout << "ending at an exitsting element" << endl;
  }
  
  // Inserting at the start
  map<size_t, size_t>::iterator start_it = step_vec.lower_bound(start);
  if (start_it->first == start) {
    // If the element exists, update it's value. 
    cout << "start exists" << endl;
    start_it->second += val;

    // Iterator points to the start element. Increment to point to 
    // next element
    ++start_it; 
  }
  else {
    // If an element does not exist, insert a new element with
    // the current value added to the previous element's value. 
    cout << "start does not exist" << endl;
    size_t prev_val;
    if (start_it == step_vec.begin()) {
      prev_val = 0;
    }
    else {
      prev_val = (--start_it)->second;
    }
    cout << "prev val: " << prev_val << endl;
    start_it = 
      step_vec.insert(pair<size_t, size_t>(start, prev_val + val)).first;

    // Iterator points to the inserted element 
    // Increment to point to next element
   ++start_it;
  }

  // Modifying values in [start + 1. end)
  cout << "inserting in the middle" << endl;
  cout << start_it->first << "\t" << end_it->first << endl;
  for (auto it = start_it; it != end_it; ++it) {
    size_t prev_val = it->second;
    cout << "prev val: " << prev_val << endl;
    it->second += val;
  }

  print_elements(); 
}


size_t
StepVector::at(const size_t pos) const {
  map<size_t, size_t>::const_iterator it = step_vec.upper_bound(pos);
  if (it == step_vec.begin())
    return 0;
  else 
    return (--it)->second;
  
}


void
StepVector::print_elements() {
  
  for (auto it = step_vec.begin(); it != step_vec.end(); ++it) {
    cout << it->first << "\t" << it->second << endl;
  }
} 
