/*
* alignment_length_distribution: compute histogram of alignment lengths
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

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <unistd.h>

#include "SamReader.hpp"
#include "SamEntry.hpp"
#include "GenomicStepVector.hpp"
#include "GenomicRegion.hpp"
#include "FeatureVector.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::string;
using std::vector;

static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all of the flags are present [default: 0]"
          << endl
      << "\t-F only include if none of the flags are present [default: 2052]"
          << endl
      << "\t-o out file prefix [required]" << endl;
  return oss.str();
}

class TestVec {
public:
  TestVec() = default;
  TestVec(const string &in) {
    vec.push_back(in);
  }
  TestVec(const TestVec &in) {
    vec.clear();
    for (size_t i = 0; i < in.size(); ++i) {
      vec.push_back(in.at(i));
    }
  }


  TestVec& operator+(const TestVec &in) {
    for (size_t i = 0; i < in.size(); ++i) {
      vec.push_back(in.at(i));
    }
    return *this;
  }

  TestVec& operator+= (const TestVec &in) {

    for (size_t i = 0; i < in.size(); ++i) {
      vec.push_back(in.at(i));
    }
    return *this;
  }

  bool operator!=(const TestVec &in) {
    bool mismatch = true;
    if (in.size() == vec.size()) {
      size_t i = 0;
      while ((i < in.size()) && (in.at(i) == vec[i]))
        ++i;
      if (i == in.size())
        mismatch = false;
    }

    return mismatch;
  }

  size_t size() const { return vec.size(); }
  string at(size_t i) const {
    return vec[i];
  }

  void print_elem() {
    for (size_t i = 0; i < vec.size(); ++i)
      cout << vec[i] << endl;
  }
private:
  vector<string> vec;
};

int
main(int argc, char* argv[]) {
  try {

    string aln_file;
    string out_prefix;
    size_t min_mapq = 0;

    uint16_t include_all = 0;
    uint16_t include_none = 0x0804;


    int opt;
    while ((opt = getopt(argc, argv, "a:q:f:F:m:o:")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'q')
        min_mapq = std::stoi(optarg);
      else if (opt == 'f')
        include_all = std::stoi(optarg);
      else if (opt == 'F')
        include_none = std::stoi(optarg);
      else if (opt == 'o')
        out_prefix = optarg;
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (aln_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    SamReader reader(aln_file);
    SamEntry entry;


    while (reader.read_sam_line(entry)) {
      if (entry.mapq >= min_mapq &&
          SamFlags::is_all_set(entry.flag, include_all) &&
          !SamFlags::is_any_set(entry.flag, include_none)) {


        SamCigar::CigarTuples tuples;
        SamCigar::string_to_tuple(entry, tuples);

        cout << entry.qname << endl;
        cout << entry.cigar << endl;
        string md_tag;
        SamTags::get_tag(entry.tags, "MD", md_tag);
        cout << md_tag << endl;
        vector<pair<size_t, string>> md_tuple;
        SamTags::md_to_tuple(md_tag, md_tuple);
        size_t ref_offset = 0;
        for (auto it = md_tuple.begin(); it != md_tuple.end(); ++it) {
          cout << it->first << "\t" << it->second << endl;
          ref_offset += it->first;
          if (it->second[0] == '^')
            ref_offset += (it->second.length() - 1);
          else if (it->second != "") {
            ++ref_offset;
            cout << "ref_offset: " << ref_offset << endl;
            size_t ref_pos = 0;
            size_t query_pos = 0;
            SamCigar::move_in_reference(tuples, ref_offset, ref_pos, query_pos);
            cout << ref_pos << "\t" << query_pos << "\t" << entry.seq[query_pos] << endl; 
          }
        }

      }
    }


    FeatureVector<string> f1{"f1"};
    FeatureVector<string> f2{"f2"};
    FeatureVector<string> f3;
    f3 = f1 + f2;
    for (size_t i = 0; i < f3.size(); ++i)
      cout << f3.at(i) << endl;

    GenomicStepVector<FeatureVector<string>> foo;
    foo.add("ch1", 5, 10, FeatureVector<string>{"foo"});
    foo.add("ch1", 7, 15, FeatureVector<string>{"bar"});
    foo.add("ch1", 1, 3, FeatureVector<string>{"fubar"});

    vector<pair<GenomicRegion, FeatureVector<string>>> out;
    foo.at(GenomicRegion{"ch1", 0, 20}, out, true);
    cout << "foo size: " << out.size() << endl;
    for (size_t i = 0; i < out.size(); ++i) {
      cout << out[i].first;
      for (size_t j = 0; j < out[i].second.size(); ++j) {
        cout << "\t" << out[i].second.at(j);
      }
      cout << endl;
    }

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
