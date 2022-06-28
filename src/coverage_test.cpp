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

using std::cout;
using std::cerr;
using std::endl;
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
    cout << "string constructor" << endl;
    vec.push_back(in);
  }
  TestVec(const TestVec &in) {
    cout << "copy constructor" << endl;
    for (size_t i = 0; i < in.size(); ++i) {
      cout << in.at(i) << endl;
      vec.push_back(in.at(i));
    }
    cout << "end copy" << endl;
  }

  TestVec& operator+(const string a) {
    vec.push_back(a);
    return *this;
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

    GenomicStepVector<size_t> coverage;

    while (reader.read_sam_line(entry)) {
      if (entry.mapq >= min_mapq &&
          SamFlags::is_all_set(entry.flag, include_all) &&
          !SamFlags::is_any_set(entry.flag, include_none)) {

        const size_t seq_len = entry.seq.length();
        size_t clip_len = 0;

        SamCigar::CigarTuples tuples;
        SamCigar::string_to_tuple(entry, tuples);

        // 'H' can only be present as the first or last operation.
        // 'H' bases are not present in SEQ, skip any 'H'.
        // 'S' can only have 'H' between them and the end.
        // 'S' bases are present in SEQ, so remove them from aligned len.
        SamCigar::CigarTuples::iterator it = tuples.begin();
        if (it->first == 'H')
          ++it;
        if (it->first == 'S')
          clip_len += it->second;

        SamCigar::CigarTuples::reverse_iterator rit = tuples.rbegin();
        if (rit->first == 'H')
          ++rit;
        if (rit->first == 'S')
          clip_len += rit->second;

        size_t aln_len = seq_len - clip_len;
        coverage.add(entry.rname, entry.pos, entry.pos + aln_len, 1);
      }
    }

    GenomicRegion region("chr1", 0, 249250621);
    vector<pair<GenomicRegion, size_t>> out;
    // coverage.at(region, out);
    coverage.at(GenomicRegion{"chr1", 0, 249250621}, out);
    for (size_t i = 0; i < 5; ++i) {
      cout << out[i].first << "\t"
           << out[i].second << endl;
    }

/*
    TestVec a{"foo"};
    a + "foo1";
    a.print_elem();

    TestVec b{"bar"};
    b + "bar2";
    b + "bar3";
    b + "bar4";
    b.print_elem();

    a = a + b;
    a + "foo3";
    a + "foo4";
    a.print_elem();


    GenomicStepVector<TestVec> gsv;
    gsv.add("ch1", 3, 10, TestVec{"a"});
    gsv.add("ch1", 3, 10, TestVec{"b"});
    vector<pair<GenomicRegion, TestVec>> out2;
    // gsv.at(GenomicRegion("ch1", 3, 10), out2);
*/
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
