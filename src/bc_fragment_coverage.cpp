/*
* Copyright (C) 2024 Rishvanth Prabakar
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
#include <algorithm>
#include <unordered_map>

#include "SamEntry.hpp"
#include "SamReader.hpp"
#include "FeatureVector.hpp"
#include "GenomicStepVector.hpp"


using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;


static void
split_string (const string &in, vector<string> &tokens,
              const char delim = ':') {

  tokens.clear();
  size_t start = 0;
  size_t end = in.find(delim);
  while (end != string::npos) {
    tokens.push_back(in.substr(start, end - start));
    start = ++end;
    end = in.find(delim, start);
  }
  tokens.push_back(in.substr(start));
}


static bool
sort_bc_index (pair<string, size_t> a, pair<string, size_t> b) {
  return (a.second < b.second);
}


static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-b barcode list [required]" << endl
      << "\t-o outfile prefix [required]" << endl
      << "\t-d name split delimeter "
          << "[default: \":\"; ignored if -t is provided]" << endl
      << "\t-c barcode field in name "
          << "[default: 7 (0 based); ignored if -t is provided]" << endl
      << "\t-q minimun mapping quality to include [default: 0]" << endl
      << "\t-f only include if all the flags are present [default: 3]" << endl
      << "\t-F only include if none of the flags are present [default: 3340]"
          << endl
      << "\t-v verbose [default: false]" << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {

    string aln_file;
    string bc_file;
    string out_prefix;

    char bc_delim = ':';
    uint8_t bc_col = 7;
    string bc_tag;

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:o:d:c:t:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'd')
        bc_delim = optarg[0];
      else if (opt == 'c')
        bc_col = std::stoi(optarg);
      else if (opt == 't')
        bc_tag = optarg;
      else if (opt == 'q')
        min_mapq = std::stoi(optarg);
      else if (opt == 'f')
        include_all = std::stoi(optarg);
      else if (opt == 'F')
        include_none = std::stoi(optarg);
      else if (opt == 'v')
        VERBOSE = true;
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (aln_file.empty() || bc_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    // process barcodes
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;

    unordered_map<string, size_t> bc_index;
    size_t bc_counter = 0;
    vector<string> bc_metadata;

    std::ifstream bc_in(bc_file);
    string line;
    while (getline(bc_in, line)) {
      vector<string> tokens;
      split_string(line, tokens, '\t');
      bc_index[tokens[0]] = bc_counter++;
      bc_metadata.push_back(tokens[0]);
    }

    if (VERBOSE) {
      cerr << "\tNumber of barcodes: " << bc_counter << endl;
    }



    // initialize sam reader
    if (VERBOSE)
      cerr << "[INITIALIZING SAM/BAM READER]" << endl;

    SamReader reader(aln_file);
    SamEntry entry1, entry2;

    // parse reference chroms from sam header
    vector<GenomicRegion> ref_chroms;
    string header;
    reader.read_sam_header(header);
    get_seq_lengths(header, ref_chroms);


    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;
    // genomic step vector to store coverage
    GenomicStepVector<FeatureVector<string>> coverage;
    size_t aln_count = 0;

    while (reader.read_pe_sam(entry1, entry2)) {

      ++aln_count;
      if (VERBOSE) {
        if (!(aln_count % 1000000)) {
          cerr << "\tprocessed " << aln_count << " alignments" << endl;
        }
      }

      if (entry1.mapq >= min_mapq && entry2.mapq >= min_mapq &&
          SamFlags::is_all_set(entry1.flag, include_all) &&
          SamFlags::is_all_set(entry2.flag, include_all) &&
          !SamFlags::is_any_set(entry1.flag, include_none) &&
          !SamFlags::is_any_set(entry2.flag, include_none)) {

        // get the cell barcode, either from a tag or the name
        string cell_bc;
        if (!bc_tag.empty()) {
          // read bc from the appropriate tag
          SamTags::get_tag(entry1.tags, bc_tag, cell_bc);
        }
        else {
          // parse the sam name to get the bc
          vector<string> tokens;
          split_string(entry1.qname, tokens, bc_delim);
          cell_bc = tokens[bc_col];
        }

        unordered_map<string, size_t>::iterator bc_it;
        bc_it = bc_index.find(cell_bc);
        if (bc_it != bc_index.end()) {

          // find the start and end location of the reads
          size_t e1_start = entry1.pos - 1;
          size_t e1_end = SamCigar::reference_end_pos(entry1);
          size_t e2_start = entry2.pos - 1;
          size_t e2_end = SamCigar::reference_end_pos(entry2);

          // find the stand and end location of the fragments
          size_t frag_start = e1_start;
          if (e2_start < e1_start)
            frag_start = e2_start;

          size_t frag_end = e1_end;
          if (e2_end > e1_end)
            frag_end = e2_end;


          // add the appriate bases for coverage
          coverage.add(entry1.rname, frag_start, frag_end,
                       FeatureVector<string>(cell_bc));


        }
      }
    }



    // write output
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;

    std::ofstream depth_file(out_prefix + "_bc_fragment_coverage.txt");
    // write header
    vector<pair<string, size_t>> bc_index_ordered;
    for (auto it = bc_index.begin(); it != bc_index.end(); ++it) {
      bc_index_ordered.push_back(std::make_pair(it->first, it->second));
    }
    std::sort(bc_index_ordered.begin(), bc_index_ordered.end(), sort_bc_index);

    depth_file << "chrom\tstart\tend";
    for (size_t i = 0; i < bc_counter; ++i) {
      depth_file << "\t" << bc_index_ordered[i].first;
    }
    depth_file << endl;

    vector<pair<GenomicRegion, FeatureVector<string>>> out;
    for (size_t i = 0; i < ref_chroms.size(); ++i) {
      coverage.at(ref_chroms[i], out);
      for (size_t j = 0; j < out.size(); ++j) {
        depth_file << out[j].first;
        // process the feature vector for each locaiton
        // Need to count how mant times each cell barcode occuers
        vector<size_t> pos_bc_count(bc_counter, 0);
        for (size_t k = 0; k < out[j].second.size(); ++k) {
          ++pos_bc_count[bc_index[out[j].second.at(k)]];
        }
        for (size_t k = 0; k < pos_bc_count.size(); ++k) {
          depth_file << "\t" << pos_bc_count[k];
        }
        depth_file << endl;
      }
    }

    depth_file.close();

  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
