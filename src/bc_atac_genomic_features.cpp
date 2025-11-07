/*
* bc_feature_matrix: count barcodes aligning to genomic features
* Copyright (C) 2023 Rishvanth Prabakar
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
#include <unistd.h>

#include "gcatlib/SamEntry.hpp"
#include "gcatlib/SamReader.hpp"
#include "gcatlib/AlignedOverlapGenomicFeature.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;

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


static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-b barcode list file [required]" << endl
      << "\t-g gtf file [required]" << endl
      << "\t-r tsr bed file [optional]" << endl
      << "\t-e enhancers bed file [optional]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-m min fragment length [default: 0]" << endl
      << "\t-M max fragment length [default: inf]" << endl
      << "\t-d name split delimeter " 
          << "[default: \":\"; ignored if -t is provided]" << endl
      << "\t-c barcode field in name " 
          << "[default: 7 (0 based); ignored if -t is provided]" << endl
      << "\t-t barcode tag in SAM/BAM file [default \"\"]" << endl 
      << "\t-q minimum mapping quality to include [default: 0]" << endl
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

    string gtf_file;
    string tsr_file;
    string enhancer_file;

    char bc_delim = ':';
    char bc_col = 7;
    string bc_tag;

    size_t min_frag_len = 0;
    size_t max_frag_len = std::numeric_limits<size_t>::max();

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:g:r:e:o:m:M:d:c:t:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'g')
        gtf_file = optarg;
      else if (opt == 'r')
        tsr_file = optarg;
      else if (opt == 'e')
        enhancer_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'm')
        min_frag_len = std::stoi(optarg);
      else if (opt == 'M')
        max_frag_len = std::stoi(optarg);
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

    if (aln_file.empty() || bc_file.empty() ||
          gtf_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    // process gtf and bed file with regions
    if (VERBOSE)
      cerr << "[PROCESSING REGIONS]" << endl;

    AlignedOverlapGenomicFeature aligned_feature;
    aligned_feature.add_gtf_features(gtf_file);

    if (!tsr_file.empty()) {
      aligned_feature.add_bed_features(tsr_file, "TSR");
    }
    if (!enhancer_file.empty()) {
      aligned_feature.add_bed_features(enhancer_file, "enhancer");
    }

    // process barcodes
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;
    aligned_feature.process_barcodes(bc_file);

    // set fragment lengths
    aligned_feature.set_min_frag_len(min_frag_len);
    aligned_feature.set_max_frag_len(max_frag_len);

    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;
    SamReader sam_reader(aln_file);
    SamEntry entry1, entry2;

    size_t aln_count = 0;

    while (sam_reader.read_pe_sam(entry1, entry2)) {
      ++aln_count;
      if (VERBOSE) {
        if (!(aln_count % 1000000)) {
          cerr << "\tprocessed " << aln_count << " pairs" << endl;
        }
      }

      // filter out low reads that do not meet criteria
      if (entry1.mapq >= min_mapq && entry2.mapq >= min_mapq &&
          SamFlags::is_all_set(entry1.flag, include_all) &&
          SamFlags::is_all_set(entry2.flag, include_all) &&
          !SamFlags::is_any_set(entry1.flag, include_none) &&
          !SamFlags::is_any_set(entry2.flag, include_none)) {

        
        // get the cell barcode, either from a tag or the name
        string cell_bc;
        if (!bc_tag.empty()) {
          // read barcode form the tag
          SamTags::get_tag(entry1.tags, bc_tag, cell_bc);
        }
        else {
          // parse the name to get the bc
          vector<string> tokens;
          split_string(entry1.qname, tokens, bc_delim);
          cell_bc = tokens[bc_col];
        }

        // figure out feature location
        aligned_feature.add(entry1, entry2, cell_bc);
      }
    }

    // write output to file
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;
    aligned_feature.feature_counts_to_file(out_prefix);
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
