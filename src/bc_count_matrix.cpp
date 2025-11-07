/*
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
#include <fstream>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>

#include "gcatlib/SamEntry.hpp"
#include "gcatlib/SamReader.hpp"
#include "gcatlib/BedReader.hpp"
#include "gcatlib/GenomicRegion.hpp"
#include "gcatlib/FeatureVector.hpp"
#include "gcatlib/GenomicStepVector.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::unordered_map;
using std::unordered_set;



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


struct RegionMetadata {
  string chrom;
  size_t start;
  size_t end;
  string region_id;
};

static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-b barcode list file [required]" << endl
      << "\t-r region bed file [required]" << endl
      << "\t-s distance to add to either side of each region [default: 0]"
        << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-m min. fragment length [default: 0]" << endl
      << "\t-M max. fragment length [default: 1024]" << endl
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
    string region_file;
    string out_prefix;

    size_t side_dist = 0;
    
    size_t min_frag_len = 0;
    size_t max_frag_len = 1024;

    char bc_delim = ':';
    char bc_col = 7;
    string bc_tag;

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:r:s:o:m:M:d:c:t:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'r')
        region_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 's')
        side_dist = std::stoi(optarg);
      else if (opt == 'm')
        min_frag_len = std::stoi(optarg);
      else if (opt == 'M')
        max_frag_len = std::stoi(optarg);
      else if (opt == 'd')
        bc_delim = optarg[0];
      else if (opt == 't')
        bc_tag = optarg;
      else if (opt == 'c')
        bc_col = std::stoi(optarg);
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

    if (aln_file.empty() || bc_file.empty() || region_file.empty() ||
        out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }


    // progess the region bed file and construct the GenomicStepVector
    if (VERBOSE)
      cerr << "[PROCESSING REGIONS]" << endl;

    GenomicStepVector<FeatureVector<string>> region;
    vector<RegionMetadata> region_metadata;

    // to index into the count matrix
    unordered_map<string, size_t> region_index;
    size_t region_counter = 0;

    BedReader region_reader(region_file);
    GenomicRegion bed_region;
    vector<string> bed_fields;
    while(region_reader.read_bed_line(bed_region, bed_fields)) {
      // store the actual feature vector
      region.add(bed_region.name, bed_region.start - side_dist,
                  bed_region.end + side_dist, 
                  FeatureVector<string>(bed_fields[0]));

      // store the metadata for final output
      RegionMetadata metadata;
      metadata.chrom = bed_region.name;
      metadata.start = bed_region.start - side_dist;
      metadata.end = bed_region.end + side_dist;
      metadata.region_id = bed_fields[0];
      region_metadata.push_back(metadata);

      // assign a uniqe index for each region to index into the count matrix
      region_index[bed_fields[0]] = region_counter++;
    }


    // process barcodes and initialize bc count vector
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;

    // to index into the count matrix
    unordered_map<string, size_t> bc_index;
    size_t bc_counter = 0;
    vector<string> bc_metadata;

    std::ifstream bc_in(bc_file);
    string line;
    while(getline(bc_in, line)) {
      vector<string> tokens;
      split_string(line, tokens, '\t');
      bc_index[tokens[0]] = bc_counter++;
      bc_metadata.push_back(tokens[0]);
    }
    bc_in.close();

    vector<vector<size_t>> region_counts(region_counter,
                                      vector<size_t>(bc_counter, 0));


    if (VERBOSE) {
      cerr << "\tNumber of regions: " << region_counter << endl;
      cerr << "\tNumber of barcodes: " << bc_counter << endl;
    }


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


        unordered_map<string, size_t>::iterator it;
        it = bc_index.find(cell_bc);
        if (it != bc_index.end()) {

          // find the start and end postion of the reads
          size_t e1_start = entry1.pos - 1;
          size_t e1_end = SamCigar::reference_end_pos(entry1);
          size_t e2_start = entry2.pos - 1;
          size_t e2_end = SamCigar::reference_end_pos(entry2);

          // find the fragment start and end positoins
          size_t frag_start = e1_start;
          if (e2_start < e1_start)
            frag_start = e2_start;

          size_t frag_end = e1_end;
          if (e2_end > e1_end)
            frag_end = e2_end;

          // check if it is in the desired fragment length 
          size_t frag_len = frag_end - frag_start;           
          if (frag_len >= min_frag_len && frag_len <= max_frag_len) {

            // check if fragment overlaps a region
            unordered_set<string> aligned_region;
            vector<pair<GenomicRegion, FeatureVector<string>>> out;
            region.at(GenomicRegion(entry1.rname, frag_start, frag_end), out);
            for (auto jt = out.begin(); jt != out.end(); ++jt) {
              for (size_t k = 0; k < jt->second.size(); ++k) {
                aligned_region.insert(jt->second.at(k));
              }
            }

            // increment count if a frament as aligned to a unique location
            if (aligned_region.size() == 1) {
              size_t col_index = bc_index[cell_bc];
              size_t row_index = region_index[*aligned_region.cbegin()];
              ++region_counts[row_index][col_index];
            }

          }

        }
      }
    }

    // write output to file
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;
    std::ofstream out(out_prefix + "_region_counts.txt");
    // write header line
    out << "chrom\tstart\tend\tregion\t";
    for (size_t i = 0; i < bc_metadata.size(); ++i) {
      out << bc_metadata[i] << "\t";
    }
    out << endl;
    // write the matrix
    for (size_t i = 0; i < region_metadata.size(); ++i) {
      out << region_metadata[i].chrom << "\t"
          << region_metadata[i].start << "\t"
          << region_metadata[i].end << "\t"
          << region_metadata[i].region_id << "\t";
      for (size_t j = 0; j < bc_metadata.size(); ++j) {
        out << region_counts[i][j] << "\t";
      }
      out << endl;
    }
    out.close();


  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
