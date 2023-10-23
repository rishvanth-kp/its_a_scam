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

#include "SamEntry.hpp"
#include "SamReader.hpp"
#include "BedReader.hpp"
#include "GenomicRegion.hpp"
#include "FeatureVector.hpp"
#include "GenomicStepVector.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::unordered_map;


struct TssMetadata {
  string chrom;
  size_t start;
  size_t end;
  string tss_id;
  char strand;
};

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
      << "\t-t TSS bed file [required]" << endl
      << "\t-b barcode list (if empty, treated as bulk sample)" << endl
      << "\t-o outfile prefix [required]" << endl
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
    string tss_file;
    string bc_file; 
    string out_prefix;

    size_t side_dist = 1000;

    char bc_delim = ':';
    uint8_t bc_col = 7;
    string bc_tag;

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;    

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:t:s:b:o:d:c:t:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 't')
        tss_file = optarg;
      else if (opt == 's')
        side_dist = std::stoi(optarg);
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

    if (aln_file.empty() || tss_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }


    // process barcodes
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;
    
    // to index the TSS enhancement martix
    unordered_map<string, size_t> bc_index;
    size_t bc_counter = 0;
    vector<string> bc_metadata;    
    bool bulk_sample = false;  

    if (bc_file.empty()) {
      // if a barocde file is not proivded, treat it as a bulk sample
      // all the reads in the sam file will be considerd as one sample
      bulk_sample = true;
      bc_index[out_prefix] = bc_counter++;
      bc_metadata.push_back(out_prefix);
    }
    else {
      // if a barcode file is provided, each barcode is a separate sample
      std::ifstream bc_in(bc_file);
      string line;
      while (getline(bc_in, line)) {
        vector<string> tokens;
        split_string(line, tokens, '\t');
        bc_index[tokens[0]] = bc_counter++;
        bc_metadata.push_back(tokens[0]);
      } 
      bc_in.close();
    }

    // create matrix to store the coverage around TSS for each sample
    vector<vector<size_t>> tss_coverage(bc_counter,
                                        vector<size_t>(2*side_dist + 1, 0));

    if (VERBOSE) {
      if (bulk_sample) {
        cerr << "\tBulk sample: " << out_prefix << endl;
      }
      else {
        cerr << "\tNumber of barocdes: " << bc_counter << endl;
      }
    }
      

    // process TSS file
    if (VERBOSE)
      cerr << "[PROCESSING TSS REGIONS]" << endl;

    GenomicStepVector<FeatureVector<string>> tss;
    vector<TssMetadata> tss_metadata;

    size_t tss_counter = 0;
  
    BedReader tss_reader(tss_file);
    GenomicRegion bed_region;
    vector<string> bed_fields;
    while (tss_reader.read_bed_line(bed_region, bed_fields)) {

      // store the tss regrion
      // storing just the TSS name. The distance from the tss can be 
      // determined from the TSS start stored in the metadata.
      tss.add(bed_region.name, bed_region.start - side_dist, 
              bed_region.end + side_dist + 1, 
              FeatureVector<string>(bed_fields[0])); 
    
 
      // store TSS metadata
      TssMetadata metadata;
      metadata.chrom = bed_region.name;
      metadata.start = bed_region.start;
      metadata.end = bed_region.end;
      metadata.tss_id = bed_fields[0];
      metadata.strand = bed_fields[1][0];

      // keep track of number of tss
      ++tss_counter;
    }   

    if (VERBOSE) 
      cerr << "\tNumber of TSS: " << tss_counter << endl;


    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;
    
    

  } 
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl; 
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
