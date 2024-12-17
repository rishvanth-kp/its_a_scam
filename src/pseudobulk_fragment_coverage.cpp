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


static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-b barcode list file [required]" << endl
      << "\t-r regions bed file [required]" << endl
      << "\t-s distance to add to either side of region [default: 0]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-d name split delimeter" 
          << "[default: \":\"]; ignored if -t is provided" << endl
      << "\t-c barcode field in name" 
          << "[default: 7 (0 based); ignored if -t is provided]" << endl
      << "\t-t barcode tag in SAM file [default \"\"]" << endl
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
  
    // parse args
    string aln_file;
    string bc_file;
    
    string regions_file;
    size_t side_dist = 0;    

    string out_prefix;

    char bc_delim = ':';
    uint8_t bc_col = 7;
    string bc_tag;

    size_t min_mapq = 0;
    size_t include_all = 0x0003;
    size_t include_none = 0x0D0C;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:r:s:o:d:c:t:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'r')
        regions_file = optarg;
      else if (opt == 's')
        side_dist = std::stoi(optarg);
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


    if (aln_file.empty() || bc_file.empty() || regions_file.empty() ||
        out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    // process barcode and pseudo-bulk info
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;

    unordered_map<string, string> bc_group;
    unordered_map<string, uint8_t> group_index;
    size_t group_counter = 0;
    size_t bc_counter = 0;

    std::ifstream bc_in(bc_file);
    string line;
    while (getline(bc_in, line)) {
      vector<string> tokens;
      split_string(line, tokens, '\t');
      
      // keep track of the group for each barcode
      bc_group[tokens[0]] = tokens[1];
      ++bc_counter;

      // if a new group is encountered, create a new index for it
      auto it = group_index.find(tokens[1]);
      if (it == group_index.end()) {
        group_index[tokens[1]] = group_counter++;
      }
    }    

    bc_in.close();
   
    if (VERBOSE) {
      cerr << "\tNumber of barcodes: " << bc_counter << endl;
      cerr << "\tNumber of groups: " << group_counter << endl;
    } 

  
    // process regions bed file
    if (VERBOSE) 
      cerr << "[PROCESSING REGIONS]" << endl;

    GenomicStepVector<uint8_t> regions;
    BedReader regions_reader(regions_file);
    GenomicRegion region_pos;

    while(regions_reader.read_bed3_line(region_pos)) {
      // store the regions in a feature vetor for look up later.
      // no need to store any info about the region, just need to know if the 
      // a given genomc location needs to be processed. 
      regions.add(region_pos.name, region_pos.start - side_dist, 
        region_pos.end + side_dist, 1); 
    }  

    // initialize sam reader
    if (VERBOSE) 
      cerr << "[INITIALIZING SAM/BAM READER]" << endl;
    
    SamReader reader(aln_file);
    SamEntry entry1, entry2;
  
    // process reference chroms from sam header
    // used later for retrieving coverage
    vector<GenomicRegion> ref_chroms;
    string header;
    reader.read_sam_header(header);
    get_seq_lengths(header, ref_chroms);
   

    // process alignments
    if (VERBOSE) 
      cerr << "[PROCESSING ALIGNMENTS]" << endl;

    // genomic step vector for storing coverage
    GenomicStepVector<FeatureVector<uint8_t>> coverage;
    // to keep track of the counts
    vector<size_t> group_frag_counts(group_counter, 0);
    vector<size_t> group_base_counts(group_counter, 0);

    size_t aln_count = 0;
    while (reader.read_pe_sam(entry1, entry2)) {
      ++aln_count;
      if (VERBOSE) {
        if (!(aln_count % 1000000)) {
          cerr << "\tprocessed " << aln_count << " alignments" << endl;   
        }
      }  

      // make sure that the fragmmt passes qc
      if (entry1.mapq >= min_mapq && entry2.mapq >= min_mapq && 
          SamFlags::is_all_set(entry1.flag, include_all) &&
          SamFlags::is_all_set(entry2.flag, include_all) &&
          !SamFlags::is_any_set(entry1.flag, include_none) &&
          !SamFlags::is_any_set(entry2.flag, include_none))  {

        // get the cell barcode, either from a tag or the name
        string cell_bc;
        if (!bc_tag.empty()) {
          // read bc from the appropriate tag
          SamTags::get_tag(entry1.tags, bc_tag, cell_bc);
        } 
        else {
          // parse the name to get the bc
          vector<string> tokens;
          split_string(entry1.qname, tokens, bc_delim);
          cell_bc = tokens[bc_col];
        }

        // check the if the barcode needs to be processed
        unordered_map<string, string>::iterator bc_it;
        bc_it = bc_group.find(cell_bc);
        if (bc_it != bc_group.end()) {
          string cell_group = bc_it->second;

          // find the start and end location of each read
          size_t e1_start = entry1.pos - 1;
          size_t e1_end = SamCigar::reference_end_pos(entry1);
          size_t e2_start = entry2.pos - 1;
          size_t e2_end = SamCigar::reference_end_pos(entry2);
  
          // fing the start and end location of the fragment
          size_t frag_start = e1_start;
          if (e2_start < e1_start)
            frag_start = e2_start;

          size_t frag_end = e1_end;
          if (e2_end >= e1_end) 
            frag_end = e2_end;
  
          // check if the fragment overlaps a region
          vector<pair<GenomicRegion, uint8_t>> out;
          regions.at(GenomicRegion(entry1.rname, frag_start, frag_end), out);
          if (out.size() > 0) { 
            // added the fragment for coverage
            coverage.add(entry1.rname, frag_start, frag_end, 
                         FeatureVector<uint8_t>(group_index[cell_group]));
            // keep track of counts
            group_base_counts[group_index[cell_group]] += (frag_end - frag_start);
            ++group_frag_counts[group_index[cell_group]];
          }
        }
        
      }

    }

    
    // write output
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;


    // write the group coverage
    std::ofstream depth_file(out_prefix + "_pseudobulk_fragment_coverage.txt");
    // write header
    vector<string> group_index_ordered(group_counter, "");
    for(auto it = group_index.begin(); it != group_index.end(); ++it) {
      group_index_ordered[it->second] = it->first;
    }

    depth_file << "chrom\tstart\tend";
    for (size_t i = 0; i < group_counter; ++i) {
      depth_file << "\t" << group_index_ordered[i];
    }  
    depth_file << endl;

    // write the rest of the file
    vector<pair<GenomicRegion, FeatureVector<uint8_t>>> out;
    for (size_t i = 0; i < ref_chroms.size(); ++i) {
      coverage.at(ref_chroms[i], out);
      for (size_t j = 0; j < out.size(); ++j) {
        depth_file << out[j].first;
        // process the feature vector for each location
        // need to count how many times each group occurs
        vector<size_t> pos_group_count(group_counter, 0);
        for (size_t k = 0; k < out[j].second.size(); ++k) {
          ++pos_group_count[out[j].second.at(k)]; 
        }
        for (size_t k = 0; k < pos_group_count.size(); ++k) {
          depth_file << "\t" << pos_group_count[k];
        }
        depth_file << endl;
      }
    }

    depth_file.close();


    // write the group counts
    std::ofstream counts_file(out_prefix + "_pseudobulk_fragment_counts.txt");
    // write header
    counts_file << "group\tfrag_count\tbase_count" << endl;

    // write counts for each group
    for (size_t i = 0; i < group_counter; ++i) {
      counts_file << group_index_ordered[i]
                  << "\t" << group_frag_counts[i]
                  << "\t" << group_base_counts[i] << endl;
    }

    counts_file.close();


    
  } 
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
