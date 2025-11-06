/*
* sc_rnaseq_feature_matrix: count barcodes aligning to genomic features
*   in scRNA-seq data
*
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
#include <unistd.h>

#include "gcatlib/SamEntry.hpp"
#include "gcatlib/SamReader.hpp"
#include "gcatlib/AlignedGenomicFeature.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;


static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-b barcode list file [required]" << endl
      << "\t-g GTF file [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-c cell barcode tag in SAM/BAM file [default: CB]" << endl
      << "\t-u UMI tag in SAM/BAM file [default: \"\"]" << endl
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all the flags are present [default: 0]" << endl
      << "\t-F only include if none of the flags are present [defauly: 2308]" 
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
  
    string bc_tag = "CB";
    string umi_tag;

    size_t min_mapq = 0;
    size_t include_all = 0;
    size_t include_none = 0x0904;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:g:o:c:u:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'g')
        gtf_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'c')
        bc_tag = optarg;
      else if (opt == 'u')
        umi_tag = optarg;
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


    // process gtf file
    if (VERBOSE)
      cerr << "[PROCESSING GTF FILE]" << endl;

    AlignedGenomicFeature aligned_feature;
    aligned_feature.preprocess_gff(gtf_file);
  

    // process barcodes
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;
    
    aligned_feature.process_barcodes(bc_file);

    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;
      
    SamReader sam_reader(aln_file);
    SamEntry entry;

    size_t aln_count = 0;
    
    while (sam_reader.read_sam_line(entry)) {
      ++aln_count;
      if (VERBOSE) {
        if (!(aln_count % 1000000)) {
          cerr << "\tprocessed " << aln_count << " pairs" << endl;
        }
      }

      // filter out low quality alignments and that do not meet
      // sam flag criteria
      if (entry.mapq >= min_mapq &&
          SamFlags::is_all_set(entry.flag, include_all) &&
          !SamFlags::is_any_set(entry.flag, include_none)) {

      
        // get barcode
        string cell_bc;
        SamTags::get_tag(entry.tags, bc_tag, cell_bc);       
        
        // added to features
        if (umi_tag.empty()) {
          aligned_feature.add(entry, cell_bc);
        }
        else {
          string umi;
          SamTags::get_tag(entry.tags, umi_tag, umi);
          aligned_feature.add(entry, cell_bc, umi);
        }      
  
      }

    }    

    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;

    // write output
    aligned_feature.feature_count_to_file(
      out_prefix + "_gex_feature_counts.txt");

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
