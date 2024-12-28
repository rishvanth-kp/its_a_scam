/*
* add_clipped_bases: get soft clipped sequences from a SAM/BAM file
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
#include <sstream>
#include <fstream>
#include <unistd.h>

#include "gcatlib/SamReader.hpp"
#include "gcatlib/SamEntry.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file [required]" << endl
      << "\t-o outfile prefix [required]" << endl
      << "\t-q minimum quality score to include [default: 0]" << endl
      << "\t-f only include if all the flags are present [default: 0]" << endl
      << "\t-F only include if none of the flags are present [default: 1280]"
        << endl;
  return oss.str();
}

int 
main (int argc, char* argv[]) {
  try {
  
    string aln_file;
    string out_prefix;
  
    size_t min_mapq = 0;
    size_t include_all = 0x0000;
    size_t include_none = 0x0500; 

    int opt;
    while ((opt = getopt(argc, argv, "a:o:q:f:F:")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'q')
        min_mapq = std::stoi(optarg);
      else if (opt == 'f')
        include_all = std::stoi(optarg);
      else if (opt == 'F')
        include_none = std::stoi(optarg);
      else 
        throw std::runtime_error(print_usage(argv[0]));
    }
 
    if (aln_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }


    SamReader sam_reader(aln_file);
    SamEntry entry;

    std::ofstream out(out_prefix + "_clipped_bases.txt");
    std::ofstream out_sa(out_prefix + "_sa_tag.txt");

    while (sam_reader.read_sam_line(entry)) {
      if (entry.mapq >= min_mapq && 
          SamFlags::is_all_set(entry.flag, include_all) &&
          !SamFlags::is_any_set(entry.flag, include_none)) {

        // Only hard-clipped bases are allowed on either side of 
        // soft-clipped bases, so it has to the first or second 
        // cigar (or last or second to last).
        
        // parse cigar string
        SamCigar::CigarTuples cigar;
        SamCigar::string_to_tuple(entry, cigar);

        // check for clipped bases at the start of alignment
        size_t start_clip_len = 0;
        if (cigar[0].first == SamCigar::Cigar::soft_clip) {
          start_clip_len = cigar[0].second; 
        }
        else if (cigar[1].first == SamCigar::Cigar::soft_clip) {
          start_clip_len = cigar[1].second; 
        } 
          
        // get the start clipped bases
        string start_clip_bases = "*";
        if (start_clip_len) {
          start_clip_bases = entry.seq.substr(0, start_clip_len); 
        }

        // check for clipped bases at the end of alignment
        size_t end_clip_len = 0;
        if (cigar[cigar.size() - 1].first == SamCigar::Cigar::soft_clip) {
          end_clip_len = cigar[cigar.size() - 1].second; 
        }
        else if (cigar[cigar.size() - 2].first == SamCigar::Cigar::soft_clip) {
          end_clip_len = cigar[cigar.size() - 2].second;
        }

        // get the end clipped bases
        string end_clip_bases = "*";
        if (end_clip_len) {
          end_clip_bases = entry.seq.substr(entry.seq.length() - end_clip_len);
        }

        // wirte output if there is atlest one clipped base
        if (start_clip_len || end_clip_len) {
          out << entry.rname << "\t"
              << entry.pos << "\t"
              << SamCigar::reference_end_pos(entry) << "\t"
              << start_clip_len << "\t"
              << start_clip_bases << "\t"
              << end_clip_len << "\t"
              << end_clip_bases << endl; 
        }
      
        // print sa tags if found
        string sa_tag;
        if (SamTags::get_tag(entry.tags, "SA", sa_tag)) {
          out_sa << entry.rname << "\t"
                 << entry.pos << "\t"
                 << SamCigar::reference_end_pos(entry) << "\t"
                 << sa_tag << endl;
        } 
                    
      }
    }

    out.close();
    out_sa.close(); 
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
