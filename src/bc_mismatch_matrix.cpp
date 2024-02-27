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
#include <string>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <algorithm>
#include <unordered_map>

#include "SamEntry.hpp"
#include "SamReader.hpp"
#include "FeatureVector.hpp"
#include "GenomicStepVector.hpp"
#include "AlignmentMismatches.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::unordered_map;


static void
count_nucs(const string &nucs,
           size_t &a_count, size_t &t_count,
           size_t &g_count, size_t &c_count) {

  a_count = 0;
  t_count = 0;
  g_count = 0;
  c_count = 0;

  for (size_t i = 0; i < nucs.length(); ++i) {
    if (nucs[i] == 'A' || nucs[i] == 'a')
      ++a_count;
    else if (nucs[i] == 'T' || nucs[i] == 't')
      ++t_count;
    else if (nucs[i] == 'G' || nucs[i] == 'g')
      ++g_count;
    else if (nucs[i] == 'C' || nucs[i] == 'c')
      ++c_count;
  }
}

static void
add_cell_to_sample(const AlignmentMismatch &cell_mm,
                   const size_t min_depth, const float min_vaf, 
                   const vector<GenomicRegion> &ref_chroms,
                   const size_t id,
                   GenomicStepVector<FeatureVector<string>> &sample_mm) {

  
  // process each chromosome
  for(size_t i = 0; i < ref_chroms.size(); ++i) {
    vector<GenomicRegion> out_region;
    vector<size_t> out_cov;
    vector<string> out_mm;
    cell_mm.at(ref_chroms[i], out_region, out_cov, out_mm);
  
    for (size_t j = 0; j < out_region.size(); ++j) {
      if (out_cov[j] >= min_depth) {
  
        // count the number of each base from the string
        size_t a_count, t_count, g_count, c_count;
        count_nucs(out_mm[j], a_count, t_count, g_count, c_count);

        string max_base;
        float max_vaf = 0;
        // A
        const float a_vaf = static_cast<float>(a_count) /
                            static_cast<float>(out_cov[j]);
        if (a_vaf > max_vaf) {
          max_base = "A";
          max_vaf = a_vaf;
        }
        // T
        const float t_vaf = static_cast<float>(t_count) /
                            static_cast<float>(out_cov[j]);
        if (t_vaf > max_vaf) {
          max_base = "T";
          max_vaf = t_vaf;
        }
        // G
        const float g_vaf = static_cast<float>(g_count) /
                            static_cast<float>(out_cov[j]);
        if (g_vaf > max_vaf) {
          max_base = "G";
          max_vaf = g_vaf;
        }
        // C
        const float c_vaf = static_cast<float>(c_count) /
                            static_cast<float>(out_cov[j]);
        if (c_vaf > max_vaf) {
          max_base = "C";
          max_vaf = c_vaf;
        }


        // add if vaf greater than threshold
        if (max_vaf >= min_vaf) {
          sample_mm.add(out_region[j], 
                        FeatureVector<string>(max_base + std::to_string(id)));
        }

/*
        // A
        const float a_vaf = static_cast<float>(a_count) /
                            static_cast<float>(out_cov[j]);
        if (a_vaf >= min_vaf) {
          sample_mm.add(out_region[j], 
                        FeatureVector<string>("A" + std::to_string(id)));
        }
        
        // T
        const float t_vaf = static_cast<float>(t_count) /
                            static_cast<float>(out_cov[j]);
        if (t_vaf >= min_vaf) {
          sample_mm.add(out_region[j], 
                        FeatureVector<string>("T" + std::to_string(id)));
        }

        // G
        const float g_vaf = static_cast<float>(g_count) /
                            static_cast<float>(out_cov[j]);
        if (g_vaf >= min_vaf) {
          sample_mm.add(out_region[j], 
                        FeatureVector<string>("G" + std::to_string(id)));
        }

        // C
        const float c_vaf = static_cast<float>(c_count) /
                            static_cast<float>(out_cov[j]);
        if (c_vaf >= min_vaf) {
          sample_mm.add(out_region[j], 
                        FeatureVector<string>("C" + std::to_string(id)));
        }
*/


      }
    }

  }
}
 

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
      << "\t-b barcode file [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-d minimun depth at base to consider [default 1]" << endl
      << "\t-l minimum VAF at base to consider [default: 0]" << endl
      << "\t-t cell barcode tag [default: CB]" << endl
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all of the flags are present [default: 0]"
        << endl
      << "\t-F only include if none of the flags are present [default: 2058]"
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

    size_t min_mapq = 0;
    size_t min_depth = 1;
    float min_vaf = 0;

    string bc_tag = "CB";

    uint16_t include_all = 0;
    uint16_t include_none = 0x804;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:b:o:d:l:t:q:f:F:v")) != -1) {
      if (opt == 'a')
        aln_file = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'd')
        min_depth = std::stoi(optarg);
      else if (opt == 'l')
        min_vaf = std::stof(optarg);
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
    
    // index to keep track of barcodes
    unordered_map<string, size_t> bc_index;
    size_t bc_counter = 0;
    // store to barcode string to write to output
    vector<string> bc_metadata;

    std::ifstream bc_in(bc_file);
    if (!bc_in) {
      throw std::runtime_error(bc_file + " does not exist.");
    }
    string line;
    while (getline(bc_in, line)) {
      vector<string> tokens;
      split_string(line, tokens, '\t');
      bc_index[tokens[0]] = bc_counter++;
      bc_metadata.push_back(tokens[0]);
    }
    bc_in.close();

    if (VERBOSE) {
      cerr << "\tNumber of barcodes: " << bc_counter << endl;
    }


    // process genomic regions from sam file
    SamReader sam_reader(aln_file);
    SamEntry e;

    if (VERBOSE)
      cerr << "[PROCESSING SAM/BAM HEADER]" << endl;
    string header;
    vector<GenomicRegion> ref_chroms;
    sam_reader.read_sam_header(header);
    get_seq_lengths(header, ref_chroms);

    
    // process alignments
    if (VERBOSE)
      cerr << "[PROCESSING ALIGNMENTS]" << endl;

    bool first = true;
    string prev_bc;
    size_t prev_bc_id;

    size_t read_counter = 0;

    GenomicStepVector<FeatureVector<string>> sample_mm;
    AlignmentMismatch* cell_mm;
      
    unordered_map<string, size_t>::iterator bc_it;

    while (sam_reader.read_sam_line(e)) {
      string e_bc;
      SamTags::get_tag(e.tags, bc_tag, e_bc);

      bc_it = bc_index.find(e_bc);

      if (bc_it != bc_index.end()) {

        if (first) {
          prev_bc = e_bc;
          prev_bc_id = bc_it->second;
          cell_mm = new AlignmentMismatch();
          first = false;
        }

        if (e.mapq >= min_mapq &&
            SamFlags::is_all_set(e.flag, include_all) &&
            !SamFlags::is_any_set(e.flag, include_none)) {

          // new cell
          if (prev_bc != e_bc) {
            if (VERBOSE) {
              cerr << prev_bc << ": " << read_counter << endl;
            }

            // add previous cell to sample_mm
            add_cell_to_sample(*cell_mm, min_depth, min_vaf, 
                               ref_chroms, prev_bc_id, sample_mm);

            // destroy 
            delete cell_mm;

            // reset every thing for next cell
            cell_mm = new AlignmentMismatch();
            prev_bc = e_bc;
            prev_bc_id = bc_it->second;
            read_counter = 1;

            // add first entry of next cell to mismatchs
            cell_mm->add(e);
          }
          // same old cell
          else {
            // add to mismatches
            cell_mm->add(e);

            ++read_counter;
          }

        }

      }

    }


    // add last cell info
    if (VERBOSE) {
      cerr << prev_bc << ": " << read_counter << endl;
    }
    add_cell_to_sample(*cell_mm, min_depth, min_vaf, ref_chroms, 
                        prev_bc_id, sample_mm);


    // destroy
    delete cell_mm; 


    // write output
    if (VERBOSE)
      cerr << "[WRITING MISMATCH OUTPUT]" << endl;

    std::ofstream out_file(out_prefix + "_mm_matrix.txt");

    // write header
    vector<pair<string, size_t>> bc_index_ordered;
    for (auto it = bc_index.begin(); it != bc_index.end(); ++it) {
      bc_index_ordered.push_back(std::make_pair(it->first, it->second));
    }
    std::sort(bc_index_ordered.begin(), bc_index_ordered.end(), sort_bc_index); 
    
    out_file << "chrom\tstart\tend";
    for (size_t i = 0; i < bc_index_ordered.size(); ++i) {
      out_file << "\t" << bc_index_ordered[i].first;
    }
    out_file << endl;

    // write the rest of the matrix
    vector<pair<GenomicRegion, FeatureVector<string>>> out_mm;
    for (size_t i = 0; i < ref_chroms.size(); ++i) {
      sample_mm.at(ref_chroms[i], out_mm);
      for (size_t j = 0; j < out_mm.size(); ++j) {
        out_file << out_mm[j].first;
        vector<string> mm_at_base(bc_index_ordered.size());
        for (size_t k = 0; k < out_mm[j].second.size(); ++k) {
          char mm_base = out_mm[j].second.at(k)[0];
          size_t mm_sample = std::stoi(out_mm[j].second.at(k).substr(1));
          mm_at_base[mm_sample] += mm_base; 
        }
        for (size_t k = 0; k < mm_at_base.size(); ++k) {
          if (mm_at_base[k].length()) {
            out_file << "\t" << mm_at_base[k];
          }
          else {
            out_file << "\t.";
          }
        }
        out_file << endl;
      }
    }
    
    out_file.close();

    if (VERBOSE)
      cerr << "[COMPUTING JACCARD INDEX]" << endl;
    
    vector<vector<size_t>> mm_counts(bc_index_ordered.size(), 
                                     vector<size_t>(bc_index_ordered.size(), 0));

    vector<pair<GenomicRegion, FeatureVector<string>>> out;
    for (size_t i = 0; i < ref_chroms.size(); ++i) {
      sample_mm.at(ref_chroms[i], out);
      for (size_t j = 0; j < out.size(); ++j) {

        // process the union and intersecon for each region
        const size_t region_len = out[j].first.end - out[j].first.start;
        for (size_t k = 0; k < out[j].second.size(); ++k) {
          // add to union
          const char mm1 = out[j].second.at(k)[0];
          const size_t sample1 = std::stoi(out[j].second.at(k).substr(1));
          mm_counts[sample1][sample1] += region_len;
          for (size_t l = k + 1; l < out[j].second.size(); ++l) {
            const char mm2 = out[j].second.at(l)[0];
            const size_t sample2 = std::stoi(out[j].second.at(l).substr(1));
            if ((mm1 == mm2) && (sample1 != sample2)) {
              mm_counts[sample1][sample2] += region_len;
            }

          }
        }

      }
    }

    // comput and write the jaccard index
    std::ofstream jaccard_out(out_prefix + "_jaccard.txt");
    for (size_t i = 0; i < bc_index_ordered.size(); ++i) {
      jaccard_out << bc_index_ordered[i].first << "\t";
    }
    jaccard_out << endl;
    for (size_t i = 0; i < bc_index_ordered.size(); ++i) {
      jaccard_out << bc_index_ordered[i].first << "\t";
      for (size_t j = 0; j < bc_index_ordered.size(); ++j) {
        size_t isect;
        if (i != j)
          isect = mm_counts[i][j] + mm_counts[j][i];
        else
          isect = mm_counts[i][j];
        // subtracting isect to account for double couting
        const size_t uni = mm_counts[i][i] + mm_counts[j][j] - isect;

        float jaccard = 0;
        if (uni) {
          jaccard = static_cast<float>(isect) / static_cast<float>(uni);
        }
   
        jaccard_out << jaccard << "\t";
      }
      jaccard_out << endl;
    }

    jaccard_out.close();

  }


  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
