/*
* mismatch_profile: report mismatch location ans statistics
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
#include <unordered_map>

#include "SamEntry.hpp"
#include "SamReader.hpp"
#include "AlignmentMismatches.hpp"
#include "GenomicStepVector.hpp"
#include "FeatureVector.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
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

static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file list [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-d minimum depth at a base to consider [default: 1]" << endl
      << "\t-v minimum VAF at a base to consider [default: 0.0]" << endl
      << "\t-q minimum mapping quality to include [default: 0]" << endl
      << "\t-f only include if all of the flags are present [default: 0]"
          << endl
      << "\t-F only include if none of the flags are present [default: 2052]"
          << endl;
    return oss.str();
}

int
main(int argc, char* argv[]) {
  try {

    string aln_file_list;
    string out_prefix;
    size_t min_mapq = 0;
    size_t min_depth = 1;
    float min_vaf = 0;

    uint16_t include_all = 0;
    uint16_t include_none = 0x804;

    int opt;
    while ((opt = getopt(argc, argv, "a:o:d:v:q:f:F:")) != -1) {
      if (opt == 'a')
        aln_file_list = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'd')
        min_depth = std::stoi(optarg);
      else if (opt == 'v')
        min_vaf = std::stof(optarg);
      else if (opt == 'q')
        min_mapq = std::stoi(optarg);
      else if (opt == 'f')
        include_all = std::stoi(optarg);
      else if (opt == 'F')
        include_none = std::stoi(optarg);
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (aln_file_list.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }



    std::ifstream aln_file(aln_file_list);
    if (!aln_file)
      throw std::runtime_error("cannot open " + aln_file_list);

    GenomicStepVector<FeatureVector<string>> sample_mm;

    size_t n_samples = 0;
    unordered_map<size_t, string> sample_map;

    vector<GenomicRegion> ref_chroms;
    string aln_file_name;
    while (getline(aln_file, aln_file_name)) {
      sample_map.insert(make_pair(n_samples, aln_file_name));

      SamReader reader(aln_file_name);
      SamEntry entry;

      // parse genomic regions from sam file
      string header;
      reader.read_sam_header(header);
      get_seq_lengths(header, ref_chroms);

      AlignmentMismatch mismatches;

      while (reader.read_sam_line(entry)) {
        if (entry.mapq >= min_mapq &&
            SamFlags::is_all_set(entry.flag, include_all) &&
            !SamFlags::is_any_set(entry.flag, include_none)) {

          mismatches.add(entry);
        }
      }

      for (size_t i = 0; i < ref_chroms.size(); ++i) {
        vector<GenomicRegion> out_region;
        vector<size_t> out_cov;
        vector<string> out_mm;
        mismatches.at(ref_chroms[i], out_region, out_cov, out_mm);
        for (size_t j = 0; j < out_region.size(); ++j) {

          if (out_cov[j] >= min_depth) {
            size_t a_count, t_count, g_count, c_count;
            count_nucs(out_mm[j], a_count, t_count, g_count, c_count);

            const float a_vaf = static_cast<float>(a_count) /
                                static_cast<float>(out_cov[j]);
            if (a_vaf >= min_vaf) {
              sample_mm.add(out_region[j],
                FeatureVector<string>{"A" + std::to_string(n_samples)});
            }

            const float t_vaf = static_cast<float>(t_count) /
                                static_cast<float>(out_cov[j]);
            if (t_vaf >= min_vaf) {
              sample_mm.add(out_region[j],
                FeatureVector<string>{"T" + std::to_string(n_samples)});
            }

            const float g_vaf = static_cast<float>(g_count) /
                                static_cast<float>(out_cov[j]);
            if (g_vaf >= min_vaf) {
              sample_mm.add(out_region[j],
                FeatureVector<string>{"G" + std::to_string(n_samples)});
            }

            const float c_vaf = static_cast<float>(c_count) /
                                static_cast<float>(out_cov[j]);
            if (c_vaf >= min_vaf) {
              sample_mm.add(out_region[j],
                FeatureVector<string>{"C" + std::to_string(n_samples)});
            }

          }
        }
      }

      ++n_samples;
    }


    vector<vector<size_t>> mm_counts(n_samples, vector<size_t>(n_samples, 0));

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
    for (size_t i = 0; i < n_samples; ++i) {
      jaccard_out << sample_map[i] << "\t";
    }
    jaccard_out << endl;
    for (size_t i = 0; i < n_samples; ++i) {
      jaccard_out << sample_map[i] << "\t";
      for (size_t j = 0; j < n_samples; ++j) {
        size_t isect;
        if (i != j)
          isect = mm_counts[i][j] + mm_counts[j][i];
        else
          isect = mm_counts[i][j];
        // subtracting isect to account for double couting
        const size_t uni = mm_counts[i][i] + mm_counts[j][j] - isect;

        const float jaccard = static_cast<float>(isect) /
          static_cast<float>(uni);
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
