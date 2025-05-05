/*
* Copyright (C) 2025 Rishvanth Prabakar
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

#include "gcatlib/SamEntry.hpp"
#include "gcatlib/SamReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;

static void
split_string (const string &in, vector<string> &tokens,
              const char delim = '\t') {

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


static char
complement (const char in) {
  if (in == 'A')
    return 'T';
  else if (in == 'T')
    return 'A';
  else if (in == 'G')
    return 'C';
  else if (in == 'C')
    return 'G';
  else if (in == 'a')
    return 't';
  else if (in == 't')
    return 'a';
  else if (in == 'g')
    return 'c';
  else if (in == 'c')
    return 'g';
  else
    return 'N';
}

static void
inplace_rev_comp (string &s) {
  std::transform(s.begin(), s.end(), s.begin(), complement);
  std::reverse(s.begin(), s.end());
}


static void
process_delimated_file (const string &in_file, const string &out_file,
                        const char bc_delim, const uint8_t bc_col,
                        const unordered_map<string, string> &bc_match) {

  std::ifstream in(in_file);
  std::ofstream out(out_file);

  size_t in_count = 0;
  string line;
  while (getline(in, line)) {
    ++in_count;
    if (!(in_count % 1000000)) {
      cout << "\tprocessed " << in_count << " entries" << endl;
    }

    // parse the line
    string bc_matched;
    vector<string> tokens;
    split_string(line, tokens, bc_delim);

    // match the barcode
    unordered_map<string, string>::const_iterator it;
    it = bc_match.find(tokens[bc_col]);
    if (it != bc_match.end()) {
      tokens[bc_col] = it->second;
    }

    // write output
    out << tokens[0];
    for (size_t i = 1; i < tokens.size(); ++i) {
      out << bc_delim << tokens[i];
    }
    out << endl;

  }

  in.close();
  out.close();

}


static string
print_usage (const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-a aligned SAM/BAM file or TSV/CSV file [required]" << endl
      << "\t-g GEX barcode list file [required]" << endl
      << "\t-b ATAC barcode list file [required]" << endl
      << "\t-o out file name [required]" << endl
      << "\t-d split delimeter"
          << "[default: \":\"]; ignored if -t is provided" << endl
      << "\t-c barcode field column"
          << "[default: 7 (0 based); ignored if -t is provided]" << endl
      << "\t-t barcode tag in SAM file [default \"\"]" << endl
      << "\t-v verbose [default: false]" << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {

    // parse args
    string in_file;
    string out_file;

    string gex_bc_file;
    string atac_bc_file;

    char bc_delim = ':';
    uint8_t bc_col = 7;
    string bc_tag;

    bool VERBOSE = false;

    int opt;
    while ((opt = getopt(argc, argv, "a:g:b:o:d:c:t:v")) != -1) {
      if (opt == 'a')
        in_file = optarg;
      else if (opt == 'g')
        gex_bc_file = optarg;
      else if (opt == 'b')
        atac_bc_file = optarg;
      else if (opt == 'o')
        out_file = optarg;
      else if (opt == 'd') {
        if (optarg[0] == 't')
          bc_delim = '\t';
        else
          bc_delim = optarg[0];
      }
      else if (opt == 'c')
        bc_col = std::stoi(optarg);
      else if (opt == 't')
        bc_tag = optarg;
      else if (opt == 'v')
        VERBOSE = true;
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (in_file.empty() || out_file.empty() ||
        gex_bc_file.empty() || atac_bc_file.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    // process barcodes
    if (VERBOSE)
      cerr << "[PROCESSING BARCODES]" << endl;

    std::ifstream gex_bc_in(gex_bc_file);
    std::ifstream atac_bc_in(atac_bc_file);

    unordered_map<string, string> bc_match;
    string gex_line, atac_line;
    while (getline(gex_bc_in, gex_line) && getline(atac_bc_in, atac_line)) {
      inplace_rev_comp(atac_line);
      bc_match[atac_line] = gex_line;
    }

    gex_bc_in.close();
    atac_bc_in.close();

    // process input
    if (VERBOSE)
      cerr << "[PROCESSING INPUT]" << endl;

    // chcek if input is a sam or a bam file
    if (in_file.length() < 4) {
      throw std::runtime_error("Invalid input file name.");
    }

    string file_type = in_file.substr(in_file.length() - 4, 4);
    if (file_type == ".sam" || file_type == ".bam") {
      if (VERBOSE) {
        cerr << "\tProcessing a SAM/BAM file" << endl;
      }


    }
    else {
      if (VERBOSE) {
        cerr << "\tProcessing a " << bc_delim << " delimated file." << endl;
        process_delimated_file (in_file, out_file, bc_delim, bc_col, bc_match);
      }


    }

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
