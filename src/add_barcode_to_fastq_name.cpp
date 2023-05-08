/*
* add_barcode_to_fastq_name: read barcode fastq file and applend it to
*   fastq name
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

#include "FastqReader.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

static void
split_name (const string &in, vector<string> &tokens) {

  const char delim = ' ';
  tokens.clear();
  // the first character is always a '@' which should be skipped.
  size_t start = 1;
  size_t end = in.find(delim);
  while (end != string::npos) {
    tokens.push_back(in.substr(start, end - start));
    start = ++end;
    end = in.find(delim, start);
  }
  tokens.push_back(in.substr(start));
}

static void
add_barcode_to_name(FastqEntry &e, const FastqEntry &bc,
                    const size_t bc_col, const size_t bc_split_pos,
                    const char bc_split_delim) {

  // split the name field at spaces.
  vector<string> name_entries;
  split_name(e.name, name_entries);

  if (!bc_split_pos) {
    // no need to split the barcode
    name_entries[bc_col] = name_entries[bc_col] + ":" + bc.seq;
  }
  else {
    // split the bc and add delimeter
    name_entries[bc_col] = name_entries[bc_col] + ":" +
      bc.seq.substr(0, bc_split_pos) + bc_split_delim +
      bc.seq.substr(bc_split_pos);
  }

  // make the entry with the barcode the first
  e.name = "@" + name_entries[bc_col] + " ";
  // then write the remaining entries
  for (size_t i = 0; i < name_entries.size(); ++i) {
    if (i != bc_col)
      e.name += name_entries[i];
    if (i != name_entries.size() - 1)
      e.name += " ";
  }

}

static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-1 fastq file 1 [required]" << endl
      << "\t-2 fastq file 2 [required for paired-end data]" << endl
      << "\t-b barcode fastq [required]" << endl
      << "\t-o out file prefix [required]" << endl
      << "\t-c space delimated column to append barcode to [default: 0]"
        << endl
      << "\t-s barcode read split position [default: 0]" << endl
      << "\t-d barcode read split delimeter [default: '+']"  << endl;
  return oss.str();
}

int
main (int argc, char* argv[]) {
  try {

    string in_file_1;
    string in_file_2;
    string bc_file;
    string out_prefix;

    size_t bc_col = 0;
    size_t bc_split_pos = 0;
    char bc_split_delim = '+';

    int opt;
    while ((opt = getopt(argc, argv, "1:2:b:o:c:s:d:")) != -1) {
      if (opt == '1')
        in_file_1 = optarg;
      else if (opt == '2')
        in_file_2 = optarg;
      else if (opt == 'b')
        bc_file = optarg;
      else if (opt == 'o')
        out_prefix = optarg;
      else if (opt == 'c')
        bc_col = std::stoi(optarg);
      else if (opt == 's')
        bc_split_pos = std::stoi(optarg);
      else if (opt == 'd')
        bc_split_delim = optarg[0];
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (in_file_1.empty() || bc_file.empty() || out_prefix.empty()) {
      throw std::runtime_error(print_usage(argv[0]));
    }

    bool se = true;
    if (!in_file_2.empty())
      se = false;

    FastqReader bc_reader(bc_file);
    FastqEntry bc;

    if (se) {
      FastqReader se_reader(in_file_1);
      FastqEntry fq1;

      std::ofstream out_fq1(out_prefix + ".fastq");

      while (se_reader.read_se_entry(fq1)) {

        // read the barcode
        bc_reader.read_se_entry(bc);

        // add the barcode sequence to name
        add_barcode_to_name(fq1, bc, bc_col, bc_split_pos, bc_split_delim);

        // write the fastq file
        out_fq1 << fq1;
      }
      out_fq1.close();
    }
    else {
      FastqReader pe_reader(in_file_1, in_file_2);
      FastqEntry fq1, fq2;

      std::ofstream out_fq1(out_prefix + "_1.fastq");
      std::ofstream out_fq2(out_prefix + "_2.fastq");

      while (pe_reader.read_pe_entry(fq1, fq2)) {

        // read the barcode
        bc_reader.read_se_entry(bc);

        // add the bacode sequence to both the names
        add_barcode_to_name(fq1, bc, bc_col, bc_split_pos, bc_split_delim);
        add_barcode_to_name(fq2, bc, bc_col, bc_split_pos, bc_split_delim);

        // write the fastq files
        out_fq1 << fq1;
        out_fq2 << fq2;
      }
      out_fq1.close();
      out_fq2.close();

    }

  }

  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
