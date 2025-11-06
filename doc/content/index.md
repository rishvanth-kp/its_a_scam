# It's a scam

**I**t's a scam **T**o **S**ell **A**nother **S**ingle **C**ell
**A**nalysis **M**ethod (It's a scam) is a collection of programs to
pre-process and quality control single cell data.

The era of designing bioinformatic pipelines using open source
tools is giving way to monolithic, and often proprietary, tools designed
for end-to-end analysis. While these tools can be convenient for going
from sequencing data to an almost complete data analysis "in a few
clicks", they are severely limiting to perform any custom analysis. We
developed its_a_scam for processing and quality control of single cell
genomic data, freeing users from the reliance on proprietary tools for
single cell data analysis.

its_a_scam is built of top on `gcatlib` and `HTSlib`, and provides tools
for wide range of single cell data processing and quality control.  It's
a scam can be used in conjugation with commonly used open source tools
to build entire single cell analysis pipelines. 

## Getting started
After [installation](https://its-a-scam.readthedocs.io/en/latest/installation/),
its_a_scam can be used for scRNA-seq and scATAC-seq count matrix
generation and quality control. Refer to following documents to get
started:

1. [10x scATAC-seq count matrix generation](
https://its-a-scam.readthedocs.io/en/latest/10x_atac_seq/).

1. scRNA-seq count matrix generation: There are several well
designed tools already available for this task. We describe the procedure 
to get started with using STARsolo for [10x data](
https://its-a-scam.readthedocs.io/en/latest/10x_gex/) and for
[Parse data](https://its-a-scam.readthedocs.io/en/latest/10x_gex/) for
Parse data.

1. [scATAC-seq quality control](
https://its-a-scam.readthedocs.io/en/latest/10x_atac_seq_qc/).

1. [scRNA-seq quality control](
https://its-a-scam.readthedocs.io/en/latest/sc_gex_qc/).



## Contact
Rishvanth Prabakar: [kaliapp@cshl.edu](mailto:kaliapp@cshl.edu). Please
let me know any issues, feedback, or ways to improve this tool. If you
would like any other metrics or programs added in, please let me know
those as well. 

## License
Copyright (C) 2025 Rishvanth Prabakar

Authors: Rish Prabakar

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
