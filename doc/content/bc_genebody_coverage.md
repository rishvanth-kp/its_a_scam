# bc_genebody_coverage

## Description
`bc_genebody_coverage` is used for creating per cell gene body coverage
plots. It takes as input a list of regions (such as transcript regions
or exonic regions of genes). It divides each region into identical
number of bins, and computes the number of reads that are aligned in
each bin. For each cell, the read counts are summed over all regions, in
effect creating a metagene gene body coverage. 
 

## Parameters
```
bc_genebody_coverage [options]
  -a aligned SAM/BAM file [required]
  -b barcode list file [required]
  -r regions bed file [required]
  -o out file prefix [required]
  -n number of divisions to for a region [default: 100]
  -c cell barcode tag in SAM/BAM file [default: CB]
  -u UMI tag in SAM/BAM file [default: ""]
  -q minimum mapping quality to include [default: 0]
  -f only include if all the flags are present [default: 0]
  -F only include if none of the flags are present [default: 2308]
  -v verbose [default: false]
```
`-a`: The input aligned paired end SAM or BAM file. The cell barcode
should be present either in a SAM tag or encoded in the query name
(QNAME) field, see options `\t`, `-d`, and `-c`. [Required].

`-b`: A list of barcodes to retain in the output file. [Required].

`-r` A bed-6 file with the regions of interest. Specifically the first 3
standard columns, the 4th column with the region name, and the 6th column
with the stand orientation info are used. See below on how to generate
this file [Required]. 

`-o`: Output file prefix. [Default = ""].

`-n`: Number of divisions to split each region into. [Default: 100].

`-c`: The sam tag containing the cell barcode. [Default: CB]

`-u`: The sam tag containing the UMI. If this tag is provided only the
first unique CB+UMI combination is used, and all other identical CB+UMI
pairs are discarded. This tag is useful to not overcount the alignments
by retraining just one alignment per UMI. By default, all alignments are
counted disregarding the UMI information. [Default: ""]

`-q`: Include only fragments for which both reads have a mapping
quality greater than or equal to the specified value. [Default: 0]

`-f`: Include only fragments for which both reads have all of the
specified flags present. For example, to include paired reads (0x1 or
0b1) and reads that are mapped in proper pair (0x2 or 0b2), specify the
sum of of these flags (1 + 2 = 3). [Default: 0]. 

`-F`: Include only fragments for which both reads have none of the
specified frags present. For example, to exclude PCR duplicates (0x400
or 0b1024) and supplementary alignments (0x800 or 0b2048), specify the
sum of these flags (1024 + 2048 = 3072). [default: 2308].

`-v`: Verbosity flag. Set the flag to print progress and statistics to
standard error. [Default: false].

## Input and output file description

### Description of input files
1. An aligned SAM/BAM file that has a cell barcode in a sam tag.  
2. A TSV file with the cell barcode in the first column. Only the
alignments containing barcodes in the file are retained in the output
SAM file.  
1. A bed-6 file with the regions of interest. This file is expected to
contain all the 6 bed columns. Additionally, all the regions that have a
same name (in the 4th columns) are expected to be contiguous, and they
should be sorted by genome positions within each regions. This can be
accomplished by sorting a bed file with:

```
cat in_unsorted.bed | sort -k1,1 -k4,4 -k2,2n > out_sorted.bed
```

### Description of output files
1. `<output_prefix>_genebody_coverage.txt`: A file with one row per cell. The
first column contains the cell ID, followed by number columns specified
in the `n` parameter (default: 100). They contain the number of reads
aligned to a specific bin summed over all regions. This file contains a
header line.  
