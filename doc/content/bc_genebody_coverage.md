# bc_genebody_coverage

## Description
`bc_genebody_coverage` is used for creating per cell gene body coverage
plots. It takes as input a list of regions (such as transcript regions
or exonic regions of genes). The bed lines that have the same name are
combined into one region (for example, the exons of a transcript). It
divides each region into the same number of bins, and computes the number
of reads that are aligned in each bin. Each alignment is counted once, in
the bin that contains its leftmost aligned position (the POS field of
the SAM entry). The bins are ordered from the 5' to the 3' end of each
region, according to its strand. For each cell, the read counts are
summed over all regions, in effect creating a metagene gene body
coverage.

## Parameters
```
bc_genebody_coverage [options]
  -a aligned SAM/BAM file [required]
  -b barcode list file [required]
  -r regions bed file [required]
  -o out file prefix [required]
  -n number of divisions for a region [default: 100]
  -c cell barcode tag in SAM/BAM file [default: CB]
  -u UMI tag in SAM/BAM file [default: ""]
  -q minimum mapping quality to include [default: 0]
  -f only include if all the flags are present [default: 0]
  -F only include if none of the flags are present [default: 2308]
  -v verbose [default: false]
```

`-a`: The input aligned SAM or BAM file. Each alignment is processed
separately. The cell barcode should be present in the SAM tag specified
by `-c`. [Required].

`-b`: A list of barcodes to include in the output file. Only the
alignments with these barcodes are counted. [Required].

`-r`: A bed-6 file with the regions of interest. Specifically, the first 3
standard columns, the 4th column with the region name, and the 6th column
with the strand orientation info are used. See below on how to generate
this file. [Required].

`-o`: Output file prefix. [Required].

`-n`: Number of divisions (bins) to split each region into. [Default: 100].

`-c`: The SAM tag containing the cell barcode. [Default: CB].

`-u`: The SAM tag containing the UMI. If this tag is provided, only the
first alignment of each CB+UMI combination is used, and all other
alignments with an identical CB+UMI pair are discarded. This tag is
useful to not overcount the alignments by retaining just one alignment
per UMI. By default, all alignments are counted disregarding the UMI
information. [Default: ""].

`-q`: Include only alignments with a mapping quality greater than or
equal to the specified value. [Default: 0].

`-f`: Include only alignments that have all of the specified flags
present. For example, to include paired reads (0x1, decimal 1) and reads
that are mapped in proper pair (0x2, decimal 2), specify the sum of
these flags (1 + 2 = 3). [Default: 0].

`-F`: Include only alignments that have none of the specified flags
present. For example, to exclude PCR duplicates (0x400, decimal 1024) and
supplementary alignments (0x800, decimal 2048), specify the sum of these
flags (1024 + 2048 = 3072). The default (4 + 256 + 2048 = 2308) excludes
alignments that are unmapped (0x4), are secondary alignments (0x100), or
are supplementary alignments (0x800). [Default: 2308].

`-v`: Verbosity flag. Set the flag to print progress and statistics to
standard error. [Default: false].

## Input and output file description

### Description of input files
1. An aligned SAM/BAM file that has a cell barcode in a SAM tag.
2. An uncompressed TSV file with the cell barcode in the first column.
Only the alignments containing barcodes in the file are counted.
3. A bed-6 file with the regions of interest. This file is expected to
contain all the 6 bed columns. Additionally, all the regions that have
the same name (in the 4th column) are expected to be contiguous, and they
should be sorted by genome positions within each region; this is not
checked. This can be accomplished by sorting a bed file with:

```
cat in_unsorted.bed | sort -k1,1 -k4,4 -k2,2n > out_sorted.bed
```

### Description of output files
1. `<output_prefix>_genebody_coverage.txt`: A file with one row per cell,
in the same order as the barcode list. The first column contains the
cell barcode, followed by the number of columns specified by `-n`
(default: 100). They contain the number of alignments that start in each
bin, summed over all regions, ordered from the 5' to the 3' end. This
file contains a header line with `barcode`, followed by the bin numbers
(starting from 0).
