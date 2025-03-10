# barcode_count  


## Description
`barcode_count` is used to count the number of fragments that belong to
each cell barcode for ATAC-seq analysis. 

## Parameters
```
barcode_count [options]
	-a aligned SAM/BAM file [required]
	-o out file prefix [required]
	-m minimun barcode count to output [default: 1000]
	-w optional cell barcode whitelist [default: ""]
	-r reverse complement the cell barcodes [default: true]
	-d name split delimeter[default: ":"; ignored if -t is provided]
	-c barcode field in name[default: 7 (0 based); ignored if -t is provided]
	-t barcode tag in SAM file [default: ""]
	-q minimum mapping quality to include [default: 0]
	-f only include if all the flags are present [default: 3]
	-F only include if none of the flags are present [default: 3340]
	-v verbose [default: false]
```

`-a`: The input aligned paired end SAM or BAM file. The cell barcode
should be present either in a SAM tag or encoded in the query name
(QNAME) field, see options `\t`, `-d`, and `-c`. [Required].

`-o`: Output file prefix to name out files. [Required].

`-m`: Minimum number of fragments (after meeting the minimum mapping
quality and include/exclude flag criteria) that are required for a cell
barcode to be included in the output. [Default: 1000].

`-w`: Path to a whitelist file (usually provided by 10x) that contains
a list of known cell barcodes. A cell barcode is included in the output
only if the barcode is present in the whitelist file. If a whitelist
file is not provided, then all barcodes that meet the above minimum
fragment counts are included. Of note, as of now, no barcode error 
correction is done. [Default: none].

`-r`:

`-d`:

`-c`:

`-t`:

`-q`: Include only fragments for which both reads have a mapping
quality greater than or equal to the specified value. [Default: 0]

`-f`: Include only fragments for which both reads have all of the
specified flags present. For example, to include paired reads (0x1 or
0b1) and reads that are mapped in proper pair (0x2 or 0b2), specify the
sum of of these flags (1 + 2 = 3). [Default: 3]. 

`-F`: Include only fragments for which both reads have none of the
specified frags present. For example, to exclude PCR duplicates (0x400
or 0b1024) and supplementary alignments (0x800 or 0b2048), specify the
sum of these flags (1024 + 2048 = 3072). [default: 3340].

`-v`: Verbosity flag. Set the flag to print progress and statistics to
standard error. [Default: false].



## Input and output file description
