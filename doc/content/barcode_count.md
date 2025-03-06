# barcode_count  


## Description
`barcode_count` is used to count the number of fragmets that belong to
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

## Input and output file description
