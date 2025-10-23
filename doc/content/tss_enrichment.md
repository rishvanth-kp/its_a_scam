# tss_enrichment 


## Description

## Parameters
```
tss_enrichment [options]
	-a aligned SAM/BAM file [required]
	-r TSS region bed file [required]
	-b barcode list (if empty, treated as bulk sample)
	-o outfile prefix [required]
	-n min fragment length [default: 10]
	-x max fragment length [default: 1000]
	-d name split delimeter [default: ":"; ignored if -t is provided]
	-c barcode field in name [default: 7 (0 based); ignored if -t is provided]
	-t barcode tag in SAM/BAM file [default ""]
	-q minimum mapping quality to include [default: 0]
	-f only include if all the flags are present [default: 3]
	-F only include if none of the flags are present [default: 3340]
	-v verbose [default: false]
```

## Input and output file description
