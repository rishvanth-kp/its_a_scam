# bc_filter_alignments


## Description

## Parameters
```
bc_filter_alignments [options]
	-a aligned SAM/BAM file [required]
	-b barcode list file [required]
	-o out file name [required]
	-d name split delimeter[default: ":"]; ignored if -t is provided
	-c barcode field in name[default: 7 (0 based); ignored if -t is provided]
	-t barcode tag in SAM file [default ""]
	-q minimum mapping quality to include [default: 0]
	-f only include if all the flags are present [default: 3]
	-F only include if none of the flags are present [default: 3340]
	-v verbose [default: false]
```

## Input and output file description
