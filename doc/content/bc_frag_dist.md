# bc_frag_dist


## Description

## Parameters
```
bc_frag_dist [options]
	-a aligned SAM/BAM file [required]
	-b barcode list file [required]
	-o out file prefix [required]
	-m max fragment length to track [default: 600]
	-d name split delimeter [default: ":"]
	-c barcode field in name [default: 7 (0 based)]
	-q minimum mapping quality to include [default: 0]
	-f only include if all the flags are present [default: 3]
	-F only include if none of the flags are present [default: 3340]
	-v verbose [default: false]
```

## Input and output file description
