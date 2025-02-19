# bc_feature_matrix


## Description

## Parameters
```
bc_feature_matrix [options]
	-a aligned SAM/BAM file [required]
	-b barcode list file [required]
	-g gtf file [optional]
	-t tsr bed file [optional]
	-e enhancers bed file [optional]
	-o out file prefix [required]
	-m min fragment length [default: 0]
	-M max fragment length [default: inf]
	-d name split delimeter [default: ":"]
	-c barcode field in name [default: 7 (0 based)]
	-q minimum mapping quality to include [default: 0]
	-f only include if all the flags are present [default: 3]
	-F only include if none of the flags are present [default: 3340]
	-v verbose [default: false]
```

## Input and output file description
