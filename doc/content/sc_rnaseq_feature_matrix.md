# sc_rnaseq_feature_matrix 


## Description

## Parameters
```
sc_rnaseq_feature_matrix [options]
	-a aligned SAM/BAM file [required]
	-b barcode list file [required]
	-g GTF file [required]
	-o out file prefix [required]
	-c cell barcode tag in SAM/BAM file [default: CB]
	-u UMI tag in SAM/BAM file [default: ""]
	-q minimum mapping quality to include [default: 0]
	-f only include if all the flags are present [default: 0]
	-F only include if none of the flags are present [defauly: 2052]
	-v verbose [default: false]
```

## Input and output file description
