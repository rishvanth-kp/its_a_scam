# bc_coverage 


## Description

## Parameters
```
bc_coverage [options]
	-a aligned SAM/BAM file [required]
	-b barcode list [required]
	-o outfile prefix [required]
	-d name split delimeter [default: ":"; ignored if -t is provided]
	-c barcode field in name [default: 7 (0 based); ignored if -t is provided]
	-q minimun mapping quality to include [default: 0]
	-f only include if all the flags are present [default: 3]
	-F only include if none of the flags are present [default: 3340]
	-v verbose [default: false]
```

## Input and output file description
