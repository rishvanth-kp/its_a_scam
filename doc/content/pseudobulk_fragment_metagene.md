# pseudobulk_fragment_metagene 


## Description

## Parameters
```
pseudobulk_fragment_metagene [options]
	-a aligned SAM/BAM file [required]
	-b barcode list file [required]
	-r regions bed file [required]
	-o out file prefix [required]
	-n number of divisions for a region [default: 100]
	-m min. fragment length [default: 0]
	-M max. fragment length [default: 1024]
	-d name split delimeter[default: ":"]; ignored if -t is provided
	-c barcode field in name[default: 7 (0 based); ignored if -t is provided]
	-t barcode tag in SAM file [default ""]
	-q minimum mapping quality to include [default: 0]
	-f only include if all the flags are present [default: 3]
	-F only include if none of the flags are present [default: 3340]
	-v verbose [default: false]
```

## Input and output file description
