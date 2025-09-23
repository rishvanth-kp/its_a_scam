# bc_in_region 


## Description

## Parameters
```
bc_in_region [options]
  -a aligned SAM/BAM file [required]
  -b barcode list file [required]
  -r regions bed file [required]
  -s distance to add to either side of each region [default: 0]
  -o out file prefix [required]
  -d name split delimeter [default: ":"]
  -c barcode field in name [default: 7 (0 based)]
  -q minimum mapping quality to include [default: 0]
  -f only include if all the flags are present [default: 3]
  -F only include if none of the flags are present [default: 3340]
  -v verbose [default: false]
```

## Input and output file description
