# barcode_match

## Description
Match the scATAC-seq cell barcode to GEX cell barcode for 10x multiome
datasets. 

## Parameters
```
barcode_match [options]
  -a aligned SAM/BAM file or TSV/CSV file [required]
  -g GEX barcode list file [required]
  -b ATAC barcode list file [required]
  -o out file name [required]
  -d split delimeter; Ignored if -t is provided. [default: ":"]
  -c barcode field column; Ignored if -t is provided. [default: 7 (0 based)]
  -t barcode tag in SAM file [default ""]
  -s suffix to add to output barcode [default: ""]
  -v verbose [default: false];
```

`-a`

`-g`

`-b`

`-o`

`-d`

`-c`

`-t`

`-s`

`-v`

## Input and output file description
