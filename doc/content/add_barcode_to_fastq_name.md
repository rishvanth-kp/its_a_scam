# add_barcode_to_fastq_name 


## Description

## Parameters
```
add_barcode_to_fastq_name [options]
        -1 fastq file 1 [required]
        -2 fastq file 2 [required for paired-end data]
        -b barcode fastq [required]
        -o out file prefix [required]
        -c space delimated column to append barcode to [default: 0]
        -f first base of barcode (0 based; ignored if -s > 0) [default: 0]
        -s barcode read split position [default: 0]
        -d barcode read split delimeter [default: '+']
```

## Input and output file description
