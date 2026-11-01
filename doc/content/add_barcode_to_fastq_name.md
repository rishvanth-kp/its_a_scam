# add_barcode_to_fastq_name

## Description
`add_barcode_to_fastq_name` is used to append the cell barcodes from a
separate index fastq file to the name line in single/paired-end fastq
files.

Most common aligners are designed for mapping single/paired-end reads to the
reference genome, and at the moment, cannot process the third fastq file
containing the cell barcode. This limitation can be easily overcome by
appending the cell barcode information from the third fastq file to the
name field of the input fastq files. The first space-delimited column of
the name field is retained as the query name (QNAME) in the aligned SAM
output, and so the cell barcode information can be parsed for each
alignment. To ensure that the cell barcode is retained, the column the
cell barcode is added to (see `-c`) is moved to the front of the name
field.

Unfortunately, `add_barcode_to_fastq_name` does not support gzipped
fastq files at the moment. The fastq files need to be decompressed
before running this program and then re-compressed after. `pigz`
is a useful program for fast compression.

## Parameters
```
add_barcode_to_fastq_name [options]
        -1 fastq file 1 [required]
        -2 fastq file 2 [required for paired-end data]
        -b barcode fastq [required]
        -o out file prefix [required]
        -c space delimited column to append barcode to [default: 0]
        -f first base of barcode (0 based; ignored if -s > 0) [default: 0]
        -s barcode read split position [default: 0]
        -d barcode read split delimiter [default: '+']
```

`-1`: Name of the first fastq file of a single/paired-end sample. [Required].

`-2`: Name of the second fastq file of a paired-end sample. If this
argument is not set, the input is treated as a single-end sample. [Required
for paired-end data].

`-b`: Name of the fastq file that contains the cell barcode. [Required].

`-o`: Output file prefix to name out files. [Required].

`-c`: The fastq name line can contain multiple columns that are
delimited with a space. The cell barcode is added to the column
specified by this parameter (0 based). The default parameter works for
fastq files that are obtained directly from an Illumina sequencer. This
parameter is useful for files downloaded from SRA, which have the
standard Illumina name in the second column (use `-c 1`). The column
must exist in the fastq name line; this is not checked. [Default: 0].

`-f`: The starting base of the cell barcode in the index fastq entry.
This parameter is useful when the cell barcode does not start in the
first base. This parameter is usually 0, unless some additional
technical bases were sequenced in the index reads. Must be less than the
length of the index read. [Default: 0].

`-s`: Split the index fastq entry at this specified position, and add
two entries to the fastq name. These entries are delimited by the `-d`
parameter. This is useful when the index read contains both the UMI and
cell barcode, and to add both these to the fastq name. Must be less than
the length of the index read. [Default: 0].

`-d`: The delimiter used when splitting the barcode. [Default: '+'].

## Input and output file description
### Description of input files
1. Single-end (specified in `-1`) or paired-end (specified in `-1` and
`-2`) fastq files that are uncompressed.
2. Index fastq file that is uncompressed. The index fastq file must
contain one entry for each read (or read pair), in the same order as the
input fastq file(s). The read names are not compared between the files,
so it is up to the user to provide them in the same order.

### Description of output files
1. Uncompressed fastq file(s) that are named `<output_prefix>.fastq`
for single-end files (or `<output_prefix>_1.fastq` and
`<output_prefix>_2.fastq` for paired-end files) that have the index
fastq sequence added to the name of the fastq files.

For example, for the following fastq name line and index read:
```
@A00:1:FC:1:1:100:200 1:N:0:ACGT
AAAACCCCGGGGTTTT
```
the name line in the output fastq file is:
```
@A00:1:FC:1:1:100:200:AAAACCCCGGGGTTTT  1:N:0:ACGT     (default)
@A00:1:FC:1:1:100:200:CCCCGGGGTTTT  1:N:0:ACGT         (-f 4)
@A00:1:FC:1:1:100:200:AAAACCCC+GGGGTTTT  1:N:0:ACGT    (-s 8)
```
