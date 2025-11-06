# bc_gex_genomic_features 

## Description

`bc_gex_genomic_features` is used to determine the fraction of aligned bases 
to different genomic regions at the single cell level.
In an scRNA-seq experiment we would expect a lot reads aligning to the
exonic regions and, depending on the protocol, to the intronic regions.
This program determines the number of aligned bases to these regions. 

## Parameters
```
bc_gex_genomic_features [options]
	-a aligned SAM/BAM file [required]
	-b barcode list file [required]
	-g GTF file [required]
	-o out file prefix [required]
	-c cell barcode tag in SAM/BAM file [default: CB]
	-u UMI tag in SAM/BAM file [default: ""]
	-q minimum mapping quality to include [default: 0]
	-f only include if all the flags are present [default: 0]
	-F only include if none of the flags are present [defauly: 2308]
	-v verbose [default: false]
```
`-a`: The input aligned paired end SAM or BAM file. The cell barcode
should be present either in a SAM tag or encoded in the query name
(QNAME) field, see options `\t`, `-d`, and `-c`. [Required].

`-b`: A list of barcodes to retain in the output file. [Required].

`-g`: GTF file that is used to parse the gene regions. [Required].

`-o`: Output file prefix. [Default = ""].

`-c`: The sam tag containing the cell barcode. [Default: CB]

`-u`: The sam tag containing the UMI. If this tag is provided only the
first unique CB+UMI combination is used, and all other identical CB+UMI
pairs are discarded. This tag is useful to not overcount the alignments
by retraining just one alignment per UMI. By default, all alignments are
counted disregarding the UMI information. [Default: ""]

`-q`: Include only fragments for which both reads have a mapping
quality greater than or equal to the specified value. [Default: 0]

`-f`: Include only fragments for which both reads have all of the
specified flags present. For example, to include paired reads (0x1 or
0b1) and reads that are mapped in proper pair (0x2 or 0b2), specify the
sum of of these flags (1 + 2 = 3). [Default: 0]. 

`-F`: Include only fragments for which both reads have none of the
specified frags present. For example, to exclude PCR duplicates (0x400
or 0b1024) and supplementary alignments (0x800 or 0b2048), specify the
sum of these flags (1024 + 2048 = 3072). [default: 2308].

`-v`: Verbosity flag. Set the flag to print progress and statistics to
standard error. [Default: false].


## Input and output file description

### Description of input files
1. An aligned SAM/BAM file that has a cell barcode in a sam tag.  
1. A TSV file with the cell barcode in the first column. Only the
alignments containing barcodes in the file are retained in the output
SAM file.  
1. A GTF file with the genomic regions of interest. The GTF file
downloaded from [Gencode](https://www.gencodegenes.org/) can be used
directly.  

### Description of output files.
`<output_prefix>_gex_feature_counts.txt`: This file contains one row per
cell. The first column contains the cell barcode, the second column
contains the total base count over all regions (that can be used for
normalizing the counts for each cell), and the subsequent columns are
the genomic regions. The contain the number of reads aligned the a
regions. This file has the header line describing the genomic regions.  
