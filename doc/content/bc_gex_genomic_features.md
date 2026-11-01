# bc_gex_genomic_features

## Description

`bc_gex_genomic_features` is used to determine the fraction of aligned bases
to different genomic regions at the single cell level.
In a scRNA-seq experiment, we would expect a lot of reads aligning to the
exonic regions and, depending on the protocol, to the intronic regions.
This program determines the number of aligned bases to these regions.

The aligned bases of each alignment (the `M`, `=`, and `X` operations of
the CIGAR string; spliced regions and soft clipped bases are not counted)
are assigned to the following genomic regions, which are obtained from
the GTF file:

- `pc_cds`, `pc_utr`, `pc_both`, `pc_none`, and `pc_intron`: Regions of
  protein coding genes that are in a CDS, in a UTR, in both (from
  different transcripts), in an exon that is neither a CDS nor a UTR, and
  in an intron.
- Regions of other genes are named by the gene type (for example,
  `lncRNA`).
- `multi_gene`: Regions that overlap more than one gene.
- `chrM_<gene_type>`: Genes on `chrM`.
- `intergenic`: Regions that are not in a gene.

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
	-F only include if none of the flags are present [default: 2308]
	-v verbose [default: false]
```

`-a`: The input aligned SAM or BAM file. Each alignment is processed
separately. The cell barcode should be present in the SAM tag specified
by `-c`. [Required].

`-b`: A list of barcodes to include in the output file. Only the
alignments with these barcodes are counted. [Required].

`-g`: GTF file that is used to parse the gene regions. [Required].

`-o`: Output file prefix. [Required].

`-c`: The SAM tag containing the cell barcode. [Default: CB].

`-u`: The SAM tag containing the UMI. If this tag is provided, only the
first alignment of each CB+UMI combination is used, and all other
alignments with an identical CB+UMI pair are discarded. This tag is
useful to not overcount the alignments by retaining just one alignment
per UMI. By default, all alignments are counted disregarding the UMI
information. [Default: ""].

`-q`: Include only alignments with a mapping quality greater than or
equal to the specified value. [Default: 0].

`-f`: Include only alignments that have all of the specified flags
present. For example, to include paired reads (0x1, decimal 1) and reads
that are mapped in proper pair (0x2, decimal 2), specify the sum of
these flags (1 + 2 = 3). [Default: 0].

`-F`: Include only alignments that have none of the specified flags
present. For example, to exclude PCR duplicates (0x400, decimal 1024) and
supplementary alignments (0x800, decimal 2048), specify the sum of these
flags (1024 + 2048 = 3072). The default (4 + 256 + 2048 = 2308) excludes
alignments that are unmapped (0x4), are secondary alignments (0x100), or
are supplementary alignments (0x800). [Default: 2308].

`-v`: Verbosity flag. Set the flag to print progress and statistics to
standard error. [Default: false].


## Input and output file description

### Description of input files
1. An aligned SAM/BAM file that has a cell barcode in a SAM tag.
1. An uncompressed TSV file with the cell barcode in the first column.
Only the alignments containing barcodes in the file are counted.
1. An uncompressed GTF file with the genomic regions of interest. The GTF
file downloaded from [Gencode](https://www.gencodegenes.org/) can be used
directly. The `gene`, `exon`, `CDS`, and `UTR` features, and the
`gene_type` attribute of the GTF file are used. The chromosome names must
match the chromosome names in the SAM/BAM file; bases that align to a
chromosome that is not in the GTF file are not counted.

### Description of output files
`<output_prefix>_gex_feature_counts.txt`: This file contains one row per
cell. The first column (`barcode`) contains the cell barcode, the second
column (`bases`) contains the total number of aligned bases counted over
all regions (that can be used for normalizing the counts for each cell),
and the subsequent columns are the genomic regions. These columns contain
the number of aligned bases in each genomic region. This file has a
header line with the names of the genomic regions. The genomic region
columns are not in a fixed order, and the rows are not in the same order
as the barcode list.
