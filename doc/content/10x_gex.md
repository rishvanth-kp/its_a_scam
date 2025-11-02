# Analysis of 10x gene-expression datasets

## Background
Single cell transcriptomic analysis is widely used to for unbiased cell
type and cell state identification in tissues. A vast majority of single
cell RNA-seq experiments are done using the 10x platform.  10x offers
several variants of scRNA-seq protocol -- most common is the the 3'
protocol, and other variants include the 5' protocol and the multiome
(single nucleus RNA-seq + scATAC seq) protocol. The analysis of data
generated from these protocols are similar with only minor differences
between them.  Since the 10x 3' protocol is the most commonly used, we
will use that to describe the analysis, and then point out the
differences when using other protocols. 

In this tutorial we will start with the fastq files and use STARsolo to
align these read to the reference genome and generate a cell-gene count
matrix.  This tutorial is a starting point to perform scRNA-seq analysis
from scratch.  STARsolo is a very well documented, reading the [STAR
documentation](
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf),
[STARsolo documentation](
https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md), [the
paper](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1) are
highly recommended, especially for using the software beyond what is
described here. 

## Analysis of 10x 3' datasets
### Reference genome and annotation preparation:
The reference genome and the gene annotations can be downloaded from
[Gencode](https://www.gencodegenes.org/). 

To align reads and generate the count matrices using STARsolo, a genome
and the annotations should be pre-processed to generate an index. This
steps needs to performed only once before the first use (for a specif
version of the the genome and annotation).  

```
STAR --runMode genomeGenerate --runThreadN {threads} 
  --genomeDir {output_dir} --genomeFastaFiles {reference.fa} 
  --sjdbGTFfile {annotation.gtf}
```

Where `{reference.fa}` and `{annotation.gtf}` are the reference genome
and the gene annotations, respectively. `{output_dir}` is the output
directory to store the index. Specify the number of `{threads}` as
optimal for your computational resource.

### Aligning reads and count matrix generation:
10x 3' protocol uses an Illumina sequencer in paired-end mode, however,
these are not the typical paired end reads (that align to close-by
regions on the reference genome). For the 10x 3' protocol one read
originates from the 3' end of a transcript (typically ~100bp), and the
other read contains the UMI and cell barcode. The length of the barcode
read varies by the version of the protocol. The 3'v3 protocol uses a 12
bp UMI and a 16bp cell barcode, whereas the 3'v2 protocol uses a 10bp
UMI and a 16bp cell barcode. The length of the UMI and cell barcode for
your version of the protocol can be obtained from the protocol
specification document. So for the v3 protocol this read should be 28 bp
and for v2 protocol should be 26bp. 

STARsolo does all the steps involved in a typical scRNA-seq workflow: 

1. Aligns the fastq files to the reference genome taking into account
that a read can span a splice junction. 
1. Cell barcode and UMI error correction to remove potential
sequencing errors.
1. Assigning a read to a gene based on the defined gene model. 
1. Cell filtering to remove empty drops that does not contain a cell and
doublet drops that contain multiple cells.

Each of these steps are controlled with user defined parameters to STAR
described [here](
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).
Parameters used specifically for single-cell and 10x data are described
[here]( https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).  
For reads generated using a standard 10x 3'v3
protocol, STAR can be run as follows:

```
STAR --genomeDir {genome_dir} --runThreadN {threads} 
  --readFilesIn {read_1.fastq.gz} {read_I.fastq.gz} --readFilesCommand zcat
  --outFileNamePrefix {output_prefix}
  --soloCBwhitelist {inclusion_list.txt}
  --soloType CB_UMI_Simple --soloUMIlen {umi_length}
  --soloBarcodeReadLength {1/0} --clipAdapterType CellRanger4 
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
  --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR
  --soloCellFilter EmptyDrops_CR
  --soloFeatures GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS
  --soloMultiMappers EM
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
  --outSAMtype BAM SortedByCoordinate --outFilterScoreMin 30 
```
Where

`--genomeDir {genome_dir}`: The directory containing the indexed
reference genome (the `output_dir` in the previous step). 

`--runThreadN {threads}`: Number of threads to use. 

`--readFilesIn {read_1.fastq.gz} {read_I.fastq.gz}`: `read_1.fastq.gz`
is the read file that contains the 3' end of the transcript and 
`read_I.fastq.gz` is the read file that contains the cell barcode and
UMI. If you are unsure which file contains what, you can usually infer
them from the first few lines of the file:
```
zcat {read_1.fastq.gz} | head -n 20
zcat {read_I.fastq.gz} | head -n 20 
```
The sequence/quality line in the transcript read will be ~100bp,
whereas it will be much smaller at 26-28bp in the barcode/UMI read. 

`--outFileNamePrefix {output_prefix}`: The prefix used for the output
file names.

`--soloCBwhitelist {inclusion_list.txt}`: A file containing a list of
known cell barcodes. This file can be found [here](
https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist)

`--soloUMIlen {umi_length}`: The length of the UMI. This varies based on
the version of the 10x3' protocol used. The v2 protocol uses a 10bp UMI
and the v3 protocol uses a 12bp UMI. 

`--soloBarcodeReadLength {1/0}`: The length of the reads in the 
`read_I.fastq.gz` file should ideally be the sum of the barcode length
(16bp) and the UMI length (12bp for v3). Setting this parameter to `1`
explicitly checks that the index reads is of this specified length.
However, sometimes these index reads can be longer than needed, in
these cases set the parameter to `0` so that the read length is not
checked (and any bases over the needed are ignored).
 

**Choice of gene model:**
Naively, a gene consists of exons interspersed with introns. Ideally,
a mature transcript that is sequenced should align to the exons. But in
reality the captured transcripts can be immature (or due to technical
reasons) and can align to both introns and exons. The choice of how a
alignment to the exon/intron is quantified is defined by the gene model
used (specified by the `soloFeatures` parameter). At the most basic
level, just the reads aligning to exons can be counted and reads
aligning to introns can be ignored (as done by the default `Gene`
option).  However, this is too restrictive since we are expected to
capture immature transcripts that could have reads aligning to introns.
Using the `GeneFull` counts all reads toward a gene if it aligned to
ether the intron of exons of a genes. Options, `GeneFull_ExonOverIntron`
and `GeneFull_Ex50pAS` prioritize exonic reads over intronic reads. The
details of the prioritization are described [here](
https://github.com/alexdobin/STAR/issues/1460). Multiple options can be
specified to `soloFeatures` to get them all in a single run of STAR. 

**Dealing with reads aligned to overlapping genes (multi-mapping reads):**
A lot of genes overlap on the genome. So if a read
aligns to such an overlapping region, which gene should it be
counted toward? The choice is controlled by the `soloMultiMappers`. 
The easiest option is to simply ignore such reads, just
count the reads that align to non-overlapping gene regions (as done by
the default `Unique` option). Other options include uniformly assigning
the read to either gene (`Uniform` option), or distributing the reads
using a log-likelihood expectation-maximization algorithm (`EM` option). 
Multiple options can be specified to `soloMultiMappers` to get them all in
a single run of STAR.

The choice of gene model and the algorithm to deal with multi-mapping
reads depends on the downstream application. 

**Cell filtering and count matrix generation:**
Cell-gene count matrix is a large matrix where the rows correspond to a
gene and the columns correspond to a cell. Each element of the matrix
represents the count assigned to a specific gene for a specific cell. 

STARsolo stores the count matrices within a folder labeled
`{output_prefix}_Solo.out`. This folder contains a sub-folder for each
`soloFeatures` specified, each of which contain a sub-folder called `raw`
and `filtered`. 

The `raw` folders contain the counts for all cell barcodes that are
specified in the inclusion list, i.e. the counts for all cells prior to
cell filtering (to remove empty drops, doublets, etc.).  This folder
contains the following files:

1. `barcodes.tsv`: contains all the cell barcodes in the 
inclusion list i.e the column names of the count matrix.  
1. `features.tsv`: contains all the gene names i.e the row names of the
count matrix.
1. `matrix.mtx`: a matrix file containing the gene counts for each
cells and gene under the `Unique` appraoch to deal with multi-mapping
reads.
1. other `*.mtx` files: one matrix file for every `soloMultiMappers`
option that containing the counts under the respective multi-mapping
appraoch. 

The `filtered` folder contains the count matrix for just the cells that
pass the cell filtration. This folder contains `barcodes.tsv`,
`features.tsv`, and `matrix.mtx` which are identical to the above
described files. The only difference is that the `baroces.tsv` contains
much fewer cells, and thus the count matrix will have fewer columns. 

To use the count matrix using the `GeneFull_Ex50pAS` 
and the `Unique` option, the matrix from the `filtered` folder can be
directed loaded into your tool of choice. For example, to load the count
matrix into `Seurat`:

```
library(Seurat)

# Read the filtered count matrix
sample.mtx <- ReadMtx(mtx = "GeneFull_Ex50pAS/filtered/matrix.mtx",
  cells = "GeneFull_Ex50pAS/filtered/barcodes.tsv",
  features = "GeneFull_Ex50pAS/filtered/features.tsv")

# Convert the count matrix into a SeuratObject
sample <- CreateSeuratObject(counts = sample.mtx)
``` 

To load the the count matrices under the other multi-mapping model, we
need to load the count matrix for all cells and then subset it to
include only the filtered cells (since only the raw count matrices are
provided for these models). As an example, to use the `GeneFull_Ex50pAS`
gene model with the `EM` option:

```
library(Seurat)

# Read the raw count matrix
sample.mtx <- ReadMtx(mtx = "GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx",
  cells = "GeneFull_Ex50pAS/raw/barcodes.tsv",
  features = "GeneFull_Ex50pAS/raw/features.tsv")

# Convert the count matrix into a SeuratObject
sample <- CreateSeuratObject(counts = sample.mtx)

# Read the filtered barcodes
filtered.bc <- read.table("GeneFull_Ex50pAS/filtered/barcodes.tsv")

# Subset the SeuratObject to contain just the filtered barcodes
sample <- sample[, filtered.bc$V1]
``` 


Generating the count matrices using other gene models or multi-mapping
options are identical to the above, but just use the appropriate count
matrices.


## Analysis of 10x multiome datasets
The 10x multiome protocol performs single-nucleus sequencing as opposed
single-cell sequencing done in the 10x 3' protocol.  

The analysis of the gene expression data from the 10x multiome protocol
is identical to the 10x 3' protocol. However, the multiome data is
expected to contain a higher fraction of reads aligning to introns since
the nucleus contains more immature transcripts. The `--soloFeatures
Gene` which counts just the exonic reads should not be used with the
multiome protocol. 

## Analysis of 10x 5' datasets
The 10x 5' protocol sequences a transcript from the 5' end, and they
align in the opposite orientation as compared to the 3' protocol. 
To ensure that STAR counts the reads that align in the reverse
orientation to the reference genome, just this additional parameter
`--soloStrand Reverse` needs to be added to the command used for the
10x3 protocol. Be sure to verify that the length of the cell barcode and
UMI match the version of the protocol that you are using. 
