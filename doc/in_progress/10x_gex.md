# Analysis of 10x gene-expression datasets

## Background
Single cell transcriptomic analysis is widely used to for unbiased cell
type and cell cell state identification in tissues. A vast majority of
single cell RNA-seq experiments are done using the 10x platform. 
10x offers several variants of scRNA-seq protocol -- most common is the
the 3' protocol, and other variants include the 5' protocol and the
multiome (single nucleus RNA-seq + scATAC seq) protocol. The analysis for
these protocols are similar with only minor differences between them.
Since the 10x 3' protocol is the most commonly sues, we will use that to
describe the analysis, and then point out the differences when using other
protocols. 

This in this tutorial we will start with the fastq files of sames that
were sequencing using one of these 10x protocols and on an Illumina
sequencer. We will use the STARsolo alignment software to align these
read and generate a cell-gene count matrix. STARsolo is a very
documented. This tutorial is a starting point to get started from the
scRNA-seq analysis from scratch. Reading the STAR documentation and the
paper is highly recommended, especially for using the software beyond
what is described here. 

## Analysis of 10x 3' datasets
**Reference genome and annotation preparation:**
The reference genome and the gene annotations can be downloaded from
Gencode. 

To align reads and generate the count matrices using STARsolo, a genome
and the annotations should be pre-processed to generate an index. This
steps needs to performed only once before the first use (for a specif
version of the the genome and annotation).  

```
STAR --runMode genomeGenerate --runThreadN {threads} 
  --genomeDir {output_dir} --genomeFastaFiles {reference.fa} 
  --sjdbGTFfile {annotation.gtf}
```


**Aligning reads and count matrix generation:**

```
STAR --genomeDir {genome_dir} --runThreadN {threads} 
  --readFilesIn {read_1.fastq.gz} {read_I.fastq.gz} --readFilesCommand zcat
  --outFileNamePrefix {output_prefix}
  --soloCBwhitelist {inclusion_list.txt}
  --soloType CB_UMI_Simple --soloUMIlen {umi_lenght}
  --soloBarcodeReadLength {1/0} --clipAdapterType CellRanger4 
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
  --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR
  --soloCellFilter EmptyDrops_CR
  --soloFeatures GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS
  --soloMultiMappers EM
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
  --outSAMtype BAM SortedByCoordinate --outFilterScoreMin 30 
```

**Choice of gene model for count matrix generation:**

**Quality control:**

## Analysis of 10x 5' datasets
