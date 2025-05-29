# Analysis of 10x ATAC-seq datasets

## Background

A thypical ATAC-seq experiment generated paried-end quenced data that
are greneated from DNA fragments that are flanked by two transposast 
cut sites. 

10x single-cell ATAC-seq and multiome sequenciing perfoms ATAC-seq at a
single cell resoluton by addion on a cell barcode sequnec the
transposate cut sites. A detatiled explation of this
procdure can be found here xyz. This cell barocde is sequenced as a 
sepatrate read file. Therefore, a 10x single-cell ATAC-seq enxprieemtn
geneated 3 read files: a piair of files that coreespiont to the
ATAC-seq fragments, and a third file that contains a line matched cell
barode for each pair of fragmetns. This cell barcode provides a cell
identyut for each fragment for deconvoluton of the each fragment to a
cell. 

The key steps in the analysis of scATAC-seq data are: 
1. Aligning the peaired-end reads to the reference genome in a way that
preserved the cell barocode information in the thiord read.
2. Removal of PCR duplicated taking into accoiunt that cell barcoide
infomation. 
3. Converting the alinmened file in to a fragments file for use by other
stanged scATAC-seq tool.

The following steps are optional and could be perforned either using
other tools or performed as needed for downsosterema analysis:
1. Counting the number of fragmetns per cell.
2. Filter out unwated alignmetns from the SAM/BAM file. 
3. Quality control.
4. Peak calling.
5. Other downstream analysis.

Each of these steps are described in detail below.

## Analysis of 10x ATAC-seq datasets
**Reference genome preparation:**


**Add cell barcode to fastq names:** 
A bulk ATAC-seq expreiment typically consitst of paired end FASTQ files
contain a pair of line matched files wich correspond to the sequgneced
regrions from the ends of the fragment between Tn5 transposase cut
sites. For a single cell ATAC-seq   


**Align the reads to the reference genome**

**Remove/mark PCR duplicates**

**Count number of fragments per cell barcode**

**Filter unwanted alignments from SAM/BAM file**

**Quality control**

**Merge samples**

**Peak calling**
