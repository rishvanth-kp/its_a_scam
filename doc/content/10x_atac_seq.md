# Analysis of 10x ATAC-seq datasets

## Background
This tutorial is intended to serve a walk through of a 10x ATAC-seq
pre-processing pipeline. The inputs are the FASTQ files and the outputs
are the cell-peak count matrices and/or fragment files.  These files can
then be loaded into standard tools for downstream analysis.  

The key steps in the analysis of scATAC-seq data are: 

1. Aligning the reads to the reference genome in a way that the cell barcode 
information is preserved.
2. Removal of PCR duplicated taking into account that cell barcode
information. 
3. Counting the number of fragments per cell.
4. Filter out unwanted alignments from the BAM file. 
5. Quality control.
6. Peak calling.
7. Cell-peak count matrix generation.

Alternatively, after step 2 the BAM file can be used to generate a
fragments file that can be used as an input to any standard scATAC-seq
tools such as [`signac`](https://stuartlab.org/signac/).  Downstream
peak-calling, QC, and count matrix can be generated within `signac`.

Each of these steps are described in detail below.

## Analysis of 10x ATAC-seq datasets
### Reference genome preparation
Reference genomes can be downloaded from the UCSC genome browser or
Gencode. Be sure to down the right version of the reference genome and
to keep the reference genome version consistent across all analysis.

`BWA`, the program used to align reads to the reference genome, requires
the genome to be pre-processed once before first use: 

``` 
bwa index {reference.fa} 
```


### Add cell barcode to fastq names
A bulk ATAC-seq experiment generates paired-end sequencing data
from DNA fragments that are flanked by two transposase cut sites.
Paired-end fastq files consists of a pair of line matched files which
correspond to the sequenced regions from the ends of the fragment.  10x
single-cell ATAC-seq and multiome sequencing performs ATAC-seq at a
single cell resolution by addition of a cell barcode sequence at the
transposase cut sites. This cell barcode is sequenced as a separate
fastq file. A 10x single-cell ATAC-seq experiment generates 3 read
files: paired-end fastq files that correspond to the ATAC-seq fragments,
and a third fastq file that contains a line matched cell barcode for
each pair of fragments. This cell barcode provides a cell identity for
each fragment.

[`bwa`](https://github.com/lh3/bwa) is designed for fast alignment of
single/paired-end reads to the reference genome, and at the moment,
cannot process the third FASTQ file containing the cell barcode.  This
limitation can be easily overcome by appending the cell barcode
information from the third FASTQ file to the name field of the paired
end FASTQ files. The name field is retained in the aligned SAM output of
BWA, and so the cell barcode information can be parsed for each
alignment.

The cell barcode can be appended to the paired-end FASTQ files using the
program `add_barcode_to_fastq_name`. Unfortunately, gzipped fastq files
are not supported at the moment. The fastq files need to be
decompressed, add the cell barcode, and re-compressed again (his
limitation will the fixed soon):

```
pigz -p {threads} -k -d read_1.fastq.gz read_2.fastq.gz read_I.fastq.gz

add_barcode_to_fastq_name -1 read_1.fastq -2 read_2_fastq 
 -b read_I_fastq.gz -o read_bc;

pigz -p {threads} read_bc_1.fastq read_bc_2.fastq

rm read_1.fastq read_2.fastq read_I.fastq
```

### Align the reads to the reference genome
The reads are aligned to the reference genome using `bwa`. The input to
BWA are the pre-built pre-built reference genome index and the fastq
files. 

The input reads usually contain primer and adapter sequence form the
library preparation. `bwa` performs a local alignment on the reads and
soft clips any of these unwanted technical bases that would not align to
the reference genome. Therefore, it is not required to explicitly trim the
adapter sequences prior to processing with `bwa`.

```
bwa mem -t {threads} {path_to_bwa_index} read_bc_1.fastq.gz read_bc_2.fastq.gz 
  1> sample.sam 2> sample_bwa_log.txt
```

The output of `bwa` is a sam file. The first few lines of a sam file
contain a header that starts with `@`. The subsequent lines contain one
entry for each read, and thus two entries for a read pair (but could
contain more than two entries when a read is split and aligned to more
than one location on the reference, these are called supplementary
alignments).  A detailed specification of a sam file format can be found
[here](https://samtools.github.io/hts-specs/SAMv1.pdf).

Sam files are convenient to view on a text editor. However, these files
consume a lot of disk space since they are not compressed. Bam files are
the compressed version of sam files, and the formats can be converted
with [`samtools`](https://www.htslib.org/) (the sam file can then be
deleted). The supplementary alignments can also be removed at the same
time.  

```
samtools view -@ {threads} -b -F 0x800 -o sample.bam sample.sam
```

### Mark PCR duplicates
Several steps in the library preparation process of going from a few
pico-gram levels of DNA in a cell to typically nano-gram amounts of DNA
required for Illumina library preparation involves PCR amplification.
This could lead to capturing multiple copies of the same DNA fragment in
the sequencing process. Since these reads arise from the same fragment,
the presence of more than one read pair does not provide any additional
information, but could lead to a bias in the downstream analysis.

PCR duplicates for bulk sequencing experiments are removed based on the
assumption that two read pairs that align to exactly the same location
on the reference genome are more likely due to a PCR duplicate rather
than two different fragments. Multiple read pairs that align exactly to
the same reference genome location are considered as PCR duplicates. A
graphical representation of a PCR duplicate can be here.

For single-cell sequencing experiments, two fragments from different
cells could align to the same reference genome location (and so meet the
above criteria for being a PCR duplicate). However, these should
not be considered as a PCR duplicates since the fragment belongs to
different cells. Therefore for removing PCR duplicates in single-cell
data, an additional criteria of read pairs having distinct cell barcodes
is applied.  

PCR duplicates are removed with `samtools` with a set of four commands.
The first command sorts the alignments in a bam file based on the read
names so that all the alignments from the same read pair are adjacent in
the bam file.

```
samtools collate -@ {threads} -o sample_collate.bam sample.bam
```

`bwa` and other mapping tools typically do not provide the correct
information about the insert size for the read pairs. 
`samtools fixmate` fixes these.

```
samtools fixmate -@ {threads} -m sample_collate.bam sample_fixmate.bam
```

The bam file is then sorted based on the genome coordinates to
facilitate removing PCR duplicates.

```
samtools sort -@ {threads} -o sample_sorted.bam sample_fixmate.bam
```

Finally, the PCR duplicates are marked and only the best alignment for
each set of duplicates are unmarked.

```
samtools markdup -@ {threads} --barcode-name sample_sorted.bam
  sample_makrdup.bam
```

The final output of these steps is a bam file in which all the duplicate
reads are marked in the SAM flag (`0x400`). 

### Barcode matching for 10x multiome data
This step in needed only for joint analysis of gene-expression and
ATAC-seq 10x multiome data. The scRNA-seq and scATAC-seq data from the
same cell have different cell barcodes, and so they cannot easily be
associated with one another. This step converts the scATAC-seq barcode
to the respective scRNA-seq barcode from the same cell so that data can
be associated. The location of these cell barcode files can be 
found [here](
https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist). 

```
barcode_match -a sample_makrdup.bam -g gex_barcode.txt 
  -b atac_barcode.txt -o sample_bc_match.sam
``` 

Only the cell barcodes that are contained in the provided barcode files
are included in the output. A cell barcode might not be contained in the
barcode file do to a sequencing error in the barcode. It is possible to
error correct the barcode prior to this step, which would increase the
number of fragments per cell. However, from experience, the advantage
gained from this step is minimal and so we have not included cell
barcode error correction as a part of this pipeline.  


### A split in the road
At this stage, the subsequent pre-processing step can be done any one
of two ways:

1. Perform peak calling and count matrix generation. To do this,
skip the next step and proceed with the analysis.    

2. Convert the barcode matched sam file to a fragments file.  The
fragments file is taken as input by scATAC-seq analysis tools such as
[`signac`](https://stuartlab.org/signac/), and all the downstream
pre-processing steps can be done there.  To do this, perform only the
next step. 

### Convert sam/bam file to fragments file
While the SAM file contains all the information in regard to an aligned
read, it typically has more information than for most applications and so
can occupy a lot more space than necessary. The information that is needed
for most downstream ATAC-seq analysis are the start and end location on
the reference genome and the cell barcode for each fragment. The
fragments file contains just this information. The fragments file is
similar to a BED file: it is a 5 column tab-separated file with the
chromosome, start location, end location, cell barcode, and number of
PCR duplicates for each fragment. 

```
bam_to_fragments -a sample_bc_match.sam -q 30 -o sample_unsort
```

This will generate a file called `sample_unsort_fragments.tsv` that
needs to be sorted by genome coordinates, compressed, and indexed. 

```
cat sample_unsort_fragments.tsv | sort -k1,1 -k2,2n > sample_fragments.tsv

bgzip sample_fragments.tsv

tabix --preset=bed sample_fragments.tsv.gz
```

The compressed `sample_fragments.tsv.gz` file can be directly be loaded
into most downstream scATAC-seq analysis programs. These programs
typically look for the indexed file `sample_fragments.tsv.gz.tbi` in the
same directory that contains the fragments file. The peak calling and
count matrix generation can be done within these programs.

Of note, the fifth column in the generated fragments file is always 1
since we have already removed PCR duplicates. This should not affect any
downstream analysis, but it might affect any QC plots that are related
to the number of number of PCR duplicates.


### Count number of fragments per cell barcode
Most cells in a 10x experiment would contain too few alignments than
needed for any meaningful analysis. Typically any cells that have fewer
than 1000 alignments are discarded (although this number can change
depending on the downstream analysis). 

After aligning a sample to the reference genome and removing PCR
duplicates, we can determine the number of aligned fragments in each
cell.

```
barcode_count -a {input} -o sample -q 30 -m 1000 
  -w atac_barcode.txt
```

Where `{input}.bam` file is `sample_makrdup.bam` for 10x ATAC-seq
protocol or the `sample_bc_match.sam` for the 10x multiome protocol. The
`actc_barcode.txt` can be obtained [here](
https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist), 
and it contains a list of all
the known cell barcodes. 

This program will generate a file called `sample_bc_counts.txt`. This
is a 2 column TSV file that contains the cell barcode in the first column
and the number of alignments fragments for that cell in the second
column. Only cell barcodes that are in the provided barcode list and
have at least a 1000 aligned fragments are included in the output. 


### Filter unwanted alignments from bam file
Now, we can get rid of all the cell barcodes that have fewer than a 1000
aligned fragments as these are very low quality cells, and are not
needed for downstream analysis.

```
bc_filter_alignments {input} -b sample_bc_counts.txt -q 30
  -o sample_bc_filtered.sam

samtools  view -@ {threads} -o sample_bc_filtered.bam sample_bc_filtered.sam
```

Where `{input}.bam` file is `sample_makrdup.bam` for 10x ATAC-seq
protocol or the `sample_bc_match.sam` for the 10x multiome protocol.

This will generate a new bam file that have just the alignments that
belong to cells contained in `sample_bc_counts.txt`. 
 

### Peak calling
Peak calling then be done using standard tools such as [MACS](
https://macs3-project.github.io/MACS/docs/callpeak.html). The input
to the peak calling program would be the bam file generated in the
above step. Peak callers generally produce a list of peaks in bed
format (or the peaks can be easily converted to a bed format).

### Cell-peak count matrix generation
We can finally generate the cell-peak matrix using

```
bc_count_matrix -a sample_bc_filtered.bam -b sample_bc_counts.txt -t
peaks.bed -q 30 -o sample
``` 

This will generate a file called `<output_prefix>_region_counts.txt`
that contains the count matrix. The rows of the file are the regions and
the columns are the cells. The first 4 columns contain the the chrom
name, start, end, and region name. The cells start from the 5th columns.
This file has a header line that has the cell barcodes.  
