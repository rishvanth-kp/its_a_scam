# 10x ATAC-seq count matrix generation: quick start 

This document is intended to get started with generating a count matrix
starting from fastq files generated from the 10x scATAC-seq or Multiome
protocol. It does not provide any explanation of the steps, go here 
for a detailed explanation.

**Add cell barcode to fastq names**: 
Decompress the input fastq files, add barcode to fastq names, compress
the output, and delete the temporary files: 

``` 
pigz -p {threads} -k -d read_1.fastq.gz
read_2.fastq.gz read_I.fastq.gz

add_barcode_to_fastq_name -1 read_1.fastq -2 read_2_fastq 
 -b read_I_fastq.gz -o read_bc

pigz -p {threads} read_bc_1.fastq read_bc_2.fastq

rm read_1.fastq read_2.fastq read_I.fastq
```

**Align the reads to the reference genome**: 
Align reads to the reference genome, convert the sam output to bam, and
delete the temporary files:

```
bwa mem -t {threads} {path_to_bwa_index} read_bc_1.fastq.gz read_bc_2.fastq.gz 
  1> sample.sam 2> sample_bwa_log.txt
samtools view -@ {threads} -b -F 0x800 -o sample.bam sample.sam
rm sample.sam
```

**Mark PCR duplicates**:
Marking PCR duplicates is done in a series of steps with `samtools`:
```
samtools collate -@ {threads} -o sample_collate.bam sample.bam

samtools fixmate -@ {threads} -m sample_collate.bam sample_fixmate.bam

samtools sort -@ {threads} -o sample_sorted.bam sample_fixmate.bam

samtools markdup -@ {threads} --barcode-name sample_sorted.bam
  sample_makrdup.bam

rm sample.bam sample_collate.bam sample_fixmate.bam sample_sorted.bam 
```


**Barcode matching for 10x multiome data**:
Skip this step for 10x scATAC-seq.  This step in needed only for joint
analysis of gene-expression and ATAC-seq 10x multiome data. Match the
ACTA-seq barcodes to the GEX barcodes, convert the output to bam, and
delete the temporary files:
```
barcode_match -a sample_makrdup.bam -g gex_barcode.txt 
  -b atac_barcode.txt -o sample_bc_match.sam

samtools view -@ {threads} -o sample_bc_match.bam sample_bc_match.sam

rm sample_bc_match.sam
``` 

**Convert sam/bam file to fragments file**:
Covert the bam file to a fragments file and index it:

```
bam_to_fragments -a sample_bc_match.bam -q 31 -o sample_unsort

cat sample_unsort_fragments.tsv | sort -k1,1 -k2,2n > sample_fragments.tsv

bgzip sample_fragments.tsv

tabix --preset=bed sample_fragments.tsv.gz

rm sample_unsort_fragments.tsv
```

The compressed `sample_fragments.tsv.gz` file can be directly be loaded
into most downstream scATAC-seq analysis programs. These programs
typically look for the indexed file `sample_fragments.tsv.gz.tbi` in the
same directory that contains the fragments file. The peak calling and
count matrix generation can be done within these programs.
