# Analysis of Parse Bioscience whole-transcriptome datasets

## Background

## Analysis of Parse 10x datasets
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

Of note, this step is exactly the same for 10x GEX analysis. So if you
have the genome already indexed with STAR, this step can be skipped.


### Alignment and gene quantification with STAR

```
STAR --genomeDir {params.indexDir} 
  --readFilesIn {input.r1} {input.rI} --readFilesCommand zcat
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
  --outSAMtype BAM SortedByCoordinate --soloCellFilter EmptyDrops_CR
  --outFilterScoreMin 30 --soloCBmatchWLtype EditDist_2
  --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR
  --soloFeatures GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS
  --soloMultiMappers EM --soloBarcodeReadLength 0
  --soloType CB_UMI_Complex --soloCBwhitelist
  {params.cbWhitelist3} {params.cbWhitelist2} {params.cbWhitelist1}
  --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85
  --soloUMIposition 0_0_0_9 --clipAdapterType CellRanger4
  --clip5pAdapterSeq AAGCAGTGGTATCAACGCAGAGTGAATGGG
  --runThreadN {threads} --outFileNamePrefix {params.outPrefix}
```

Where:

`--genomeDir {genome_dir}`: The directory containing the indexed
reference genome (the `output_dir` in the previous step). 

`--runThreadN {threads}`: Number of threads to use. 

`--readFilesIn {read_1.fastq.gz} {read_I.fastq.gz}`: `read_1.fastq.gz`
is the read file that contains the transcript and `read_I.fastq.gz` is
the read file that contains the cell barcode and UMI.

`--outFileNamePrefix {output_prefix}`: The prefix used for the output
file names.

`--soloCBwhitelist {params.cbWhitelist3} {params.cbWhitelist2} {params.cbWhitelist1}`:

 
### 3' gene quantification

### Random hexamer collapsing and whole transcriptome quantification

### Merging sub-libraries

### Sample demultiplexing 
