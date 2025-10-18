# Analysis of Parse Bioscience whole-transcriptome datasets

## Background


## A brief description of the Parse whole-transcriptome protocol



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

### Description of reads generated using the Parse protocol
The Parse protocol generates paired-end reads, however, these are not
the typical paired-end reads. One read contains the transcript and the
other read contains the cell barcode and UMI.  The transcript read can
either originate form the 3' end captured by a poly-T primer, or an
internal region captured by a random hexamer primer.  The other reads
contain several the three cell barcodes and the UMI that are at fixed
locations on the read. For example, for the WTv2 protocol: the UMI is
at location 0-9 (0-based, inclusive intervals), the barcode 3 is at
location 10-17, barcode 2 at location 48-55, and barcode 1 at
location 78-85. These based in between the barcodes and UMIs are
filler bases and can be ignored. These positions of the UMI and barcode
locations in the read varies by protocol, and can be count in the
protocol user manual. 

Furthermore, a sample is split into 8 different sublibraries. Each
sublibrary gets a different Illumina barcodes, and are thus sequenced as
8 indipanes libraries each yielding a pair of fastq files. Note that
this Illumina barcode implicitly acts at the 4th cell barcode. 

### Alignment and gene quantification with STAR

**Aligning one sublibary**:

STARsolo does the following steps:

1. Aligns the fastq files to the reference genome taking into account
that a read can span a splice junction.
1. Cell barcode and UMI error correction to remove potential sequencing
errors. STARsolo takes into account that the a cell barcode consists of
3 distinct sequences.  
1. Assigning a read to a gene based on the defined gene model. 

```
STAR --genomeDir {genome_dir} --runThreadN {threads}
  --readFilesIn {read_1.fastq.gz} {read_I.fastq.gz} 
  --readFilesCommand zcat
  --outFileNamePrefix {output_prefix}
  --soloType CB_UMI_Complex --soloCBwhitelist
  {cb_2_3.txt} {cb_2_3.txt} {cb_1.txt}
  --soloCBposition {cb_3_pos} {cb_2_pos} {cb_1_pod}
  --soloUMIposition 0_0_0_9 --clipAdapterType CellRanger4
  --clip5pAdapterSeq AAGCAGTGGTATCAACGCAGAGTGAATGGG
  --soloCBmatchWLtype EditDist_2 --soloCellFilter EmptyDrops_CR
  --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR
  --soloFeatures GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS
  --soloMultiMappers EM --soloBarcodeReadLength 0
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
  --outSAMtype BAM SortedByCoordinate --outFilterScoreMin 30   
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

`--soloCBwhitelist {cb_2_3.txt} {cb_2_3.txt} {cb_1.txt}`:

`--soloCBposition {cb_3_pos} {cb_2_pos} {cb_1_pod}`: 
For the WTv2 protocol the positions of the barcode are 
`0_10_0_17 0_48_0_55 0_78_0_85`. The format of these parameters are
described under the [STARsolo documentation](
https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#complex-barcodes).
Refer the the user manual for the other versions of the protocol. 

`--soloUMIposition {umi_pos}`: For the WTv2 protocol the UMI position
is `0_0_0_9`. The format of this parameter is identical to the above. 
Refer the the user manual for the other versions of the protocol. 

Some of the parameter choices used here were derived from 
[here](https://github.com/alexdobin/STAR/issues/1517) and 
[here](https://github.com/alexdobin/STAR/issues/1945) 

The choice of gene model and the way to deal with multi-mapping reads
are identical to the 10x protocol, and are described 
[here](). 

**Aligning all sublibraries:**
STAR needs to be run indecently for each sublibrary, and so the the
above step needs to be done 8 times.

### Cell-gene count matrix generation
 
**3' gene quantification**:

**Random hexamer collapsing and whole transcriptome quantification**:

### Merging sub-libraries


