# Analysis of Parse Biosciences whole-transcriptome datasets

## Background
Single-cell RNA-seq relies on barcoding each cell with a unique
barcode, so that the reads generated can then be de-multiplexed and
assigned to a cell. The [combinatorial
barcoding](https://www.science.org/doi/10.1126/science.aam8940) or the
[SPLIT-seq](https://www.science.org/doi/10.1126/science.aam8999)
protocol is an approach to accomplish this without the need for
additional instrumentation to isolate individual cells. Instead it uses
each cell as an isolated compartment to tag all transcripts within a
cell with a unique combination of barcodes.

The cells are mildly fixed and permeabilized so that the transcripts are
contained within the cell. The cells are split and deposited into wells
(typically in 96 well plates) that contain a unique barcode in each
well. Since the number of cells is much larger than the number of
wells, each well would contain multiple cells that get the same barcode.
The cells from all the wells are pooled and then split again into
multiple wells each with a unique barcode. Every cell now gets a second
barcode. The probability of two cells going to the same well in both
the barcoding rounds is very low, and thus each cell would have a
unique combination of cell barcodes. If the starting number of cells
is high, then the process of pooling the cells and splitting them to
add another cell barcode can be repeated multiple times to ensure that
the probability of multiple cells getting the same combination of barcodes
is low.

After sequencing, each read pair contains both the cell barcodes and
the transcripts. Each read can be assigned to a cell using this cell
barcode information. Conceptually, aligning a read from the SPLIT-seq
protocol is identical to the 10x protocol, but with much longer cell
barcodes split over multiple positions on a read.

## A brief description of the Parse whole-transcriptome protocol
The Parse Biosciences whole-transcriptome protocol is an implementation
of the combinatorial barcoding approach for scRNA-seq. This protocol
uses 3 sets of cell barcodes with UMIs, and allows for multiplexing
cells from multiple samples.

Briefly, cells from one or more samples are deposited into distinct
wells in the first barcode plate. Each well contains an oligo that is
attached to a transcript on one end (to serve as the reverse
transcriptase primer) and has an overhang on the other end with the cell
barcode sequence.  Each well contains two kinds of oligos: (1) poly-T
oligos that attach to the poly-A tails at the 3' end of the transcript,
and (2) random hexamer oligos that attach to random locations on the
transcripts. Each of these oligos contains a distinct cell barcode, so a
well contains two different barcodes: one for the 3' primers and the
other for the internal random hexamer primers. The advantage of having
two different barcodes is that we can computationally separate out the
two kinds of reads and generate either a 3' count matrix or a
whole-transcriptome count matrix, as needed for the downstream
application. Moreover, since each well contains a unique set of
barcodes, the samples in distinct wells can be separated out by this
barcode. The RNA molecules are then reverse-transcribed.

The cells are then pooled and split to add the second set of barcodes,
and the process is repeated to add the third set. The cells are pooled
again and split into 8 different tubes. The cDNA amplification and
library preparation are then performed. Since each of these 8 tubes is
processed and sequenced independently, this acts as an implicit 4th cell
barcode.


## Analysis of Parse datasets
### Reference genome and annotation preparation:
The reference genome and the gene annotations can be downloaded from
[Gencode](https://www.gencodegenes.org/).

To align reads and generate the count matrices using STARsolo, a genome
and the annotations should be pre-processed to generate an index. This
step needs to be performed only once before the first use (for a
specific version of the genome and annotation).

```
STAR --runMode genomeGenerate --runThreadN {threads} \
  --genomeDir {output_dir} --genomeFastaFiles {reference.fa} \
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
either originate from the 3' end captured by a poly-T primer, or an
internal region captured by a random hexamer primer.  The other read
contains the three cell barcodes and the UMI that are at fixed locations
on the read. For example, for the WTv2 protocol: the UMI is at location
0-9 (0-based, inclusive intervals), the barcode 3 is at location 10-17,
barcode 2 at location 48-55, and barcode 1 at location 78-85. These
bases in between the barcodes and UMIs are filler bases and can be
ignored. These positions of the UMI and barcode locations in the read
vary by protocol, and can be found in the protocol user manual.

Furthermore, a sample is split into 8 different sublibraries. Each
sublibrary gets a different Illumina barcode, and they are thus
sequenced as 8 independent libraries each yielding a pair of fastq
files. Note that this Illumina barcode implicitly acts as the 4th cell
barcode.

### Alignment and gene quantification with STAR

**Aligning one sublibrary**:

STARsolo does the following steps:

1. Aligning the fastq files to the reference genome, taking into account
that a read can span a splice junction.
1. Correcting the cell barcodes and UMIs to remove potential sequencing
errors. STARsolo takes into account that a cell barcode consists of
3 distinct sequences.
1. Assigning a read to a gene based on the defined gene model.

```
STAR --genomeDir {genome_dir} --runThreadN {threads} \
  --readFilesIn {read_1.fastq.gz} {read_I.fastq.gz} \
  --readFilesCommand zcat \
  --outFileNamePrefix {output_prefix} \
  --soloType CB_UMI_Complex \
  --soloCBwhitelist {cb_2_3.txt} {cb_2_3.txt} {cb_1.txt} \
  --soloCBposition {cb_3_pos} {cb_2_pos} {cb_1_pos} \
  --soloUMIposition {umi_pos} --clipAdapterType CellRanger4 \
  --clip5pAdapterSeq AAGCAGTGGTATCAACGCAGAGTGAATGGG \
  --soloCBmatchWLtype EditDist_2 --soloCellFilter EmptyDrops_CR \
  --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloFeatures GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS \
  --soloMultiMappers EM --soloBarcodeReadLength 0 \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
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

`--soloCBwhitelist {cb_2_3.txt} {cb_2_3.txt} {cb_1.txt}`: One-column
files that contain the lists of valid barcodes for barcode 3, barcode 2,
and barcode 1, respectively (barcodes 2 and 3 use the same list). These
files can be made from the `sequence` column of the CSV cell barcode
files from Parse; `cb_1.txt` must contain both the poly-T and the random
hexamer barcodes. Keep the barcodes in this order (3, 2, 1) here and in
`--soloCBposition`, since `parseRandhexCollapse.R` (see below) finds
each barcode by its position in the cell barcode.

`--soloCBposition {cb_3_pos} {cb_2_pos} {cb_1_pos}`:
For the WTv2 protocol the positions of the barcode are
`0_10_0_17 0_48_0_55 0_78_0_85`. The format of these parameters is
described under the [STARsolo documentation](
https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#complex-barcodes).
Refer to the user manual for the other versions of the protocol.

`--soloUMIposition {umi_pos}`: For the WTv2 protocol the UMI position
is `0_0_0_9`. The format of this parameter is identical to the above.
Refer to the user manual for the other versions of the protocol.

Some of the parameter choices used here were derived from
[here](https://github.com/alexdobin/STAR/issues/1517) and
[here](https://github.com/alexdobin/STAR/issues/1945).

The choice of gene model and the way to deal with multi-mapping reads
are identical to the 10x protocol, and are described
[here](10x_gex.md).

**Aligning all sublibraries:**
STAR needs to be run independently for each sublibrary, and so the
above step needs to be done 8 times.

### Cell-gene count matrix generation
Now that the reads are aligned and the count matrices are generated for
all barcodes, we need to do the following:

1. Get rid of cells that contain too few read counts.
2. Demultiplex the samples based on the first cell barcode.
3. Generate a 3' cell-gene count matrix by counting only the reads that are
aligned to the 3' end (i.e., the poly-T primer reads).
4. Generate a whole-transcriptome cell-gene count matrix by merging
the counts from the poly-T reads and the random hexamer reads.



The first cell barcode is involved in several different functions:

1. It acts as the sample barcode, where each sample is assigned one or more
barcodes. The sample barcodes are based on the wells to which the sample
was added in the first round of barcoding.
2. Within each sample it functions as a standard cell barcode.
3. Each cell contains two different barcode sequences: one for fragments
starting at the 3' end (using a poly-T primer) and another for fragments
that start internal to the transcript (using a random hexamer primer).

Taking these into account, generating the count matrix is a bit of an
intricate process. We can use the
[`parseRandhexCollapse.R`](https://github.com/rishvanth-kp/its_a_scam/blob/main/scripts/parseRandhexCollapse.R)
script in the `scripts` directory to generate Seurat objects with the
count matrices. The script requires the R packages `tidyverse`,
`optparse`, and `Seurat`.

```
parseRandhexCollapse.R -r {count_matrix_dir} -m {matrix.mtx} \
  -n {min_count} -s {sample_ids.csv} -l {sublibrary_id} \
  -a {cb_1.csv} -b {cb_2_3.csv} -o {outfile_prefix}
```

where:

`-r {count_matrix_dir}`: Path to a STARsolo `raw` folder that contains
the `barcodes.tsv`, `features.tsv`, and matrix files (for example,
`{output_prefix}Solo.out/GeneFull_Ex50pAS/raw`).

`-m {matrix.mtx}`: The name of the matrix file to use, for example,
`UniqueAndMult-EM.mtx` for the `EM` option. [Default: matrix.mtx].

`-n {min_count}`: Keep only cells that have at least `min_count`
counts. This is used only if `-f` is not set. [Default: 500].

`-f {filtered_barcodes}`: Optional. A file with the cell barcodes to
keep, one per line and without a header (for example, the `barcodes.tsv`
file in the STARsolo `filtered` folder). If this is set, `-n` is
ignored.

`-s {sample_ids.csv}`: A 2-column CSV file (without a header) that
contains the well ID in the first column and the sample name in the
second column. This file is used for sample demultiplexing based on the
wells to which samples were loaded in the first cell barcode plate. The
well IDs must match the `well` column of the `-a` file.

`-l {sublibrary_id}`: The ID for the sublibrary. Typically a number
between 1 and 8 indicating the sublibrary tube ID that acts as the
implicit 4th barcode.

`-a {cb_1.csv}`: The CSV file from Parse with the first cell barcodes.
It must contain the columns `sequence`, `bci` (barcode index), `well`,
and `stype` (`T` for poly-T and `R` for random hexamer barcodes). The
`T` and `R` barcodes must be in the same order, so that the n-th `T`
barcode and the n-th `R` barcode are from the same well; this is not
checked.

`-b {cb_2_3.csv}`: The CSV file from Parse with the second and third
cell barcodes. It must contain the columns `sequence` and `bci`.

`-o {outfile_prefix}`: Output file prefix.

The script writes two Seurat objects:

1. `<outfile_prefix>_3prime.Rds`: The 3' count matrix, with only the
poly-T reads.
2. `<outfile_prefix>_collapse.Rds`: The whole-transcriptome count
matrix, in which the random hexamer counts of each cell are added to its
poly-T counts (if the random hexamer barcode of the cell is also kept by
`-n` or `-f`).

The cells are named `<bci1>_<bci2>_<bci3>__<sublibrary_id>`, where
`bci1`, `bci2`, and `bci3` are the two-digit barcode indices from the
Parse CSV files, and the `sample` metadata column contains the sample
name from `-s`.

This step needs to be done independently for each sublibrary, and so
eight count matrices will be generated.

### Merging sublibraries
We can then merge all sublibraries into one count matrix by simply
appending all the matrices by columns. With Seurat, this can be done
with the `merge` function.
