# List of things to get done

## Documentation
1. Add a TSS enrichment section (`bc_tss_enrichment` and
`atacTssEnrichment.R`) to the scATAC-seq QC procedure.

## Enhancements
1. Update programs to read and write from compressed fastq files.
1. Update programs to read and write from bam files.
1. Profile code to identify bottlenecks.
1. Add multi-threading capabilities.

## Code clean-up
1. Add `-t` (barcode from a SAM tag) to the help text of `bc_coverage` and
`bc_fragment_coverage`.
1. Fix the R script help text: `-f` and `-s` in `atacFragLength.R`, `-g` in
`gexGenebodyCoverage.R`, and `-t` in `atacTssEnrichment.R` (all copied from
`atacFlagstat.R`), and the line break in the `-n` help of
`parseRandhexCollapse.R`.
1. Rename the `bc_tss_enrichment` output `<prefix>_tss_enrichmet.txt` to
`<prefix>_tss_enrichment.txt`. This may break existing scripts; update the
documentation and `atacTssEnrichment.R` usage too.
1. Review the differences between `bc_flagstat` and `samtools flagstat`:
counts are over all alignments (samtools uses primary alignments only),
QC-failed reads are not split out, and an alignment with both 0x100 and
0x800 is counted as both secondary and supplementary.
1. Fix help-text typos (`delimeter`, `minimun`, `defauly`, missing space
before `[default`) and comment typos across the programs.

## Bug fixes
Bugs and enhancements are tracked as
[GitHub issues](https://github.com/rishvanth-kp/its_a_scam/issues).
