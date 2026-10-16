# oist/LuscombeU_Selenoproteins_Pipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.1 - [2026-10-16]

### `Fixed`

- `SELENOPROFILES_ASSESS`: its environment and image now include pyfaidx, which `selenoprofiles assess` imports but the selenoprofiles4 package does not declare.
- `SELENOPROFILES_ASSESS`: a rerun in the same work directory (for example after a Slurm requeue) no longer fails on decompressing its inputs; a name that is a staged input is never overwritten.

## v0.1.0 - [2026-10-16]

Initial release of oist/LuscombeU_Selenoproteins_Pipeline, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Optional `transcripts` samplesheet column: spliced alignment of transcripts or ESTs (nf-core `minimap2/align`) and transcript support of each Sec codon (`transcripts/<sample>.sec_transcripts.tsv`).
- MultiQC section "Transcript support": Sec codons inside an aligned exon of at least one transcript, per genome.

### `Fixed`

- Parameters given on the command line with Nextflow 26.04 and later (syntax parser v2), which are strings: nf-schema casts them to their schema types during validation, and the boolean options are converted before use (`--skip_secis false` no longer skips).

### `Dependencies`

- nf-schema 2.8.0

### `Deprecated`
