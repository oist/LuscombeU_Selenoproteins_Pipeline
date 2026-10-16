# oist/LuscombeU_Selenoproteins_Pipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0dev - [unreleased]

### `Added`

- `--selenoprofiles_overrides`: a directory or `.tar.gz` archive of Selenoprofiles profiles that replace the standard profiles of the same name (new module `SELENOPROFILES_OVERRIDE`), for lineage-specific profiles.

### `Changed`

- Resource labels set from the execution trace of nine appendicularian genomes (56 to 644 Mb), with room for mammalian genomes: `SELENOPROFILES_SELENOPROFILES` from `process_high` to `process_medium`, `SECMARKER` from `process_medium` to `process_low`, and the nf-core `minimap2/align` module from `process_high` to `process_medium` (recorded as a module patch).

### `Fixed`

- Selenoprofiles ran without intron prediction (`-no_splice`) when the sample had no `domain`, which fragmented eukaryotic predictions (seen in the module test); it is now used only for bacteria and archaea.

- Tasks rerun in the same work directory (for example after a Slurm requeue) no longer fail on decompressing their input (the shell runs with `noclobber`, `bash -C`): the local modules remove a leftover decompressed file before writing it, and stop with an error, without writing, if that name is a symbolic link (a staged input), so that an input file can never be overwritten.

### `Dependencies`

### `Deprecated`

## v0.1.2 - [2026-10-16]

### `Fixed`

- `SELENOPROFILES_ASSESS`: transcripts are identified by the `Parent` of their CDS features (`-cg Parent`). The default (`ID`) failed on annotations whose CDS features have no `ID`, such as the ANISEED OSKA2016 annotation; on RefSeq annotations both give the same classes.

### `Changed`

- Test profile: the `hs_chr22` sample now has an annotation (SELENOO, with CDS features carrying only `Parent`) and a transcript covering its Sec codon (`assets/test/`), so that `-profile test` runs every process, including `SELENOPROFILES_ASSESS`, `SEC_ANNOTATION_CHECK` and a non-zero transcript support.

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
