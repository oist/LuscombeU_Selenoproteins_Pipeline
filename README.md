# oist/LuscombeU_Selenoproteins_Pipeline

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/oist/LuscombeU_Selenoproteins_Pipeline)
[![GitHub Actions CI Status](https://github.com/oist/LuscombeU_Selenoproteins_Pipeline/actions/workflows/nf-test.yml/badge.svg)](https://github.com/oist/LuscombeU_Selenoproteins_Pipeline/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/oist/LuscombeU_Selenoproteins_Pipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/oist/LuscombeU_Selenoproteins_Pipeline/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/oist/LuscombeU_Selenoproteins_Pipeline)

## Introduction

**oist/LuscombeU_Selenoproteins_Pipeline** annotates selenoprotein genes, selenocysteine (Sec) machinery genes and selenocysteine tRNA (tRNA-Sec) genes in genome assemblies. It takes one FASTA file per genome, optionally with a GFF3 annotation, and reports gene predictions, the agreement between independent methods, the status of the given annotation at each selenoprotein gene, and a per-genome Sec trait call in a MultiQC report.

![Pipeline overview](docs/images/selenoproteins_metro.svg)

1. Assembly statistics ([`seqkit stats`](https://bioinf.shenwei.me/seqkit/))
2. Selenoprotein and Sec machinery genes ([`Selenoprofiles`](https://github.com/mariottigenomicslab/selenoprofiles4))
3. Cross-check of Sec codons with protein-to-genome alignments of the profile sequences ([`miniprot`](https://github.com/lh3/miniprot))
4. SECIS elements in the 3' flanks of eukaryotic predictions, as evidence ([`Infernal`](http://eddylab.org/infernal/) with an Rfam model)
5. Assessment of a given annotation at each prediction (`selenoprofiles assess`)
6. tRNA-Sec genes ([`Secmarker`](https://github.com/didacs/secmarker), Infernal with Rfam RF01852, optionally [`tRNAscan-SE`](https://github.com/UCSC-LoweLab/tRNAscan-SE)), merged into loci
7. Summary tables and Sec trait call ([`MultiQC`](http://multiqc.info/))

Selenoprofiles finds members of known selenoprotein families only. See [usage](docs/usage.md#interpretation-limits) for the limits of each step.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet with one genome per row:

```csv title="samplesheet.csv"
sample,fasta,domain,profiles,annotation
human,GRCh38.fa.gz,eukaryota,,GRCh38.gff3.gz
ecoli,ecoli_k12.fna,bacteria,,
```

`domain` is `eukaryota`, `bacteria` or `archaea`; `profiles` and `annotation` are optional. Then run:

```bash
nextflow run oist/LuscombeU_Selenoproteins_Pipeline -r main \
   -profile <conda/mamba, or wave,singularity / wave,docker> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

Selenoprofiles 4.5.8 and Secmarker 1.0.0 have no Biocontainers image yet (Bioconda recipes pending). Use
`-profile conda` or `mamba`, or add `wave` to a container profile so that Wave builds images from the module
conda environments. With a container profile alone, these two tools must be installed on the host.

Test the installation with `-profile test,conda` (the test downloads small genomes from nf-core/test-datasets).

Without `--selenoprofiles_data`, the pipeline downloads the Selenoprofiles data from the authors' server. See [usage](docs/usage.md) and [output](docs/output.md) for details.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

oist/LuscombeU_Selenoproteins_Pipeline was written by Charles Plessy (Luscombe Unit, OIST) with Claude (Anthropic).

We thank the following people for their extensive assistance in the development of this pipeline:


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use oist/LuscombeU_Selenoproteins_Pipeline for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
