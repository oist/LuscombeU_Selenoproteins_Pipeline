# oist/LuscombeU_Selenoproteins_Pipeline: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

## Samplesheet input

The samplesheet is a comma-separated file with a header row and one genome per row. Pass it with `--input`.

```bash
--input '[path to samplesheet file]'
```

| Column       | Required | Description                                                                                                                                                                                                                                               |
| ------------ | -------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`     | yes      | Genome name, unique in the samplesheet. Letters, digits, `_`, `.` and `-` only.                                                                                                                                                                           |
| `fasta`      | yes      | Genome assembly in FASTA format, optionally gzip-compressed (`.fa`, `.fasta`, `.fna`, with or without `.gz`).                                                                                                                                             |
| `domain`     | yes      | `eukaryota`, `bacteria` or `archaea`. Selects the default profiles, disables intron prediction for prokaryotes, restricts the SECIS search to eukaryotes and sets the tRNAscan-SE model.                                                                  |
| `profiles`   | no       | Comma-separated Selenoprofiles profiles or profile sets. Defaults: `eukarya,machinery` for eukaryotes, `prokarya` for bacteria and archaea. The sets are listed in `assets/profile_sets.tsv`.                                                             |
| `annotation` | no       | Reference annotation of the genome in GFF3 format, optionally gzip-compressed. When given, `selenoprofiles assess` reports how it represents each predicted selenoprotein. The CDS `ID` attribute must identify the transcript, as in Ensembl GFF3 files. |

```csv title="samplesheet.csv"
sample,fasta,domain,profiles,annotation
human,GRCh38.fa.gz,eukaryota,,GRCh38.gff3.gz
ecoli,ecoli_k12.fna,bacteria,,
```

An [example samplesheet](../assets/test/samplesheet_test.csv) is used by the test profile.

## Reference data

Selenoprofiles needs its profile alignments. Profiles with tag or gene-ontology filters also need the gene-ontology data and the tag database that `selenoprofiles -download` fetches from the Selenoprofiles authors' server. Without `--selenoprofiles_data`, the pipeline runs that download; the compute environment then needs network access to that server. To reuse a download, pass the directory with `--selenoprofiles_data`. A directory or `.tar.gz` archive holding only the profile alignments is also accepted, for profiles that need no tag or gene-ontology data.

The SECIS and tRNA-Sec searches with cmsearch use the covariance models given by `--secis_cm` and `--trnasec_cm`, by default the Rfam models RF00031 (SECIS_1) and RF01852 (tRNA-Sec). The models need gathering thresholds, as Rfam models have.

## Annotation checks

When the samplesheet gives an annotation (GFF3 with mRNA or transcript, exon and CDS features), the pipeline:

- runs `selenoprofiles assess`;
- classifies how the annotation treats each predicted Sec codon, with the classes of Ticó et al. 2026 (PLoS Comput Biol 22:e1013885): Well annotated (split into Sec declared or not, e.g. translated as X), Stop codon, Skipped, Downstream, Upstream, Out of frame, Spliced and Absent, at gene and transcript level;
- with `--readthrough_screen`, screens every annotated transcript for a TGA that may encode Sec: upstream of the CDS start with an ATG further upstream, through a TGA stop codon, or inside the CDS. Candidates are ranked by the best SECIS hit in their 3′ UTR (`--screen_cm`, default `--secis_cm`); transcripts with a strong SECIS hit and no TGA event are listed too, as their 3′ UTR may contain an unannotated gene.

For _Oikopleura dioica_, a SECIS model built from its own elements is provided in `assets/secis/oikopleura_dioica_secis.cm`. A second model, `assets/secis/appendicularia_secis.cm`, was trained on SECIS elements of six appendicularian genomes (both _O. dioica_ cryptic species, _O. albicans_, _O. vanhoeffeni_, _O. longicauda_, _Fritillaria borealis_); it is more sensitive across species but scores some _O. dioica_ intergenic regions up to about 30 bits, so prefer the _O. dioica_ model for genome-wide scans in _O. dioica_. Neither model has a gathering threshold, so use it with `--secis_min_bits` (for example 15) for the SECIS step:

```bash
--secis_cm ${projectDir}/assets/secis/oikopleura_dioica_secis.cm --secis_min_bits 15 --readthrough_screen
```

## Sec codon layer

Each genome gets `compare/<sample>.sec_codons.gff3`: one `stop_codon_redefined_as_selenocysteine` feature per Sec codon found by Selenoprofiles and/or miniprot, with the family, the support and the SECIS score. Following Ticó et al. 2026 (Fig 5), give this layer to a gene annotation pipeline as high-scoring evidence rather than using the Selenoprofiles gene models, which are less reliable away from the Sec codon.

## Expected selenoproteome

With `--expectation` and `--expectation_lineage`, the Sec sites of each genome are counted per family and compared with the expectation. For _Oikopleura dioica_:

```bash
--expectation ${projectDir}/assets/expectations/oikopleura_dioica.csv --expectation_lineage Oikopleura
# diploid assemblies of males (XY): --expectation_copies 2 --expectation_hemizygous TXNRD,SELENOT
```

Missing families are not always absent genes: SELENOT, SEPHS2 and SELENOS are missed by Selenoprofiles 4.5.8 in _O. dioica_ and need the annotation checks or a manual search.

## Interpretation limits

- Selenoprofiles is homology-based: it finds members of known selenoprotein families only.
- The Sec trait call (`present`, `absent`, `inconsistent`) combines selenoprotein genes with supported tRNA-Sec loci: found by two or more tools and, when Secmarker reports the locus, with the G73 discriminator base of functional tRNA-Sec genes. An `absent` call can also result from an incomplete assembly; check the assembly statistics in the report.
- SECIS hits are reported as evidence and do not filter predictions. Rfam has no archaeal SECIS model, and the pipeline does not search bacterial SECIS elements, which lie inside the coding sequence.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run oist/LuscombeU_Selenoproteins_Pipeline --input ./samplesheet.csv --outdir ./results  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run oist/LuscombeU_Selenoproteins_Pipeline -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull oist/LuscombeU_Selenoproteins_Pipeline
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [oist/LuscombeU_Selenoproteins_Pipeline releases page](https://github.com/oist/LuscombeU_Selenoproteins_Pipeline/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
