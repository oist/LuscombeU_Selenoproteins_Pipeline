# oist/LuscombeU_Selenoproteins_Pipeline: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

- [Assembly statistics](#assembly-statistics): seqkit stats.
- [Selenoprotein genes](#selenoprotein-genes): Selenoprofiles predictions, miniprot cross-check, SECIS evidence and annotation assessment.
- [tRNA-Sec](#trna-sec): Secmarker, cmsearch and tRNAscan-SE predictions merged into loci.
- [Summary](#summary): tables used by the MultiQC report.
- [MultiQC](#multiqc) and [pipeline information](#pipeline-information).

### Assembly statistics

- `seqkit/`: `<sample>.tsv`, sequence statistics of each genome (`seqkit stats --all`).

### Selenoprotein genes

- `selenoprofiles/`
  - `<sample>.selenoprofiles.gtf`: predictions; the Sec codons are features of type `Selenocysteine`. The last field of `gene_id` is the Selenoprofiles label (`selenocysteine`, `cysteine`, `pseudo` and others).
  - `<sample>.selenoprofiles.faa`, `<sample>.selenoprofiles.cds.fa`: predicted proteins and coding sequences.
  - `<sample>.three_prime.fa`: 6000 nt downstream of each prediction.
  - `<sample>_alignments/`: alignment of the predictions of each family with its profile.
  - `<sample>_results.tar.gz`: complete Selenoprofiles result folder.
  - `<sample>.transcripts.tsv`, `<sample>.genes.tsv`: annotation status from `selenoprofiles assess`, when an annotation is given.
- `miniprot/`: `<sample>.gff`, alignments of the Sec-containing profile sequences (`summary/selenoprotein_queries.faa`).
- `sec/`: `<sample>.miniprot_sec.tsv`, TGA codons aligned to a query selenocysteine (`aligned`) or found right after an alignment that stops before it (`next_codon`).
- `infernal/`: `<sample>.secis.tbl.gz`, cmsearch hits of the SECIS model in the 3' flanks.
- `compare/`
  - `<sample>.sec_sites.tsv`: union of Selenoprofiles and miniprot Sec sites, with the support of each and the best SECIS score of the prediction.
  - `<sample>.genes_summary.tsv`: counts per label, Sec site agreement (Jaccard index) and number of selenoproteins with a SECIS hit.

Selenoprofiles assigns each prediction to the family of the profile that found it, by sequence similarity to the
profile and, for some families, to tagged proteins of a reference database. It does not test orthology: a prediction
may be a paralog of the named gene, or a member of a related family. This matters most for genes with close
relatives, such as PSTK (related to the kinase KTI12), eEFSec (related to other translation elongation factors) or
SECp43 (related to other RNA-binding proteins), and for lineage-specific profiles given with
`--selenoprofiles_overrides`. Orthology must be established outside the pipeline, for example by phylogeny with the
related families.

### Sec codon layer and expectation

- `compare/<sample>.sec_codons.gff3`: Sec codons as GFF3 features for annotation pipelines.
- `expectation/`, with `--expectation`: `<sample>.expectation.tsv`, expected and observed Sec sites per family; `<sample>.expectation_summary.tsv`, missing and extra families.

### Annotation checks

- `annotation/`: `<sample>.sec_annotation.tsv`, for each predicted Sec codon, the class in the given annotation (Ticó et al. 2026), the residue it gets, the transcript and protein concerned, and a note when the codon lies in another transcript's UTR; `<sample>.sec_annotation_transcripts.tsv`, the class for every overlapping transcript; `<sample>.sec_annotation_summary.tsv`, counts at gene and transcript level.
- `readthrough/`, with `--readthrough_screen`: `<sample>.readthrough.tsv`, TGA events per transcript (`start_extension`, `upstream_tga`, `stop_readthrough`, `internal_tga`, `sec_declared`); `<sample>.3utr.fa`, the 3′ UTRs searched for SECIS elements; `<sample>.candidates.tsv`, ranked candidates with the SECIS bit score, E-value and position.
- `infernal/`: `<sample>.screen.tbl.gz`, SECIS hits in the 3′ UTRs.

### Transcript support

- MultiQC section "Transcript support", per genome: Sec codons, Sec codons inside an aligned exon of at least one transcript, supporting alignments, median 3′ extension.
- `transcripts/`, when the samplesheet gives transcripts: `<sample>.sec_transcripts.tsv`, for each Sec codon of `compare/<sample>.sec_sites.tsv`, the number of transcripts whose primary spliced alignment (minimap2 `-x splice`) contains the whole codon inside one aligned exon (`transcripts`), and the longest extension of these transcripts 3′ of the codon on the strand of the gene (`max_3prime_extent`, in nucleotides), which shows whether they reach the SECIS region.

### tRNA-Sec

- `secmarker/`: `<sample>.trnasec.gff`, `.fa`, `.ss`.
- `infernal/`: `<sample>.trnasec.tbl.gz`, cmsearch hits of the tRNA-Sec model.
- `trnascanse/`: tRNAscan-SE output, when `--trnasec_tools` includes `trnascanse`.
- `compare/`: `<sample>.trnasec_loci.tsv`, merged loci with the score of each tool, and `<sample>.trnasec_summary.tsv`.

### Summary

- `summary/`: MultiQC custom-content tables (`*_mqc.tsv`) and the Sec trait call per genome (`sec_trait_mqc.tsv`).
  - `machinery_mqc.tsv`: per genome, the tRNA-Sec loci supported by two or more tools and the Selenoprofiles status of
    the machinery genes (SEPHS2, PSTK, SecS, eEFSec, SBP2) and of the associated factor SECp43, which is not required for
    selenoprotein synthesis in all cells. A gene found here is a homolog of the profile, not
    a proven ortholog (see above); a gene broken by the assembly may be missed.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
