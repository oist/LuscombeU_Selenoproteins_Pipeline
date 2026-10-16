/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SEQKIT_STATS } from '../modules/nf-core/seqkit/stats/main'
include { UNTAR as UNTAR_SELENOPROFILES_DATA } from '../modules/nf-core/untar/main'
include { SELENOPROFILES_DOWNLOAD } from '../modules/local/selenoprofiles/download/main'
include { SELENOPROFILES_SELENOPROFILES } from '../modules/local/selenoprofiles/selenoprofiles/main'
include { SELENOPROFILES_ASSESS } from '../modules/local/selenoprofiles/assess/main'
include { PROFILE_QUERIES } from '../modules/local/profile_queries/main'
include { MINIPROT_INDEX } from '../modules/nf-core/miniprot/index/main'
include { MINIPROT_ALIGN } from '../modules/nf-core/miniprot/align/main'
include { SEC_CHECK } from '../modules/local/sec_check/main'
include { INFERNAL_CMSEARCH as INFERNAL_CMSEARCH_SECIS } from '../modules/nf-core/infernal/cmsearch/main'
include { INFERNAL_CMSEARCH as INFERNAL_CMSEARCH_TRNASEC } from '../modules/nf-core/infernal/cmsearch/main'
include { SECMARKER } from '../modules/local/secmarker/main'
include { TRNASCANSE } from '../modules/nf-core/trnascanse/main'
include { COMPARE_GENES } from '../modules/local/compare_genes/main'
include { COMPARE_TRNASEC } from '../modules/local/compare_trnasec/main'
include { SEC_SUMMARY } from '../modules/local/sec_summary/main'
include { SEC_ANNOTATION_CHECK } from '../modules/local/sec_annotation_check/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TRANSCRIPTS } from '../modules/nf-core/minimap2/align/main'
include { SEC_TRANSCRIPT_SUPPORT } from '../modules/local/sec_transcript_support/main'
include { READTHROUGH_SCREEN } from '../modules/local/readthrough_screen/main'
include { RANK_READTHROUGH } from '../modules/local/rank_readthrough/main'
include { EXPECTATION_CHECK } from '../modules/local/expectation_check/main'
include { INFERNAL_CMSEARCH as INFERNAL_CMSEARCH_SCREEN } from '../modules/nf-core/infernal/cmsearch/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_selenoproteins_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SELENOPROTEINS {
    take:
    ch_samplesheet // channel: [ val(meta), path(fasta), path(annotation) or [], path(transcripts) or [] ]
    selenoprofiles_data
    trnasec_tools
    skip_miniprot
    skip_secis
    secis_cm
    trnasec_cm
    readthrough_screen
    screen_cm
    expectation
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def tools = trnasec_tools.tokenize(',')*.trim()

    def ch_genomes = ch_samplesheet.map { meta, fasta, _annotation, _transcripts -> [meta, fasta] }

    SEQKIT_STATS(ch_genomes)
    ch_multiqc_files = ch_multiqc_files.mix(SEQKIT_STATS.out.stats.map { _meta, stats -> stats })

    //
    // Selenoprotein genes
    //
    def ch_selenoprofiles_data
    if (!selenoprofiles_data) {
        SELENOPROFILES_DOWNLOAD()
        ch_selenoprofiles_data = SELENOPROFILES_DOWNLOAD.out.db
    }
    else if (file(selenoprofiles_data, checkIfExists: true).isDirectory()) {
        ch_selenoprofiles_data = channel.value(file(selenoprofiles_data, checkIfExists: true))
    }
    else {
        UNTAR_SELENOPROFILES_DATA(channel.value([[id: 'selenoprofiles_data'], file(selenoprofiles_data, checkIfExists: true)]))
        ch_selenoprofiles_data = UNTAR_SELENOPROFILES_DATA.out.untar.map { _meta, dir -> dir }
    }

    SELENOPROFILES_SELENOPROFILES(
        ch_genomes.map { meta, fasta -> [meta, fasta, meta.profiles] },
        ch_selenoprofiles_data,
    )

    def ch_miniprot_sites = ch_genomes.map { meta, _fasta -> [meta, []] }
    if (!skip_miniprot) {
        PROFILE_QUERIES(ch_selenoprofiles_data)
        MINIPROT_INDEX(ch_genomes)
        def ch_align = ch_genomes
            .join(MINIPROT_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)
            .combine(PROFILE_QUERIES.out.queries)
            .multiMap { meta, _fasta, index, queries ->
                pep: [meta, queries]
                ref: [meta, index]
            }
        MINIPROT_ALIGN(ch_align.pep, ch_align.ref)
        SEC_CHECK(
            MINIPROT_ALIGN.out.gff.join(ch_genomes, failOnMismatch: true, failOnDuplicate: true),
            PROFILE_QUERIES.out.queries,
        )
        ch_miniprot_sites = SEC_CHECK.out.sites
    }

    def ch_secis = ch_genomes.map { meta, _fasta -> [meta, []] }
    if (!skip_secis) {
        def ch_flanks = SELENOPROFILES_SELENOPROFILES.out.three_prime.branch { meta, flanks ->
            search: meta.domain == 'eukaryota' && flanks.size() > 0
            skip: true
        }
        INFERNAL_CMSEARCH_SECIS(
            ch_flanks.search.map { meta, flanks -> [meta, file(secis_cm, checkIfExists: true), flanks] },
            false,
            true,
        )
        ch_secis = INFERNAL_CMSEARCH_SECIS.out.target_summary.mix(ch_flanks.skip.map { meta, _flanks -> [meta, []] })
    }

    COMPARE_GENES(
        SELENOPROFILES_SELENOPROFILES.out.gtf.join(ch_miniprot_sites, failOnMismatch: true, failOnDuplicate: true).join(ch_secis, failOnMismatch: true, failOnDuplicate: true).map { meta, gtf, sites, secis -> [meta, gtf, sites, secis, meta.profiles] },
        file("${projectDir}/assets/profile_sets.tsv", checkIfExists: true),
    )

    //
    // Transcript support: spliced alignment of the samplesheet transcripts, Sec codons inside aligned exons
    //
    def ch_transcripts = ch_samplesheet
        .filter { _meta, _fasta, _annotation, transcripts -> transcripts }
        .multiMap { meta, fasta, _annotation, transcripts ->
            reads: [meta, transcripts]
            reference: [meta, fasta]
        }
    MINIMAP2_ALIGN_TRANSCRIPTS(ch_transcripts.reads, ch_transcripts.reference, false, '', true, false)
    SEC_TRANSCRIPT_SUPPORT(COMPARE_GENES.out.sites.join(MINIMAP2_ALIGN_TRANSCRIPTS.out.paf, failOnDuplicate: true))

    def ch_annotated = SELENOPROFILES_SELENOPROFILES.out.gtf
        .join(ch_samplesheet, failOnMismatch: true, failOnDuplicate: true)
        .filter { _meta, _gtf, _fasta, annotation, _transcripts -> annotation }
        .map { meta, gtf, fasta, annotation, _transcripts -> [meta, gtf, annotation, fasta] }
    SELENOPROFILES_ASSESS(ch_annotated)
    SEC_ANNOTATION_CHECK(ch_annotated.map { meta, gtf, annotation, _fasta -> [meta, gtf, annotation] })

    //
    // Transcript screen: TGA codons at the edges of annotated CDSs, ranked by SECIS in the 3' UTR
    //
    def ch_candidates = channel.empty()
    if (readthrough_screen) {
        READTHROUGH_SCREEN(ch_annotated.map { meta, _gtf, annotation, fasta -> [meta, annotation, fasta] })
        INFERNAL_CMSEARCH_SCREEN(
            READTHROUGH_SCREEN.out.utr.map { meta, utr -> [meta, file(screen_cm ?: secis_cm, checkIfExists: true), utr] },
            false,
            true,
        )
        RANK_READTHROUGH(READTHROUGH_SCREEN.out.events.join(INFERNAL_CMSEARCH_SCREEN.out.target_summary, failOnMismatch: true, failOnDuplicate: true))
        ch_candidates = RANK_READTHROUGH.out.candidates
    }

    //
    // tRNA-Sec
    //
    def ch_empty = ch_genomes.map { meta, _fasta -> [meta, []] }
    def ch_secmarker = ch_empty
    def ch_cmsearch = ch_empty
    def ch_trnascanse = ch_empty
    if ('secmarker' in tools) {
        SECMARKER(ch_genomes)
        ch_secmarker = SECMARKER.out.gff
    }
    if ('cmsearch' in tools) {
        INFERNAL_CMSEARCH_TRNASEC(
            ch_genomes.map { meta, fasta -> [meta, file(trnasec_cm, checkIfExists: true), fasta] },
            false,
            true,
        )
        ch_cmsearch = INFERNAL_CMSEARCH_TRNASEC.out.target_summary
    }
    if ('trnascanse' in tools) {
        TRNASCANSE(ch_genomes, false, false, false, false)
        ch_trnascanse = TRNASCANSE.out.tsv
    }
    COMPARE_TRNASEC(
        ch_secmarker.join(ch_cmsearch, failOnMismatch: true, failOnDuplicate: true).join(ch_trnascanse, failOnMismatch: true, failOnDuplicate: true)
    )

    //
    // Summary tables for MultiQC, in samplesheet order
    //
    def ch_sample_order = ch_samplesheet
        .map { meta, _fasta, _annotation, _transcripts -> meta.id }
        .collectFile(name: 'samples.txt', newLine: true, sort: false)
    //
    // Expected selenoproteome (families), for lineages with an expectation table
    //
    def ch_expectation = channel.empty()
    if (expectation) {
        EXPECTATION_CHECK(COMPARE_GENES.out.sites, file(expectation, checkIfExists: true))
        ch_expectation = EXPECTATION_CHECK.out.summary
    }

    SEC_SUMMARY(
        ch_sample_order,
        COMPARE_GENES.out.summary.map { _meta, summary -> summary }.collect(),
        COMPARE_TRNASEC.out.summary.map { _meta, summary -> summary }.collect(),
        SELENOPROFILES_ASSESS.out.genes.map { _meta, genes -> genes }.collect().ifEmpty([]),
        SEC_ANNOTATION_CHECK.out.summary.map { _meta, summary -> summary }.collect().ifEmpty([]),
        ch_candidates.map { _meta, table -> table }.collect().ifEmpty([]),
        ch_expectation.map { _meta, table -> table }.collect().ifEmpty([]),
        SEC_TRANSCRIPT_SUPPORT.out.support.map { _meta, table -> table }.collect().ifEmpty([]),
    )
    ch_multiqc_files = ch_multiqc_files.mix(SEC_SUMMARY.out.multiqc.flatten())

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'selenoproteins_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'selenoproteins'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions = ch_versions // channel: [ path(versions.yml) ]
}
