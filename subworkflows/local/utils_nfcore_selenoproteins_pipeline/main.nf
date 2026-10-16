//
// Subworkflow with functionality specific to the oist/LuscombeU_Selenoproteins_Pipeline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { samplesheetToList } from 'plugin/nf-schema'
include { paramsHelp } from 'plugin/nf-schema'
include { completionEmail } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: Display version and exit
    validate_params // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir //  string: The output directory where the results will be saved
    input //  string: Path to input samplesheet
    help // boolean: Display help message and exit
    help_full // boolean: Show the full help message
    show_hidden // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    def default_profiles = [eukaryota: 'eukarya,machinery', bacteria: 'prokarya', archaea: 'prokarya']
    channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fasta, annotation, transcripts ->
            [meta + [profiles: meta.profiles ?: default_profiles[meta.domain]], fasta, annotation ?: [], transcripts ?: []]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email //  string: email address
    email_on_fail //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def tools = params.trnasec_tools.tokenize(',')*.trim()
    def citation_text = [
        "Tools used in the workflow included:",
        "SeqKit (Shen et al. 2016),",
        "Selenoprofiles (Mariotti and Guigó 2010) with BLAST (Altschul et al. 1997), Exonerate (Slater and Birney 2005), GeneWise (Birney et al. 2004) and MAFFT (Katoh and Standley 2013),",
        params.skip_miniprot ? "" : "miniprot (Li 2023),",
        'secmarker' in tools ? "Secmarker (Santesmasses et al. 2017)," : "",
        'trnascanse' in tools ? "tRNAscan-SE (Chan et al. 2021)," : "",
        "Infernal (Nawrocki and Eddy 2013) with Rfam models (Kalvari et al. 2021)",
        "and MultiQC (Ewels et al. 2016).",
    ].findAll { text -> text }.join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def tools = params.trnasec_tools.tokenize(',')*.trim()
    def reference_text = [
        "<li>Shen W, Le S, Li Y, Hu F (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS One 11(10):e0163962. doi: 10.1371/journal.pone.0163962</li>",
        "<li>Mariotti M, Guigó R (2010). Selenoprofiles: profile-based scanning of eukaryotic genome sequences for selenoprotein genes. Bioinformatics 26(21):2656-2663. doi: 10.1093/bioinformatics/btq516</li>",
        "<li>Altschul SF et al. (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res 25(17):3389-3402. doi: 10.1093/nar/25.17.3389</li>",
        "<li>Slater GS, Birney E (2005). Automated generation of heuristics for biological sequence comparison. BMC Bioinformatics 6:31. doi: 10.1186/1471-2105-6-31</li>",
        "<li>Birney E, Clamp M, Durbin R (2004). GeneWise and Genomewise. Genome Res 14(5):988-995. doi: 10.1101/gr.1865504</li>",
        "<li>Katoh K, Standley DM (2013). MAFFT multiple sequence alignment software version 7. Mol Biol Evol 30(4):772-780. doi: 10.1093/molbev/mst010</li>",
        params.skip_miniprot ? "" : "<li>Li H (2023). Protein-to-genome alignment with miniprot. Bioinformatics 39(1):btad014. doi: 10.1093/bioinformatics/btad014</li>",
        'secmarker' in tools ? "<li>Santesmasses D, Mariotti M, Guigó R (2017). Computational identification of the selenocysteine tRNA (tRNASec) in genomes. PLoS Comput Biol 13(2):e1005383. doi: 10.1371/journal.pcbi.1005383</li>" : "",
        'trnascanse' in tools ? "<li>Chan PP, Lin BY, Mak AJ, Lowe TM (2021). tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. Nucleic Acids Res 49(16):9077-9096. doi: 10.1093/nar/gkab688</li>" : "",
        "<li>Nawrocki EP, Eddy SR (2013). Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29(22):2933-2935. doi: 10.1093/bioinformatics/btt509</li>",
        "<li>Kalvari I et al. (2021). Rfam 14: expanded coverage of metagenomic, viral and microRNA families. Nucleic Acids Res 49(D1):D192-D200. doi: 10.1093/nar/gkaa1047</li>",
        "<li>Ewels P, Magnusson M, Lundin S, Käller M (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32(19):3047-3048. doi: 10.1093/bioinformatics/btw354</li>",
    ].findAll { text -> text }.join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
