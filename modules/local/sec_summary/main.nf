process SEC_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Python standard library only: the MultiQC image provides it
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    path samples
    path gene_summaries, stageAs: 'genes/*'
    path trnasec_summaries, stageAs: 'trnasec/*'
    path assess_tables, stageAs: 'assess/*'
    path sec_annotation_tables, stageAs: 'secannot/*'
    path candidate_tables, stageAs: 'candidates/*'
    path expectation_tables, stageAs: 'expectation/*'
    path transcript_tables, stageAs: 'transcripts/*'

    output:
    path "*_mqc.tsv", emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def assess = assess_tables ? "--assess ${[assess_tables].flatten().join(' ')}" : ''
    def secannot = sec_annotation_tables ? "--sec-annotation ${[sec_annotation_tables].flatten().join(' ')}" : ''
    def candidates = candidate_tables ? "--candidates ${[candidate_tables].flatten().join(' ')}" : ''
    def expectation = expectation_tables ? "--expectation ${[expectation_tables].flatten().join(' ')}" : ''
    def transcripts = transcript_tables ? "--transcripts ${[transcript_tables].flatten().join(' ')}" : ''
    """
    sec_summary.py \\
        --samples ${samples} \\
        --genes ${[gene_summaries].flatten().join(' ')} \\
        --trnasec ${[trnasec_summaries].flatten().join(' ')} \\
        ${assess} \\
        ${secannot} \\
        ${candidates} \\
        ${expectation} \\
        ${transcripts}
    """

    stub:
    """
    touch selenoprotein_genes_mqc.tsv sec_site_agreement_mqc.tsv trnasec_mqc.tsv sec_trait_mqc.tsv
    """
}
