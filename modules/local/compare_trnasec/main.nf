process COMPARE_TRNASEC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Python standard library only: the MultiQC image provides it
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    tuple val(meta), path(secmarker_gff), path(cmsearch_table), path(trnascanse_table)

    output:
    tuple val(meta), path("*.trnasec_loci.tsv"), emit: loci
    tuple val(meta), path("*.trnasec_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def secmarker = secmarker_gff ? "--secmarker ${secmarker_gff}" : ''
    def cmsearch = cmsearch_table ? "--cmsearch ${cmsearch_table}" : ''
    def trnascanse = trnascanse_table ? "--trnascanse ${trnascanse_table}" : ''
    """
    compare_trnasec.py \\
        --sample ${meta.id} \\
        ${secmarker} \\
        ${cmsearch} \\
        ${trnascanse} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.trnasec_loci.tsv
    printf 'sample\\tloci\\tloci_multiple_tools\\tloci_supported\\n${meta.id}\\t0\\t0\\t0\\n' > ${prefix}.trnasec_summary.tsv
    """
}
