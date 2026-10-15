process EXPECTATION_CHECK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Python standard library only: the MultiQC image provides it
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    tuple val(meta), path(sites)
    path expectation

    output:
    tuple val(meta), path("*.expectation.tsv"), emit: families
    tuple val(meta), path("*.expectation_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    expectation_check.py \\
        --sample ${meta.id} \\
        --sites ${sites} \\
        --expectation ${expectation} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.expectation.tsv
    printf 'sample\\texpected\\tobserved\\tmissing\\textra\\n${meta.id}\\t0\\t0\\tnone\\tnone\\n' > ${prefix}.expectation_summary.tsv
    """
}
