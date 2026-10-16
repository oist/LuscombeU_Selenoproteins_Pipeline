process MACHINERY_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Python standard library only: the MultiQC image provides it
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    path samples
    path gtfs, stageAs: 'gtf/*'
    path trnasec_summaries, stageAs: 'trnasec/*'

    output:
    path "machinery_mqc.tsv", emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    machinery_summary.py \\
        --samples ${samples} \\
        --gtf ${[gtfs].flatten().join(' ')} \\
        --trnasec ${[trnasec_summaries].flatten().join(' ')}
    """

    stub:
    """
    touch machinery_mqc.tsv
    """
}
