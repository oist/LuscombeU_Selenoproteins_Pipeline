process SEC_TRANSCRIPT_SUPPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Python standard library only: the MultiQC image provides it
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    tuple val(meta), path(sites), path(paf)

    output:
    tuple val(meta), path("*.sec_transcripts.tsv"), emit: support
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sec_transcript_support.py \\
        --sites ${sites} \\
        --paf ${paf} \\
        --output ${prefix}.sec_transcripts.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'contig\\tstart\\tend\\tstrand\\tfamily\\tsupport\\ttranscripts\\tmax_3prime_extent\\n' > ${prefix}.sec_transcripts.tsv
    """
}
