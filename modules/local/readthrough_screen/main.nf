process READTHROUGH_SCREEN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Python standard library only: the MultiQC image provides it
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    tuple val(meta), path(annotation), path(fasta)

    output:
    tuple val(meta), path("*.readthrough.tsv"), emit: events
    tuple val(meta), path("*.3utr.fa"), emit: utr
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    readthrough_screen.py \\
        --gff ${annotation} \\
        --genome ${fasta} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'transcript\\tgene\\tprotein\\tproduct\\tcontig\\tstrand\\tevent\\textra_codons\\ttga_codons\\ttga_positions\\tutr_length\\n' > ${prefix}.readthrough.tsv
    touch ${prefix}.3utr.fa
    """
}
