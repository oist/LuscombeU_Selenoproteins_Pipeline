process SELENOPROFILES_ASSESS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/58/58fad6038205d880917370e4e521d0d0b2e07c6bbb06d858facbf8892321fd8e/data'
        : 'community.wave.seqera.io/library/blast-legacy_exonerate_mafft_obonet_pruned:6c6cfa335a4446c7' }"

    input:
    tuple val(meta), path(selenoprofiles_gtf), path(annotation), path(fasta)

    output:
    tuple val(meta), path("*.transcripts.tsv"), emit: transcripts
    tuple val(meta), path("*.genes.tsv"), emit: genes
    tuple val("${task.process}"), val('selenoprofiles'), eval("python -c 'from selenoprofiles4._version import __version__; print(__version__)'"), topic: versions, emit: versions_selenoprofiles

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_gz = fasta.getExtension() == "gz"
    def target = fasta_gz ? fasta.getBaseName() : fasta
    def annotation_gz = annotation.getExtension() == "gz"
    def gff = annotation_gz ? annotation.getBaseName() : annotation
    """
    export HOME=\$PWD
    selenoprofiles -setup > /dev/null
    ${fasta_gz ? "gunzip -c ${fasta} > ${target}" : ''}
    ${annotation_gz ? "gunzip -c ${annotation} > ${gff}" : ''}

    # assess fails on an empty prediction set, which is a valid result
    if grep -q -P '\\tSelenocysteine\\t' ${selenoprofiles_gtf}; then
        selenoprofiles assess \\
            -s ${selenoprofiles_gtf} \\
            -e ${gff} \\
            -f ${target} \\
            -o ${prefix}.transcripts.tsv \\
            -agg ${prefix}.genes.tsv \\
            ${args}
    else
        printf 'transcript_id\\ttranscript_id_ens\\tType_annotation\\n' > ${prefix}.transcripts.tsv
        cp ${prefix}.transcripts.tsv ${prefix}.genes.tsv
    fi
    rm -f ${fasta_gz ? target : ''} ${annotation_gz ? gff : ''} *.fai
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'transcript_id\\ttranscript_id_ens\\tType_annotation\\n' > ${prefix}.transcripts.tsv
    cp ${prefix}.transcripts.tsv ${prefix}.genes.tsv
    """
}
