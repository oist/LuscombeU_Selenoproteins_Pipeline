process SELENOPROFILES_ASSESS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15df41af9baed69a73306cce8a4260d7d028930475b25af8c280470a7360df61/data'
        : 'community.wave.seqera.io/library/coreutils_findutils_gawk_grep_pruned:e66263d7188b105f' }"

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
    ${fasta_gz ? "if [ -L ${target} ]; then echo 'ERROR: ${target} is a staged input; not overwritten' >&2; exit 1; fi; rm -f ${target}; gunzip -c ${fasta} > ${target}" : ''}
    ${annotation_gz ? "if [ -L ${gff} ]; then echo 'ERROR: ${gff} is a staged input; not overwritten' >&2; exit 1; fi; rm -f ${gff}; gunzip -c ${annotation} > ${gff}" : ''}

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
        rm -f ${prefix}.transcripts.tsv ${prefix}.genes.tsv
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
