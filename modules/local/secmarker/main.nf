process SECMARKER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/86/8609cefb2a5e54e83895de216adb1cdfa80d9b3a1d0ddf391d2b03d09d76840a/data'
        : 'community.wave.seqera.io/library/coreutils_findutils_gawk_grep_pruned:03d4a0e8e2a05f27' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.trnasec.gff"), emit: gff
    tuple val(meta), path("*.trnasec.fa"), emit: fasta
    tuple val(meta), path("*.trnasec.ss"), emit: structure
    tuple val("${task.process}"), val('secmarker'), eval("python -c 'from secmarker._version import __version__; print(__version__)'"), topic: versions, emit: versions_secmarker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz"
    def target = is_compressed ? fasta.getBaseName() : fasta
    """
    ${is_compressed ? "gunzip -c ${fasta} >| ${target}" : ''}

    set +e
    python -m secmarker.secmarker_lib \\
        -t ${target} \\
        -o ${prefix} \\
        -cpu ${task.cpus} \\
        ${args}
    status=\$?
    set -e
    # secmarker 1.0.0 fails with a RuntimeError when Infernal reports no hit
    if [ \$status -ne 0 ]; then
        if [ ! -s ${prefix}/infernal.out ] || grep -q -E '^ +\\([0-9]+\\) ' ${prefix}/infernal.out; then
            exit \$status
        fi
    fi

    # secmarker writes no files when it finds no tRNA-Sec
    for ext in gff fa ss; do
        if [ -f ${prefix}/trnasec.\$ext ]; then
            mv ${prefix}/trnasec.\$ext ${prefix}.trnasec.\$ext
        else
            touch ${prefix}.trnasec.\$ext
        fi
    done
    rm -rf ${prefix} ${is_compressed ? target : ''}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.trnasec.gff
    touch ${prefix}.trnasec.fa
    touch ${prefix}.trnasec.ss
    """
}
