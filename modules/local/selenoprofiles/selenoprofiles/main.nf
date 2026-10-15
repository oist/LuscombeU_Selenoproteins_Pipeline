process SELENOPROFILES_SELENOPROFILES {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/58/58fad6038205d880917370e4e521d0d0b2e07c6bbb06d858facbf8892321fd8e/data'
        : 'community.wave.seqera.io/library/blast-legacy_exonerate_mafft_obonet_pruned:6c6cfa335a4446c7' }"

    input:
    tuple val(meta), path(fasta), val(profiles)
    path selenoprofiles_data

    output:
    tuple val(meta), path("*.selenoprofiles.gtf"), emit: gtf
    tuple val(meta), path("*.selenoprofiles.faa"), emit: faa
    tuple val(meta), path("*.selenoprofiles.cds.fa"), emit: cds
    tuple val(meta), path("*.three_prime.fa"), emit: three_prime
    tuple val(meta), path("*_alignments"), emit: alignments
    tuple val(meta), path("*_results.tar.gz"), emit: results
    tuple val("${task.process}"), val('selenoprofiles'), eval("python -c 'from selenoprofiles4._version import __version__; print(__version__)'"), topic: versions, emit: versions_selenoprofiles

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz"
    def target = is_compressed ? fasta.getBaseName() : fasta
    """
    # selenoprofiles reads its configuration from \$HOME and writes BLAST indices next to the target
    export HOME=\$PWD
    selenoprofiles -setup > /dev/null
    ${is_compressed ? "gunzip -c ${fasta} > ${target}" : ''}

    if [ -d ${selenoprofiles_data}/selenoprotein_profiles ]; then
        data_option="-selenoprofiles_data_dir ${selenoprofiles_data}"
    else
        data_option="-profiles_folder ${selenoprofiles_data}"
    fi

    mkdir -p tmp ${prefix}_alignments
    selenoprofiles \\
        -o ${prefix}_results \\
        -t ${target} \\
        -s ${prefix} \\
        -p ${profiles} \\
        \$data_option \\
        -temp tmp \\
        -ncpus ${task.cpus} \\
        -output_gtf_file ${prefix}.selenoprofiles.gtf \\
        -output_fasta_file ${prefix}.selenoprofiles.faa \\
        -output_cds_file ${prefix}.selenoprofiles.cds.fa \\
        -output_three_prime_file ${prefix}.three_prime.fa \\
        ${args}

    # the container has no find(1): copy the alignments with a bash glob
    shopt -s globstar nullglob
    for ali in ${prefix}_results/**/*.ali; do cp "\$ali" ${prefix}_alignments/; done
    tar -czf ${prefix}_results.tar.gz ${prefix}_results
    rm -rf tmp ${prefix}_results ${target}.n?? formatdb.log ${is_compressed ? target : ''}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_alignments
    touch ${prefix}.selenoprofiles.gtf
    touch ${prefix}.selenoprofiles.faa
    touch ${prefix}.selenoprofiles.cds.fa
    touch ${prefix}.three_prime.fa
    tar -cf - -T /dev/null | gzip -n > ${prefix}_results.tar.gz
    """
}
