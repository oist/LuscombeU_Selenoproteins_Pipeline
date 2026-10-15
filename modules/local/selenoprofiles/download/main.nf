process SELENOPROFILES_DOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/58/58fad6038205d880917370e4e521d0d0b2e07c6bbb06d858facbf8892321fd8e/data'
        : 'community.wave.seqera.io/library/blast-legacy_exonerate_mafft_obonet_pruned:6c6cfa335a4446c7' }"

    output:
    path "selenoprofiles_data", emit: db
    tuple val("${task.process}"), val('selenoprofiles'), eval("python -c 'from selenoprofiles4._version import __version__; print(__version__)'"), topic: versions, emit: versions_selenoprofiles

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export HOME=\$PWD
    selenoprofiles -setup > /dev/null
    # -download asks for one confirmation
    echo | selenoprofiles -download -selenoprofiles_data_dir selenoprofiles_data ${args}
    rm -f selenoprofiles_data/*.tar.gz selenoprofiles_data/*.fasta.gz
    """

    stub:
    """
    mkdir -p selenoprofiles_data/selenoprotein_profiles selenoprofiles_data/gene_ontology
    touch selenoprofiles_data/uniref50_seleno_sim.fasta
    """
}
