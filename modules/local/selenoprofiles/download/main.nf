process SELENOPROFILES_DOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fcea51d99062f03670c049c64abad3e89eb20a434b9f66dc69375a4a519eeac7/data'
        : 'community.wave.seqera.io/library/coreutils_findutils_gawk_grep_pruned:e68344212d05a99d' }"

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
