process COMPARE_GENES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Python standard library only: the MultiQC image provides it
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    tuple val(meta), path(selenoprofiles_gtf), path(miniprot_sites), path(secis_tables), val(profiles)
    path profile_sets

    output:
    tuple val(meta), path("*.sec_sites.tsv"), emit: sites
    tuple val(meta), path("*.sec_codons.gff3"), emit: layer
    tuple val(meta), path("*.genes_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def miniprot = miniprot_sites ? "--miniprot-sites ${miniprot_sites}" : ''
    def secis = secis_tables ? "--secis ${[secis_tables].flatten().join(' ')}" : ''
    """
    compare_genes.py \\
        --sample ${meta.id} \\
        --selenoprofiles-gtf ${selenoprofiles_gtf} \\
        ${miniprot} \\
        ${secis} \\
        --profiles ${profiles} \\
        --profile-sets ${profile_sets} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sec_sites.tsv
    touch ${prefix}.sec_codons.gff3
    printf 'sample\\tpredictions\\tselenocysteine\\n${meta.id}\\t0\\t0\\n' > ${prefix}.genes_summary.tsv
    """
}
