process SELENOPROFILES_OVERRIDE {
    tag "selenoprofiles_data"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fcea51d99062f03670c049c64abad3e89eb20a434b9f66dc69375a4a519eeac7/data'
        : 'community.wave.seqera.io/library/coreutils_findutils_gawk_grep_pruned:e68344212d05a99d' }"

    input:
    path standard, stageAs: 'standard_data'
    path overrides, stageAs: 'overrides'

    output:
    path 'selenoprofiles_data', emit: db

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # A copy of the Selenoprofiles data made of symbolic links, in which each profile of the overrides
    # (<name>.fa with its .fa.config and optional .fa.profile_data) replaces the standard profile of the same name.
    mkdir selenoprofiles_data
    cp -rs "\$(readlink -f standard_data)/." selenoprofiles_data/
    if [ -d selenoprofiles_data/selenoprotein_profiles ]; then
        target=selenoprofiles_data/selenoprotein_profiles
    else
        target=selenoprofiles_data
    fi
    n=0
    for f in "\$(readlink -f overrides)"/*.fa; do
        [ -e "\$f" ] || continue
        name=\$(basename "\$f" .fa)
        rm -f "\$target/\$name.fa" "\$target/\$name.fa.config" "\$target/\$name.fa.profile_data"
        cp "\$f" "\$target/"
        for ext in config profile_data; do
            if [ -e "\$f.\$ext" ]; then cp "\$f.\$ext" "\$target/"; fi
        done
        n=\$((n + 1))
    done
    if [ "\$n" -eq 0 ]; then echo "ERROR: no profile (*.fa) in the overrides" >&2; exit 1; fi
    echo "\$n profiles replaced" >&2
    """

    stub:
    """
    mkdir selenoprofiles_data
    """
}
