// CALL_SITES: bin/nma_call.R; treatment or factor analysis (§6.3, §6.4, §5.5). Seeded; workers = task.cpus.
// Implemented in WP4 (R/nanomodamp: treatment_test / factor_test); the stub block keeps `-stub` runs working.
process CALL_SITES {
    tag "$name"
    label 'process_medium'

    // placeholder until the nanomodamp-r container (R 4.3 + nanomodamp + blme/car/lme4/furrr/ggplot2) is published
    container 'docker.io/library/ubuntu:22.04'

    input:
    tuple val(name), path(merged), path(analyses)

    output:
    tuple val(name), path("${name}"), emit: dir
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def all_sites = params.plot_all_sites ? '--plot_all_sites' : ''
    """
    nma_call.R \\
        --counts ${merged} \\
        --analyses ${analyses} \\
        --name ${name} \\
        --outdir . \\
        --cpus ${task.cpus} \\
        ${all_sites}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        nanomodamp: \$(Rscript -e 'cat(as.character(packageVersion("nanomodamp")))')
        blme: \$(Rscript -e 'cat(as.character(packageVersion("blme")))')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${name}/data_summary ${name}/data_raw ${name}/plots
    touch ${name}/${name}_log.txt ${name}/${name}_resolved_config.yaml
    touch versions.yml
    """
}
