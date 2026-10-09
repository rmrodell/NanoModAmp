// CALL_SITES: bin/nma_call.R; treatment or factor analysis (§6.3, §6.4, §5.5). Seeded; workers = task.cpus.
// WP0 stub: the real command is implemented in WP4; the stub block lets `-stub` runs test the wiring.
process CALL_SITES {
    tag "$name"
    label 'process_medium'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP4

    input:
    tuple val(name), path(merged), path(analyses)

    output:
    tuple val(name), path("${name}"), emit: dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "CALL_SITES is not implemented yet (WP4); use -stub" >&2
    exit 1
    """

    stub:
    """
    mkdir -p ${name}/data_summary ${name}/data_raw ${name}/plots
    touch ${name}/${name}_log.txt ${name}/${name}_resolved_config.yaml
    """
}
