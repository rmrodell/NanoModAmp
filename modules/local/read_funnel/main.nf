// READ_FUNNEL: One row per executed step (§5.2); steps skipped by the run's parameters are omitted.
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process READ_FUNNEL {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), val(steps), path(files, stageAs: 'stage??/*')

    output:
    tuple val(meta), path("${prefix}.read_funnel.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "READ_FUNNEL is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "sample_id\tstep\tfile\trecords\tbytes\n" > ${prefix}.read_funnel.tsv
    for s in ${steps.join(' ')}; do printf "${meta.id}\t\$s\tNA\t0\t0\n" >> ${prefix}.read_funnel.tsv; done
    """
}
