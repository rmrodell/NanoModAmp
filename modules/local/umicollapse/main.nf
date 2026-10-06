// UMICOLLAPSE: Skipped when umi=false (D29).
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process UMICOLLAPSE {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(bam), path(bai)
    val umicollapse_args

    output:
    tuple val(meta), path("${prefix}.dedup.bam"), emit: bam
    tuple val(meta), path("${prefix}.umicollapse.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "UMICOLLAPSE is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dedup.bam ${prefix}.umicollapse.log
    """
}
