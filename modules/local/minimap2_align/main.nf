// MINIMAP2_ALIGN
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(reads)
    path fasta
    val minimap2_args

    output:
    tuple val(meta), path("${prefix}.mapped.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "MINIMAP2_ALIGN is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mapped.bam
    """
}
