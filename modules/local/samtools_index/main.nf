// SAMTOOLS_INDEX
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam}.bai"), emit: bam_bai

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "SAMTOOLS_INDEX is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${bam}.bai
    """
}
