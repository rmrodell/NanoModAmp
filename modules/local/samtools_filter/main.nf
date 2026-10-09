// SAMTOOLS_FILTER: -F {sam_exclude_flags} (funnel step 'primary'), then -q {min_mapq} (step 'mapq_filtered').
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process SAMTOOLS_FILTER {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(bam)
    val exclude_flags
    val min_mapq

    output:
    tuple val(meta), path("${prefix}.primary.bam"), emit: primary
    tuple val(meta), path("${prefix}.filtered.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "SAMTOOLS_FILTER is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.primary.bam ${prefix}.filtered.bam
    """
}
