// SEQTK_RC
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process SEQTK_RC {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.antisense_rc.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "SEQTK_RC is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.antisense_rc.fastq.gz
    """
}
