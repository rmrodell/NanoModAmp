// CAT_FASTQ: Single .fastq.gz, or a directory whose *.fastq.gz files are concatenated in lexical order (§5.1).
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process CAT_FASTQ {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}.raw.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "CAT_FASTQ is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.raw.fastq.gz
    """
}
