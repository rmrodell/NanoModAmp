// CUTADAPT_LINKED: Linked-adapter pass, both adapters required (D4). Included as CUTADAPT_SENSE, CUTADAPT_ANTISENSE and CUTADAPT_POOL.
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process CUTADAPT_LINKED {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(reads)
    val adapter_5p
    val adapter_3p
    val min_length
    val min_overlap
    val step

    output:
    tuple val(meta), path("${prefix}.${step}.fastq.gz"), emit: reads
    tuple val(meta), path("${prefix}.${step}.cutadapt.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "CUTADAPT_LINKED is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.${step}.fastq.gz
    touch ${prefix}.${step}.cutadapt.log
    """
}
