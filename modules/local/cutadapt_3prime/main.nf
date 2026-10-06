// CUTADAPT_3PRIME: legacy 3′-only ONT trim (D31, opt-in): cutadapt -m {trim1_min_length} -O {trim1_min_overlap} -a <S3>,
// untrimmed reads kept, no antisense pass. Used only with ont_adapter_mode = 'three_prime_only'.
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process CUTADAPT_3PRIME {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(reads)
    val adapter_3p
    val min_length
    val min_overlap

    output:
    tuple val(meta), path("${prefix}.trim1_3prime.fastq.gz"), emit: reads
    tuple val(meta), path("${prefix}.trim1_3prime.cutadapt.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "CUTADAPT_3PRIME is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.trim1_3prime.fastq.gz
    touch ${prefix}.trim1_3prime.cutadapt.log
    """
}
