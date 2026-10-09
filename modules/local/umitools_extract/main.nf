// UMITOOLS_EXTRACT: --extract-method=string --3prime; UMI appended with '_' as UMICollapse expects. Skipped when umi=false (D29).
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process UMITOOLS_EXTRACT {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(reads)
    val umi_pattern

    output:
    tuple val(meta), path("${prefix}.umi_extracted.fastq.gz"), emit: reads
    tuple val(meta), path("${prefix}.umi_extract.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "UMITOOLS_EXTRACT is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.umi_extracted.fastq.gz
    touch ${prefix}.umi_extract.log
    """
}
