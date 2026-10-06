// SAMTOOLS_FAIDX
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process SAMTOOLS_FAIDX {
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    path fasta

    output:
    path("${fasta}.fai"), emit: fai

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "SAMTOOLS_FAIDX is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    """
    touch ${fasta}.fai
    """
}
