// MERGE_ORIENT: Sense reads followed by RC antisense reads; counts read IDs found in both (§6.1 step 3).
// WP0 stub: the real command is implemented in WP2; the stub block lets `-stub` runs test the wiring.
process MERGE_ORIENT {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP2

    input:
    tuple val(meta), path(sense), path(antisense_rc)

    output:
    tuple val(meta), path("${prefix}.merged.fastq.gz"), emit: reads
    tuple val(meta), path("${prefix}.both_orientations.tsv"), emit: both_orientations

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "MERGE_ORIENT is not implemented yet (WP2); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.merged.fastq.gz
    printf "sample_id\tn_reads_in_both_orientations\n${meta.id}\t0\n" > ${prefix}.both_orientations.tsv
    """
}
