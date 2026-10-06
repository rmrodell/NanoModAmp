// COUNT: bin/nma_count.R; nanomodamp::count_sites() (§6.2, §5.3).
// WP0 stub: the real command is implemented in WP3; the stub block lets `-stub` runs test the wiring.
process COUNT {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP3

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path bed

    output:
    tuple val(meta), path("${prefix}.counts.tsv"), emit: counts
    tuple val(meta), path("${prefix}.failed_regions.tsv"), emit: failed_regions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "COUNT is not implemented yet (WP3); use -stub" >&2
    exit 1
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "chr\tpos\tgene\ttotalReads\tA.count\tC.count\tG.count\tT.count\tDeletion.count\tInsertion.count\tref\tkmer\tstrand\tdelrate\n" > ${prefix}.counts.tsv
    printf "chr\tstart\tend\tgene\terror\n" > ${prefix}.failed_regions.tsv
    """
}
