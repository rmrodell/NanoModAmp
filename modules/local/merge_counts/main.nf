// MERGE_COUNTS: bin/nma_merge.R; sample_id + sample-sheet metadata + count columns, joined by name (§5.3, R-06).
// WP0 stub: the real command is implemented in WP3; the stub block lets `-stub` runs test the wiring.
process MERGE_COUNTS {
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP3

    input:
    path counts
    path samplesheet

    output:
    path("counts_merged.tsv"), emit: merged

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "MERGE_COUNTS is not implemented yet (WP3); use -stub" >&2
    exit 1
    """

    stub:
    """
    printf "sample_id\tchr\tpos\tgene\ttotalReads\tA.count\tC.count\tG.count\tT.count\tDeletion.count\tInsertion.count\tref\tkmer\tstrand\tdelrate\n" > counts_merged.tsv
    """
}
