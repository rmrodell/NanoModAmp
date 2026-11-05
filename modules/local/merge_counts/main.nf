// MERGE_COUNTS: bin/nma_merge.R; sample_id + sample-sheet metadata + count columns, joined by name (§5.3, R-06).

process MERGE_COUNTS {
    label 'process_single'

    container 'ghcr.io/rmrodell/nanomodamp-r:0.1.0' // pinned R 4.3.2 / Bioconductor 3.18 image (containers/nanomodamp-r)

    input:
    path counts
    path samplesheet

    output:
    path("counts_merged.tsv"), emit: merged

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    nma_merge.R --samplesheet ${samplesheet} --counts ${counts} --out counts_merged.tsv
    """

    stub:
    """
    printf "sample_id\tchr\tpos\tgene\ttotalReads\tA.count\tC.count\tG.count\tT.count\tDeletion.count\tInsertion.count\tref\tkmer\tstrand\tdelrate\n" > counts_merged.tsv
    """
}
