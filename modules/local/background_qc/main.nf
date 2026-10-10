// BACKGROUND_QC: bin/nma_bgqc.R; only when count_all_bases (§6.6).

process BACKGROUND_QC {
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // TODO(WP3): nanomodamp-r image from containers/nanomodamp-r (not yet built/pushed)

    input:
    path merged

    output:
    path("background"), emit: dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    nma_bgqc.R ${merged} background
    """

    stub:
    """
    mkdir -p background && touch background/background_qc.tsv background/background_qc_by_treat.tsv background/background_qc.pdf background/background_qc.png
    """
}
