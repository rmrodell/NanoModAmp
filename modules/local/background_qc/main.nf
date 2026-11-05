// BACKGROUND_QC: bin/nma_bgqc.R; only when count_all_bases (§6.6).

process BACKGROUND_QC {
    label 'process_single'

    container 'ghcr.io/rmrodell/nanomodamp-r:0.1.0' // pinned R 4.3.2 / Bioconductor 3.18 image (containers/nanomodamp-r)

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
