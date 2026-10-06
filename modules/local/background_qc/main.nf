// BACKGROUND_QC: bin/nma_bgqc.R; only when count_all_bases (§6.6).
// WP0 stub: the real command is implemented in WP3; the stub block lets `-stub` runs test the wiring.
process BACKGROUND_QC {
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP3

    input:
    path merged

    output:
    path("background"), emit: dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "BACKGROUND_QC is not implemented yet (WP3); use -stub" >&2
    exit 1
    """

    stub:
    """
    mkdir -p background && touch background/background_qc.tsv background/background_qc.pdf background/background_qc.png
    """
}
