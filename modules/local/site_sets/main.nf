// SITE_SETS: bin/nma_sitesets.R; unions/intersections/differences (§6.5).
// Implemented in WP4 (R/nanomodamp: site_sets / site_quartiles); the stub block keeps `-stub` runs working.
process SITE_SETS {
    label 'process_single'

    container 'ghcr.io/rmrodell/nanomodamp-r:0.1.0' // pinned R 4.3.2 / Bioconductor 3.18 image (containers/nanomodamp-r)

    input:
    path calling_dirs, stageAs: 'calling/*'
    path analyses

    output:
    path("site_sets"), emit: dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    nma_sitesets.R \\
        --analyses ${analyses} \\
        --calling_dir calling \\
        --outdir site_sets
    """

    stub:
    """
    mkdir -p site_sets
    """
}
