// SITE_SETS: bin/nma_sitesets.R; unions/intersections/differences (§6.5).
// WP0 stub: the real command is implemented in WP4; the stub block lets `-stub` runs test the wiring.
process SITE_SETS {
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP4

    input:
    path calling_dirs, stageAs: 'calling/*'
    path analyses

    output:
    path("site_sets"), emit: dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "SITE_SETS is not implemented yet (WP4); use -stub" >&2
    exit 1
    """

    stub:
    """
    mkdir -p site_sets
    """
}
