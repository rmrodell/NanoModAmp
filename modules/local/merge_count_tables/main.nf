// MERGE_COUNT_TABLES: --input_counts (D32). Merges counts_merged.tsv tables from earlier runs by column
// name (R-06): metadata columns are the union in first-seen order (missing -> NA, with a warning);
// duplicate sample_id or differing count/site columns fail. Writes counts_merged.tsv (§5.3 layout)
// and merge_sources.tsv (sample_id -> source table).
// WP0 stub: the real command is implemented in WP3 (bin/nma_merge.R); the stub block lets `-stub` runs test the wiring.
process MERGE_COUNT_TABLES {
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // placeholder, pinned in WP3

    input:
    path tables, stageAs: 'tables/??/*'
    val sources                         // original paths, same order as tables

    output:
    path("counts_merged.tsv"), emit: merged
    path("merge_sources.tsv"), emit: sources

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "MERGE_COUNT_TABLES is not implemented yet (WP3); use -stub" >&2
    exit 1
    """

    stub:
    """
    head -n 1 \$(ls tables/*/* | head -n 1) > counts_merged.tsv
    printf "sample_id\\tsource_table\\n" > merge_sources.tsv
    srcs=(${sources.join(' ')})
    i=0
    for t in tables/*/*; do
        tail -n +2 "\$t" >> counts_merged.tsv
        tail -n +2 "\$t" | cut -f1 | sort -u | awk -v t="\${srcs[\$i]}" 'BEGIN{OFS="\\t"} {print \$1, t}' >> merge_sources.tsv
        i=\$((i + 1))
    done
    """
}
