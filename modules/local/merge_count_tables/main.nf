// MERGE_COUNT_TABLES: --input_counts (D32). Merges counts_merged.tsv tables from earlier runs by column
// name (R-06): metadata columns are the union in first-seen order (missing -> NA, with a warning);
// duplicate sample_id or differing count/site columns fail. Writes counts_merged.tsv (§5.3 layout)
// and merge_sources.tsv (sample_id -> source table).
// bin/nma_merge.R --tables (nanomodamp::merge_count_tables()).
process MERGE_COUNT_TABLES {
    label 'process_single'

    container 'docker.io/library/ubuntu:22.04' // TODO(WP3): nanomodamp-r image from containers/nanomodamp-r (not yet built/pushed)

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
    nma_merge.R --tables ${tables.join(' ')} --sources ${sources.collect { "'${it}'" }.join(' ')} \\
        --out counts_merged.tsv --sources-out merge_sources.tsv
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
