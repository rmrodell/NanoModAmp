/*
    CALL: counts_merged.tsv + analyses YAML -> one CALL_SITES task per analysis, then site sets
    (plan §6.3–6.5, §5.4, §5.5)
*/
include { CALL_SITES } from '../../modules/local/call_sites'
include { SITE_SETS  } from '../../modules/local/site_sets'

workflow CALL {
    take:
    ch_merged        // channel: counts_merged.tsv
    analyses_file    // path: analyses YAML
    analyses_config  // map: parsed analyses YAML

    main:
    def ch_analyses = channel.value(analyses_file)
    def ch_tasks = channel.fromList(analyses_config.analyses*.name)
        .combine(ch_merged)
        .combine(ch_analyses)
    CALL_SITES(ch_tasks)

    def ch_site_sets = channel.empty()
    if (analyses_config.site_sets) {
        SITE_SETS(CALL_SITES.out.dir.map { _name, dir -> dir }.collect(), ch_analyses)
        ch_site_sets = SITE_SETS.out.dir
    }

    emit:
    calling   = CALL_SITES.out.dir   // channel: [ name, calling/<name>/ ]
    site_sets = ch_site_sets         // channel: site_sets/
}
