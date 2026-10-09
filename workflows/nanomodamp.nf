/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PREPROCESS             } from '../subworkflows/local/preprocess'
include { COUNT_SITES            } from '../subworkflows/local/count'
include { CALL                   } from '../subworkflows/local/call'
include { MERGE_COUNT_TABLES     } from '../modules/local/merge_count_tables'
include { resolveInputCounts     } from '../subworkflows/local/utils_nfcore_nanomodamp_pipeline'
include { resolveMinimap2Args    } from '../subworkflows/local/utils_nfcore_nanomodamp_pipeline'
include { loadAnalysesConfig     } from '../subworkflows/local/utils_nfcore_nanomodamp_pipeline'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanomodamp_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOMODAMP {

    take:
    ch_samplesheet // channel: [ meta, fastq file or directory ] read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    def ch_counts_merged
    if (params.input_counts) {
        //
        // D32: merge counts_merged.tsv tables from earlier runs; no preprocessing or counting
        //
        def count_tables = resolveInputCounts(params.input_counts)
        MERGE_COUNT_TABLES(channel.fromList(count_tables).collect(), count_tables*.toString())
        ch_counts_merged = MERGE_COUNT_TABLES.out.merged
    } else {
        def ch_fasta = channel.value(file(params.fasta, checkIfExists: true))
        def ch_bed   = channel.value(file(params.bed, checkIfExists: true))

        //
        // SUBWORKFLOW: preprocessing per sample (§6.1; mode from orientation_adapters / ont_adapter_mode / umi)
        //
        PREPROCESS(ch_samplesheet, ch_fasta, resolveMinimap2Args())
        PREPROCESS.out.funnel
            .map { _meta, tsv -> tsv }
            .collectFile(name: 'read_funnel.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/metrics")

        //
        // SUBWORKFLOW: counting and merge (§6.2, §5.3)
        //
        COUNT_SITES(PREPROCESS.out.bam_bai, ch_fasta, ch_bed, channel.value(file(params.input, checkIfExists: true)))
        ch_counts_merged = COUNT_SITES.out.merged
    }

    //
    // SUBWORKFLOW: site calling (§6.3–6.5), only when an analyses YAML is given
    //
    if (params.analyses) {
        CALL(ch_counts_merged, file(params.analyses, checkIfExists: true), loadAnalysesConfig(params.analyses))
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'nanomodamp_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'nanomodamp'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
