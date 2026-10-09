//
// Subworkflow with functionality specific to the rmrodell/nanomodamp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    //
    // Check parameter combinations the schema cannot express (D4, D29, D31, R-15)
    //
    validateInputParameters()

    //
    // Create channel from the sample sheet. The sheet was validated against
    // assets/schema_input.json above; it is read here with all columns so that
    // metadata columns are carried in meta (§5.1).
    //
    // With --input_counts (D32) there is no sample sheet: preprocessing and counting are skipped
    channel
        .fromList(input ? readSamplesheet(input) : [])
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet // channel: [ meta, fastq file or directory ]
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
//
// Reverse complement of a DNA sequence
//
def revcomp(String seq) {
    return seq.toUpperCase().reverse().tr('ACGT', 'TGCA')
}

//
// Parameter rules that nextflow_schema.json cannot express
//
def validateInputParameters() {
    def errors = []
    // D32: exactly one of --input (sample sheet) and --input_counts (merged count tables)
    if (params.input && params.input_counts) {
        errors << "--input and --input_counts are mutually exclusive: use --input to process FASTQs, or --input_counts to merge count tables from earlier runs and call sites."
    }
    if (!params.input && !params.input_counts) {
        errors << "Provide --input (sample sheet) or --input_counts (counts_merged.tsv tables from earlier runs)."
    }
    if (params.input_counts) {
        resolveInputCounts(params.input_counts)
    }
    if (params.input) {
        ['library_type', 'fasta', 'bed'].each { p ->
            if (!params[p]) {
                errors << "--${p} is required with --input."
            }
        }
        if (params.bed) {
            errors.addAll(checkBed(params.bed, params.bed_coordinates))
        }
    }
    if (revcomp(params.ont_adapter_sense_3p) != params.ont_adapter_antisense_5p.toUpperCase()) {
        errors << "--ont_adapter_antisense_5p must be the reverse complement of --ont_adapter_sense_3p."
    }
    if (revcomp(params.ont_adapter_sense_5p) != params.ont_adapter_antisense_3p.toUpperCase()) {
        errors << "--ont_adapter_antisense_3p must be the reverse complement of --ont_adapter_sense_5p."
    }
    // D29: pool orientation is an opt-in exception for legacy libraries without ONT adapters or UMIs
    if (params.orientation_adapters == 'pool') {
        if (params.library_type != 'mpra') {
            errors << "--orientation_adapters pool requires --library_type mpra (pool adapters exist only in MPRA libraries)."
        }
        if (params.umi) {
            errors << "--orientation_adapters pool cannot be combined with --umi true: the UMI sits outside the pool adapters and would be trimmed away. Set --umi false for libraries without UMIs (see docs/usage.md)."
        }
        if (revcomp(params.pool_adapter_3p) != params.pool_adapter_antisense_5p.toUpperCase()) {
            errors << "--pool_adapter_antisense_5p must be the reverse complement of --pool_adapter_3p."
        }
        if (revcomp(params.pool_adapter_5p) != params.pool_adapter_antisense_3p.toUpperCase()) {
            errors << "--pool_adapter_antisense_3p must be the reverse complement of --pool_adapter_5p."
        }
    }
    // D31: 3′-only ONT trimming is an opt-in exception for legacy libraries without the 5′ ONT adapter
    if (params.ont_adapter_mode == 'three_prime_only' && params.orientation_adapters != 'ont') {
        errors << "--ont_adapter_mode three_prime_only is valid only with --orientation_adapters ont."
    }
    if (params.analyses) {
        loadAnalysesConfig(params.analyses)
    }
    if (errors) {
        error("Invalid parameters:\n  - " + errors.join("\n  - "))
    }
    if (params.orientation_adapters != 'ont' || !params.umi || params.ont_adapter_mode != 'linked') {
        log.warn "Non-default preprocessing mode (orientation_adapters=${params.orientation_adapters}, ont_adapter_mode=${params.ont_adapter_mode}, umi=${params.umi}). " +
            "These settings are only for legacy libraries built without ONT adapters, the 5′ ONT adapter, or UMIs; see docs/usage.md."
    }
}

//
// BED checks (§3, D11, D19, D27/R-15): >= 6 columns; under bed0 a row with start = end is an error
// (the legacy single-site style needs --bed_coordinates one_based_start; never auto-detected)
//
def checkBed(bed, convention) {
    def errors = []
    def f = file(bed)
    if (!f.exists()) {
        return errors   // reported by schema validation
    }
    def zero_width = []
    def minus = 0
    f.eachLine { line, n ->
        if (!line.trim() || line.startsWith('#') || line.startsWith('track') || line.startsWith('browser')) {
            return
        }
        def c = line.split('\t')
        if (c.size() < 6) {
            errors << "--bed ${f.name} line ${n}: needs at least 6 tab-separated columns (chr, start, end, name, score, strand)"
            return
        }
        if (convention == 'bed0' && c[1] == c[2]) {
            zero_width << "${c[0]}:${c[1]}"
        }
        if (c[5] != '+') {
            minus++
        }
    }
    if (zero_width) {
        errors << "--bed ${f.name}: ${zero_width.size()} row(s) have start = end (e.g. ${zero_width.take(3).join(', ')}), which is empty in standard BED (--bed_coordinates bed0). " +
            "If this is a legacy single-site BED with 1-based start = end, rerun with --bed_coordinates one_based_start; otherwise convert it to start = site - 1, end = site (R-15)."
    }
    if (minus) {
        log.warn "--bed ${f.name}: ${minus} region(s) not on the '+' strand; Nano-BID-Amp amplicons are expected on '+' (D19)"
    }
    return errors
}

//
// --input_counts (D32): comma-separated counts_merged.tsv files, or directories of *.tsv files
//
def resolveInputCounts(input_counts) {
    def tables = []
    input_counts.toString().tokenize(',')*.trim().findAll { p -> p }.each { p ->
        def f = file(p)
        if (!f.exists()) {
            error("--input_counts: not found: ${p}")
        }
        if (f.isDirectory()) {
            def tsvs = f.listFiles().findAll { g -> g.name.endsWith('.tsv') }.sort { g -> g.name }
            if (!tsvs) {
                error("--input_counts: directory contains no .tsv files: ${p}")
            }
            tables.addAll(tsvs)
        } else {
            tables << f
        }
    }
    tables.each { t ->
        def header = t.withReader { r -> r.readLine() }?.split('\t') as List
        def missing = ['sample_id', 'chr', 'pos', 'totalReads', 'delrate'].findAll { c -> !(c in header) }
        if (missing) {
            error("--input_counts: ${t} is not a counts_merged.tsv table (missing columns ${missing}; docs/contracts/counts.md)")
        }
    }
    if (tables.size() != tables.unique(false).size()) {
        error("--input_counts: the same table is given more than once")
    }
    return tables
}

//
// minimap2 arguments: explicit --minimap2_args, else by library type (D5)
//
def resolveMinimap2Args() {
    if (params.minimap2_args) {
        return params.minimap2_args
    }
    return params.library_type == 'mpra' ? '-ax sr' : '-ax splice -uf'
}

//
// Read the sample sheet with all columns; resolve relative fastq paths against the
// sheet's directory; check files, empty directories and duplicates (§5.1)
//
def readSamplesheet(input) {
    def sheet = file(input, checkIfExists: true)
    def rows = []
    def seen = [] as Set
    sheet.splitCsv(header: true, strip: true).eachWithIndex { row, i ->
        if (row.values().every { v -> v == null || v.toString().trim() == '' }) {
            log.warn "Sample sheet ${sheet.name}: ignoring blank line ${i + 2}"
            return
        }
        def id = row.sample_id
        if (id in seen) {
            error("Sample sheet ${sheet.name}: duplicate sample_id '${id}'")
        }
        seen << id
        def fq = row.fastq.toString()
        def path = (fq.startsWith('/') || fq.contains('://')) ? file(fq) : file(sheet.parent.resolve(fq).toString())
        if (!path.exists()) {
            error("Sample sheet ${sheet.name}: fastq for '${id}' not found: ${path}")
        }
        if (path.isDirectory() && !path.listFiles().any { f -> f.name.endsWith('.fastq.gz') }) {
            error("Sample sheet ${sheet.name}: fastq directory for '${id}' contains no *.fastq.gz files: ${path}")
        }
        if (!path.isDirectory() && !path.name.endsWith('.fastq.gz')) {
            error("Sample sheet ${sheet.name}: fastq for '${id}' must be a .fastq.gz file or a directory: ${path}")
        }
        def meta = [id: id] + row.findAll { k, _v -> k != 'fastq' }
        rows << [meta, path]
    }
    if (!rows) {
        error("Sample sheet ${sheet.name} has no samples")
    }
    return rows
}

//
// Load the analyses YAML and check what the pipeline needs to schedule tasks (§5.4).
// Full validation against assets/schema_analyses.json happens in site calling (WP4).
//
def loadAnalysesConfig(analyses) {
    def cfg = new org.yaml.snakeyaml.Yaml().load(file(analyses, checkIfExists: true).text)
    def errors = []
    if (!(cfg instanceof Map) || !(cfg.analyses instanceof List) || !cfg.analyses) {
        error("--analyses ${analyses}: must contain a non-empty 'analyses' list (docs/contracts/analyses.md)")
    }
    def names = cfg.analyses.collect { a -> a.name }
    if (names.any { n -> !(n ==~ /[A-Za-z0-9._-]+/) }) {
        errors << "every analysis needs a name matching [A-Za-z0-9._-]+"
    }
    if (names.size() != names.unique(false).size()) {
        errors << "analysis names must be unique"
    }
    cfg.analyses.each { a ->
        if (!(a.type in ['treatment', 'factor'])) {
            errors << "analysis '${a.name}': type must be 'treatment' or 'factor'"
        }
        if (a.type == 'factor' && (!a.factor || !(a.levels instanceof List) || a.levels.size() != 2)) {
            errors << "factor analysis '${a.name}': needs 'factor' and 'levels: [baseline, experimental]'"
        }
    }
    (cfg.site_sets ?: []).each { ss ->
        def missing = (ss.of ?: []).findAll { n -> !(n in names) }
        if (missing) {
            errors << "site set '${ss.name}': unknown analyses ${missing}"
        }
    }
    if (errors) {
        error("--analyses ${analyses}:\n  - " + errors.join("\n  - "))
    }
    return cfg
}

def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "cutadapt, seqtk, UMI-tools, minimap2, SAMtools, UMICollapse, Rsamtools, blme and MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
