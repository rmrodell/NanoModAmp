// READ_FUNNEL: one row per executed step (§5.2): sample_id, step, file, records, bytes.
// FASTQ records = lines / 4 of the decompressed file (R-08: .gz counted correctly);
// BAM records = all alignment records (samtools view -c), so 'mapped' includes secondary and
// supplementary records and 'primary' shows how many the -F filter kept.
// Steps skipped by the run's parameters are omitted, not written as zero.
process READ_FUNNEL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), val(steps), path(files, stageAs: 'stage??/*')

    output:
    tuple val(meta), path("${prefix}.read_funnel.tsv"), emit: tsv
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_list = files instanceof List ? files : [files]
    def rows = [steps, file_list].transpose().collect { step, f -> "row ${step} ${f}" }.join('\n    ')
    """
    set -euo pipefail
    row() {
        local step=\$1 f=\$2 n
        case "\$f" in
            *.bam)              n=\$(samtools view -c "\$f") ;;
            *.fastq.gz|*.fq.gz) n=\$(bgzip -dc "\$f" | awk 'END { print NR / 4 }') ;;
            *)                  n=\$(awk 'END { print NR / 4 }' "\$f") ;;
        esac
        printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "${meta.id}" "\$step" "\$(basename "\$f")" "\$n" "\$(stat -L -c %s "\$f")" >> ${prefix}.read_funnel.tsv
    }
    printf "sample_id\\tstep\\tfile\\trecords\\tbytes\\n" > ${prefix}.read_funnel.tsv
    ${rows}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "sample_id\\tstep\\tfile\\trecords\\tbytes\\n" > ${prefix}.read_funnel.tsv
    for s in ${steps.join(' ')}; do printf "${meta.id}\\t\$s\\tNA\\t0\\t0\\n" >> ${prefix}.read_funnel.tsv; done
    """
}
