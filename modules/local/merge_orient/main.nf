// MERGE_ORIENT: sense reads followed by the reverse-complemented antisense reads (§6.1 step 3).
// Reads found in both orientation outputs are kept, as in the paper, and counted; above 0.1% of
// merged reads a warning is written (plan §6.1: raise as an open question).
process MERGE_ORIENT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(sense), path(antisense_rc)

    output:
    tuple val(meta), path("${prefix}.merged.fastq.gz")      , emit: reads
    tuple val(meta), path("${prefix}.both_orientations.tsv"), emit: both_orientations
    tuple val("${task.process}"), val('coreutils'), eval("cat --version | sed '1!d;s/.* //'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail
    export LC_ALL=C
    cat ${sense} ${antisense_rc} > ${prefix}.merged.fastq.gz

    ids() { bgzip -dc "\$1" | awk 'NR % 4 == 1 { print substr(\$1, 2) }' | sort -u; }
    both=\$(comm -12 <(ids ${sense}) <(ids ${antisense_rc}) | wc -l)
    merged=\$(bgzip -dc ${prefix}.merged.fastq.gz | awk 'END { print NR / 4 }')
    frac=\$(awk -v b="\$both" -v m="\$merged" 'BEGIN { printf "%.6f", (m > 0 ? b / m : 0) }')
    printf "sample_id\\tn_reads_in_both_orientations\\tn_merged_reads\\tfraction\\n%s\\t%s\\t%s\\t%s\\n" \\
        "${meta.id}" "\$both" "\$merged" "\$frac" > ${prefix}.both_orientations.tsv
    if awk -v f="\$frac" 'BEGIN { exit !(f > 0.001) }'; then
        echo "WARNING: ${meta.id}: \$both reads (\$frac of merged reads) were found in both orientations (> 0.1%; plan §6.1)" >&2
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | bgzip -c > ${prefix}.merged.fastq.gz
    printf "sample_id\\tn_reads_in_both_orientations\\tn_merged_reads\\tfraction\\n${meta.id}\\t0\\t0\\t0\\n" > ${prefix}.both_orientations.tsv
    """
}
