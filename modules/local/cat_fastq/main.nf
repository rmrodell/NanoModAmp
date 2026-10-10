// CAT_FASTQ: a single .fastq.gz, or a directory whose *.fastq.gz files are concatenated in
// lexical (C locale) order (§5.1, §6.1 step 1). gzip members concatenate into a valid gzip file.
process CAT_FASTQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}.raw.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('coreutils'), eval("cat --version | sed '1!d;s/.* //'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail
    export LC_ALL=C
    shopt -s nullglob
    src=\$(ls -d input/*)
    if [ -d "\$src" ]; then
        files=( "\$src"/*.fastq.gz )   # glob expansion is sorted; LC_ALL=C gives byte (lexical) order
        if [ \${#files[@]} -eq 0 ]; then
            echo "ERROR: FASTQ directory for sample '${meta.id}' contains no *.fastq.gz files: \$(basename "\$src")" >&2
            exit 1
        fi
        cat "\${files[@]}" > ${prefix}.raw.fastq.gz
    else
        if [ ! -s "\$src" ]; then
            echo "ERROR: FASTQ file for sample '${meta.id}' is missing or empty: \$(basename "\$src")" >&2
            exit 1
        fi
        ln -s "\$src" ${prefix}.raw.fastq.gz
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | bgzip -c > ${prefix}.raw.fastq.gz
    """
}
