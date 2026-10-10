// SAMTOOLS_FAIDX: always builds a fresh .fai (R-31). The --fasta is copied (decompressed if .gz)
// to ref/<name>.fa and indexed there, so a .fai the user left next to --fasta is never used.
process SAMTOOLS_FAIDX {
    label 'process_single'

    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    path fasta

    output:
    path("ref/${name}"),     emit: fasta
    path("ref/${name}.fai"), emit: fai
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    name = fasta.name.replaceAll(/\.gz$/, '')
    def unpack = fasta.name.endsWith('.gz') ? "gzip -cd ${fasta}" : "cat ${fasta}"
    """
    mkdir -p ref
    ${unpack} > ref/${name}
    samtools faidx ref/${name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    name = fasta.name.replaceAll(/\.gz$/, '')
    """
    mkdir -p ref
    touch ref/${name} ref/${name}.fai
    printf '"%s":\\n    samtools: stub\\n' "${task.process}" > versions.yml
    """
}
