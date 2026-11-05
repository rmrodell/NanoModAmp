// COUNT: bin/nma_count.R; nanomodamp::count_sites() (§6.2, §5.3).
// Counting options come from params through ext.args (conf/modules.config).
process COUNT {
    tag "$meta.id"
    label 'process_single'

    container 'ghcr.io/rmrodell/nanomodamp-r:0.1.0' // pinned R 4.3.2 / Bioconductor 3.18 image (containers/nanomodamp-r)

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path bed

    output:
    tuple val(meta), path("${prefix}.counts.tsv"), emit: counts
    tuple val(meta), path("${prefix}.failed_regions.tsv"), emit: failed_regions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    nma_count.R \\
        --bam ${bam} \\
        --fasta ${fasta} \\
        --bed ${bed} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "chr\tpos\tgene\ttotalReads\tA.count\tC.count\tG.count\tT.count\tDeletion.count\tInsertion.count\tref\tkmer\tstrand\tdelrate\n" > ${prefix}.counts.tsv
    printf "chr\tstart\tend\tgene\terror\n" > ${prefix}.failed_regions.tsv
    """
}
