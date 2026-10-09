/*
    COUNT: final BAMs -> per-sample counts -> counts_merged.tsv (plan §6.2, §5.3; D9, D10, D12)
*/
include { SAMTOOLS_FAIDX } from '../../../modules/local/samtools_faidx'
include { COUNT          } from '../../../modules/local/count'
include { MERGE_COUNTS   } from '../../../modules/local/merge_counts'
include { BACKGROUND_QC  } from '../../../modules/local/background_qc'

workflow COUNT_SITES {
    take:
    ch_bam_bai      // channel: [ meta, bam, bai ]
    ch_fasta        // channel: path(fasta), value
    ch_bed          // channel: path(bed), value
    ch_samplesheet  // channel: path(samplesheet), value

    main:
    // .fai is built if missing (§3)
    def fai = file("${params.fasta}.fai")
    def ch_fai
    if (fai.exists()) {
        ch_fai = channel.value(fai)
    } else {
        SAMTOOLS_FAIDX(ch_fasta)
        ch_fai = SAMTOOLS_FAIDX.out.fai
    }

    COUNT(ch_bam_bai, ch_fasta, ch_fai, ch_bed)
    MERGE_COUNTS(COUNT.out.counts.map { _meta, tsv -> tsv }.collect(), ch_samplesheet)

    def ch_background = channel.empty()
    if (params.count_all_bases) {
        BACKGROUND_QC(MERGE_COUNTS.out.merged)
        ch_background = BACKGROUND_QC.out.dir
    }

    emit:
    counts         = COUNT.out.counts            // channel: [ meta, counts.tsv ]
    failed_regions = COUNT.out.failed_regions    // channel: [ meta, failed_regions.tsv ]
    merged         = MERGE_COUNTS.out.merged     // channel: counts_merged.tsv
    background     = ch_background               // channel: background/ (count_all_bases only)
}
