/*
    PREPROCESS: raw FASTQ -> final BAM per sample (plan §6.1)

    Mode is set by run parameters (D29, D31), never by library_type:
      orientation_adapters = 'ont' (default) | 'pool' (opt-in, legacy libraries without ONT adapters)
      ont_adapter_mode     = 'linked' (default) | 'three_prime_only' (opt-in, legacy libraries without the 5′ ONT adapter)
      umi                  = true  (default) | false  (opt-in, legacy libraries without UMIs)
    The MPRA pool trim runs only for library_type 'mpra' with ONT orientation.
*/
include { CAT_FASTQ                                  } from '../../modules/local/cat_fastq'
include { CUTADAPT_3PRIME                            } from '../../modules/local/cutadapt_3prime'
include { CUTADAPT_LINKED as CUTADAPT_SENSE          } from '../../modules/local/cutadapt_linked'
include { CUTADAPT_LINKED as CUTADAPT_ANTISENSE      } from '../../modules/local/cutadapt_linked'
include { CUTADAPT_LINKED as CUTADAPT_POOL           } from '../../modules/local/cutadapt_linked'
include { SEQTK_RC                                   } from '../../modules/local/seqtk_rc'
include { MERGE_ORIENT                               } from '../../modules/local/merge_orient'
include { UMITOOLS_EXTRACT                           } from '../../modules/local/umitools_extract'
include { MINIMAP2_ALIGN                             } from '../../modules/local/minimap2_align'
include { SAMTOOLS_SORT                              } from '../../modules/local/samtools_sort'
include { SAMTOOLS_FILTER                            } from '../../modules/local/samtools_filter'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED  } from '../../modules/local/samtools_index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP     } from '../../modules/local/samtools_index'
include { UMICOLLAPSE                                } from '../../modules/local/umicollapse'
include { READ_FUNNEL                                } from '../../modules/local/read_funnel'

workflow PREPROCESS {
    take:
    ch_reads          // channel: [ meta, fastq file or directory ]
    ch_fasta          // channel: path(fasta), value
    minimap2_args     // string: resolved minimap2 arguments

    main:
    // Funnel step order (§5.2); steps a run skips are simply absent
    def funnel_steps = ['raw', 'trim1_3prime', 'trim1_sense', 'trim1_antisense', 'antisense_rc', 'merged',
                        'umi_extracted', 'trim2_pool', 'mapped', 'primary', 'mapq_filtered', 'dedup']
    def use_pool = params.orientation_adapters == 'pool'
    def (s5, s3, a5, a3) = use_pool
        ? [params.pool_adapter_5p, params.pool_adapter_3p, params.pool_adapter_antisense_5p, params.pool_adapter_antisense_3p]
        : [params.ont_adapter_sense_5p, params.ont_adapter_sense_3p, params.ont_adapter_antisense_5p, params.ont_adapter_antisense_3p]

    CAT_FASTQ(ch_reads)
    def ch_funnel = CAT_FASTQ.out.reads.map { meta, f -> [meta, 'raw', f] }

    def ch_trimmed
    def ch_both_orientations = channel.empty()
    if (params.ont_adapter_mode == 'three_prime_only') {
        // D31: legacy single 3′-adapter pass, untrimmed reads kept, no antisense pass or RC
        CUTADAPT_3PRIME(CAT_FASTQ.out.reads, params.ont_adapter_sense_3p, params.trim1_min_length, params.trim1_min_overlap)
        ch_trimmed = CUTADAPT_3PRIME.out.reads
        ch_funnel  = ch_funnel.mix(ch_trimmed.map { meta, f -> [meta, 'trim1_3prime', f] })
    } else {
        CUTADAPT_SENSE    (CAT_FASTQ.out.reads, s5, s3, params.trim1_min_length, params.trim1_min_overlap, 'trim1_sense')
        CUTADAPT_ANTISENSE(CAT_FASTQ.out.reads, a5, a3, params.trim1_min_length, params.trim1_min_overlap, 'trim1_antisense')
        SEQTK_RC(CUTADAPT_ANTISENSE.out.reads)
        MERGE_ORIENT(CUTADAPT_SENSE.out.reads.join(SEQTK_RC.out.reads))
        ch_funnel = ch_funnel
            .mix(CUTADAPT_SENSE.out.reads.map     { meta, f -> [meta, 'trim1_sense', f] })
            .mix(CUTADAPT_ANTISENSE.out.reads.map { meta, f -> [meta, 'trim1_antisense', f] })
            .mix(SEQTK_RC.out.reads.map           { meta, f -> [meta, 'antisense_rc', f] })
            .mix(MERGE_ORIENT.out.reads.map       { meta, f -> [meta, 'merged', f] })
        ch_trimmed = MERGE_ORIENT.out.reads
        ch_both_orientations = MERGE_ORIENT.out.both_orientations
    }

    if (params.umi) {
        UMITOOLS_EXTRACT(ch_trimmed, params.umi_pattern)
        ch_trimmed = UMITOOLS_EXTRACT.out.reads
        ch_funnel  = ch_funnel.mix(ch_trimmed.map { meta, f -> [meta, 'umi_extracted', f] })
    }

    if (params.library_type == 'mpra' && !use_pool) {
        CUTADAPT_POOL(ch_trimmed, params.pool_adapter_5p, params.pool_adapter_3p,
                      params.trim2_min_length, params.trim2_min_overlap, 'trim2_pool')
        ch_trimmed = CUTADAPT_POOL.out.reads
        ch_funnel  = ch_funnel.mix(ch_trimmed.map { meta, f -> [meta, 'trim2_pool', f] })
    }

    MINIMAP2_ALIGN(ch_trimmed, ch_fasta, minimap2_args)
    SAMTOOLS_SORT(MINIMAP2_ALIGN.out.bam)
    SAMTOOLS_FILTER(SAMTOOLS_SORT.out.bam, params.sam_exclude_flags, params.min_mapq)
    SAMTOOLS_INDEX_FILTERED(SAMTOOLS_FILTER.out.bam)
    ch_funnel = ch_funnel
        .mix(SAMTOOLS_SORT.out.bam.map       { meta, f -> [meta, 'mapped', f] })
        .mix(SAMTOOLS_FILTER.out.primary.map { meta, f -> [meta, 'primary', f] })
        .mix(SAMTOOLS_FILTER.out.bam.map     { meta, f -> [meta, 'mapq_filtered', f] })
    def ch_bam_bai = SAMTOOLS_INDEX_FILTERED.out.bam_bai

    if (params.umi) {
        UMICOLLAPSE(ch_bam_bai, params.umicollapse_args)
        SAMTOOLS_INDEX_DEDUP(UMICOLLAPSE.out.bam)
        ch_bam_bai = SAMTOOLS_INDEX_DEDUP.out.bam_bai
        ch_funnel  = ch_funnel.mix(UMICOLLAPSE.out.bam.map { meta, f -> [meta, 'dedup', f] })
    }

    // One funnel table per sample, steps in contract order
    def ch_funnel_in = ch_funnel
        .map { meta, step, f -> [meta.id, meta, step, f] }
        .groupTuple()
        .map { _id, metas, steps, files ->
            def order = steps.withIndex().sort { a, b -> funnel_steps.indexOf(a[0]) <=> funnel_steps.indexOf(b[0]) }*.get(1)
            [metas[0], order.collect { i -> steps[i] }, order.collect { i -> files[i] }]
        }
    READ_FUNNEL(ch_funnel_in)

    emit:
    bam_bai = ch_bam_bai                  // channel: [ meta, bam, bai ]
    funnel  = READ_FUNNEL.out.tsv         // channel: [ meta, tsv ]
    both_orientations = ch_both_orientations // channel: [ meta, tsv ]; empty with ont_adapter_mode = 'three_prime_only'
}
