/*
    PREPROCESS: raw FASTQ -> final BAM per sample (plan §6.1)

    Mode is set by run parameters (D29, D31), never by library_type:
      orientation_adapters = 'ont' (default) | 'pool' (opt-in, legacy libraries without ONT adapters)
      ont_adapter_mode     = 'linked' (default) | 'three_prime_only' (opt-in, legacy libraries without the 5′ ONT adapter)
      umi                  = true  (default) | false  (opt-in, legacy libraries without UMIs)
    The MPRA pool trim runs only for library_type 'mpra' with ONT orientation.

    Tool arguments (adapters, lengths, filters) are set in conf/modules.config from params.
*/
include { CAT_FASTQ                                    } from '../../../modules/local/cat_fastq'
include { CUTADAPT as CUTADAPT_3PRIME                  } from '../../../modules/nf-core/cutadapt'
include { CUTADAPT as CUTADAPT_SENSE                   } from '../../../modules/nf-core/cutadapt'
include { CUTADAPT as CUTADAPT_ANTISENSE               } from '../../../modules/nf-core/cutadapt'
include { CUTADAPT as CUTADAPT_POOL                    } from '../../../modules/nf-core/cutadapt'
include { SEQTK_SEQ as SEQTK_RC                        } from '../../../modules/nf-core/seqtk/seq'
include { MERGE_ORIENT                                 } from '../../../modules/local/merge_orient'
include { UMITOOLS_EXTRACT                             } from '../../../modules/nf-core/umitools/extract'
include { MINIMAP2_ALIGN                               } from '../../../modules/nf-core/minimap2/align'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PRIMARY       } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_MAPQ          } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED    } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP       } from '../../../modules/nf-core/samtools/index'
include { UMICOLLAPSE                                  } from '../../../modules/nf-core/umicollapse'
include { READ_FUNNEL                                  } from '../../../modules/local/read_funnel'

workflow PREPROCESS {
    take:
    ch_reads          // channel: [ meta, fastq file or directory ]
    ch_fasta          // channel: path(fasta), value
    minimap2_args     // string: resolved minimap2 arguments (passed to MINIMAP2_ALIGN via meta)

    main:
    // Funnel step order (§5.2); steps a run skips are simply absent
    def funnel_steps = ['raw', 'trim1_3prime', 'trim1_sense', 'trim1_antisense', 'antisense_rc', 'merged',
                        'umi_extracted', 'trim2_pool', 'mapped', 'primary', 'mapq_filtered', 'dedup']
    def use_pool = params.orientation_adapters == 'pool'

    // single_end is needed by the nf-core modules; it is removed again before emitting
    CAT_FASTQ(ch_reads.map { meta, reads -> [meta + [single_end: true], reads] })
    def ch_funnel = CAT_FASTQ.out.reads.map { meta, f -> [meta, 'raw', f] }

    def ch_trimmed
    def ch_both_orientations = channel.empty()
    if (params.ont_adapter_mode == 'three_prime_only') {
        // D31: legacy single 3′-adapter pass, untrimmed reads kept, no antisense pass or RC
        CUTADAPT_3PRIME(CAT_FASTQ.out.reads)
        ch_trimmed = CUTADAPT_3PRIME.out.reads
        ch_funnel  = ch_funnel.mix(ch_trimmed.map { meta, f -> [meta, 'trim1_3prime', f] })
    } else {
        // Linked passes, both adapters required (D4); ONT or pool adapters by orientation_adapters
        CUTADAPT_SENSE(CAT_FASTQ.out.reads)
        CUTADAPT_ANTISENSE(CAT_FASTQ.out.reads)
        SEQTK_RC(CUTADAPT_ANTISENSE.out.reads)
        MERGE_ORIENT(CUTADAPT_SENSE.out.reads.join(SEQTK_RC.out.fastx))
        ch_funnel = ch_funnel
            .mix(CUTADAPT_SENSE.out.reads.map     { meta, f -> [meta, 'trim1_sense', f] })
            .mix(CUTADAPT_ANTISENSE.out.reads.map { meta, f -> [meta, 'trim1_antisense', f] })
            .mix(SEQTK_RC.out.fastx.map           { meta, f -> [meta, 'antisense_rc', f] })
            .mix(MERGE_ORIENT.out.reads.map       { meta, f -> [meta, 'merged', f] })
        ch_trimmed = MERGE_ORIENT.out.reads
        ch_both_orientations = MERGE_ORIENT.out.both_orientations
    }

    if (params.umi) {
        UMITOOLS_EXTRACT(ch_trimmed)
        ch_trimmed = UMITOOLS_EXTRACT.out.reads
        ch_funnel  = ch_funnel.mix(ch_trimmed.map { meta, f -> [meta, 'umi_extracted', f] })
    }

    if (params.library_type == 'mpra' && !use_pool) {
        CUTADAPT_POOL(ch_trimmed)
        ch_trimmed = CUTADAPT_POOL.out.reads
        ch_funnel  = ch_funnel.mix(ch_trimmed.map { meta, f -> [meta, 'trim2_pool', f] })
    }

    // minimap2 (args by library type, D5; ext.args reads meta.minimap2_args), sorted inside the module (§6.1 step 6)
    MINIMAP2_ALIGN(
        ch_trimmed.map { meta, f -> [meta + [minimap2_args: minimap2_args], f] },
        ch_fasta.map { fasta -> [[id: fasta.baseName], fasta] },
        true, '', false, false
    )
    def ch_mapped = MINIMAP2_ALIGN.out.bam.map { meta, bam -> [meta.findAll { k, _v -> k != 'minimap2_args' }, bam] }

    // -F {sam_exclude_flags} (step 'primary'), then -q {min_mapq} (step 'mapq_filtered') (§6.1 step 7)
    SAMTOOLS_VIEW_PRIMARY(ch_mapped.map { meta, bam -> [meta, bam, []] }, [[], [], []], [[], []], [[], []], '')
    SAMTOOLS_VIEW_MAPQ(SAMTOOLS_VIEW_PRIMARY.out.bam.map { meta, bam -> [meta, bam, []] }, [[], [], []], [[], []], [[], []], '')
    SAMTOOLS_INDEX_FILTERED(SAMTOOLS_VIEW_MAPQ.out.bam)
    ch_funnel = ch_funnel
        .mix(ch_mapped.map                     { meta, f -> [meta, 'mapped', f] })
        .mix(SAMTOOLS_VIEW_PRIMARY.out.bam.map { meta, f -> [meta, 'primary', f] })
        .mix(SAMTOOLS_VIEW_MAPQ.out.bam.map    { meta, f -> [meta, 'mapq_filtered', f] })
    def ch_bam_bai = SAMTOOLS_VIEW_MAPQ.out.bam.join(SAMTOOLS_INDEX_FILTERED.out.index)

    if (params.umi) {
        // UMICollapse with default settings (D7)
        UMICOLLAPSE(ch_bam_bai, 'bam')
        SAMTOOLS_INDEX_DEDUP(UMICOLLAPSE.out.bam)
        ch_bam_bai = UMICOLLAPSE.out.bam.join(SAMTOOLS_INDEX_DEDUP.out.index)
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
    bam_bai           = ch_bam_bai.map { meta, bam, bai -> [meta.findAll { k, _v -> k != 'single_end' }, bam, bai] } // channel: [ meta, bam, bai ]
    funnel            = READ_FUNNEL.out.tsv.map { meta, tsv -> [meta.findAll { k, _v -> k != 'single_end' }, tsv] }    // channel: [ meta, tsv ]
    both_orientations = ch_both_orientations.map { meta, tsv -> [meta.findAll { k, _v -> k != 'single_end' }, tsv] }  // channel: [ meta, tsv ]; empty with three_prime_only
    reads             = ch_trimmed.map { meta, f -> [meta.findAll { k, _v -> k != 'single_end' }, f] }                  // channel: [ meta, fastq.gz ]: reads as mapped
}
