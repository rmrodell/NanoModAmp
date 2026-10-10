#!/usr/bin/env nextflow
/*
    Local harness for the PREPROCESS subworkflow (WP2), for machines without nf-test.
    Run from the repository root, one invocation per mode, then check with check_outputs.py:

      nextflow run tests/fixtures/preprocess/harness.nf -profile apptainer \
          --samples mpra_ont,S_dir --library_type mpra --outdir <dir>

    The same fixtures and expectations are used by
    subworkflows/local/preprocess/tests/main.nf.test in CI.
*/
include { PREPROCESS } from '../../../subworkflows/local/preprocess'

params.samples = null   // comma-separated fixture names under fastq/ (file <name>.fastq.gz or directory <name>/)

workflow {
    def fixture_dir = "${projectDir}/fastq"
    def ch_reads = channel.fromList(params.samples.tokenize(',')).map { name ->
        def f = file("${fixture_dir}/${name}.fastq.gz")
        [[id: name], f.exists() ? f : file("${fixture_dir}/${name}", type: 'dir', checkIfExists: true)]
    }
    def mm2 = params.minimap2_args ?: (params.library_type == 'mpra' ? '-ax sr' : '-ax splice -uf')
    PREPROCESS(ch_reads, channel.value(file("${projectDir}/ref.fa", checkIfExists: true)), mm2)

    def out = file(params.outdir)
    PREPROCESS.out.reads.subscribe { meta, f -> f.copyTo(out.resolve("reads/${meta.id}.fastq.gz")) }
    PREPROCESS.out.funnel.subscribe { meta, f -> f.copyTo(out.resolve("funnel/${meta.id}.read_funnel.tsv")) }
    PREPROCESS.out.both_orientations.subscribe { meta, f -> f.copyTo(out.resolve("both/${meta.id}.tsv")) }
    PREPROCESS.out.bam_bai.subscribe { meta, bam, bai ->
        bam.copyTo(out.resolve("final_bam/${meta.id}.bam"))
        bai.copyTo(out.resolve("final_bam/${meta.id}.bam.bai"))
    }
}
