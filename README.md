# rmrodell/nanomodamp


[![GitHub Actions CI Status](https://github.com/rmrodell/NanoModAmp/actions/workflows/ci.yml/badge.svg)](https://github.com/rmrodell/NanoModAmp/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/rmrodell/NanoModAmp/actions/workflows/linting.yml/badge.svg)](https://github.com/rmrodell/NanoModAmp/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.7-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**rmrodell/nanomodamp** processes Nano-BID-Amp data (bisulfite-induced deletions read by Nanopore
amplicon sequencing) from demultiplexed FASTQ to per-site deletion counts and pseudouridine site
calls with standard plots. It reproduces the analyses of *RNA sequence, structure, and cell type
specific features drive pseudouridylation by PUS7* (Figures 2 and 3), with every difference
documented in the [Change Register](CHANGE_REGISTER.md).

> [!WARNING]
> **Under construction (WP0 scaffold).** Processes are stubs; only `-stub` runs work. See
> [PROGRESS.md](PROGRESS.md) and [plan.md](plan.md).

Steps (plan §2):

1. Concatenate FASTQs per sample.
2. Orient reads with linked ONT adapters (both required), reverse-complement antisense reads, merge.
3. Extract 3′ UMIs (umi_tools); MPRA only: trim the pool adapters (cutadapt).
4. Map (minimap2; `-ax sr` for MPRA, `-ax splice -uf` for endogenous), keep primary alignments with MAPQ ≥ 30.
5. Deduplicate (UMICollapse); report a read funnel.
6. Count bases, deletions and insertions per site (Rsamtools) and merge with sample-sheet metadata.
7. Call sites: treatment (BS vs input) and generic factor tests (blme::bglmer), site sets, plots (PDF + PNG), MultiQC.

Libraries built without ONT adapters, the 5′ ONT adapter, or UMIs (such as some of the paper's
runs) need opt-in settings; see [docs/usage.md](docs/usage.md#libraries-without-ont-adapters-or-umis-legacy-libraries).
All other data, including MPRA in cellulo data, use the defaults.

## Usage

```csv title="samplesheet.csv"
sample_id,fastq,treat,rep,celltype,vector
HepG2_WT_1_BS,fastq/HepG2_WT_1_BS.fastq.gz,BS,1,HepG2,WT
HepG2_WT_1_input,fastq/HepG2_WT_1_input.fastq.gz,input,1,HepG2,WT
```

```bash
nextflow run rmrodell/nanomodamp \
   -profile <docker/singularity/apptainer/conda/institute> \
   --input samplesheet.csv --library_type mpra \
   --fasta pool.fa --bed sites.bed \
   --analyses analyses.yaml --outdir <OUTDIR>
```

Test your installation with `-profile test` (synthetic data) and `-profile test_golden`
(paper data subset). See [usage](docs/usage.md), [output](docs/output.md),
[methods](docs/methods.md) and the [contracts](docs/contracts/README.md).

## Credits

rmrodell/nanomodamp was originally written by Rebecca Rodell.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use rmrodell/nanomodamp for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
