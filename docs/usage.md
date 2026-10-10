# rmrodell/nanomodamp: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

rmrodell/nanomodamp goes from demultiplexed Nano-BID-Amp FASTQ files to per-site deletion
counts and, optionally, pseudouridine site calls with standard plots. One pipeline run
processes one **sequencing run of one library type** (`endogenous` or `mpra`).

## Sample sheet

Pass a comma-separated sample sheet with a header row using `--input`:

```csv title="samplesheet.csv"
sample_id,fastq,treat,rep,celltype,vector
HepG2_WT_1_BS,fastq/HepG2_WT_1_BS.fastq.gz,BS,1,HepG2,WT
HepG2_WT_1_input,fastq/HepG2_WT_1_input.fastq.gz,input,1,HepG2,WT
HepG2_WT_2_BS,fastq/barcode07/,BS,2,HepG2,WT
```

| Column      | Description |
| ----------- | ----------- |
| `sample_id` | Unique sample name; letters, digits, `.`, `_` and `-` only. |
| `fastq`     | A `.fastq.gz` file, **or** a directory whose `*.fastq.gz` files are concatenated in lexical order. Relative paths are resolved against the sample sheet's directory. |
| `treat`     | Exactly `input` or `BS` (D14: no implicit recoding; label your no-enzyme controls `input`). |
| `rep`       | Batch-paired replicate: rep *k* of every condition comes from the same batch (D17). |
| any other   | Metadata (e.g. `celltype`, `vector`). Carried into `counts/counts_merged.tsv` and usable in analysis subsets and factors. |

Duplicate `sample_id`s, missing files, empty FASTQ directories and `treat` values other than
`input`/`BS` stop the run at startup. Blank lines are ignored with a warning. The full
contract is in [docs/contracts/samplesheet.md](contracts/samplesheet.md).

## Library type, reference and targets

| Parameter | Meaning |
| --------- | ------- |
| `--library_type` | `endogenous` (transcript amplicons) or `mpra` (oligo pool). Sets the minimap2 default (`mpra`: `-ax sr`; `endogenous`: `-ax splice -uf`, D5) and enables the MPRA pool-adapter trim. |
| `--fasta` | Transcript sequences (endogenous) or the oligo pool (MPRA), plain or `.gz`. The `.fai` index is always rebuilt by the pipeline; a `.fai` next to the FASTA is ignored, because a stale index (e.g. left over after the FASTA was rewritten) silently gives wrong reference bases (R-31). |
| `--bed` | Amplicons or single sites, at least 6 columns, `+` strand. |
| `--bed_coordinates` | `bed0` (default): standard 0-based, half-open BED; a single site at 1-based position *p* is `start = p-1, end = p`. `one_based_start`: legacy reading of the paper's BEDs, where start and end are both 1-based and inclusive (single site `start = end = p`). Under `bed0`, a row with start = end stops the run with a message pointing to `one_based_start`; the convention is never guessed (D27, R-15). |

## Analyses

Site calling runs when `--analyses analyses.yaml` is given; otherwise the pipeline stops after
counting. See [docs/contracts/analyses.md](contracts/analyses.md) and
`assets/analyses_example.yaml`.

## Libraries without ONT adapters or UMIs (legacy libraries)

**Use the defaults for all data**, including endogenous, MPRA in vitro and **MPRA in cellulo**
libraries: `--orientation_adapters ont`, `--ont_adapter_mode linked` and `--umi true`. Reads are
oriented with both ONT adapters required (D4), UMIs are extracted, and duplicates are removed with
UMICollapse.

Three opt-in settings exist **only** for legacy libraries that were built without some of these
elements. They are never chosen automatically (not from `--library_type`, and not from any sample
name or label such as "in cellulo"). Set them only when you know the library lacks the adapter or
UMI.

| Setting | Use only when the library has… | What changes |
| ------- | ------------------------------ | ------------ |
| `--orientation_adapters pool --umi false` (D29) | no ONT adapters and no UMIs (MPRA only) | Reads are oriented with the MPRA pool adapters (`--pool_adapter_5p/3p` sense, `--pool_adapter_antisense_5p/3p` antisense, `--trim1_min_length/--trim1_min_overlap`); the separate pool trim, UMI extraction and UMICollapse are skipped; the filtered BAM is the final BAM. |
| `--umi false` | ONT adapters but no UMIs | UMI extraction and UMICollapse are skipped. |
| `--ont_adapter_mode three_prime_only` (D31) | no 5′ ONT adapter on the reads | One cutadapt pass trims the 3′ ONT adapter (`-a <ont_adapter_sense_3p>`, `-m/-O` from `trim1_*`) and **keeps untrimmed reads**; there is no antisense pass or reverse complement. Only valid with `--orientation_adapters ont`. |

**Example: the golden test package.** The paper's in-cellulo MPRA samples from sequencing run
20241114 were made from a library without ONT adapters or UMIs, and the endogenous samples from
run 20250418 lack the 5′ ONT adapter. `conf/test_golden.config` is the only place the pipeline
sets these modes, and only for those samples:

```bash
# paper run 20241114 (in-cellulo samples, library without ONT adapters/UMIs)
nextflow run rmrodell/nanomodamp -profile test_golden,docker --golden_experiment mpra_incell --outdir results_incell
# equivalent explicit parameters for such a legacy library
nextflow run rmrodell/nanomodamp -profile docker --input samplesheet.csv --library_type mpra \
    --fasta pool.fa --bed sites.bed \
    --orientation_adapters pool --umi false --outdir results
```

New in-cellulo data sequenced with the standard protocol (ONT adapters and UMIs) use the
defaults; do **not** add these options.

**Read funnel.** `metrics/read_funnel.tsv` lists only the steps that ran:

| Mode | Steps |
| ---- | ----- |
| default (`ont`, `linked`, `umi`) | `raw, trim1_sense, trim1_antisense, antisense_rc, merged, umi_extracted, [trim2_pool], mapped, primary, mapq_filtered, dedup` |
| `pool` + `umi=false` | `raw, trim1_sense, trim1_antisense, antisense_rc, merged, mapped, primary, mapq_filtered` |
| `three_prime_only` | `raw, trim1_3prime, umi_extracted, [trim2_pool], mapped, primary, mapq_filtered, dedup` |

`trim2_pool` appears for `--library_type mpra` with ONT orientation. Skipped steps are omitted,
not written as zero. The pipeline log also prints a warning whenever a non-default mode is used.

**Validation errors** (the run stops at startup):

- `--orientation_adapters pool` with `--umi true`: the UMI sits outside the pool adapters and would be trimmed away.
- `--orientation_adapters pool` with `--library_type endogenous`: pool adapters exist only in MPRA libraries.
- `--ont_adapter_mode three_prime_only` with `--orientation_adapters pool`.
- An antisense adapter that is not the reverse complement of the matching sense adapter (ONT, and pool when used).

## Merging runs that need different preprocessing (D32)

Preprocessing parameters apply to a whole pipeline run. When samples that belong in one analysis
come from sequencing runs that need **different** preprocessing (for example a legacy library and
a standard one), process each run separately up to counting, then merge the count tables and call
sites once:

```bash
# 1. and 2.: preprocess and count each run (no --analyses)
nextflow run rmrodell/nanomodamp -profile docker --input run1.csv --library_type endogenous \
    --fasta tx.fa --bed sites.bed --ont_adapter_mode three_prime_only --outdir run1
nextflow run rmrodell/nanomodamp -profile docker --input run2.csv --library_type endogenous \
    --fasta tx.fa --bed sites.bed --outdir run2
# 3.: merge the count tables and call sites (no --input)
nextflow run rmrodell/nanomodamp -profile docker \
    --input_counts run1/counts/counts_merged.tsv,run2/counts/counts_merged.tsv \
    --analyses analyses.yaml --outdir merged
```

`--input_counts` takes comma-separated `counts_merged.tsv` files or a directory of them. It cannot
be combined with `--input`. Preprocessing and counting are skipped; the tables are merged by column
name and written to `counts/counts_merged.tsv`, with `counts/merge_sources.tsv` recording which
table each sample came from. Metadata columns are the union of the tables (missing values become
NA, with a warning). The run stops if a `sample_id` appears in more than one table or if the count
or site columns differ. Use the same reference, BED and `--bed_coordinates` for all runs you merge,
and keep `rep` values batch-paired across runs (D17); give replicates from different runs distinct
`rep` values unless they truly come from the same batch.

**Golden example.** The paper's endogenous samples come from run 20250418 (no 5′ ONT adapter,
`three_prime_only`) and run 20251022 (defaults). `test_golden` runs them as three invocations:

```bash
nextflow run rmrodell/nanomodamp -profile test_golden,docker --golden_experiment endogenous_20250418 --outdir golden_endogenous_20250418
nextflow run rmrodell/nanomodamp -profile test_golden,docker --golden_experiment endogenous_20251022 --outdir golden_endogenous_20251022
nextflow run rmrodell/nanomodamp -profile test_golden,docker --golden_experiment endogenous --outdir golden_endogenous
```

The third reads `golden_endogenous_*/counts/counts_merged.tsv` relative to the launch directory;
pass `--input_counts` to point elsewhere.

## Running the pipeline

```bash
nextflow run rmrodell/nanomodamp -profile docker --input ./samplesheet.csv --library_type mpra \
    --fasta pool.fa --bed sites.bed \
    --analyses analyses.yaml --outdir ./results
```

Test profiles: `test` (small synthetic dataset; until WP1 it is stub data and must be run with
`-stub`), `test_golden` (paper data subset, `--golden_experiment endogenous_20250418 |
endogenous_20251022 | endogenous | mpra_invitro | mpra_incell`; endogenous needs three runs, see above) and `test_full` (documented, not run in CI).

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

Parameters can also be given in a YAML/JSON file with `-params-file params.yaml`. Do not use
`-c <file>` to specify parameters, as this will result in errors; custom config files given
with `-c` must only be used for tuning process resource specifications and similar settings.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull rmrodell/nanomodamp
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [rmrodell/nanomodamp releases page](https://github.com/rmrodell/nanomodamp/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
