# nf-core/rnastructurome: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnastructurome/usage](https://nf-co.re/rnastructurome/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. At minimum it must be a comma-separated file with the columns `sample`, `fastq_1`, `cell_line`, `condition`, and `replicate`, plus a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,cell_line,condition,replicate
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,HEK293T,treated,1
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,HEK293T,treated,1
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,HEK293T,treated,1
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the columns `sample`, `fastq_1`, `cell_line`, `condition`, and `replicate` to be present.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,cell_line,condition,replicate
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,HEK293T,treated,1
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,HEK293T,treated,2
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,HEK293T,treated,3
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,HEK293T,untreated,1
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,HEK293T,untreated,2
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,,HEK293T,untreated,3
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,,HEK293T,untreated,3
```

| Column      | Description                                                                                                                                                                            |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`    | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`   | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                           |
| `fastq_2`   | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                           |
| `cell_line` | Cell-line identifier used to pair samples for `rf-norm`.                                                                                                                               |
| `condition` | Sample condition for `rf-norm`. Allowed values are `treated`, `untreated`, and `denatured`.                                                                                           |
| `replicate` | Replicate identifier used to pair samples for `rf-norm`.                                                                                                                               |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

For `rf-norm`, samples are grouped by identical `cell_line` and `replicate` values.

- `treated` may be analysed on its own.
- `untreated` requires a matching `treated` sample with the same `cell_line` and `replicate`.
- `denatured` requires matching `treated` and `untreated` samples with the same `cell_line` and `replicate`.

`rf-norm` defaults are selected automatically from the probing principle and available controls:

- `RT-stop` with matching `untreated`: Ding scoring (`-sm 1`) with Box-plot normalisation (`-nm 3`)
- `RT-stop` without `untreated`: Rouskin scoring (`-sm 2`) with 90% Winsorizing (`-nm 2`)
- `MaP` with matching `untreated` and optional `denatured`: Siegfried scoring (`-sm 3`) with Box-plot normalisation (`-nm 3`)
- `MaP` without `untreated`: Zubradt scoring (`-sm 4`) with Box-plot normalisation (`-nm 3`)

### Cutadapt behavior by principle

Cutadapt settings depend on the sample `principle` (`RT-stop` or `MaP`):

| Principle | 5' quality trimming (`--cutadapt-5quality`) | 3' quality trimming (`--cutadapt-3quality`) | Notes                                                           |
| --------- | ------------------------------------------- | ------------------------------------------- | --------------------------------------------------------------- |
| `RT-stop` | Forced to `0`                               | Default `20`                                | 5' quality trimming is intentionally disabled for RT-stop data. |
| `MaP`     | Default `20`                                | Default `20`                                | Can be overridden with CLI parameters.                          |

Adapter trimming is optional and uses this precedence:

1. Per-sample values from the samplesheet (`adapter_5p`, `adapter_3p`)
2. Global parameters (`--cutadapt_adapter_5p`, `--cutadapt_adapter_3p`)

If no adapters are specified, adapter trimming is skipped and only quality/length trimming is applied.

Additional clipping controls are configurable (optional):

- `--cutadapt-len` (default: `25`): minimum read length kept after clipping.
- `--cutadapt-min-align` (default: `1`): minimum adapter overlap in nt to trigger adapter trimming.
- `--cutadapt-trim-N` (default: enabled): trims terminal `N` bases. Set `--cutadapt-trim-N false` to disable.

### Optional UMI extraction

Enable UMI extraction before cutadapt by providing a UMI pattern.

When enabled, each sample must have a UMI pattern available. Pattern precedence is:

1. Per-sample `umi_pattern` in the samplesheet
2. Global `--umi_pattern` fallback

The pipeline passes this pattern to `umi_tools extract` as `--bc-pattern`.

Pattern handling:

- If `umi_pattern` contains only `N`, `C`, `X`, the pipeline uses `umi_tools` string mode (`--bc-pattern`).
- If `umi_pattern` contains IUPAC degenerate bases (for example `D`), the pipeline automatically switches to regex mode (`--extract-method=regex`) and converts the pattern accordingly (e.g. `D -> [AGT]`).

### RNAframework rf-count options

The pipeline exposes these `rf-count` flags:

- `--rfcount_img` (default: `true`): enables statistics plots (`-g`).
- `--rfcount_trim_5prime` (default: `0`): number of 5' bases trimmed in `rf-count` (`-t5`).
- `--rfcount_mask_file` (optional): path to mask file (`--mask-file`).
- `--rfcount_primary_only` (default: `false`): primary alignments only (`--primary-only`).
- `--rfcount_paired_only` (default: `false`): paired-end reads where both mates map (`--paired-only`).
- `--rfcount_properly_paired` (default: `false`): paired-end reads mapped in proper pairs (`--properly-paired`).
- `--rfcount_map_sort_by_read_name` (default: `false`): in mutation mode, pre-sort read pairs by read name (`--sort-by-read-name`).
- `--rfcount_map_discard_shorter` (default: `1`): in mutation mode, discard reads shorter than this length (`--discard-shorter`).
- `--rfcount_map_min_quality` (default: `20`): in mutation mode, minimum base quality for mutation calls (`--min-quality`).
- `--rfcount_map_collapse_consecutive` (default: `true`): in mutation mode, collapse consecutive mutations (`--collapse-consecutive`).
- `--rfcount_map_max_collapse_distance` (default: `2`): max distance for mutation collapsing (`--max-collapse-distance`).
- `--rnaframework_container` (default: `docker.io/rnastructurome/rnaframework:2.9.6-r1`): container image used for local RNAframework modules.
- `--rnaframework_r_path` (default: `/usr/bin/R`): path to `R` inside the RNAframework container for `rf-count -g`.

Paired-end default behavior:

- If a sample is paired-end and neither `--rfcount_paired_only` nor `--rfcount_properly_paired` is set, the pipeline uses `--properly-paired` by default.
- `--rfcount_paired_only` and `--rfcount_properly_paired` are mutually exclusive; set only one.

Reference FASTA behavior:

- `rf-count -f` is always populated from the pipeline FASTA (`--fasta` or genome-config FASTA); no separate `rf-count` FASTA flag is required.
- Mutation-mode options above are applied only for MaP samples (`principle=MaP`), where `rf-count` runs with `-m`.
- In MaP mutation mode (`-m`), `--rfcount_trim_5prime` (`-t5`) has no effect (RNAframework behavior).

Build the default local RNAframework+R image before running with Docker:

```bash
docker build --platform linux/amd64 -t docker.io/rnastructurome/rnaframework:2.9.6-r1 docker/rnaframework-r
```

### RNAframework rf-norm options

The pipeline exposes these `rf-norm` flags:

- `--rfnorm_remap_reactivities` (default: `false`): remaps normalized reactivities to the 0-1 range (`--remap-reactivities`).
- `--rfnorm_reactive_bases` (optional): reactive bases used for normalization windows (`--reactive-bases`), e.g. `AC` for DMS.
- `--rfnorm_norm_window` (optional): normalization window size (`--norm-window`).
- `--rfnorm_window_offset` (optional): normalization window offset (`--window-offset`).
- `--rfnorm_dynamic_window` (optional): dynamically resize normalization windows to include at least this many reactive bases (`--dynamic-window`).
- `--rfnorm_norm_independent` (default: `false`): normalize each reactive base independently (`--norm-independent`).
- `--rfnorm_norm_factor` (optional): use a fixed normalization factor for all transcripts (`--norm-factor`). For 90% Winsorizing, provide two comma-separated values for the 5th and 95th percentiles.
- `--rfnorm_raw` (default: `false`): score raw reactivities without normalization (`--raw`).
- `--rfnorm_pseudocount` (optional): Ding scoring pseudocount (`--pseudocount`).
- `--rfnorm_max_score` (optional): Ding scoring maximum score (`--max-score`).
- `--rfnorm_ignore_lower_than_untreated` (default: `false`): set reactivities lower than untreated to zero for Ding/Siegfried methods (`--ignore-lower-than-untreated`).
- `--rfnorm_max_untreated_mut` (optional): maximum untreated mutation rate for Siegfried scoring (`--max-untreated-mut`).
- `--rfnorm_max_mutation_rate` (optional): maximum mutation rate for MaP scoring methods (`--max-mutation-rate`).
- `--rfnorm_mean_coverage` (default: `0`): discard transcripts below this mean coverage (`--mean-coverage`).
- `--rfnorm_median_coverage` (default: `0`): discard transcripts below this median coverage (`--median-coverage`).
- `--rfnorm_nan` (default: `10`): positions below this coverage are reported as NaN (`--nan`).
- `--rfnorm_img` (default: `false`): enables rf-norm plots of raw reactivity data across samples (`--img`).
- `--rnaframework_r_path` (default: `/usr/bin/R`): path to `R` inside the RNAframework container for both `rf-count -g` and `rf-norm -g`.

Notes:

- If `--rfnorm_reactive_bases` is not provided, the pipeline sets `AC` automatically for `method=DMS`. All other methods fall back to the RNAFramework default (`N`, all bases).
- `--rfnorm_img` automatically appends `-R ${params.rnaframework_r_path}` so `R` must be available in the RNAframework container.
- `--rfnorm_dynamic_window` is most useful together with `--rfnorm_reactive_bases` for base-specific chemistries such as DMS.
- The pipeline still selects `rf-norm` scoring (`-sm`) and normalization (`-nm`) defaults automatically from `principle` plus the available treated / untreated / denatured controls.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/rnastructurome --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/rnastructurome -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/rnastructurome
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rnastructurome releases page](https://github.com/nf-core/rnastructurome/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

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
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
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

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

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
