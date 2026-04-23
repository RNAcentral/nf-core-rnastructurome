<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-rnastructurome_logo_dark.png">
    <img alt="nf-core/rnastructurome" src="docs/images/nf-core-rnastructurome_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/rnastructurome)
[![GitHub Actions CI Status](https://github.com/nf-core/rnastructurome/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/rnastructurome/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/rnastructurome/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/rnastructurome/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/rnastructurome/results)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/rnastructurome)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rnastructurome-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/rnastructurome)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/rnastructurome** is a bioinformatics pipeline for the analysis of chemical-based high-throughput RNA structure probing data. It accepts FASTQ files from SHAPE or DMS experiments using either the **RT-stop** or **mutational profiling (MaP)** principle, and processes them from raw reads through alignment and deduplication to per-base reactivity scores and RNA secondary structure predictions using the [RNAFramework](https://rnaframework.readthedocs.io) toolkit.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/rnastructurome \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --fasta transcripts.fa \
   --gtf annotation.gtf.gz \
   --outdir <OUTDIR>
```

For `rf-norm`, the samplesheet must include `cell_line`, `condition`, and `replicate` columns so samples are paired correctly during normalisation.

The minimum required samplesheet columns are `sample`, `fastq_1`, `cell_line`, `condition`, and `replicate`.

- Samples are grouped by identical `cell_line` and `replicate`.
- `treated` may be analysed on its own.
- `untreated` requires a matching `treated` sample with the same `cell_line` and `replicate`.
- `denatured` requires matching `treated` and `untreated` samples with the same `cell_line` and `replicate`.

`rf-norm` defaults are selected automatically from the probing principle and available controls:

- `RT-stop` with matching `untreated`: Ding scoring (`-sm 1`) with Box-plot normalisation (`-nm 3`)
- `RT-stop` without `untreated`: Rouskin scoring (`-sm 2`) with 90% Winsorizing (`-nm 2`)
- `MaP` with matching `untreated` and optional `denatured`: Siegfried scoring (`-sm 3`) with Box-plot normalisation (`-nm 3`)
- `MaP` without `untreated`: Zubradt scoring (`-sm 4`) with Box-plot normalisation (`-nm 3`)

Additional `rf-norm` parameters exposed by the pipeline:

- `--rfnorm_remap_reactivities`: remap normalized reactivities to the 0-1 range.
- `--rfnorm_reactive_bases <string>`: set the reactive bases used for normalization windows, e.g. `AC` for DMS.
- `--rfnorm_norm_window <int>`: set the normalization window size.
- `--rfnorm_window_offset <int>`: set the normalization window offset.
- `--rfnorm_dynamic_window <int>`: use dynamic normalization windows with at least this many reactive bases.
- `--rfnorm_norm_independent`: normalize each reactive base independently.
- `--rfnorm_norm_factor <float[,float]>`: supply a fixed normalization factor for all transcripts. For 90% Winsorizing, provide two comma-separated values.
- `--rfnorm_raw`: score raw reactivities without applying normalization.
- `--rfnorm_pseudocount <float>`: set the Ding pseudocount.
- `--rfnorm_max_score <float>`: set the Ding maximum score.
- `--rfnorm_ignore_lower_than_untreated`: set reactivities lower than untreated to zero for Ding/Siegfried methods.
- `--rfnorm_max_untreated_mut <float>`: set the Siegfried untreated mutation cutoff.
- `--rfnorm_max_mutation_rate <float>`: set the MaP mutation-rate cutoff.
- `--rfnorm_mean_coverage <float>`: discard transcripts below this mean coverage.
- `--rfnorm_median_coverage <float>`: discard transcripts below this median coverage.
- `--rfnorm_nan <int>`: report positions below this coverage as NaN. Default: `10`.
- `--rfnorm_img`: generate rf-norm plots. This automatically uses `--rnaframework_r_path` to locate `R` inside the RNAframework container.

If `--rfnorm_reactive_bases` is not provided, the pipeline sets `AC` automatically for samples with `method=DMS`. All other methods fall back to the RNAFramework default (`N`, all bases).

Example:

```csv
sample,fastq_1,cell_line,condition,replicate
HEK293T_treated_rep1,HEK293T_treated_rep1.fastq.gz,HEK293T,treated,1
HEK293T_untreated_rep1,HEK293T_untreated_rep1.fastq.gz,HEK293T,untreated,1
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/rnastructurome/usage) and the [parameter documentation](https://nf-co.re/rnastructurome/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/rnastructurome/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/rnastructurome/output).

## Credits

nf-core/rnastructurome was originally written by RNAcentral.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#rnastructurome` channel](https://nfcore.slack.com/channels/rnastructurome) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
