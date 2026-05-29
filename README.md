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

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/rnastructurome)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rnastructurome-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/rnastructurome)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/rnastructurome** is a bioinformatics pipeline for the analysis of chemical-based high-throughput RNA structure probing data. It accepts FASTQ files from SHAPE or DMS experiments using either the **RT-stop** or **mutational profiling (MaP)** principle, and processes them from raw reads through alignment and deduplication to per-base reactivity scores, RNA secondary structure predictions, and publication-quality 2D structure diagrams.

The pipeline handles reference resolution automatically: supply a transcript FASTA and GTF directly, configure them via `params.genomes`, or let the pipeline download them from Ensembl by organism name. Samples are grouped by cell line and replicate so that treated, untreated, and denatured controls are paired correctly for normalisation.

Default pipeline steps:

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter and quality trimming ([`Cutadapt`](https://cutadapt.readthedocs.io/)) with principle-aware settings (5′ quality trimming disabled for RT-stop)
3. Optional UMI extraction ([`UMI-tools`](https://umi-tools.readthedocs.io/))
4. Alignment to transcript reference ([`Bowtie`](http://bowtie-bio.sourceforge.net/) for RT-stop; [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/) for MaP)
5. BAM sorting, indexing, and QC ([`SAMtools`](https://www.htslib.org/))
6. Optional UMI-aware deduplication ([`UMI-tools dedup`](https://umi-tools.readthedocs.io/)) or duplicate marking ([`SAMtools markdup`](https://www.htslib.org/))
7. Per-base reactivity counting ([`rf-count`](https://rnaframework.readthedocs.io/en/latest/rf-count/))
8. Reactivity normalisation with automatic control-pairing and scoring-method selection ([`rf-norm`](https://rnaframework.readthedocs.io/en/latest/rf-norm/))
9. RNA secondary structure prediction ([`rf-fold`](https://rnaframework.readthedocs.io/en/latest/rf-fold/))
10. Reactivity WIG tracks and genomic BigWig files ([`rf-wiggle`](https://rnaframework.readthedocs.io/en/latest/rf-wiggle/))
11. 2D structure diagram drawing: template-matched diagrams via [`R2DT`](https://github.com/RNAcentral/R2DT) where a template exists, [`ViennaRNA`](https://www.tbi.univie.ac.at/RNA/) RNAplot fallback for the remainder, both coloured by reactivity
12. RDAT export combining per-transcript reactivity and structure ([`rnaframework_to_rdat`](bin/rnaframework_to_rdat.py))
13. Aggregated QC report ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,cell_line,condition,replicate
HEK293T_treated_rep1,HEK293T_treated_rep1.fastq.gz,,HEK293T,treated,1
HEK293T_untreated_rep1,HEK293T_untreated_rep1.fastq.gz,,HEK293T,untreated,1
```

Each row is one sample. `fastq_2` is optional (leave empty for single-end). `cell_line`, `condition`, and `replicate` are used to pair treated/untreated/denatured controls for `rf-norm`.

Supported `condition` values: `treated`, `untreated`, `denatured`.

Now run the pipeline:

```bash
nextflow run nf-core/rnastructurome \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --fasta transcripts.fa \
   --gtf annotation.gtf.gz \
   --outdir <OUTDIR>
```

If you omit `--fasta` and `--gtf`, add an `organism` column to your samplesheet (e.g. `Homo sapiens`) and the pipeline will download the reference from Ensembl automatically.

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

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
