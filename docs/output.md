# nf-core/rnastructurome: Output

## Introduction

This document describes the output produced by the pipeline. Most of the early steps (QC, trimming, alignment) follow standard RNA-seq conventions; this page focuses on the RNA Framework modules and downstream outputs that are specific to this pipeline.

## Output layout

All paths are relative to the top-level output directory specified with `--outdir`.

```text
<outdir>/
├── fastqc/
├── cat/                               (if multi-lane FASTQ merging occurs)
├── cutadapt/
├── star/
├── bowtie/                            (if --transcriptome)
├── samtools/
├── umi_tools/                         (if UMI pattern provided)
├── count/
│   └── <sample>/*.{rc,rci}
│       └── plots/*.pdf
├── norm/
│   ├── <group>/
│   │   ├── xml/*.xml
│   │   ├── plots/*.pdf
│   │   └── wiggle/*.wig
│   └── merged_bw/*.bw
├── jackknife/                         (if --jackknife_reference provided)
├── fold/
│   ├── <group>/
│   │   ├── dotbracket/*.db
│   │   ├── 2D-structures/*.svg
│   │   ├── summaries/*.pdf
│   │   ├── dotplot/*.dp
│   │   ├── bp/*.bp
│   │   └── rdat/*.rdat
│   ├── merged_bp/*.bp
│   └── shannon_bw/*.bw
├── eval/                              (if --eval_reference provided)
├── reference/                         (only when reference downloaded from Ensembl)
├── multiqc/
└── pipeline_info/
```

## Standard steps (QC, trimming, alignment)

Read QC (FastQC), adapter trimming (Cutadapt), alignment (STAR), and SAMtools post-processing follow standard nf-core/rnaseq conventions. See the [nf-core/rnaseq output documentation](https://nf-co.re/rnaseq/3.14.0/docs/output/) for a description of those output files.

---

## RNA Framework

This pipeline uses modules from [RNA Framework](https://rnaframework-docs.readthedocs.io/en/latest/).

### rf-count

<details markdown="1">
<summary>Output files</summary>

- `count/<sample>/`
  - `<sample>.rc`: Transcript-level raw per-position count file used as input to `rf-norm`.
  - `<sample>.rc.rci`: Index sidecar for the transcript-level count file.
  - `<sample>.rfcount_summary.tsv`: Per-sample count summary from `rf-count` on the transcriptome route.
  - `<sample>.rfcount_genome_summary.tsv`: Per-sample count summary from `rf-count-genome` on the genome route.
  - `<sample>.plus.rc` / `<sample>.minus.rc`: Genome-coordinate strand-specific raw count files from `rf-count-genome`, when produced on the genome route.
  - `index.rci`: RNA Framework count index file, when produced by `rf-count` or `rf-count-genome`.
  - `plots/*.pdf`: Per-base count distribution plots.
  - `error.out` / `samtools.log`: Diagnostic logs produced by RNA Framework, when present.

</details>

[`rf-count`](https://rnaframework-docs.readthedocs.io/en/latest/rf-count/) calculates per-base RT-stop or mutation counts and read coverage from aligned reads. On the transcriptome route, the pipeline runs `rf-count` directly against transcript-coordinate BAM files. On the genome route, it runs [`rf-count-genome`](https://rnaframework-docs.readthedocs.io/en/latest/rf-count-genome/) first, then converts genome-coordinate count files to transcript-level `.rc` files with `rf-rctools extract`.

The transcript-level `<sample>.rc` and `<sample>.rc.rci` files are passed to `rf-norm` for reactivity normalisation.

### rf-norm

<details markdown="1">
<summary>Output files</summary>

- `norm/<group>/`
  - `xml/<transcript>.xml`: Normalised reactivity profiles in RNA Framework XML format, used as input for `rf-jackknife` and `rf-fold`.
  - `rfnorm.log`: Raw `rf-norm` console output.
  - `plots/<transcript>.pdf`: Per-transcript normalisation plots.
  - `wiggle/<transcript>.wig`: Per-transcript reactivity wiggle tracks produced by `rf-wiggle`.
- `norm/merged_bw/`
  - `<cell_line>_reactivity.bw`: Genomic-coordinate reactivity BigWig. When multiple replicates are present the per-replicate transcript-level tracks are averaged position-by-position, remapped to genomic coordinates using the GTF, and converted to BigWig format.

</details>

[`rf-norm`](https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/) normalises raw RT-stop or MaP counts into per-nucleotide reactivity scores. Output is grouped by `cell_line` + `replicate` (e.g. `HEK293T_1`). The XML files are the primary output passed to downstream structure-prediction steps.

### rf-jackknife

<details markdown="1">
<summary>Output files</summary>

- `jackknife/<group>/`
  - `FMI.csv`: FMI (Fowlkes–Mallows Index) scores as a semicolon-delimited matrix — rows are slope values, columns are intercept values, cells are FMI scores. The optimal slope/intercept pair is the cell with the highest value.
  - `rfjackknife.log`: Raw `rf-jackknife` console output including the optimal slope/intercept summary.
  - `FMI.pdf`: Heatmap of FMI scores across the grid _(only when `--rfjackknife_img` is enabled)_.

</details>

[`rf-jackknife`](https://rnaframework-docs.readthedocs.io/en/latest/rf-jackknife/) identifies optimal slope and intercept parameters for reactivity-guided structure folding by performing a grid search and comparing predicted structures to a set of reference structures using the Fowlkes–Mallows Index.

This step only runs when `--jackknife_reference` is provided. The optimal slope/intercept values from the CSV should be passed as `--rffold_slope` and `--rffold_intercept` for subsequent pipeline runs.

### rf-fold

<details markdown="1">
<summary>Output files</summary>

- `fold/<group>/`
  - `dotbracket/<transcript>.db`: Predicted secondary structure in dot-bracket notation.
  - `2D-structures/<transcript>.svg`: 2D structure diagram.
  - `summaries/<transcript>.pdf`: Fold summary plot.
  - `rdat/<transcript>.rdat`: RDAT-format file combining normalised reactivity (XML) and predicted structure (dot-bracket).
  - `rffold.log`: Raw `rf-fold` console output.
  - `dotplot/<transcript>.dp`: Dot-plot of base-pairing probabilities _(only when `--rffold_dotplot` is enabled)_.
  - `bp/<transcript>.bp`: Base-pair list derived from dot-plots _(only when `--rffold_dotplot` is enabled)_.
  - `conversion_warnings.log`: Warnings from dot-plot to base-pair conversion, if any.
- `fold/merged_bp/`
  - `<cell_line>_merged.bp`: All per-transcript base-pair entries merged into a single file per cell line, suitable for genome browser visualisation.
- `fold/shannon_bw/` _(only when `--rffold_shannon_entropy` is enabled, which is the default)_
  - `<cell_line>_shannon.bw`: Genomic-coordinate per-position Shannon entropy BigWig across the transcript ensemble for a cell line.

</details>

[`rf-fold`](https://rnaframework-docs.readthedocs.io/en/latest/rf-fold/) predicts RNA secondary structures from normalised reactivity profiles using ViennaRNA. Replicates for the same cell line are folded together.

:::note
- CT output is additionally produced when `--rffold_ct` is enabled.
- Genomic BigWigs convert Ensembl chromosome names (e.g. `1`, `X`, `MT`) to UCSC-style names (`chr1`, `chrX`, `chrM`) for compatibility with genome browsers.
- If `rf-fold` fails to produce structures for any transcript the task exits with an error rather than silently continuing.
:::

---

## Reference files

<details markdown="1">
<summary>Output files</summary>

- `reference/<organism>/`
  - `ensembl_source_urls.txt`: URLs of files downloaded from Ensembl.
  - `<organism>.fa`: Downloaded reference FASTA (genome FASTA for the default STAR route; transcript FASTA for `--transcriptome`).

</details>

Reference files are only published when they are downloaded automatically from Ensembl (i.e. `--genome_fasta` / `--transcriptome_fasta` / `--fasta` was not provided). When a local FASTA is supplied, nothing is written to this directory.

---

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: MultiQC report aggregating QC metrics across all samples.
  - `multiqc_data/`: Directory containing parsed statistics from each tool.
  - `multiqc_plots/`: Per-plot exports in PDF, PNG, and SVG formats.

</details>

[MultiQC](http://multiqc.info) aggregates results from all pipeline tools into a single report. The report includes:

- FastQC read quality metrics
- Cutadapt trimming statistics
- Bowtie/Bowtie2/STAR alignment rates
- SAMtools flagstat, idxstats, and stats
- A custom **Count Progression** table showing mapped reads before/after deduplication and the number of transcripts covered by `rf-count`

---

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_<timestamp>.html`: Nextflow execution report with per-process resource usage.
  - `execution_timeline_<timestamp>.html`: Timeline view of task execution.
  - `execution_trace_<timestamp>.txt`: Tab-separated trace file with runtime and resource metrics per task.
  - `pipeline_dag_<timestamp>.html`: Directed acyclic graph of the pipeline.
  - `params_<timestamp>.json`: Snapshot of all parameters used for the run.
  - `nf_core_rnastructurome_software_mqc_versions.yml`: Software versions for all tools used.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. The files listed above are generated by default when the pipeline finishes. All files are timestamped so successive runs in the same output directory accumulate without overwriting prior results.
