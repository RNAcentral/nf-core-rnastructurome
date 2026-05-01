# nf-core/rnastructurome: Output

## Introduction

This document describes the files produced by the pipeline and where to find them.
All paths below are relative to `--outdir`.

## Output layout at a glance

```text
<outdir>/
├── fastqc/
├── cat/                               (if multi-lane FASTQ merging occurs)
├── cutadapt/
├── bowtie/                            (RT-stop: alignment BAMs + Bowtie index)
│   └── bowtie/                        (Bowtie index files)
├── samtools/                          (sorted/deduplicated BAMs + QC)
├── umi_tools/                         (if UMI options are enabled)
├── count/
│   └── <sample>/
│       └── plots/
├── norm/
│   ├── <group>/                       (one per cell_line + replicate combination)
│   │   ├── xml/
│   │   ├── plots/
│   │   └── wiggle/
│   └── merged_bw/                     (per-cell_line genomic reactivity BigWigs)
├── fold/
│   ├── <group>/                       (one per cell_line)
│   │   ├── dotbracket/
│   │   ├── 2D-structures/
│   │   ├── summaries/
│   │   ├── dotplot/                   (if --rffold_dotplot)
│   │   ├── bp/                        (if --rffold_dotplot)
│   │   └── rdat/
│   ├── merged_bp/                     (merged base-pair files per cell_line)
│   └── shannon_bw/                    (if --rffold_shannon_entropy)
├── reference/                         (only when reference downloaded from Ensembl)
├── multiqc/
└── pipeline_info/
```

## Pipeline overview

Major stages:

1. Read QC (`FastQC`)
2. Read processing (`cutadapt`, optional `umi_tools`)
3. Alignment (`bowtie`/`bowtie2`, `samtools`)
4. RNAframework reactivity processing (`rf-count`, `rf-norm`, `rf-fold`)
5. Aggregated reporting (`MultiQC` + `pipeline_info`)

## FastQC

### Files

- `fastqc/*_fastqc.html`
- `fastqc/*_fastqc.zip`

### Meaning

FastQC is run on raw reads and reports per-base quality, sequence content, adapters, and overrepresented sequences.

## FASTQ merging

### Files

- `cat/<sample>.merged.fastq.gz`

### Meaning

When a sample has reads split across multiple FASTQ files (e.g. multi-lane sequencing), they are concatenated before trimming. This directory is only present when merging occurs.

## Cutadapt

### Files

- `cutadapt/<sample>.trim.fastq.gz`
- `cutadapt/<sample>.cutadapt.log`

### Meaning

Adapter-trimmed reads and the associated trimming log. Trimming parameters are derived from sample metadata and pipeline parameters.

## Alignment

### `bowtie/` (RT-stop experiments)

- `bowtie/<sample>.bam` — aligned reads
- `bowtie/<sample>.out` — Bowtie alignment summary log
- `bowtie/bowtie/<index>.*` — Bowtie index files built from the reference FASTA

For MaP experiments, Bowtie2 is used and produces an equivalent structure.

### `samtools/`

- `samtools/<sample>.sorted.bam` + `.sorted.bam.bai` — coordinate-sorted, duplicate-marked BAM and index
- `samtools/<sample>.bam` + `.bam.bai` — pre-sort BAM and index
- `samtools/<sample>.flagstat` — mapping summary
- `samtools/<sample>.idxstats` — per-reference mapping counts
- `samtools/<sample>.stats` — full samtools stats report
- `samtools/<reference>.sorted.fa.fai` — FASTA index

### Meaning

The sorted, duplicate-marked BAMs are consumed by `rf-count`. QC metrics are captured by MultiQC.

## RNAframework

### `count/` (`rf-count`)

Per-sample directories, e.g. `count/<sample>/`.

- `<sample>.rc` — raw per-position counts (RT-stop or MaP depending on sample principle)
- `index.rci` — index sidecar for the `.rc` file
- `<sample>.rfcount_summary.tsv` — per-sample summary metrics parsed from the `rf-count` log
- `plots/base_stats.pdf` — count plots (always generated)

### `norm/` (`rf-norm`, `rf-wiggle`)

Per-group directories under `norm/<group>/`, where `<group>` is the `cell_line_replicate` combination (e.g. `HEK293T_1`).

- `xml/<transcript>.xml` — normalized reactivity profiles used for downstream folding
- `rfnorm.log` — raw `rf-norm` console output
- `plots/<transcript>.pdf` — normalization plots (always generated)
- `wiggle/<transcript>.wig` — per-transcript reactivity wiggle tracks (one file per transcript, produced by `rf-wiggle`)

Merged and averaged BigWig under `norm/merged_bw/`:

- `merged_bw/<cell_line>_reactivity.bw` — genomic-coordinate reactivity BigWig for a cell line; when multiple replicates are available the per-replicate transcript tracks are averaged position-by-position, remapped to genomic coordinates with the GTF, and converted to BigWig

### `fold/` (`rf-fold`)

Per-group directories under `fold/<group>/`, where `<group>` is the cell line name (e.g. `HEK293T`). Replicates are folded together when present.

- `dotbracket/<transcript>.db` — predicted secondary structures in dot-bracket format
- `2D-structures/<transcript>.svg` — 2D structure diagram images
- `summaries/<transcript>.pdf` — fold summary plots
- `dotplot/<transcript>.dp` — dot-plot files (only when `--rffold_dotplot` is enabled)
- `bp/<transcript>.bp` — base-pair list derived from dot-plots (only when `--rffold_dotplot` is enabled)
- `rdat/<transcript>.rdat` — RDAT-format file combining reactivity (XML) and structure (dot-bracket)
- `rffold.log` — raw `rf-fold` console output
- `conversion_warnings.log` — warnings from dot-plot to base-pair conversion, if any

Merged base-pair files under `fold/merged_bp/`:

- `merged_bp/<cell_line>_merged.bp` — all per-transcript base-pair entries merged into a single file per cell line (for genome browser visualisation)

Shannon entropy BigWigs under `fold/shannon_bw/` (only when `--rffold_shannon_entropy` is enabled, which is the default):

- `shannon_bw/<cell_line>_shannon.bw` — genomic-coordinate per-position Shannon entropy across the transcript ensemble for a cell line

Notes:

- CT output is additionally produced when `--rffold_ct` is enabled.
- Genomic BigWigs convert common Ensembl chromosome names (`1`, `X`, `MT`) to UCSC-style names (`chr1`, `chrX`, `chrM`) for compatibility with IGV genomes such as `hg38`.
- If `rf-fold` reports errors or produces no structures, the task fails hard (pipeline stops instead of silently continuing).

## Reference

### Files

- `reference/<organism>/ensembl_source_urls.txt` — URLs of files downloaded from Ensembl
- `reference/<organism>/<organism>.sorted.fa` — sorted reference transcript FASTA

### Meaning

Only produced when the reference FASTA is downloaded from Ensembl (i.e. `--fasta` was not provided). When a local FASTA is supplied via `--fasta`, no reference files are published to the output directory.

## MultiQC

### Files

- `multiqc/multiqc_report.html`
- `multiqc/multiqc_data/`
- `multiqc/multiqc_plots/` — per-plot exports in PDF, PNG, and SVG

### Meaning

`multiqc_report.html` is the main summary report. It aggregates:

- FastQC read quality results
- Cutadapt trimming statistics
- Bowtie/Bowtie2 alignment rates
- Samtools flagstat, idxstats, and stats
- A custom `Count Progression` table with mapped reads before/after deduplication and `rf-count` covered transcripts

## Pipeline information

### Files

- `pipeline_info/execution_report_<timestamp>.html`
- `pipeline_info/execution_timeline_<timestamp>.html`
- `pipeline_info/execution_trace_<timestamp>.txt`
- `pipeline_info/pipeline_dag_<timestamp>.html`
- `pipeline_info/params_<timestamp>.json`
- `pipeline_info/nf_core_rnastructurome_software_mqc_versions.yml`

### Meaning

These files are for run provenance, traceability, and resource/runtime troubleshooting. Files are timestamped so multiple runs in the same output directory accumulate without overwriting.

## Troubleshooting missing outputs

1. Check process status in `.nextflow.log` and `pipeline_info/execution_trace_<timestamp>.txt`.
2. For fold problems, inspect `fold/<group>/rffold.log`.
3. Re-run with `-resume` after fixes to avoid recomputing successful tasks.
