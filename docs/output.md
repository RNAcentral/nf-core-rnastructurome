# nf-core/rnastructurome: Output

## Introduction

This document describes the output produced by the pipeline. Most of the early steps (QC, trimming, alignment) follow standard RNA-seq conventions and their outputs are described in the [nf-core/rnaseq](https://nf-co.re/rnaseq/3.14.0/docs/output/) documentation; this page focuses on the RNA Framework modules and downstream outputs specific to this pipeline.

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
│   ├── rfcount_summary_all_samples.tsv
│   └── <sample>/*.{rc,rci}
├── norm/
│   ├── <reference>.norm_factors.txt      (if cross-experiment normalisation enabled)
│   ├── <reference>.rfnormfactor.log
│   ├── <group>_<replicate>/
│   │   ├── xml/*.xml
│   │   ├── wiggle/*.wig
│   │   └── rfnorm.log
│   ├── genome_bw/*.bw
│   └── transcript_bw/*.bw
├── jackknife/                         (if --jackknife_reference provided)
│   └── <group>/
│       ├── FMI.csv
│       └── rfjackknife.log
├── fold/
│   ├── <group>/
│   │   ├── structures/r2dt/*.svg
│   │   ├── structures/viennarna/*.svg
│   │   ├── rdat/*.rdat
│   │   ├── rffold.log
│   │   └── extracted_structures/      (if --structextract)
│   │       ├── dotbracket/*.db        (extracted motifs, dot-bracket)
│   │       └── images/*_ss.svg        (per-motif diagrams, reactivity-coloured)
│   ├── genome_bp/*.bp
│   ├── transcript_bp/*.bp
│   ├── shannon_genome_bw/*.bw
│   └── shannon_transcript_bw/*.bw
├── correlate/                         (if --correlate_replicates and >1 replicate)
│   ├── <group>.{pearson,spearman}.rfcorrelate.log
│   └── <group>.{pearson,spearman}_correlate/
│       ├── matrix.csv
│       └── pairwise/*.tsv
├── eval/                              (if --rfeval_reference provided)
│   └── <group>/
│       ├── *.csv
│       └── rfeval.log
├── reference/                         (only when reference downloaded from Ensembl)
│   ├── *.sorted.fa
│   ├── *.annotation.gtf.gz
│   ├── ensembl_fasta_source_url.txt
│   └── ensembl_gtf_source_url.txt
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
<summary>Output files — transcriptome route</summary>

- `count/rfcount_summary_all_samples.tsv`: Aggregated `rf-count` summary across all samples in the run (the per-sample summaries are merged into this single TSV rather than published individually).
- `count/<sample>/`
  - `<sample>.rc`: Transcript-level raw per-position count file used as input to `rf-norm`.
  - `index.rci`: Index sidecar for the transcript-level count file.

</details>

<details markdown="1">
<summary>Output files — genome route</summary>

- `count/<sample>/`
  - `<sample>.plus.rc` / `<sample>.minus.rc`: Genome-coordinate strand-specific raw count files.
  - `<sample>.rc`: Transcript-level count file extracted from genome coordinates via `rf-rctools extract`, used as input to `rf-norm`.
  - `<sample>.rc.rci`: Index sidecar for the transcript-level count file.
  - The aggregated `count/rfcount_summary_all_samples.tsv` (above) also covers genome-route samples.

</details>

[`rf-count`](https://rnaframework-docs.readthedocs.io/en/latest/rf-count/) calculates per-base RT-stop or mutation counts and read coverage from aligned reads. On the transcriptome route, the pipeline runs `rf-count` directly against transcript-coordinate BAM files. On the genome route, it runs [`rf-count-genome`](https://rnaframework-docs.readthedocs.io/en/latest/rf-count-genome/) on the genome first, then converts genome-coordinate count files to transcript-level `.rc` files with `rf-rctools extract`.

The transcript-level `<sample>.rc` and `<sample>.rc.rci` files are passed to `rf-norm` for reactivity normalisation.

### rf-normfactor (optional)

<details markdown="1">
<summary>Output files </summary>

- `norm/<reference>.norm_factors.txt`: Transcriptome-wide, cross-experiment normalisation factors for the reference, applied to every group's `rf-norm` via `-nf`.
- `norm/<reference>.rfnormfactor.log`: Raw `rf-normfactor` console output.

</details>

[`rf-normfactor`](https://rnaframework-docs.readthedocs.io/en/latest/rf-normfactor/) derives a single set of normalisation factors across all of a reference's samples so that reactivities are on a common scale for cross-sample comparison. It runs once per reference and only when cross-experiment normalisation is enabled for that reference (see `--rfnorm_use_normfactor` in the [usage docs](usage.md); auto-enabled for a reference with paired untreated controls or more than one treated sample). When it does not run, `rf-norm` normalises each group independently and no factor file is produced.

### rf-norm

<details markdown="1">
<summary>Output files </summary>

- `norm/<group>/`
  - `xml/<transcript>.xml`: Normalised reactivity profiles in RNA Framework XML format, used as input for `rf-fold` and optionally `rf-jackknife`.
  - `rfnorm.log`: Raw `rf-norm` console output.
  - `wiggle/<transcript>.wig`: Per-transcript reactivity wiggle tracks produced by `rf-wiggle`.
- `norm/genome_bw/`
  - `<sample_group>_reactivity_genome.bw`: Genomic-coordinate reactivity BigWig. When multiple replicates are present the per-replicate transcript-level tracks are averaged position-by-position, remapped to genomic coordinates using the GTF, and converted to BigWig format.
- `norm/transcript_bw/`
  - `<sample_group>_reactivity_transcript.bw`: Transcript-coordinate reactivity BigWig. Same averaged reactivity data as `genome_bw/` but in transcript coordinates, suitable for visualisation alongside transcript-level annotations.

</details>

[`rf-norm`](https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/) normalises raw RT-stop or MaP counts into per-nucleotide reactivity scores. Output is grouped by `sample_group` + `replicate` (e.g. `HEK293T_1`). The XML files are the primary output passed to downstream structure-prediction steps. The BigWigs provide reactivity tracks (averaged across replicates where applicable) in both genomic and transcript coordinates for genome browser visualisation.

### rf-correlate (replicate QC)

<details markdown="1">
<summary>Output files </summary>

- `correlate/<group>.pearson_correlate/` and `correlate/<group>.spearman_correlate/`
  - `matrix.csv`: Overall pairwise correlation matrix between the group's replicates.
  - `pairwise/<repA>_vs_<repB>.tsv`: Per-transcript correlation coefficients and p-values for each replicate pair.
  - `*.pdf`: Correlation heatmap (only with `--correlate_img`).
- `correlate/<group>.{pearson,spearman}.rfcorrelate.log`: Raw `rf-correlate` console output for each method.

</details>

[`rf-correlate`](https://rnaframework-docs.readthedocs.io/en/latest/rf-correlate/) quantifies replicate reproducibility by correlating reactivity profiles between the replicates of a sample group (transcriptome-wide and per-transcript). It runs for groups with more than one replicate (controlled by `--correlate_replicates`) and computes both **Pearson** (reactivity-capped) and **Spearman** (rank-based) correlations, published to separate `.pearson`/`.spearman` directories. The overall pairwise correlations are also summarised in the **MultiQC report** as a per-sample-group table (number of replicates, mean Pearson and mean Spearman correlation), giving an at-a-glance reproducibility check.

### rf-fold

<details markdown="1">
<summary>Output files </summary>

- `fold/<group>/`
  - `structures/r2dt/<transcript>.svg`: Template-based 2D structure diagram drawn by R2DT, when a matching template is available.
  - `structures/viennarna/<transcript>.svg`: RNAplot 2D structure diagram drawn by ViennaRNA for transcripts not drawn by R2DT.
  - `structures/r2dt.log`: log of transcripts drawn/skipped by R2DT.
  - `rdat/<transcript>.rdat`: RDAT-format file (compatible with RMDB) which summarises results (sequence, dot-bracket notation, reactivities) and key parameters for how it was produced.
  - `rffold.log`: Raw `rf-fold` console output.
  - `conversion_warnings.log`: Warnings from dot-plot to base-pair conversion, if any.
- `fold/genome_bp/`
  - `<sample_group>_genome.bp`: Base-pair arcs merged into a single file per sample group in genome coordinates, suitable for arc diagram visualisation in a genome browser such as IGV.
- `fold/transcript_bp/`
  - `<sample_group>_transcript.bp`: Same base-pair arcs in transcript coordinates (transcript ID as chromosome, 1-based transcript positions). Suitable for visualisation against a transcript-level reference in IGV.
- `fold/shannon_genome_bw/`
  - `<sample_group>_shannon_genome.bw`: Genomic-coordinate per-position Shannon entropy BigWig across the transcript ensemble for a sample group.
- `fold/shannon_transcript_bw/`
  - `<sample_group>_shannon_transcript.bw`: Transcript-coordinate Shannon entropy BigWig, equivalent to `shannon_genome_bw/` but in transcript coordinates.

</details>

[`rf-fold`](https://rnaframework-docs.readthedocs.io/en/latest/rf-fold/) predicts RNA secondary structures from normalised reactivity profiles. Replicates for the same sample group are combined and folded together using ViennaRNA RNAfold in a windowed manner, which improves accuracy for long transcripts by folding overlapping sequence windows and merging the results. Key outputs per transcript include 2D structure diagrams, [RDAT](https://rmdb.stanford.edu/deposit/specs/) files summarising the predicted structure and reactivity values. Genome-wide Shannon entropy profiles and base-pair arcs are additionally provided as BigWig and `.bp` tracks for visualisation in a genome browser, in both genome and transcript coordinates.

- **[R2DT](https://github.com/RNAcentral/R2DT)** — used when a matching template exists in the R2DT library (rRNA, snRNA, tRNA, etc.). Produces layouts comparable across organisms. Published to `fold/<group>/structures/r2dt/`.
- **[ViennaRNA](https://www.tbi.univie.ac.at/RNA/)** — fallback for transcripts without an R2DT template. `RNAplot` draws an energy-minimised 2D diagram from the dot-bracket structure. Published to `fold/<group>/structures/viennarna/`.

---

### rf-structextract (optional)

<details markdown="1">
<summary>Output files </summary>

- `fold/<group>/extracted_structures/dotbracket/*.db`: Extracted structural motifs for the sample group in dot-bracket notation (one multi-record file per transcript, each record named `<transcript>_<start>-<end>`), filtered to high-confidence, low-reactivity / low-Shannon elements meeting the configured criteria.
- `fold/<group>/extracted_structures/images/*_ss.svg`: One 2D diagram per extracted motif, drawn with ViennaRNA RNAplot and coloured by SHAPE reactivity in the same style as the rf-fold structure plots. Written unless `--structextract_plot false`.

</details>

[`rf-structextract`](https://rnaframework-docs.readthedocs.io/en/latest/rf-structextract/) extracts well-defined structural elements from the rf-fold output by combining the predicted structures with per-base reactivity and Shannon entropy. Only runs when `--structextract` is set (see [usage docs](usage.md) for the selection-criteria parameters). It identifies substructures whose bases are consistently below the transcript median for both reactivity and Shannon entropy — the signature of stably folded regions — and that pass the configured length, pairing, and (optionally) thermodynamic-significance filters. Outputs are nested under the sample group's `fold/<group>/` directory.

---

### rf-jackknife (optional)

<details markdown="1">
<summary>Output files </summary>

- `jackknife/<group>/`
  - `FMI.csv`: FMI (Fowlkes–Mallows Index) scores as a semicolon-delimited matrix — rows are slope values, columns are intercept values, cells are FMI scores. The optimal slope/intercept pair is the cell with the highest value.
  - `rfjackknife.log`: Raw `rf-jackknife` console output including the optimal slope/intercept summary.

</details>

[`rf-jackknife`](https://rnaframework-docs.readthedocs.io/en/latest/rf-jackknife/) identifies optimal slope and intercept parameters for reactivity-guided structure folding by performing a grid search and comparing predicted structures to a set of reference structures using the Fowlkes–Mallows Index.

This step only runs when `--jackknife_reference` is provided.


### rf-eval (optional)

<details markdown="1">
<summary>Output files</summary>

- `eval/<group>/`
  - `*.csv`: Per-transcript evaluation metrics — the Unpaired Coefficient, DSCI, and AUROC — comparing reactivity profiles to the reference.
  - `rfeval.log`: Raw `rf-eval` console output.

</details>

[`rf-eval`](https://rnaframework-docs.readthedocs.io/en/latest/rf-eval/) evaluates how well normalised reactivity profiles agree with a set of reference structures. It reports three metrics per transcript:

- **Unpaired Coefficient** — fraction of highly reactive bases that are unpaired
- **DSCI** — probability that a randomly selected unpaired base has higher reactivity than a paired base
- **AUROC** — area under the ROC curve treating reactivity as a classifier of unpaired bases

This step only runs when `--rfeval_reference` is provided.

---

## Reference files

<details markdown="1">
<summary>Output files</summary>

- `reference/`
  - `ensembl_fasta_source_url.txt`: URL of the FASTA file downloaded from Ensembl, for provenance.
  - `ensembl_gtf_source_url.txt`: URL of the GTF file downloaded from Ensembl, for provenance.
  - `<organism>.sorted.fa`: Chromosome-sorted reference FASTA, ready for reuse in subsequent runs via `--genome_fasta` / `--transcriptome_fasta`.
  - `<organism>.annotation.gtf.gz`: Ensembl gene annotation GTF, ready for reuse via `--gtf`.

</details>

Reference files are only published when they are downloaded automatically from Ensembl (i.e. `--genome_fasta` / `--transcriptome_fasta` / `--fasta` was not provided). When a local FASTA is supplied, nothing is written to this directory. The sorted FASTA can be passed directly to subsequent runs to skip the download step.

---

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: MultiQC report aggregating QC metrics across all samples.

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
