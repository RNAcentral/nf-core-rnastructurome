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
│   ├── <reference>.norm_factors.txt      (if cross-experiment normalisation enabled)
│   ├── <group>/
│   │   ├── xml/*.xml
│   │   ├── plots/*.pdf
│   │   └── wiggle/*.wig
│   ├── genome_bw/*.bw
│   └── transcript_bw/*.bw
├── jackknife/                         (if --jackknife_reference provided)
├── fold/
│   ├── <group>/
│   │   ├── structures/r2dt/*.svg
│   │   ├── structures/viennarna/*.svg
│   │   ├── summaries/*.pdf
│   │   └── rdat/*.rdat
│   ├── genome_bp/*.bp
│   ├── transcript_bp/*.bp
│   ├── shannon_genome_bw/*.bw
│   └── shannon_transcript_bw/*.bw
├── eval/                              (if --eval_reference provided)
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

- `count/<sample>/`
  - `<sample>.rc`: Transcript-level raw per-position count file used as input to `rf-norm`.
  - `index.rci`: Index sidecar for the transcript-level count file.
  - `<sample>.rfcount_summary.tsv`: Per-sample count summary.
  - `plots/*.pdf`: Per-base count distribution plots.

</details>

<details markdown="1">
<summary>Output files — genome route</summary>

- `count/<sample>/`
  - `<sample>.plus.rc` / `<sample>.minus.rc`: Genome-coordinate strand-specific raw count files.
  - `<sample>.rc`: Transcript-level count file extracted from genome coordinates via `rf-rctools extract`, used as input to `rf-norm`.
  - `<sample>.rc.rci`: Index sidecar for the transcript-level count file.
  - `<sample>.rfcount_genome_summary.tsv`: Per-sample count summary.
  - `plots/*.pdf`: Per-base count distribution plots.

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
  - `plots/<transcript>.pdf`: Per-transcript normalisation plots.
  - `wiggle/<transcript>.wig`: Per-transcript reactivity wiggle tracks produced by `rf-wiggle`.
- `norm/genome_bw/`
  - `<sample_group>_reactivity_genome.bw`: Genomic-coordinate reactivity BigWig. When multiple replicates are present the per-replicate transcript-level tracks are averaged position-by-position, remapped to genomic coordinates using the GTF, and converted to BigWig format.
- `norm/transcript_bw/`
  - `<sample_group>_reactivity_transcript.bw`: Transcript-coordinate reactivity BigWig. Same averaged reactivity data as `genome_bw/` but in transcript coordinates, suitable for visualisation alongside transcript-level annotations.

</details>

[`rf-norm`](https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/) normalises raw RT-stop or MaP counts into per-nucleotide reactivity scores. Output is grouped by `sample_group` + `replicate` (e.g. `HEK293T_1`). The XML files are the primary output passed to downstream structure-prediction steps. Per-transcript reactivity plots are available as PDFs; the BigWigs provide reactivity tracks (averaged across replicates where applicable) in both genomic and transcript coordinates for genome browser visualisation.

### rf-fold

<details markdown="1">
<summary>Output files </summary>

- `fold/<group>/`
  - `structures/r2dt/<transcript>.svg`: Template-based 2D structure diagram drawn by R2DT, when a matching template is available.
  - `structures/viennarna/<transcript>.svg`: RNAplot 2D structure diagram drawn by ViennaRNA for transcripts not drawn by R2DT.
  - `structures/r2dt.log`: log of transcripts drawn/skipped by R2DT.
  - `summaries/<transcript>.pdf`: Per-transcript fold summary plot.
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

[`rf-fold`](https://rnaframework-docs.readthedocs.io/en/latest/rf-fold/) predicts RNA secondary structures from normalised reactivity profiles. Replicates for the same sample group are combined and folded together using ViennaRNA RNAfold in a windowed manner, which improves accuracy for long transcripts by folding overlapping sequence windows and merging the results. Key outputs per transcript include 2D structure diagrams, [RDAT] (https://rmdb.stanford.edu/deposit/specs/) files summarising the predicted structure and reactivity values, and per-transcript summary PDF plots. Genome-wide Shannon entropy profiles and base-pair arcs are additionally provided as BigWig and `.bp` tracks for visualisation in a genome browser, in both genome and transcript coordinates.

- **[R2DT](https://github.com/RNAcentral/R2DT)** — used when a matching template exists in the R2DT library (rRNA, snRNA, tRNA, etc.). Produces layouts comparable across organisms. Published to `fold/<group>/structures/r2dt/`.
- **[ViennaRNA](https://www.tbi.univie.ac.at/RNA/)** — fallback for transcripts without an R2DT template. `RNAplot` draws an energy-minimised 2D diagram from the dot-bracket structure. Published to `fold/<group>/structures/viennarna/`.

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
  - `*.csv`: Per-transcript evaluation metrics including sensitivity, PPV, and FMI comparing predicted structures to the reference.
  - `rfeval.log`: Raw `rf-eval` console output.

</details>

[`rf-eval`](https://rnaframework-docs.readthedocs.io/en/latest/rf-eval/) evaluates how well normalised reactivity profiles agree with a set of reference structures. It reports three metrics per transcript:

- **Unpaired Coefficient** — fraction of highly reactive bases that are unpaired
- **DSCI** — probability that a randomly selected unpaired base has higher reactivity than a paired base
- **AUROC** — area under the ROC curve treating reactivity as a classifier of unpaired bases

This step only runs when `--eval_reference` is provided.

 Enabled by `--rfeval_reference`. 

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
