# nf-core/rnastructurome: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnastructurome/usage](https://nf-co.re/rnastructurome/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files supplied with `-c` can be used for infrastructure settings such as resources, executors, containers, or module arguments, but should not be used for ordinary pipeline parameters.

## Samplesheet input

The easiest way to run this pipeline is to create a full samplesheet that contains all of the information about each sample in a comma-separated file, including the header row shown below and pass it as `--input '[path to samplesheet file]'`:

```csv title="full_samplesheet.csv"
sample,sample_id,fastq_1,fastq_2,method,principle,cell_line,condition,replicate,organism,pH,adapter_3p,adapter_5p,umi_pattern
HEK293T_treated_r1,GSM000001,/data/treated_r1.fastq.gz,,SHAPE,RT-stop,HEK293T,treated,1,Homo sapiens,7.5,,,
HEK293T_untreated_r1,GSM000002,/data/untreated_r1.fastq.gz,,SHAPE,RT-stop,HEK293T,untreated,1,Homo sapiens,7.5,,,
```

An [example samplesheet](../assets/samplesheet.csv) is provided.

However, you can also provide a more minimal version if for example you don't need to specify some of the options, like in the example above you could decide to not provide the columns with information about the adapters or umi pattern if you are happy to use the default options. 
You can also provide a very minimal samplesheet with just the information required about each invidivual sample and pass the uniform values across all samples as parameters. For example:

```csv title="minimal_samplesheet.csv"
sample,fastq_1,cell_line,condition,replicate
HEK293T_treated_r1,/data/treated_r1.fastq.gz,HEK293T,treated,1
HEK293T_untreated_r1,/data/untreated_r1.fastq.gz,HEK293T,untreated,1
```
then pass essential but uniform options like this:

```bash 
--input /path/to/samplesheet.csv --method SHAPE --principle RT-stop --organism Homo sapiens
```

### Column reference

| Column        | Required | Description                                                                                           |
| ------------- | -------- | ----------------------------------------------------------------------------------------------------- |
| `sample`      | yes      | Sample identifier. Re-use the same name to concatenate re-sequenced samples (multi-lane support).     |
| `sample_id`   | no       | Unique identifier for an individual sequencing run.                                                   |
| `fastq_1`     | yes      | Read 1 FASTQ path (`.fastq.gz` / `.fq.gz`).                                                           |
| `fastq_2`     | no       | Read 2 FASTQ path for paired-end data.                                                                |
| `cell_line`   | yes      | Group key used for control pairing in `rf-norm`.                                                      |
| `condition`   | yes      | One of `treated`, `untreated`, `denatured`.                                                           |
| `replicate`   | yes      | Replicate key used for control pairing in `rf-norm`.                                                  |
| `method`      | no       | Probing chemistry: `SHAPE` or `DMS`. Controls chemistry-specific defaults. Falls back to `--method`.  |
| `principle`   | no       | Readout principle: `rt-stop` or `map`. Controls alignment, rf-count, and rf-norm defaults. Falls back to `--principle`. |
| `organism`    | no       | Used for automatic genome/transcriptome download from Ensembl/NCBI (e.g. `Homo sapiens`). Falls back to `--organism`. |
| `pH`          | no       | DMS reaction pH. `pH >= 8.0` sets reactive bases to `ACGU`; otherwise defaults to `AC`.              |
| `adapter_3p`  | no       | 3′ adapter sequence passed to Cutadapt.                                                               |
| `adapter_5p`  | no       | 5′ adapter sequence passed to Cutadapt.                                                               |
| `umi_pattern` | no       | UMI pattern passed to umi_tools. Supplying this enables UMI extraction before trimming.               |


## Other considerations and parameters

### Using your own genome/transcriptome files

This pipeline expects an `organism` parameter to be passed either in the samplesheet or as a global parameter, which will trigger an automatic download of FASTA and GTF files from Ensembl. In the cases of bacteria and viruses where Ensembl does not have a reference genome/transcriptome, it will look this up in NCBI and download the files from there instead. However, you can pass your own files if preferred. For this you need to use the `--fasta` and `--gtf` flags. If these flags are used, the pipeline will prioritise them over the Ensembl download even when `organism` is passed. Example of usage:

```bash
--input /path/to/samplesheet.csv --fasta /path/to/genome_or_transcriptome.fa.gz --gtf /path/to/annotation.gtf.gz
```

### Using genome versus transcriptome

By default the pipeline aligns to the genome with STAR and extracts transcript-level counts using the GTF. To align directly to the transcriptome instead, pass `--transcriptome` and the rest happens automatically. See the [Alignment routes](#alignment-routes) section for more detail on when each route is appropriate.

### Using optional modules

Two optional RNAframework modules are available:

1. **rf-jackknife** — tunes folding parameters against reference RNA structures. Pass `--jackknife_reference` with a path to a `.db` file of known structures. When provided, rf-jackknife runs between rf-norm and rf-fold and calibrates the slope and intercept parameters passed to rf-fold.
2. **rf-eval** — evaluates the agreement between your reactivity data and a set of reference RNA structures, computing metrics such as AUROC and DSCI that can be used as quality control. Pass `--rfeval_reference` with a path to a `.db` file of known structures to enable it.

## More details for key steps

### Adapter trimming

Trimming is done with Cutadapt. Defaults differ by principle because RT-stop experiments encode signal at the read 3′ end and disabling 5′ quality trimming preserves those positions.

| Principle | `--cutadapt_5quality` | `--cutadapt_3quality` |
|-----------|-----------------------|-----------------------|
| `RT-stop` | forced to `0` | default `20` |
| `MaP` | default `20` | default `20` |

Adapter precedence: per-sample columns (`adapter_5p`, `adapter_3p`) → global flags (`--cutadapt_adapter_5p`, `--cutadapt_adapter_3p`) → fallback `AGATCGGAAGAGC`.

Set `--cutadapt_quality_only` to skip adapter trimming and perform quality/length filtering only.

### UMI handling (optional)

By default the pipeline uses SAMtools markdup for deduplicating reads but UMI-tagged samples can be run through this pipeline by enabling `umi_pattern` (per-sample column or `--umi_pattern` if the pattern is the same for all samples).

- Patterns containing only `N/C/X` use `umi_tools --bc-pattern` mode directly.
- Patterns that contain IUPAC (e.g. containing `D`) are converted to `--extract-method=regex` mode automatically. Each IUPAC character is expanded to a regex character class (`D` → `[AGT]`, `N` → `[ACGT]`, etc.) and the pattern is anchored to the read start as a named capture group. For example, `DNNN` becomes `^(?P<umi_1>[AGT][ACGT][ACGT][ACGT])`.

### Reference resolution

The `organism` value (per-sample column or `--organism`) is used to resolve and download the reference automatically. Latin binomials (`Homo sapiens`) are preferred, but a small set of shorthands are also accepted: `human`, `mouse`, `rat`, and `yeast`.

For organisms not in Ensembl, the pipeline falls back to NCBI using pre-configured RefSeq accessions. The following are already built-in: Influenza A, SARS-CoV-2, Dengue virus, Zika virus, HIV-1, Rotavirus A, and Escherichia coli (K-12 MG1655). For any other organism, the pipeline will attempt an NCBI auto-search.

The pipeline always downloads the most recent available assembly - the current Ensembl release for eukaryotes, or the pinned RefSeq accession for NCBI organisms. To keep the reference version fixed across runs, run the pipeline once and then reuse the local files published under `reference/` by passing them with `--fasta` and `--gtf`.

Resolution order for each file type:

| File | 1st | 2nd | 3rd | 4th |
|------|-----|-----|-----|-----|
| FASTA | `--fasta` | `params.genomes[key]` | Ensembl download | NCBI fallback |
| GTF | `--gtf` | `params.genomes[key].gtf` | Ensembl download | Synthetic GTF (NCBI) |

### Alignment routes

**Default: genome alignment with STAR**

The pipeline downloads a soft-masked genome FASTA (`dna_sm.toplevel.fa.gz`) and GTF from Ensembl, builds a STAR index, and aligns with splice-junction awareness. Genome-coordinate counts from `rf-count-genome` are then converted to transcript-level RC files by `rf-rctools extract` using the GTF before passing to `rf-norm`. This is the recommended route for most experiments: STAR handles spliced reads correctly, the soft-masked genome reduces spurious multi-mappers from repetitive elements, and aligning to the genome avoids the need to choose between competing transcript isoforms. 

`--star_multimap_nmax` (default `10`) sets the maximum number of genome locations a read is allowed to map to — reads exceeding this are discarded. The default of 10 is appropriate for most protein-coding genes. If you are probing highly repetitive RNAs such as rRNA or snRNA you may want to increase this (e.g. `--star_multimap_nmax 50`), though be aware that allowing too many multimappers can introduce noise if reads are assigned ambiguously. Additional STAR flags can be passed via `ext.args` in a custom config.

**Optional: transcriptome alignment with Bowtie**

Add `--transcriptome` to align directly to a transcript FASTA instead of the genome. In this option, the pipeline downloads Ensembl cDNA + ncRNA FASTA files; for bacteria and viruses not in Ensembl it falls back to NCBI and builds a transcript FASTA from the genome assembly automatically. The aligner used depends on the probing principle: Bowtie for RT-stop and Bowtie2 for MaP. This route works well for bacteria and viruses where genome annotation is sparse or absent, or when you prefer to map directly to a curated set of transcripts.

`--bowtie_k` (default `1`) sets the maximum number of alignments to report per read — the equivalent of `--star_multimap_nmax` for the transcriptome route. Setting it to `1` means only uniquely mapping reads are kept. Increase it if you want to retain reads that map to multiple transcripts (e.g. paralogs or transcript isoforms), though the same caveats about ambiguous assignment apply. Use `--bowtie_all` instead to report all valid alignments which is useful if you want various isoforms to be reported. Additional flags can be passed via `ext.args` in a custom config.

### rf-count

`rf-count` (transcriptome route) or `rf-count-genome` (genome route) quantifies chemical probing signal from the aligned BAM files. What exactly is counted depends on the principle: for RT-stop experiments it tallies read 3′ ends that accumulate at modified bases; for MaP it measures per-position mutation rates. A per-transcript coverage and reactivity plot is generated for every sample.

On the genome route, strandedness is handled automatically — RT-stop libraries are always treated as second-strand (this is fixed by experimental design), while for MaP the pipeline infers strandedness per sample using `RSeQC infer_experiment`. If inference is ambiguous you can override it with `--rfcount_strandedness first|second|unstranded`. After counting, `rf-rctools extract` automatically converts the genome-coordinate RC files to transcript-level RC files using the GTF before passing to `rf-norm`.

A few parameters worth knowing about:

`--rfcount_trim_5prime` trims a fixed number of bases from the 5′ end of each read before counting. This is useful for RT-stop experiments where the first few bases after the adapter can carry sequence-context bias that inflates apparent stop rates.

`--rfcount_mask_file` accepts a BED file of regions to exclude from counting entirely — handy for masking rRNA or other highly-expressed contaminating transcripts that would otherwise dominate the output.

`--rfcount_primary_only` restricts counting to primary alignments. If you have increased `--star_multimap_nmax` or `--bowtie_k` to allow multimappers, you may want to pair this with `--rfcount_primary_only` to avoid double-counting reads that aligned to multiple loci.

For the full list of available options see the [rf-count documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-count/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

### rf-norm

`rf-norm` normalises per-position counts into reactivity scores. Samples are grouped by `cell_line + replicate`, and treated samples are normalised against their matched controls within the same group. The scoring and normalisation method are selected automatically based on the probing principle and which conditions are present: RT-stop with an untreated control uses Ding scoring with box-plot normalisation; without an untreated control it falls back to Rouskin scoring with Winsorizing. MaP with an untreated control uses Siegfried scoring; without, Zubradt. 

You can override the normalisation method alone (without changing scoring) with `--rfnorm_norm_method`:
- `1` is 2-8% normalisation (takes the top 10% of reactivities, discards the very highest 2%, and uses the mean of the remaining 8% as the scaling factor)
- `2` is 90% Winsorizing (clips any value above the 95th or below the 5th percentile to those boundaries, then divides all values by the 95th percentile)
- `3` is box-plot normalisation (removes outliers beyond 1.5× IQR then divides by the mean of the next top 10%), which is the default for most conditions
- `4` is Mitchell normalisation (MaP only, uses the higher of the mean 90th–95th percentile reactivity or the 75th percentile of non-zero reactivities as the scaling factor)

Treated samples are usually paired with untreated by matching `cell_line + replicate` exactly. However, in cases where the authors did not create an exact matching untreated sample for some specific treatments but have one for other samples in the same dataset, the pipeline falls back to an untreated sample sharing the same cell_line and replicate — for example, `MDA-MB-231_DMSO_treated_r1` will match exactly to `MDA-MB-231_DMSO_untreated_r1`, but `MDA-MB-231_MTX_treated_r1` will also pair with `MDA-MB-231_DMSO_untreated_r1` if no exact matching untreated sample exists. A warning is emitted when a fallback is used; if more than one candidate matches, the pipeline errors. Disable this behaviour with `--fuzzy_untreated_pairing false`, in which case unmatched groups proceed without a negative control.

For DMS experiments, a few defaults change automatically. `--rfnorm_reactive_bases` is set to `AC` (or `ACGU` when `pH ≥ 8`). `--rfnorm_dynamic_window` defaults to `50` when `pH < 8`; it controls the size of the sliding window used to compute local normalisation factors along the transcript — a smaller window is better suited to DMS at physiological pH where reactivity can vary sharply over short stretches. `--rfnorm_nan` defaults to `100` rather than `1000`; it sets the minimum number of reads required at a position for a reactivity value to be reported — positions with fewer reads are set to NaN instead of reporting a potentially unreliable values. These can all be overridden explicitly if needed.

`--rfnorm_mean_coverage` and `--rfnorm_median_coverage` set minimum coverage thresholds below which transcript positions are masked as NaN — useful for filtering out poorly-covered transcripts. `--rfnorm_chunk_size` (default `5000` on genome route, disabled on transcriptome route) splits the treated RC file into chunks of N transcripts and runs a separate rf-norm job per chunk in parallel, which can substantially reduce wall time for large genome-route experiments.

For the full list of available options see the [rf-norm documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

## rf-jackknife (optional)

[`rf-jackknife`](https://rnaframework-docs.readthedocs.io/en/latest/rf-jackknife/) calibrates rf-fold slope and intercept parameters against a set of reference structures. Enabled by `--jackknife_reference`. It runs between rf-norm and rf-fold, iterating over a slope/intercept grid and scoring each combination with the FMI (Fowlkes–Mallows Index). `rf-fold` is gated on jackknife completing for each group.

Two E. coli rRNA calibration references are bundled under `assets/ecoli_rrna_calibration/`:

- `ecoli_k12_rrna_reference_collab.db` — pseudoknots marked with `[]`; use with default `--rfjackknife_keep_pseudoknots true`
- `ecoli_k12_rrna_reference_crw.db` — CRW structures; all pairs in `()` notation; `--rfjackknife_keep_pseudoknots` has no effect

| Flag | Default | Description |
|------|---------|-------------|
| `--jackknife_reference` | — | Path to reference `.db` file — **required to enable** |
| `--rfjackknife_pool_all` | `true` | Pool all groups into one jackknife run; optimal slope/intercept passed to rf-fold automatically. Set `false` to run per group (slope/intercept must then be set manually via `--rffold_slope`/`--rffold_intercept`) |
| `--rfjackknife_slope` | `0,5` | Slope range (`min,max`) |
| `--rfjackknife_intercept` | `-3,0` | Intercept range (`min,max`) |
| `--rfjackknife_slope_step` | `0.2` | Slope grid increment |
| `--rfjackknife_intercept_step` | `0.2` | Intercept grid increment |
| `--rfjackknife_keep_pseudoknots` | `true` | Retain pseudoknotted base-pairs in reference (`-kp`) |
| `--rfjackknife_keep_lonelypairs` | `true` | Retain lonely base-pairs (`-kl`) |
| `--rfjackknife_mfmi` | `false` | Use modified FMI (`-m`) |
| `--rfjackknife_relaxed` | `false` | Relaxed FMI criteria (Deigan et al. 2009) (`-x`) |
| `--rfjackknife_img` | `false` | Generate FMI heatmap PDF (requires R) |
| `--rfjackknife_rf_fold_params` | `'-md 600'` | Extra flags passed to rf-fold inside jackknife (`-rp`) |

Output: `jackknife/<group>/` — FMI CSV per slope/intercept combination; optional heatmap PDF.

For the full list of available options see the [rf-jackknife documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-jackknife/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

## rf-fold

[`rf-fold`](https://rnaframework-docs.readthedocs.io/en/latest/rf-fold/) predicts RNA secondary structures from normalised reactivity profiles. The pipeline groups XMLs by `cell_line` (merging replicates) and folds them together. Dot-bracket output is the default; CT format is optional.

| Flag | Default | Description |
|------|---------|-------------|
| `--rffold_slope` | `4.6` | Reactivity-to-constraint slope (`-sl`); auto-overridden by jackknife when `--rfjackknife_pool_all true` |
| `--rffold_intercept` | `-2.2` | Reactivity-to-constraint intercept (`-in`); auto-overridden by jackknife |
| `--rffold_window` | `true` | Enable windowed folding (`-w`) |
| `--rffold_ct` | `false` | Write CT format alongside dot-bracket (`-ct`) |
| `--rffold_unconstrained` | `false` | Fold without reactivity constraints (`-u`) |
| `--rffold_vienna_no_lonely_pairs` | `false` | No-lonely-pairs mode (`-nlp`) |
| `--rffold_vienna_constrained` | `false` | Hard constraints mode (`-hc`) |
| `--rffold_vienna_max_bp_span` | — | Maximum base-pair span (`-md`) |
| `--rffold_only_common` | `false` | Keep only structures common across replicates (`-oc`) |
| `--rffold_fold_constraint_file` | — | External constraint file (`-c`) |
| `--rffold_dotplot` | `true` | Generate dot-plot output (`-dp`) |
| `--rffold_shannon_entropy` | `true` | Compute Shannon entropy (`-sh`); produces `fold/shannon_bw/` BigWigs |
| `--rffold_vienna_rnaplot` | `RNAplot` | RNAplot layout engine (`-vrp`) |

For the full list of available options see the [rf-fold documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-fold/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

## rf-eval (optional)

[`rf-eval`](https://rnaframework-docs.readthedocs.io/en/latest/rf-eval/) evaluates how well normalised reactivity profiles agree with a set of reference structures. Enabled by `--rfeval_reference`. It runs per `cell_line + replicate` group and reports three metrics per transcript:

- **Unpaired Coefficient** — fraction of highly reactive bases that are unpaired
- **DSCI** — probability that a randomly selected unpaired base has higher reactivity than a paired base
- **AUROC** — area under the ROC curve treating reactivity as a classifier of unpaired bases

| Flag | Default | Description |
|------|---------|-------------|
| `--rfeval_reference` | — | Path to reference `.db` file — **required to enable** |
| `--rfeval_reactivity_cutoff` | `0.7` | Reactivity threshold for unpaired classification (`-c`) |
| `--rfeval_img` | `false` | Generate metric plots — ROC curves, histograms (requires R) |
| `--rfeval_ignore_terminal` | `false` | Exclude terminal base-pairs from calculations (`-it`) |
| `--rfeval_keep_pseudoknots` | `false` | Retain pseudoknotted base-pairs (`-kp`) |
| `--rfeval_keep_lonelypairs` | `false` | Retain lonely base-pairs (`-kl`) |

Output: `eval/<group>/` — per-transcript CSV with Unpaired Coefficient, DSCI, and AUROC; optional PDF plots.

For the full list of available options see the [rf-eval documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-eval/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

## Visualisation

### 2D structure diagrams

After `rf-fold`, the pipeline draws reactivity-coloured 2D structure diagrams for each transcript. Nucleotides are coloured by normalised reactivity averaged across replicates.

- **[R2DT](https://github.com/RNAcentral/R2DT)** — used when a matching template exists in the R2DT library (rRNA, snRNA, tRNA, etc.). Produces layouts comparable across organisms. Published to `fold/<group>/structures/r2dt/`.
- **[ViennaRNA](https://www.tbi.univie.ac.at/RNA/)** — fallback for transcripts without an R2DT template. `RNAplot` draws an energy-minimised 2D diagram from the dot-bracket structure. Published to `fold/<group>/structures/viennarna/`.

Transcripts drawn by R2DT are recorded in `r2dt_drawn_ids.txt` so ViennaRNA does not duplicate them.

### Browser tracks

After normalisation, `rf-wiggle` converts each normalised XML to per-transcript WIG tracks. These are merged across transcripts and, when multiple replicates share the same cell line, averaged position-by-position before genomic remapping and BigWig conversion via the GTF.

Outputs:

- `norm/merged_bw/<cell_line>_reactivity.bw` — genomic reactivity BigWig
- `norm/transcript_bw/` — transcript-coordinate reactivity BigWigs
- `fold/shannon_bw/<cell_line>_shannon.bw` — genomic Shannon entropy BigWig (when `--rffold_shannon_entropy` is enabled, which is the default)
- `fold/transcript_shannon_bw/` — transcript-coordinate Shannon entropy BigWigs

## RDAT export

The pipeline produces [RDAT](https://rmdb.stanford.edu/tools/rdat_format/)-format files bundling normalised reactivity (from rf-norm XML) with the predicted structure (from rf-fold dot-bracket output). One `.rdat` file is written per transcript with both a reactivity profile and a predicted structure. Files are published under `fold/<group>/rdat/`.

## Outputs at a glance

Main output areas under `--outdir`:

- `fastqc/`, `cutadapt/`, `star/` (or `bowtie*/` with `--transcriptome`), `samtools*/` — preprocessing and alignment
- `count/` — count tables and plots
- `norm/<group>/` — normalised XML, normalisation plots, and per-transcript wiggle tracks
- `norm/merged_bw/` — per-cell-line genomic reactivity BigWigs
- `norm/transcript_bw/` — per-cell-line transcript-coordinate reactivity BigWigs
- `jackknife/<group>/` — FMI CSV and optional heatmap (only when `--jackknife_reference` is set)
- `fold/<group>/` — secondary structures, fold reports, optional CT, optional dotplots, and SVG diagrams
- `fold/merged_bp/` — merged base-pair files per cell line
- `fold/transcript_merged_bp/` — transcript-coordinate merged base-pair files
- `fold/shannon_bw/` — genomic Shannon entropy BigWigs
- `fold/transcript_shannon_bw/` — transcript-coordinate Shannon entropy BigWigs
- `eval/<group>/` — per-transcript evaluation metrics (only when `--rfeval_reference` is set)
- `multiqc/` — aggregated QC report

For full output details, see [output documentation](output.md).

## Running the pipeline

Typical usage (genome route, auto-downloads reference from Ensembl):

```bash
nextflow run nf-core/rnastructurome \
  --input samplesheet.csv \
  --outdir results \
  -profile docker
```

With local reference files:

```bash
nextflow run nf-core/rnastructurome \
  --input samplesheet.csv \
  --outdir results \
  --fasta genome.fa.gz \
  --gtf annotation.gtf.gz \
  -profile docker
```

Transcriptome route with jackknife calibration:

```bash
nextflow run nf-core/rnastructurome \
  --input samplesheet.csv \
  --outdir results \
  --transcriptome \
  --jackknife_reference assets/ecoli_rrna_calibration/ecoli_k12_rrna_reference_collab.db \
  -profile docker
```

For repeated runs with the same parameters, use a params file:

```bash
nextflow run nf-core/rnastructurome -profile docker -params-file params.yaml
```

## Core Nextflow arguments

These options are part of Nextflow itself and use a single hyphen. Pipeline parameters use a double hyphen.

### `-profile`

Use this parameter to choose configuration profiles. Common profiles: `docker`, `singularity`, `apptainer`, `podman`, `conda`. Multiple profiles can be combined (later profiles override earlier ones):

```bash
-profile test,docker
```

If `-profile` is not specified, the pipeline runs locally and expects all software on `PATH`.

### `-resume`

Resume from cached work where inputs and process configuration match prior runs:

```bash
nextflow run nf-core/rnastructurome -profile docker \
  --input samplesheet.csv --outdir results -resume
```

Use `nextflow log` to list previous run names; resume a specific run with `-resume <run-name>`.

### `-c`

Load a Nextflow config file for process resources, executors, or module-specific `ext.args`. Do not use `-c` for pipeline parameters (`input`, `outdir`, `fasta`, etc.) — use CLI flags or `-params-file` for those.

## Custom configuration

### Resource requests

Processes that fail with retryable resource-related exit codes are automatically retried with increased resources. Override per-process resources with a custom config:

```groovy title="custom_resources.config"
process {
    withName: 'STAR_GENOMEGENERATE' {
        cpus   = 16
        memory = '120.GB'
        time   = '24.h'
    }
}
```

### Custom tool arguments

Some modules expose `ext.args` for advanced flags not covered by pipeline parameters:

```groovy title="custom_args.config"
process {
    withName: 'STAR_ALIGN_RTSTOP' {
        ext.args = '--outFilterMultimapNmax 20'
    }
}
```

### nf-core/configs

For reusable institutional profiles, see [nf-core/configs](https://github.com/nf-core/configs) and the [nf-core configuration docs](https://nf-co.re/docs/usage/configuration).

## Updating and reproducibility

Update cached pipeline code:

```bash
nextflow pull nf-core/rnastructurome
```

Pin a release for reproducibility:

```bash
nextflow run nf-core/rnastructurome -r <VERSION> ...
```

The release, parameters, software versions, and execution trace are recorded under `pipeline_info/`. Use Nextflow `-bg` or `screen`/`tmux` for long sessions. Cap JVM memory if needed: `NXF_OPTS='-Xms1g -Xmx4g'`.

## Troubleshooting

### General debugging steps

1. Check `.nextflow.log` for the first error message.
2. Check `pipeline_info/execution_trace_<timestamp>.txt` to identify which process failed and its exit status.
3. Use `-resume` after fixing issues — all previously successful tasks will be cached and skipped.

### Missing outputs

- Confirm the process completed in the execution trace (`status = COMPLETED`).
- For `rf-fold` problems, inspect `fold/<group>/rffold.log` and check that input XMLs from `rf-norm` are non-empty.

### rf-jackknife: "0 imported" error

If `rf-jackknife` exits with `Error: No reference structure passed checks`, the most common cause is a sequence mismatch between the reference `.db` file and the mapping FASTA. To diagnose:

```bash
python3 -c "
import re
with open('path/to/transcript.xml') as f:
    c = f.read()
m = re.search(r'<sequence>(.*?)</sequence>', c, re.DOTALL)
print(m.group(1).strip().replace('\n','').replace('\t','').replace(' ',''))
" > /tmp/xml_seq.txt

sed -n '2p' path/to/reference.db > /tmp/db_seq.txt
diff /tmp/xml_seq.txt /tmp/db_seq.txt
```

Fix: rebuild the reference DB using sequences from the same FASTA the pipeline maps to, retaining the original dot-bracket annotation.

If any flag in `--rfjackknife_rf_fold_params` is unrecognised, jackknife exits with `Error: Invalid RF Fold parameters`. Note that `-x` (`--relaxed`) is a rf-jackknife flag — do not include it in `--rfjackknife_rf_fold_params`.

### rf-jackknife: ID mismatch

Transcript IDs in the `.db` file must match those in the XML files. RNA Framework derives IDs from the first whitespace-delimited token of each FASTA header — so `>16S_rRNA U00096.3:...` becomes `16S_rRNA`.

### Strandedness inference failures

If `RSeQC infer_experiment` is ambiguous for MaP samples on the genome route, supply `--rfcount_strandedness first|second|unstranded` to override.

### Container issues

If Singularity processes fail with locale or `TERM` errors, ensure `singularity.autoMounts = true` is set in your config.
