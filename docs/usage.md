# nf-core/rnastructurome: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnastructurome/usage](https://nf-co.re/rnastructurome/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files supplied with `-c` can be used for infrastructure settings such as resources, executors, containers, or module arguments, but should not be used for ordinary pipeline parameters.

## Required inputs

1. A samplesheet (`--input`) — see [Samplesheet input](#samplesheet-input) below.
2. Either an `organism` column in the samplesheet (triggers automatic Ensembl/NCBI download) or explicit `--genome_fasta` / `--gtf` flags.
3. A reference structure file (`--jackknife_reference`) — optional path to a `.db` file of known structures. When provided, `rf-jackknife` runs between `rf-norm` and `rf-fold` to calibrate slope/intercept parameters. When omitted, `rf-fold` runs directly using whatever slope/intercept params are configured.

## Alignment routes

### Default: genome alignment with STAR

The pipeline downloads a soft-masked genome assembly (`dna_sm.toplevel.fa.gz`) from Ensembl and the corresponding GTF, then builds a STAR index and aligns with splice-junction awareness. This always fetches the current Ensembl release (e.g. GRCh38 for human, not the older GRCh37/hg19).

For organisms not present in Ensembl (bacteria, viruses) the pipeline falls back to NCBI automatically using pre-configured accessions (see `conf/viral_genomes.config`, `conf/bacterial_genomes.config`) or an auto-search.

After alignment, genome-level counts from `rf-count-genome` are converted to transcript-level RC files by `rf-rctools extract` using the GTF before passing to `rf-norm`.

### Optional: transcriptome alignment with Bowtie

Add `--transcriptome` to use the older route: Ensembl cDNA + ncRNA FASTA → sorted transcript FASTA → Bowtie (RT-stop) / Bowtie2 (MaP).

```bash
nextflow run main.nf -profile docker \
  --input samplesheet.csv \
  --transcriptome \
  --outdir results
```

## Reference genome options

The minimum reference inputs for the default genome route are a genome FASTA and a GTF annotation. The pipeline can either download these automatically from Ensembl/NCBI using the `organism` column, or use files you provide explicitly.

### Explicit reference file specification

| Flag | Default | Description |
|------|---------|-------------|
| `--genome_fasta` | *(auto)* | Path to a local genome FASTA (skips Ensembl/NCBI download in STAR route). `--fasta` is accepted as a legacy alias. |
| `--transcriptome_fasta` | *(auto)* | Path to a local transcript FASTA (skips download in `--transcriptome` route). `--fasta` is also accepted here. |
| `--gtf` | *(auto)* | Path to a local GTF (skips Ensembl GTF download). |
| `--transcriptome` | `false` | Use transcriptome alignment route (Bowtie/Bowtie2) instead of genome (STAR). |

For example:

```bash
nextflow run nf-core/rnastructurome -profile docker \
  --input samplesheet.csv \
  --genome_fasta genome.fa.gz \
  --gtf annotation.gtf.gz \
  --outdir results
```

### Automatic reference resolution

The pipeline resolves a `reference_key` from the per-sample `organism` column or the global `--organism` flag.

- Latin binomials (`Homo sapiens`) are normalised to Ensembl format (`homo_sapiens`) automatically.
- Values already in `genus_species` format are lower-cased and used as-is.
- Shorthand values (`human`) must be resolved via `params.genomes` or `--ensembl_species_map`.

### Genome FASTA (STAR route, default)

1. `--genome_fasta` / `--fasta` if provided — used directly.
2. `params.genomes[reference_key].genome_fasta` if configured.
3. Ensembl — downloads `*.dna_sm.toplevel.fa.gz` for the species (tries main Ensembl, then EnsemblGenomes metazoa/fungi/plants/protists divisions).
4. NCBI fallback — for bacteria and viruses not present in Ensembl, using accessions from `conf/bacterial_genomes.config` or `conf/viral_genomes.config`, or via auto-search.

### Transcript FASTA (`--transcriptome` route)

1. `--transcriptome_fasta` / `--fasta` if provided — used directly.
2. `params.genomes[reference_key].transcript_fasta` (also `transcriptome` or `cdna`) if configured.
3. Ensembl — downloads `cdna.all.fa.gz` and, when present, `ncrna.fa.gz`, then merges them.  If `ncrna.fa.gz` is absent the pipeline continues with cDNA only.
4. NCBI fallback — as above.

### GTF annotation (both routes)

1. `--gtf` if provided — used directly.
2. `params.genomes[reference_key].gtf` if configured.
3. Ensembl — downloads the species GTF (prefers non-`abinitio` `.gtf.gz`).
4. For NCBI references — a synthetic single-exon GTF is generated automatically from the NCBI FASTA by `ncbi_gtf.py`.

Ensembl source can be tuned with:

- `--ensembl_release` (`current` by default; `latest` is treated the same, and values such as `114` or `release-114` are also accepted)
- `--ensembl_base_url` (default `https://ftp.ensembl.org/pub`)

### Reference reuse

When references are downloaded automatically, the sorted FASTA, GTF, and source URL provenance are published under `reference/`. Reuse these files in subsequent runs with `--genome_fasta` / `--transcriptome_fasta` and `--gtf` to avoid repeated downloads and to keep the reference version fixed.

## Samplesheet input

Create a samplesheet with information about the samples you would like to analyse, then pass it with:

```bash
--input '[path to samplesheet file]'
```

The samplesheet must be comma-separated, include a header row, and contain at least:

- `sample`
- `fastq_1`
- `cell_line`
- `condition`
- `replicate`

`fastq_2` is optional (leave empty for single-end).

### Multiple runs of the same sample

Use the same `sample` value when the same biological sample has been sequenced more than once, for example across multiple lanes. The pipeline concatenates all FASTQs for that sample before downstream analysis.

All rows for a repeated `sample` must have the same layout: either all single-end or all paired-end.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,cell_line,condition,replicate,organism,principle
HEK293T_treated_r1,lane1_R1.fastq.gz,,HEK293T,treated,1,Homo sapiens,RT-stop
HEK293T_treated_r1,lane2_R1.fastq.gz,,HEK293T,treated,1,Homo sapiens,RT-stop
HEK293T_untreated_r1,lane1_untreated_R1.fastq.gz,,HEK293T,untreated,1,Homo sapiens,RT-stop
```

### Minimal samplesheet

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,cell_line,condition,replicate,organism,principle
HEK293T_treated_r1,/data/treated_r1.fastq.gz,,HEK293T,treated,1,Homo sapiens,RT-stop
HEK293T_untreated_r1,/data/untreated_r1.fastq.gz,,HEK293T,untreated,1,Homo sapiens,RT-stop
```

### Full samplesheet

The samplesheet can include optional per-sample metadata columns. A more complete example with both single-end and paired-end rows is shown below.

```csv title="samplesheet.csv"
sample,sample_id,fastq_1,fastq_2,method,principle,cell_line,condition,replicate,organism,pH,adapter_3p,adapter_5p,umi_pattern
HEK293T_treated_r1,GSM000001,/data/treated_r1.fastq.gz,,SHAPE,RT-stop,HEK293T,treated,1,Homo sapiens,7.5,,,
HEK293T_untreated_r1,GSM000002,/data/untreated_r1.fastq.gz,,SHAPE,RT-stop,HEK293T,untreated,1,Homo sapiens,7.5,,,
HEK293T_map_r1,GSM000003,/data/map_r1_1.fastq.gz,/data/map_r1_2.fastq.gz,DMS,MaP,HEK293T,treated,1,Homo sapiens,8.0,AGATCGGAAGAGC,AGATCGGAAGAGC,NNNNNN
```

### Column reference

| Column      | Description                                                                                      |
| ----------- | ------------------------------------------------------------------------------------------------ |
| `sample`    | Sample identifier. Re-used IDs are concatenated before downstream analysis (multi-lane support). |
| `fastq_1`   | Read 1 FASTQ path (`.fastq.gz` / `.fq.gz`).                                                      |
| `fastq_2`   | Read 2 FASTQ path for paired-end data (optional).                                                |
| `cell_line` | Group key used for control pairing in `rf-norm`.                                                 |
| `condition` | One of `treated`, `untreated`, `denatured`.                                                      |
| `replicate` | Replicate key used for control pairing in `rf-norm`.                                             |

Optional per-sample columns supported by the pipeline include:

- `sample_id` (external accession or display identifier; falls back to global `--sample_id` when set)
- `method` (`SHAPE` or `DMS`; controls chemistry-specific defaults)
- `organism` (used for automatic genome/transcriptome download from Ensembl/NCBI; e.g. `Homo sapiens`)
- `principle` (`RT-stop` or `MaP`) — **required**; controls chemistry-specific alignment parameters, rf-count mutation counting, and rf-norm scoring/normalisation defaults
- `pH` (used for DMS reactive-base defaults; `pH >= 8.0` treats all bases as potentially reactive)
- `adapter_5p`
- `adapter_3p`
- `umi_pattern`

An [example samplesheet](../assets/samplesheet.csv) is provided.

## Condition pairing logic (`rf-norm`)

Samples are grouped by identical `cell_line + replicate`.

- `treated` can be analyzed alone.
- `untreated` requires a matching `treated` sample.
- `denatured` requires matching `treated` and `untreated` samples.

Default `rf-norm` scoring/normalization is selected automatically based on probing principle and available controls:

- RT-stop + untreated: Ding (`-sm 1`) + Box-plot (`-nm 3`)
- RT-stop only: Rouskin (`-sm 2`) + 90% Winsorizing (`-nm 2`)
- MaP + untreated (and optional denatured): Siegfried (`-sm 3`) + Box-plot (`-nm 3`)
- MaP only: Zubradt (`-sm 4`) + Box-plot (`-nm 3`)

Advanced override:

- `--rfnorm_norm_method 2|3` forces the `rf-norm` normalisation mode while leaving scoring-method auto-selection unchanged.

### Fuzzy untreated pairing

By default, if a treated group has no exact `cell_line + replicate` untreated match, the pipeline falls back to any untreated sample that shares the same **cell_line base token** (the portion before the first underscore) at the same replicate. For example, `MDA-MB-231_MTX_treated_r1` will automatically pair with an `MDA-MB-231_untreated_r1` control even though the `cell_line` values differ (`MDA-MB-231_MTX` vs `MDA-MB-231`).

The pipeline warns when a fallback is used:

```
[WARN] No exact untreated match for 'MDA-MB-231_MTX_r1' — falling back to 'MDA-MB-231_r1' (shared cell_line base token at same replicate).
```

If more than one untreated group matches the base token at the same replicate, the pipeline errors rather than choosing arbitrarily. To disable fuzzy matching entirely and require exact `cell_line + replicate` pairing, set:

```bash
--fuzzy_untreated_pairing false
```

When disabled, treated groups with no exact untreated match proceed without a negative control (scoring method 2 or 4 instead of 1 or 3).

## Adapter trimming options

The pipeline uses Cutadapt for quality and adapter trimming. Trimming defaults depend on the probing principle because RT-stop experiments encode signal at the read end.

| Principle | `--cutadapt_5quality` | `--cutadapt_3quality` | Notes                                         |
| --------- | --------------------- | --------------------- | --------------------------------------------- |
| `RT-stop` | forced to `0`         | default `20`          | 5' quality trimming is intentionally disabled |
| `MaP`     | default `20`          | default `20`          | both are user-configurable                    |

Adapter precedence:

1. Per-sample (`adapter_5p`, `adapter_3p`)
2. Global (`--cutadapt_adapter_5p`, `--cutadapt_adapter_3p`)

If no adapters are provided, the pipeline falls back to `AGATCGGAAGAGC` for both the 5' and 3' adapter. Set `--cutadapt_quality_only` to skip adapter trimming entirely and perform quality/length filtering only.

## Unique Molecular Identifiers (UMI)

Enable with `umi_pattern` (sample-level or global `--umi_pattern`).

- Patterns containing only `N/C/X` use direct `umi_tools --bc-pattern` mode.
- IUPAC patterns (for example with `D`) are converted automatically to regex mode.

UMI extraction runs before Cutadapt so that molecular tags are available during deduplication.

### Duplicate handling

- UMI-tagged samples use `umi_tools dedup`.
- Non-UMI samples use `samtools markdup` by default.

## RNAframework modules

### `rf-count`

- Plot generation is always enabled (`-g -R`).
- FASTA is always passed from pipeline reference resolution.
- MaP-specific mutation options are applied only for `principle=MaP`.

Key parameters:

- `--rfcount_trim_5prime` (`-t5`)
- `--rfcount_mask_file` (`--mask-file`)
- `--rfcount_primary_only` (`--primary-only`)
- `--rfcount_paired_only` (`--paired-only`)
- `--rfcount_properly_paired` (`--properly-paired`)
- `--rfcount_map_sort_by_read_name` (`--sort-by-read-name`, MaP)
- `--rfcount_map_discard_shorter` (`--discard-shorter`, MaP)
- `--rfcount_map_min_quality` (`--min-quality`, MaP)
- `--rfcount_map_collapse_consecutive` (`--collapse-consecutive`, MaP)
- `--rfcount_map_max_collapse_distance` (`--max-collapse-distance`, MaP)

Paired-end default:

- If neither paired filter is set, pipeline defaults to `--properly-paired`.

#### Genome route: strandedness

On the genome route, `rf-count-genome` requires a library strandedness value (`--library-strandedness`) to correctly assign reads to the plus or minus strand. Getting this wrong puts counts on the wrong strand and produces incorrect reactivity profiles.

- **RT-stop**: always treated as second-strand (`second`) — this is fixed by experimental design and requires no user input.
- **MaP**: strandedness is inferred automatically per sample using `RSeQC infer_experiment`. If inference is ambiguous (very small BAM, missing BED annotation), you can override it with `--rfcount_strandedness` (`first`, `second`, or `unstranded`).

#### Genome route: `rf-rctools extract`

In the default STAR/genome route, `rf-count-genome` produces genome-coordinate `.rc` files. Before `rf-norm`, `rf-rctools extract` converts these to transcript-level RC files using the GTF. This step runs automatically and requires no user input.

### `rf-norm`

- Plot generation is always enabled (`--img -R`).

Key parameters:

- `--rfnorm_remap_reactivities`
- `--rfnorm_reactive_bases`
- `--rfnorm_norm_window`
- `--rfnorm_window_offset`
- `--rfnorm_dynamic_window`
- `--rfnorm_norm_independent`
- `--rfnorm_norm_factor`
- `--rfnorm_norm_method`
- `--rfnorm_raw`
- `--rfnorm_pseudocount`
- `--rfnorm_max_score`
- `--rfnorm_ignore_lower_than_untreated`
- `--rfnorm_max_untreated_mut`
- `--rfnorm_max_mutation_rate`
- `--rfnorm_mean_coverage`
- `--rfnorm_median_coverage`
- `--rfnorm_nan`

Notes:

- For `method=DMS`, pipeline defaults `--rfnorm_reactive_bases` to `AC`, or `ACGU` when `pH >= 8`. It defaults `--dynamic-window` to `50` only when `pH < 8`.
- `--rfnorm_nan` defaults to `100` for DMS and `1000` for all other methods. Override with an explicit value to apply the same threshold regardless of method.
- `R` must exist at `--rnaframework_r_path` in the active runtime environment. The default is `/usr/bin/R`.

### `rf-wiggle`

After normalization, `rf-wiggle` converts each normalized XML to per-transcript WIG tracks. These are merged across transcripts and, when multiple replicates share the same cell line, averaged position-by-position before genomic remapping and BigWig conversion via the GTF, producing a single `norm/merged_bw/<cell_line>_reactivity.bw`.

### `rf-jackknife`

[`rf-jackknife`](https://rnaframework-docs.readthedocs.io/en/latest/rf-jackknife/) runs between `rf-norm` and `rf-fold` **only when `--jackknife_reference` is provided**. It iteratively folds a set of known reference structures across a grid of slope/intercept values and scores each combination with the FMI (Fowlkes–Mallows Index), allowing you to identify optimal normalization parameters.

`rf-fold` only starts for a fold group after `rf-jackknife` has completed successfully for that group (when it is enabled).

Required parameter:

| Flag | Description |
|------|-------------|
| `--jackknife_reference` | Path to a reference `.db` structure file — **required to enable rf-jackknife** |

Optional tuning parameters:

| Flag | Default | Description |
|------|---------|-------------|
| `--rfjackknife_slope` | `0,5` | Slope range to test (`min,max`) |
| `--rfjackknife_intercept` | `-3,0` | Intercept range to test (`min,max`) |
| `--rfjackknife_slope_step` | `0.2` | Slope grid increment |
| `--rfjackknife_intercept_step` | `0.2` | Intercept grid increment |
| `--rfjackknife_keep_pseudoknots` | `true` | Retain pseudoknotted base-pairs (`[]` notation) in the reference structure before comparison (`-kp`). Enable when the reference `.db` file marks pseudoknots explicitly. Has no effect on references that use only `()` notation (e.g. CRW). |
| `--rfjackknife_keep_lonelypairs` | `true` | Retain lonely base-pairs (1 bp helices) in the reference structure (`-kl`) |
| `--rfjackknife_mfmi` | `false` | Use modified FMI instead of standard FMI (`-m`) |
| `--rfjackknife_relaxed` | `false` | Use relaxed FMI criteria (Deigan et al. 2009) (`-x`) |
| `--rfjackknife_img` | `false` | Generate FMI heatmap PDF (requires R) |
| `--rfjackknife_rf_fold_params` | `'-md 600'` | Additional parameters passed to rf-fold inside jackknife (`-rp`) |
| `--rfjackknife_pool_all` | `true` | Pool XMLs from all cell-line groups into a single jackknife run (output: `jackknife/all_groups/`). The optimal slope/intercept is automatically passed to `rf-fold` in the same run. Set to `false` to run one jackknife per group instead (slope/intercept must then be supplied manually via `--rffold_slope`/`--rffold_intercept`). |

Two E. coli rRNA calibration references are bundled under `assets/ecoli_rrna_calibration/`:

- `ecoli_k12_rrna_reference_collab.db` — structures from the RNA Framework author's current reference set; pseudoknots marked with `[]`, use with default `--rfjackknife_keep_pseudoknots true`
- `ecoli_k12_rrna_reference_crw.db` — structures from the Comparative RNA Website (CRW); all pairs in `()` notation, no pseudoknot distinction; `--rfjackknife_keep_pseudoknots` has no effect

Both use sequences extracted from `ecoli_k12_rrna.fa` (the same FASTA the calibration samples are mapped to) so no sequence mismatch errors occur.

Output: `jackknife/<group>/` — CSV of FMI scores per slope/intercept combination; optional heatmap PDF.

### `rf-fold`

Default behavior:

- Dot-bracket output is default (CT is optional).
- Windowed folding is enabled by default (`-w`); disable with `--rffold_window false`.
- Graphical fold reports are always enabled (`-g -R`).
- RNAplot overlays are enabled via `-vrp` (default `RNAplot`).
- Dotplot generation is controlled by `--rffold_dotplot` (default `true`).

Supported pipeline options and mapped flags:

- `--rffold_ct` -> `-ct`
- `--rffold_window` -> `-w`
- `--rffold_unconstrained` -> `-i`
- `--rffold_vienna_no_lonely_pairs` -> `-nlp`
- `--rffold_vienna_constrained` -> `-hc`
- `--rffold_vienna_max_bp_span` -> `-md`
- `--rffold_only_common` -> `-oc`
- `--rffold_fold_constraint_file` -> `-c`
- `--rffold_dotplot` -> `-dp`
- `--rffold_shannon_entropy` -> `-sh`
- `--rffold_vienna_rnaplot` -> `-vrp`
- `--rffold_slope` -> `-sl` (slope for reactivity-to-constraint conversion; default `4.6`; overridden automatically when `--rfjackknife_pool_all` is used)
- `--rffold_intercept` -> `-in` (intercept for reactivity-to-constraint conversion; default `-2.2`; overridden automatically when `--rfjackknife_pool_all` is used)

Replicate handling:

- The pipeline groups normalized XMLs by `cell_line` and folds them together, so replicates are folded as a combined set.
- Single-XML groups run as standard one-sample `rf-fold`.

BigWig outputs from `rf-fold`:

- Shannon entropy BigWig (`fold/shannon_bw/<cell_line>_shannon.bw`) is produced when `--rffold_shannon_entropy` is enabled (default: `true`). Per-transcript WIG files from `rf-fold` are merged across all transcripts, remapped through the GTF to genomic coordinates, and converted to BigWig format.

### RNAframework runtime settings

- `--rnaframework_container`: RNAframework image for local modules when using container-based profiles
- `--rnaframework_r_path`: `R` binary used by RNAframework plotting (default: `/usr/bin/R`)

## 2D structure visualisation

After `rf-fold`, the pipeline draws coloured 2D structure diagrams for every transcript that passes coverage filters.

### R2DT (template-matched diagrams)

[R2DT](https://github.com/RNAcentral/R2DT) draws structures using a curated library of templates derived from known RNA families (rRNA, snRNA, tRNA, etc.). When a template exists for a transcript, R2DT produces a layout that is directly comparable across organisms and studies.

- Diagrams are generated for all transcripts in the fold output for which R2DT finds a matching template.
- Nucleotides are coloured by normalised SHAPE/DMS reactivity averaged across replicates.
- Transcripts drawn by R2DT are recorded in `r2dt_drawn_ids.txt` so ViennaRNA does not duplicate them.

### ViennaRNA (fallback diagrams)

For transcripts without an R2DT template, the pipeline uses [`ViennaRNA`](https://www.tbi.univie.ac.at/RNA/) RNAplot to draw an energy-minimised 2D diagram from the dot-bracket structure produced by `rf-fold`.

- All transcripts not covered by R2DT receive a ViennaRNA diagram.
- Nucleotides are coloured with the same reactivity scale as R2DT diagrams.

R2DT diagrams are published to `fold/<group>/structures/r2dt/`; ViennaRNA fallback diagrams are published to `fold/<group>/structures/viennarna/`.

## Optional analysis modules

These modules are not run automatically. They are invoked by passing the appropriate flag and require user-supplied reference data.

### rf-eval — benchmark reactivities against known structures

[`rf-eval`](https://rnaframework-docs.readthedocs.io/en/latest/rf-eval/) evaluates how well the pipeline's normalised reactivity profiles agree with a set of experimentally validated or published secondary structures. It reports three metrics per transcript:

- **Unpaired Coefficient** — fraction of highly reactive bases that are unpaired
- **DSCI** — probability that a randomly selected unpaired base has higher reactivity than a paired base
- **AUROC** — area under the ROC curve treating reactivity as a classifier of unpaired bases

Key parameters:

| Flag | Default | Description |
|------|---------|-------------|
| `--rfeval_reactivity_cutoff` | `0.7` | Reactivity threshold for unpaired classification |
| `--rfeval_img` | `false` | Generate metric plots — distributions, ROC curves, histograms (requires R) |
| `--rfeval_ignore_terminal` | `false` | Exclude terminal base pairs from calculations |

Output: `eval/<group>/` — CSV with per-transcript Unpaired Coefficient, DSCI, and AUROC; optional PDF plots.

## RDAT export

The pipeline produces [RDAT](https://rmdb.stanford.edu/tools/rdat_format/)-format files that bundle per-transcript reactivity profiles (from the normalised XML) with the predicted secondary structure (from `rf-fold` dot-bracket output). RDAT is a community standard for depositing structure probing data in the RNA Mapping Database (RMDB).

- One `.rdat` file is written per transcript that has both a normalised reactivity profile and a predicted structure.
- Files are published under `fold/<group>/rdat/`.
- Use `--rnaframework_container` to point to the RNAFramework container image if running without a network connection to pull images automatically.

## Outputs at a glance

Main output areas under `--outdir`:

- `fastqc/`, `cutadapt/`, `star/` (or `bowtie*/` with `--transcriptome`), `samtools*/` — preprocessing and alignment
- `count/` — count tables and always-on plots
- `norm/<group>/` — normalized XML, normalization plots, and per-transcript wiggle tracks
- `norm/merged_bw/` — per-cell-line genomic reactivity BigWigs (`<cell_line>_reactivity.bw`); reactivity values are averaged across replicates before genomic remapping and conversion when multiple replicates are available
- `norm/transcript_bw/` — per-cell-line transcript-coordinate reactivity BigWigs
- `jackknife/<group>/` — FMI CSV and optional heatmap (only when `--jackknife_reference` is set)
- `fold/<group>/` — inferred secondary structures (dot-bracket by default), fold reports, optional CT, optional dotplots, and `structures/r2dt/` or `structures/viennarna/` SVG diagrams
- `fold/merged_bp/` — merged base-pair files per cell line (produced from dotplots when `--rffold_dotplot` is enabled)
- `fold/transcript_merged_bp/` — transcript-coordinate merged base-pair files
- `fold/shannon_bw/` — per-cell-line genomic Shannon entropy BigWigs (`<cell_line>_shannon.bw`; produced when `--rffold_shannon_entropy` is enabled, which is the default)
- `fold/transcript_shannon_bw/` — transcript-coordinate Shannon entropy BigWigs
- `eval/<group>/` — per-transcript evaluation metrics CSV and optional plots (when rf-eval is enabled)
- `multiqc/` — final aggregated QC report

For full output details, see [output documentation](output.md).

## Running the pipeline

Typical usage (genome route, auto-downloads reference from Ensembl):

```bash
nextflow run nf-core/rnastructurome \
  --input ./samplesheet.csv \
  --outdir ./results \
  -profile docker
```

With local genome files:

```bash
nextflow run nf-core/rnastructurome \
  --input ./samplesheet.csv \
  --outdir ./results \
  --jackknife_reference ./known_structures.db \
  --genome_fasta ./genome.fa.gz \
  --gtf ./annotation.gtf.gz \
  -profile docker
```

Transcriptome route (Bowtie/Bowtie2):

```bash
nextflow run nf-core/rnastructurome \
  --input ./samplesheet.csv \
  --outdir ./results \
  --transcriptome \
  -profile docker
```

The pipeline creates:

```bash
work                # Nextflow working directory
<OUTDIR>            # Published results
.nextflow.log       # Nextflow run log
```

If you repeatedly use the same parameters, put them in a params file and run with `-params-file`:

```bash
nextflow run nf-core/rnastructurome -profile docker -params-file params.yaml
```

## Core Nextflow arguments

These options are part of Nextflow itself and use a single hyphen. Pipeline parameters use a double hyphen.

### `-profile`

Use this parameter to choose configuration profiles for software packaging and execution backends.

Common software profiles include:

- `docker`
- `singularity`
- `apptainer`
- `podman`
- `conda`

Cluster or institution-specific profiles can also be loaded, either from this repository or from `nf-core/configs`.

Multiple profiles can be combined:

```bash
-profile test,docker
```

Order matters: later profiles override earlier ones. If `-profile` is not specified, the pipeline runs locally and expects all software to be available on `PATH`, which is not recommended for reproducible analyses.

### `-resume`

Resume from cached work where inputs and process configuration match prior runs.

```bash
nextflow run nf-core/rnastructurome -profile docker \
  --input samplesheet.csv \
  --outdir results \
  -resume
```

You can also resume a specific run name: `-resume <run-name>`. Use `nextflow log` to list previous run names.

### `-c`

Use `-c` to load a Nextflow config file for process resources, executors, containers, publish behaviour, or module-specific `ext.args`.

Do not use `-c` for routine pipeline parameters such as `input`, `outdir`, `genome_fasta`, or `gtf`; use CLI flags or `-params-file` for those.

## Updating the pipeline

Update cached pipeline code:

```bash
nextflow pull nf-core/rnastructurome
```

When running from a local checkout, update the checkout with Git instead.

## Reproducibility

Pin a release for reproducibility:

```bash
nextflow run nf-core/rnastructurome -r <VERSION> ...
```

The release, parameters, software versions, and execution trace are recorded under `pipeline_info/`. For published or shared analyses, keep the `params_<timestamp>.json`, trace file, MultiQC report, and exact samplesheet alongside the final outputs.

## Custom configuration

### Resource requests

Default CPU, memory, and time requests are defined in the pipeline config. Processes that fail with retryable resource-related exit codes are automatically retried with increased resources.

Use a custom Nextflow config for local resource tuning:

```groovy title="custom_resources.config"
process {
    withName: 'STAR_GENOMEGENERATE' {
        cpus   = 16
        memory = '120.GB'
        time   = '24.h'
    }
}
```

Run with:

```bash
nextflow run nf-core/rnastructurome -profile docker -c custom_resources.config ...
```

### Custom containers

Container and conda environments are defined per process. Override them only when you need a patched or site-specific image, and keep overrides scoped to the affected process.

```groovy title="custom_container.config"
process {
    withName: 'RNAFRAMEWORK_RFCOUNT_GENOME' {
        container = 'docker.io/dincarnato/rnaframework:2.9.6'
    }
}
```

### Custom tool arguments

Some modules expose `ext.args` for advanced tool flags that are not regular pipeline parameters. Use this sparingly and scope it to the specific process.

```groovy title="custom_args.config"
process {
    withName: 'STAR_ALIGN_RTSTOP' {
        ext.args = '--outFilterMultimapNmax 20'
    }
}
```

### nf-core/configs

The pipeline can load institutional profiles from `nf-core/configs`. If your organisation needs a reusable shared profile, test it locally with `-c` first, then consider contributing it to `nf-core/configs`.

For more details, see:

- [nf-core configuration docs](https://nf-co.re/docs/usage/configuration)
- [nf-core/configs](https://github.com/nf-core/configs)
- [Nextflow config docs](https://www.nextflow.io/docs/latest/config.html)

## Running in the background

Use Nextflow `-bg`, or run inside `screen`/`tmux` for long sessions.

## Nextflow memory requirements

If needed, cap Nextflow JVM memory:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Troubleshooting

### General debugging steps

1. Check `.nextflow.log` for the first error message — Nextflow prints the root cause before cascading failures.
2. Check `pipeline_info/execution_trace_<timestamp>.txt` to identify which process failed and its exit status.
3. Use `-resume` after fixing issues — all previously successful tasks will be cached and skipped.

### Missing outputs

If expected output files are absent:

- Confirm the process completed successfully in the execution trace (`status = COMPLETED`).
- For processes that ran in a cluster scratch directory (`scratch = true` in your profile), intermediate inputs are not persisted to the work directory — only declared outputs are copied back.
- For `rf-fold` problems, inspect `fold/<group>/rffold.log` and check that the input XMLs from `rf-norm` are non-empty.

### rf-jackknife: "0 imported" error

If `rf-jackknife` exits with `Error: No reference structure passed checks`, the import step failed. Common causes:

**Sequence mismatch between reference DB and mapping reference**

The reference `.db` file supplied via `--jackknife_reference` contains sequences that must match the sequences in the genome/transcriptome FASTA used for mapping. If the reference DB was built from a different genome assembly or annotation source (e.g. CRW vs. Ensembl), the per-nucleotide sequence comparison will fail even when the transcript IDs match.

To diagnose, extract the sequence from a failed sample's XML and diff it against the corresponding entry in the `.db` file:

```bash
# Get the sequence from a failed rf-norm XML
python3 -c "
import re
with open('path/to/transcript.xml') as f:
    c = f.read()
m = re.search(r'<sequence>(.*?)</sequence>', c, re.DOTALL)
print(m.group(1).strip().replace('\n','').replace('\t','').replace(' ',''))
" > /tmp/xml_seq.txt

# Get the sequence from the reference DB (e.g. line 2 for the first entry)
sed -n '2p' path/to/reference.db > /tmp/db_seq.txt

diff /tmp/xml_seq.txt /tmp/db_seq.txt
```

The fix is to rebuild the reference DB using sequences extracted from the same FASTA the pipeline maps to, while retaining the original dot-bracket structure annotation (valid as long as sequence lengths are identical).

**Invalid rf-fold parameters passed via `--rfjackknife_rf_fold_params`**

rf-jackknife validates the rf-fold parameters before starting. If any flag is unrecognised by the version of RNA Framework in the container, you will see `Error: Invalid RF Fold parameters`. Check which flags are supported by running:

```bash
singularity exec <container.img> rf-fold --help 2>&1 | grep -E "^\s+\-"
```

Note that `-x` (`--relaxed`) is a **rf-jackknife** flag, not an rf-fold flag — it should not be included in `--rfjackknife_rf_fold_params`.

### rf-jackknife: ID mismatch warning

If the transcript IDs in your reference `.db` file do not match the IDs in the XML files, no structures will be imported. RNA Framework derives transcript IDs from the first whitespace-delimited token of each FASTA header — so `>16S_rRNA U00096.3:...` becomes `16S_rRNA`. Ensure the entries in your `.db` file use the same short IDs.

### Strandedness inference failures (genome route)

On the genome route, the pipeline runs `RSeQC infer_experiment` to detect strandedness for MaP samples. If the BAM is very small or the BED annotation is missing, inference can be ambiguous. Supply `--rfcount_strandedness` explicitly to set the fallback value used by `rf-count-genome`:

```bash
--rfcount_strandedness first   # or second, unstranded
```

### Container / singularity issues

If running with Singularity and processes fail with locale or `TERM` errors, ensure `singularity.autoMounts = true` is set in your config. The pipeline exports `TERM` inside each task script to work around missing terminal environment variables inside containers.
