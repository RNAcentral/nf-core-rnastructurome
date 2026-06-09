# nf-core/rnastructurome: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnastructurome/usage](https://nf-co.re/rnastructurome/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

`nf-core/rnastructurome` processes structure-probing RNA-seq data from FASTQ to normalized reactivities and secondary structure outputs.

High-level workflow:

1. Input QC and trimming (`FastQC`, `cutadapt`, optional `umi_tools extract`)
2. Reference download — genome FASTA + GTF from Ensembl (NCBI fallback for bacteria/viruses)
3. Alignment (`STAR`) and BAM post-processing (`samtools`, optional `umi_tools dedup`)
4. Reactivity counting (`rf-count`) for RT-stop and MaP
5. Reactivity normalization (`rf-norm`) using available controls
6. Structure inference (`rf-fold`) from normalized XMLs
7. Aggregated reporting (`MultiQC`)

The probing principle (`RT-stop` or `MaP`) is read from the `principle` column in the samplesheet CSV and controls chemistry-specific trimming, alignment and normalization parameters automatically.

## Quick start

```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheet.csv \
  --outdir results
```

The pipeline downloads the genome FASTA and GTF automatically from Ensembl for each organism in the samplesheet.

Supply your own genome files to skip the download:

```bash
nextflow run main.nf -profile docker \
  --input samplesheet.csv \
  --genome_fasta genome.fa.gz \
  --gtf annotation.gtf.gz \
  --outdir results
```

Resume after an interruption:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results -resume
```

## Required inputs

1. A samplesheet (`--input`) — see [Samplesheet input](#samplesheet-input) below.
2. Either an `organism` column in the samplesheet (triggers automatic Ensembl/NCBI download) or explicit `--genome_fasta` / `--gtf` flags.

## Alignment routes

### Default: genome alignment with STAR

The pipeline downloads a soft-masked genome assembly (`dna_sm.toplevel.fa.gz`) from Ensembl and the corresponding GTF, then builds a STAR index and aligns with splice-junction awareness. This always fetches the current Ensembl release (e.g. GRCh38 for human, not the older GRCh37/hg19).

For organisms not present in Ensembl (bacteria, viruses) the pipeline falls back to NCBI automatically using pre-configured accessions (see `conf/viral_genomes.config`, `conf/bacterial_genomes.config`) or an auto-search.

### Optional: transcriptome alignment with Bowtie

Add `--transcriptome` to use the older route: Ensembl cDNA + ncRNA FASTA → sorted transcript FASTA → Bowtie (RT-stop) / Bowtie2 (MaP).

```bash
nextflow run main.nf -profile docker \
  --input samplesheet.csv \
  --transcriptome \
  --outdir results
```

## Reference input flags

| Flag | Default | Description |
|------|---------|-------------|
| `--genome_fasta` | *(auto)* | Path to a local genome FASTA (skips Ensembl/NCBI download in STAR route). `--fasta` is accepted as a legacy alias. |
| `--transcriptome_fasta` | *(auto)* | Path to a local transcript FASTA (skips download in `--transcriptome` route). `--fasta` is also accepted here. |
| `--gtf` | *(auto)* | Path to a local GTF (skips Ensembl GTF download). |
| `--transcriptome` | `false` | Use transcriptome alignment route (Bowtie/Bowtie2) instead of genome (STAR). |

## Reference resolution

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

## Samplesheet input

The samplesheet must be CSV and include:

- `sample`
- `fastq_1`
- `cell_line`
- `condition`
- `replicate`

`fastq_2` is optional (leave empty for single-end).

### Minimal example

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,cell_line,condition,replicate
HEK293T_treated_r1,/data/treated_r1.fastq.gz,,HEK293T,treated,1
HEK293T_untreated_r1,/data/untreated_r1.fastq.gz,,HEK293T,untreated,1
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

- `organism` (used for automatic genome/transcriptome download from Ensembl/NCBI; e.g. `Homo sapiens`)
- `principle` (`RT-stop` or `MaP`) — **required**; controls chemistry-specific alignment parameters, rf-count mutation counting, and rf-norm scoring/normalisation defaults
- `adapter_5p`
- `adapter_3p`
- `umi_pattern`
- `principle` (`RT-stop` or `MaP`, if not already set upstream)

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

## Preprocessing behavior

### Cutadapt by principle

| Principle | `--cutadapt-5quality` | `--cutadapt-3quality` | Notes                                         |
| --------- | --------------------- | --------------------- | --------------------------------------------- |
| `RT-stop` | forced to `0`         | default `20`          | 5' quality trimming is intentionally disabled |
| `MaP`     | default `20`          | default `20`          | both are user-configurable                    |

Adapter precedence:

1. Per-sample (`adapter_5p`, `adapter_3p`)
2. Global (`--cutadapt_adapter_5p`, `--cutadapt_adapter_3p`)

If no adapters are provided, the pipeline falls back to `AGATCGGAAGAGC` for both the 5' and 3' adapter. Set `--cutadapt_quality_only` to skip adapter trimming entirely and perform quality/length filtering only.

### Optional UMI extraction

Enable with `umi_pattern` (sample-level or global `--umi_pattern`).

- Patterns containing only `N/C/X` use direct `umi_tools --bc-pattern` mode.
- IUPAC patterns (for example with `D`) are converted automatically to regex mode.

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
- `R` must exist at `--rnaframework_r_path` in the active runtime environment. The default is `/usr/bin/R`.
- After normalization, per-transcript WIG files are produced by `rf-wiggle` and merged across transcripts. If multiple replicates share the same cell line, their merged tracks are averaged position-by-position before BigWig conversion, producing a single `norm/merged_bw/<cell_line>_reactivity.bw`.

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

Both diagram types are published to `fold/<group>/2D-structures/`.

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
- `fold/<group>/` — inferred secondary structures (dot-bracket by default), fold reports, optional CT, optional dotplots
- `fold/merged_bp/` — merged base-pair files per cell line (produced from dotplots when `--rffold_dotplot` is enabled)
- `fold/shannon_bw/` — per-cell-line genomic Shannon entropy BigWigs (`<cell_line>_shannon.bw`; produced when `--rffold_shannon_entropy` is enabled, which is the default)
- `multiqc/` — final aggregated QC report

For full output details, see [output documentation](output.md).

## Running the pipeline

Typical usage (genome route, auto-downloads reference from Ensembl):

```bash
nextflow run nf-core/rnastructurome --input ./samplesheet.csv --outdir ./results -profile docker
```

With local genome files:

```bash
nextflow run nf-core/rnastructurome --input ./samplesheet.csv --outdir ./results \
  --genome_fasta ./genome.fa.gz --gtf ./annotation.gtf.gz -profile docker
```

Transcriptome route (Bowtie/Bowtie2):

```bash
nextflow run nf-core/rnastructurome --input ./samplesheet.csv --outdir ./results \
  --transcriptome -profile docker
```

The pipeline creates:

```bash
work                # Nextflow working directory
<OUTDIR>            # Published results
.nextflow.log       # Nextflow run log
```

For repeated runs, use a params file:

```bash
nextflow run nf-core/rnastructurome -profile docker -params-file params.yaml
```

```yaml title="params.yaml"
input: "./samplesheet.csv"
outdir: "./results/"
fasta: "./transcripts.fa"
gtf: "./annotation.gtf.gz"
```

> [!WARNING]
> Use `-params-file` for pipeline parameters. Do not use `-c` for ordinary pipeline params.

## Core Nextflow arguments

### `-profile`

Select execution and software backend presets. Common options include `docker`, `singularity`, `podman`, `apptainer`, and `conda`.

Multiple profiles can be combined, for example:

```bash
-profile test,docker
```

Order matters: later profiles override earlier ones.

### `-resume`

Resume from cached work where inputs and process configuration match prior runs.

### `-c`

Use only for Nextflow config overrides (resources/execution behavior), not routine pipeline params.

## Reproducibility and updates

Update cached pipeline code:

```bash
nextflow pull nf-core/rnastructurome
```

Pin a release for reproducibility:

```bash
nextflow run nf-core/rnastructurome -r 1.3.1 ...
```

## Custom configuration

For resources, tool args, custom containers, and institutional configs, see:

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
