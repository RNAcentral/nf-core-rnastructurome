# nf-core/rnastructurome: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnastructurome/usage](https://nf-co.re/rnastructurome/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

`nf-core/rnastructurome` processes structure-probing RNA-seq data from FASTQ to normalized reactivities and secondary structure outputs.

High-level workflow:

1. Input QC and trimming (`FastQC`, `cutadapt`, optional `umi_tools extract`)
2. Alignment (`bowtie` / `bowtie2`) and BAM post-processing (`samtools`, optional `umi_tools dedup`)
3. Reactivity counting (`rf-count`) for RT-stop and MaP
4. Reactivity normalization (`rf-norm`) using available controls
5. Structure inference (`rf-fold`) from normalized XMLs
6. Aggregated reporting (`MultiQC`)

## Quick start

```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheet.csv \
  --fasta transcripts.fa \
  --outdir results
```

If resuming after an interruption or fix:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --fasta transcripts.fa --outdir results -resume
```

## Required inputs

1. A samplesheet (`--input`)
2. A transcript reference FASTA (`--fasta`) or a genome configuration that resolves to transcript FASTA

Reference resolution behavior:

- If `--fasta` is set, that file is used directly and should be a transcript FASTA.
- If `--fasta` is not set, the pipeline first tries transcript FASTA keys in `params.genomes[reference_key]`, where `reference_key` is resolved from `organism`:
  - `transcript_fasta`
  - `transcriptome`
  - `cdna`
- If no transcript FASTA path is configured, it falls back to Ensembl auto-download by species:
  - Uses `params.genomes[reference_key].ensembl_species` when defined, then `--ensembl_species_map` aliases, otherwise treats the reference key itself as Ensembl species if it matches `genus_species` (e.g. `homo_sapiens`, `saccharomyces_cerevisiae`).
  - Latin binomials such as `Homo sapiens` or `Saccharomyces cerevisiae` are normalized automatically to Ensembl species format (`homo_sapiens`, `saccharomyces_cerevisiae`).
  - Downloads `cdna.all.fa.gz` and, when available, `ncrna.fa.gz` from Ensembl FTP, then merges available files into one transcript FASTA used for mapping and RNAframework.
- Ensembl source can be tuned with:
  - `--ensembl_release` (`current` by default; accepts values like `114` or `release-114`)
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

| Column | Description |
| --- | --- |
| `sample` | Sample identifier. Re-used IDs are concatenated before downstream analysis (multi-lane support). |
| `fastq_1` | Read 1 FASTQ path (`.fastq.gz` / `.fq.gz`). |
| `fastq_2` | Read 2 FASTQ path for paired-end data (optional). |
| `cell_line` | Group key used for control pairing in `rf-norm`. |
| `condition` | One of `treated`, `untreated`, `denatured`. |
| `replicate` | Replicate key used for control pairing in `rf-norm`. |

Optional per-sample columns supported by the pipeline include:

- `organism` (preferred for transcriptome auto-resolution; e.g. `Homo sapiens`)
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

## Preprocessing behavior

### Cutadapt by principle

| Principle | `--cutadapt-5quality` | `--cutadapt-3quality` | Notes |
| --- | --- | --- | --- |
| `RT-stop` | forced to `0` | default `20` | 5' quality trimming is intentionally disabled |
| `MaP` | default `20` | default `20` | both are user-configurable |

Adapter precedence:

1. Per-sample (`adapter_5p`, `adapter_3p`)
2. Global (`--cutadapt_adapter_5p`, `--cutadapt_adapter_3p`)

If no adapters are provided, the pipeline falls back to `AGATCGGAAGAGC` for both the 5' and 3' adapter. Set `--cutadapt_quality_only` to skip adapter trimming entirely and perform quality/length filtering only.

### Optional UMI extraction

Enable with `umi_pattern` (sample-level or global `--umi_pattern`).

- Patterns containing only `N/C/X` use direct `umi_tools --bc-pattern` mode.
- IUPAC patterns (for example with `D`) are converted automatically to regex mode.

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

- For `method=DMS`, pipeline defaults `--rfnorm_reactive_bases` to `AC` if unset.
- `R` must exist at `--rnaframework_r_path` inside the RNAframework container.

### `rf-fold`

Default behavior:

- Dot-bracket output is default (CT is optional).
- Graphical fold reports are always enabled (`-g -R`).
- RNAplot overlays are enabled via `-vrp` (default `/usr/bin/RNAplot`).
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

- The pipeline groups normalized XMLs by `cell_line + principle + method + rf-norm scoring/normalization mode` and folds them together.
- Single-XML groups run as standard one-sample `rf-fold`.

### RNAframework runtime settings

- `--rnaframework_container`: RNAframework image for local modules
- `--rnaframework_r_path`: `R` binary inside that container

## Outputs at a glance

Main output areas under `--outdir`:

- `fastqc/`, `cutadapt/`, `bowtie*/`, `samtools*/` for preprocessing and alignment
- `count/`: count tables and always-on plots
- `norm/`: normalized XML and always-on normalization plots
- `fold/`: inferred secondary structures (dot-bracket by default), fold reports, optional CT, optional dotplots
- `multiqc/`: final aggregated QC report

For full output details, see [output documentation](output.md).

## Running the pipeline

Typical usage:

```bash
nextflow run nf-core/rnastructurome --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
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
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
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
