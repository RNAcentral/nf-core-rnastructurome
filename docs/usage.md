# nf-core/rnastructurome: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnastructurome/usage](https://nf-co.re/rnastructurome/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files supplied with `-c` can be used for infrastructure settings such as resources, executors, containers, or module arguments, but should not be used for ordinary pipeline parameters.

## Samplesheet input

The easiest way to run this pipeline is to create a full samplesheet that contains all of the information about each sample in a comma-separated file, including the header row shown below and pass it as `--input '[path to samplesheet file]'`:

```csv title="full_samplesheet.csv"
sample,sample_id,fastq_1,fastq_2,method,principle,sample_group,condition,replicate,organism,pH,adapter_3p,adapter_5p,umi_pattern
HEK293T_treated_r1,GSM000001,/data/treated_r1.fastq.gz,,SHAPE,RT-stop,HEK293T,treated,1,Homo sapiens,7.5,,,
HEK293T_untreated_r1,GSM000002,/data/untreated_r1.fastq.gz,,SHAPE,RT-stop,HEK293T,untreated,1,Homo sapiens,7.5,,,
```

An [example samplesheet](../assets/samplesheet.csv) is provided.

However, you can also provide a more minimal version if for example you don't need to specify some of the options, like in the example above you could decide to not provide the columns with information about the adapters or umi pattern if you are happy to use the default options. 
You can also provide a very minimal samplesheet with just the information required about each invidivual sample and pass the uniform values across all samples as parameters. For example:

```csv title="minimal_samplesheet.csv"
sample,fastq_1,sample_group,condition,replicate
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
| `sample_group` | yes      | Group key used for control pairing in `rf-norm`.                                                      |
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

`rf-count` (transcriptome route) or `rf-count-genome` (genome route) quantifies chemical probing signal from the aligned BAM files. What exactly is counted depends on the principle: for RT-stop experiments it tallies read 3′ ends that accumulate at modified bases; for MaP it measures per-position mutation rates.

Per-transcript coverage plots are disabled by default (`--rfcount_img false`) because they are slow to generate at transcriptome scale. Enable with `--rfcount_img true` if you want them.

On the genome route, strandedness is handled automatically — RT-stop libraries are always treated as second-strand (this is fixed by experimental design), while for MaP the pipeline infers strandedness per sample using `RSeQC infer_experiment`. If inference is ambiguous you can override it with `--rfcount_strandedness first|second|unstranded`. After counting, `rf-rctools extract` automatically converts the genome-coordinate RC files to transcript-level RC files using the GTF before passing to `rf-norm`.

A few parameters worth knowing about:

`--rfcount_trim_5prime` trims a fixed number of bases from the 5′ end of each read before counting. This is useful for RT-stop experiments where the first few bases after the adapter can carry sequence-context bias that inflates apparent stop rates.

`--rfcount_mask_file` accepts a BED file of regions to exclude from counting entirely — handy for masking rRNA or other highly-expressed contaminating transcripts that would otherwise dominate the output.

`--rfcount_primary_only` restricts counting to primary alignments. If you have increased `--star_multimap_nmax` or `--bowtie_k` to allow multimappers, you may want to pair this with `--rfcount_primary_only` to avoid double-counting reads that aligned to multiple loci.

For the full list of available options see the [rf-count documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-count/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

### rf-norm

`rf-norm` normalises per-position counts into reactivity scores. Samples are grouped by `sample_group + replicate`, and treated samples are normalised against their matched controls within the same group. The scoring and normalisation method are selected automatically based on the probing principle and which conditions are present: RT-stop with an untreated control uses Ding scoring with box-plot normalisation; without an untreated control it falls back to Rouskin scoring with Winsorizing. MaP with an untreated control uses Siegfried scoring; without, Zubradt. 

You can override the normalisation method alone (without changing scoring) with `--rfnorm_norm_method`:
- `1` is 2-8% normalisation (takes the top 10% of reactivities, discards the very highest 2%, and uses the mean of the remaining 8% as the scaling factor)
- `2` is 90% Winsorizing (clips any value above the 95th or below the 5th percentile to those boundaries, then divides all values by the 95th percentile)
- `3` is box-plot normalisation (removes outliers beyond 1.5× IQR then divides by the mean of the next top 10%), which is the default for most conditions
- `4` is Mitchell normalisation (MaP only, uses the higher of the mean 90th–95th percentile reactivity or the 75th percentile of non-zero reactivities as the scaling factor)

Treated samples are usually paired with untreated by matching `sample_group + replicate` exactly. However, in cases where the authors did not create an exact matching untreated sample for some specific treatments but have one for other samples in the same dataset, the pipeline falls back to an untreated sample sharing the same `sample_group` base token and replicate — for example, `MDA-MB-231_DMSO_treated_r1` will match exactly to `MDA-MB-231_DMSO_untreated_r1`, but `MDA-MB-231_MTX_treated_r1` will also pair with `MDA-MB-231_DMSO_untreated_r1` if no exact matching untreated sample exists. A warning is emitted when a fallback is used; if more than one candidate matches, the pipeline errors. Disable this behaviour with `--fuzzy_untreated_pairing false`, in which case unmatched groups proceed without a negative control.

For DMS experiments, a few defaults change automatically. `--rfnorm_reactive_bases` is set to `AC` (or `ACGU` when `pH ≥ 8`). `--rfnorm_dynamic_window` defaults to `50` when `pH < 8`; it controls the size of the sliding window used to compute local normalisation factors along the transcript — a smaller window is better suited to DMS at physiological pH where reactivity can vary sharply over short stretches. `--rfnorm_nan` defaults to `100` rather than `1000`; it sets the minimum number of reads required at a position for a reactivity value to be reported — positions with fewer reads are set to NaN instead of reporting a potentially unreliable values. These can all be overridden explicitly if needed.

Per-transcript reactivity plots are disabled by default (`--rfnorm_img false`) because generating them via R for thousands of transcripts is the main source of rf-norm runtime on large datasets. Enable with `--rfnorm_img true` if you want them.

For the full list of available options see the [rf-norm documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

### rf-normfactor (cross-experiment normalisation)

By default `rf-norm` normalises each `sample_group + replicate` group independently. This is fine for looking at one sample, but it means reactivities are not on a common scale across samples — so comparing reactivity between replicates or conditions can be misleading. `rf-normfactor` solves this by deriving a single set of transcriptome-wide normalisation factors across all of a reference's samples at once, which `rf-norm` then applies (via `-nf`) to every group. This is the approach used for cross-sample analyses in recent transcriptome-wide SHAPE-MaP studies.

Whether this runs is decided **per reference** by `--rfnorm_use_normfactor`:

- unset (default) — **auto**: enabled for a reference that has paired treated/untreated controls **or** more than one treated sample. A reference with a single treated sample and no control keeps independent per-sample box-plot normalisation.
- `true` — force on for every reference.
- `false` — force off; always use per-sample normalisation.

The factor is computed once per reference across all that reference's treated samples (each paired with its own resolved untreated/denatured control), so within a run different references are normalised independently but all samples mapping to the same reference share one scale. Cross-experiment normalisation with Siegfried scoring requires every treated sample on a reference to have a matched untreated control — if some do but others do not, the pipeline errors rather than mispair; add the missing controls or set `--rfnorm_use_normfactor false`.

The factor calculation reuses the same scoring and normalisation settings as the downstream `rf-norm` (scoring method, `--rfnorm_reactive_bases`, `--rfnorm_pseudocount`, `--rfnorm_max_score`, `--rfnorm_ignore_lower_than_untreated`, `--rfnorm_max_untreated_mut`, `--rfnorm_max_mutation_rate`, `--rfnorm_median_coverage`), so factors are computed on the same footing as the reactivities. The minimum per-base coverage used when calculating factors is set with `--rfnorm_normfactor_min_coverage` (default `1000`, matching the `-mc` value used in published transcriptome-wide protocols).

> **Note**: because auto is the default, multi-sample runs (or any run with paired untreated controls) now use cross-experiment normalisation rather than independent per-sample normalisation, which changes the scale of reported reactivities relative to earlier behaviour. Set `--rfnorm_use_normfactor false` to restore per-sample normalisation. Auto is also disabled while `--rfnorm_chunk_size` is in use.

For the full list of available options see the [rf-normfactor documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-normfactor/).

### rf-correlate (replicate reproducibility QC)

When a sample group has more than one replicate, `rf-correlate` measures how reproducible the replicates are by computing pairwise correlations between their reactivity profiles (transcriptome-wide and per-transcript). This runs by default — set `--correlate_replicates false` to disable it — and is a no-op for single-replicate groups. The overall pairwise correlations are summarised into a **MultiQC table** (one row per sample group, reporting the number of replicates and the mean and minimum pairwise correlation), so a low number flags a discordant replicate before it dilutes the folded consensus.

By default it uses Pearson correlation; switch to Spearman with `--correlate_spearman true`. Other options: `--correlate_min_values` sets the minimum number of covered positions required to correlate a transcript (a value between 0 and 1 is read as a fraction of transcript length), `--correlate_ignore_sequence` tolerates sequence differences (e.g. SNVs) between compared transcripts, and `--correlate_img true` additionally writes the rf-correlate correlation heatmap PDF. The full per-transcript pairwise TSVs and the correlation matrix are published under `correlate/`.

For the full list of available options see the [rf-correlate documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-correlate/).

### rf-fold

`rf-fold` predicts RNA secondary structures from normalised reactivity profiles using ViennaRNA. The pipeline groups XMLs by `sample_group`, merging replicates, and folds them together. Dot-bracket output is the default; CT format can be enabled with `--rffold_ct`.

#### rf-structextract (optional)

Once structures are folded, `rf-structextract` can pull out the high-confidence structural elements — substructures whose bases show consistently low reactivity **and** low Shannon entropy (the signature of a well-defined, stably folded region) and that meet thermodynamic and geometric criteria. It runs after rf-fold on each `sample_group`, using the fold output (structures + Shannon entropy) together with the group's rf-norm reactivity profiles. Enable it with `--structextract true`; it is off by default.

The selection criteria are exposed as parameters, matching the rf-structextract defaults: window size for the median reactivity/Shannon scan (`--structextract_win_size`, 50 nt), the minimum transcript length evaluated (`--structextract_min_transcript_len`, 500 nt), the minimum fraction of bases that must sit below the transcript median (`--structextract_min_below_median`, 0.7), the minimum paired-base fraction (`--structextract_min_paired_frac`, 0.45), and motif length bounds (`--structextract_min_motif_len` 50, `--structextract_max_motif_len` unset). You can restrict output to multiway-junction elements with `--structextract_multiway_only`, or skip the reactivity or Shannon test individually with `--structextract_ignore_react` / `--structextract_ignore_shannon`. To additionally keep only motifs whose folding free energy is significantly lower than expected by chance, set `--structextract_eval_energy true` (tuned with `--structextract_pvalue`, `--structextract_n_shufflings`, and `--structextract_dinucl_shuffle`).

Results are written under the sample group's fold directory at `fold/<group>/extracted_structures/`, with the dot-bracket motifs in `dotbracket/` and one 2D diagram per motif in `images/`. The diagrams are drawn with ViennaRNA RNAplot and coloured by SHAPE reactivity in the same style as the rf-fold structure plots (the motif's reactivity is sliced from its parent transcript); disable them with `--structextract_plot false`.

For the full list of available options see the [rf-structextract documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-structextract/).

After folding, the pipeline draws reactivity-coloured 2D structure diagrams using [R2DT](https://github.com/RNAcentral/R2DT) where a template exists (rRNA, snRNA, tRNA, and other well-characterised RNA families), falling back to ViennaRNA's `RNAplot` for everything else. Note that R2DT is only available when running with a container profile (`docker`, `singularity`, `apptainer`) — under `conda` or `mamba` it is not currently available, and all diagrams will fall back to ViennaRNA's `RNAplot` instead.

Some important parameters to consider are `--rffold_slope` (default `4.6`) and `--rffold_intercept` (default `-2.2`), which control how reactivity values are converted into folding constraints. These defaults were determined using this pipeline on E. coli total RNA SHAPE-MaP data from [Borovska et al. 2026](https://www.nature.com/articles/s41587-025-02739-0) and should serve as a good starting point for most experiments. You can use `rf-jackkife` with your own data to determine the optimal values for your dataset.

`--rffold_shannon_entropy` (default `true`) computes per-position Shannon entropy alongside the predicted structure, which gives a measure of folding confidence. `--rffold_only_common` keeps only transcripts covered in at least N XML experiments — when not set explicitly, the pipeline enables this automatically for fold groups with more than one replicate, setting N to the number of replicates in that group so only transcripts present in all replicates are retained. `--rffold_unconstrained` folds without using reactivity data at all, useful as a baseline comparison.

ViennaRNA `RNAplot` structure diagrams from rf-fold are disabled by default (`--rffold_img false`). Enable with `--rffold_img true` to generate them. Note that R2DT template-based diagrams (see above) are independent of this flag and are always attempted when running with a container profile.

For the full list of available options see the [rf-fold documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-fold/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

### rf-jackknife
This module can be run for an individual sample or by pooling the samples together to get a consensus slope/intercept.

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
  --jackknife_reference ecoli_rrna_calibration/ecoli_k12_rrna_reference_collab.db \
  -profile docker
```

For repeated runs with the same parameters, use a params file:

```bash
nextflow run nf-core/rnastructurome -profile docker -params-file params.yaml
```

Example `params.yaml`:

```yaml
input: /path/to/samplesheet.csv
outdir: results
organism: Homo sapiens
method: SHAPE
principle: RT-stop

# trim 2 nt from 5' end to remove sequence-context bias
rfcount_trim_5prime: 2

# use 2-8% normalisation instead of the default box-plot
rfnorm_norm_method: 1
```

## Reproducibility

Always specify the pipeline version when running on your data using `-r` (one hyphen):

```bash
nextflow run nf-core/rnastructurome -r 1.0.0 -profile docker --input samplesheet.csv --outdir results
```

This pins the exact version of the pipeline code and all software containers, so re-running with the same `-r` tag months later will produce identical results. The version is recorded in the MultiQC report and in the RDAT files so you always know what was used.

For complete reproducibility, save your parameters to a `params.yaml` file (see [Running the pipeline](#running-the-pipeline)) and share it alongside your data. When doing so, avoid including cluster-specific paths or institutional profile names — use relative paths or published dataset identifiers instead so that others can reproduce the run in their own environment.

## Core Nextflow arguments

These options are part of Nextflow itself and use a single hyphen. Pipeline parameters use a double hyphen.

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments or analyses.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods — see below. We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility; when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [nf-core/configs](https://github.com/nf-core/configs) at run time, making institutional cluster profiles available automatically. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Multiple profiles can be combined; they are loaded in sequence so later profiles can overwrite earlier ones:

```bash
-profile test,docker
```

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on `PATH`. This is not recommended as it can lead to different results on different machines.

| Profile | Description |
|---|---|
| `test` | Minimal test using the STAR genome-alignment route. Uses human mitochondrial chromosome (MT-RNR1) test data — no other parameters needed. |
| `test_transcriptome` | Minimal test using the Bowtie2 transcriptome route. Uses a single-transcript FASTA (ENST00000389680 / MT-RNR1) to exercise the `--transcriptome` path — no other parameters needed. |
| `test_prokaryote` | Prokaryote transcriptome-route test using E. coli 16S rRNA DMS-MaP data and a 16S reference structure for jackknife calibration — no other parameters needed. Significantly faster with `conda` than with container profiles since it includes rf-jackknife. |
| `docker` | Use [Docker](https://docs.docker.com/engine/installation/) containers. |
| `singularity` | Use [Singularity](https://www.sylabs.io/guides/3.0/user-guide/) containers. |
| `podman` | Use [Podman](https://podman.io/) containers. |
| `shifter` | Use [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/) containers. |
| `charliecloud` | Use [Charliecloud](https://hpc.github.io/charliecloud/) containers. |
| `apptainer` | Use [Apptainer](https://apptainer.org/) containers. |
| `wave` | Enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ≥ 24.03.0-edge). |
| `conda` | Use [Conda](https://conda.io/miniconda.html). Please only use Conda as a last resort when containers are not possible. Note that R2DT structure diagrams are not available under conda/mamba — ViennaRNA RNAplot is used as fallback. |
| `arm64` | Applies overrides supplying ARM-compatible containers and Conda environments. See [Running on Linux ARM architectures](#running-on-linux-arm-architectures). |

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file for process resources, executors, or module-specific `ext.args`. See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18), it will automatically be resubmitted with a higher resource request (2× original, then 3× original). If it still fails after the third attempt then the pipeline execution is stopped.

Computationally intensive steps in this pipeline include STAR alignment, rf-count, rf-norm (especially in chunked genome mode), and rf-fold. These are labelled `process_high` or `process_medium` and will benefit most from tuning.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) sections of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the [nf-core/configs](https://github.com/nf-core/configs) git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the nf-core/configs repository with the addition of your config file, associated documentation file (see examples in [nf-core/configs/docs](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the [main Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time. Some HPC setups also allow you to run Nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

### Nextflow memory requirements

In some cases, the Nextflow Java virtual machine can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Troubleshooting

### General debugging steps

1. Check `.nextflow.log` for the first error message.
2. Check `pipeline_info/execution_trace_<timestamp>.txt` to identify which process failed and its exit status.
3. Use `-resume` after fixing issues — all previously successful tasks will be cached and skipped.
