# `subset_fastq_from_bam.py`

Utility to build a small FASTQ subset for `nf-core/rnastructurome` from transcript-aligned BAM files.

It is intended for the case where:

- the original FASTQ files are large
- you already have BAMs on Codon
- you want a small reproducible test subset without relying on host `samtools`

The script runs `samtools view` from one of three backends:

- host `samtools` on `PATH`
- `singularity exec` / `apptainer exec`
- Docker

It then finds three transcripts spanning low, medium, and high alignment support, extracts read names for those transcripts, and writes subset FASTQ files plus a replacement samplesheet.

## What a sample is

Each sample is a row in your existing pipeline samplesheet.

The script matches BAMs to samples using:

1. `sample_id` if that column exists
2. otherwise `sample`

For example, with:

```csv
sample,sample_id,fastq_1,fastq_2,cell_line,condition,replicate,organism
treated_rep1,treated_rep1,/data/treated.fastq.gz,,HEK293T,treated,1,Homo sapiens
untreated_rep1,untreated_rep1,/data/untreated.fastq.gz,,HEK293T,untreated,1,Homo sapiens
```

- `treated_rep1` is one sample
- `untreated_rep1` is a second sample
- each sample must have its own BAM mapping via `--bam sample_key=/path/to/sample.bam`

The script does not infer sample identity from BAM filenames alone. The sample key in `--bam` must match the samplesheet row key.

## Requirements

- one of:
  - host `samtools`
  - Singularity / Apptainer
  - Docker
- transcript-aligned BAMs, one per samplesheet row
- original FASTQ paths still present in the samplesheet

Default Singularity image:

```text
 /hps/nobackup/agb/rnacentral/chemprob/nf-core-rnastructurome/work/singularity/img/depot.galaxyproject.org-singularity-samtools-1.22.1--h96c455f_0.img
```

Default Docker image:

```text
docker.io/rnastructurome/rnaframework:2.9.6-r1
```

## Basic usage

```bash
python bin/subset_fastq_from_bam.py \
  --samplesheet /path/to/samplesheet.csv \
  --bam treated_rep1=/path/to/treated_rep1.bam \
  --bam untreated_rep1=/path/to/untreated_rep1.bam \
  --bam denatured_rep1=/path/to/denatured_rep1.bam \
  --container-engine singularity \
  --singularity-image /path/to/samtools.img \
  --output-dir /path/to/subset_out
```

## How transcript selection works

Default behavior:

1. Count primary BAM alignments per transcript across all supplied BAMs
2. Drop transcripts below `--min-alignments` (default: `10`)
3. Choose three transcripts at low, medium, and high quantiles (default: `0.2,0.5,0.8`)
4. Collect read names from each BAM for those transcripts
5. Filter each sample FASTQ to only those reads

You can bypass auto-selection and force exact transcripts:

```bash
python bin/subset_fastq_from_bam.py \
  --samplesheet /path/to/samplesheet.csv \
  --bam treated_rep1=/path/to/treated_rep1.bam \
  --bam untreated_rep1=/path/to/untreated_rep1.bam \
  --bam denatured_rep1=/path/to/denatured_rep1.bam \
  --transcript ENST00000331789 \
  --transcript ENST00000456328 \
  --transcript ENST00000619216 \
  --output-dir /path/to/subset_out
```

## Main arguments

- `--samplesheet`: existing pipeline samplesheet CSV
- `--bam SAMPLE=/path/to/sample.bam`: BAM for one samplesheet row; repeat once per sample
- `--output-dir`: destination for subset FASTQs and reports
- `--min-alignments`: minimum per-transcript alignment count for auto-selection
- `--quantiles`: low, medium, high selection quantiles
- `--transcript`: manually specify exactly three transcripts instead of auto-selection
- `--sample-id-column`: preferred samplesheet key column; default is `sample_id`, with fallback to `sample`
- `--container-engine`: `auto`, `host`, `singularity`, or `docker`
- `--container-image`: Docker image used to provide `samtools`
- `--singularity-image`: Singularity / Apptainer image used to provide `samtools`
- `--container-platform`: Docker platform, default `linux/amd64`

## Output files

The output directory contains:

- subset FASTQ files for each sample
- `selected_transcripts.tsv`: the three chosen transcripts and aggregate alignment counts
- `subset_summary.tsv`: per-sample selected read counts and transcript breakdown
- `<samplesheet>.subset.csv`: a new samplesheet pointing to the subset FASTQs

Example output names:

```text
subset_out/
  treated_rep1.R1.treated.subset.fastq.gz
  untreated_rep1.R1.untreated.subset.fastq.gz
  selected_transcripts.tsv
  subset_summary.tsv
  samplesheet.subset.csv
```

## Notes and limitations

- BAM reference names are assumed to be transcript identifiers, not chromosome names
- duplicate, supplementary, secondary, and unmapped records are ignored during counting
- paired-end FASTQ is supported if `fastq_2` is populated in the samplesheet
- the script filters the original FASTQ files named in the samplesheet; it does not reconstruct reads from BAM
