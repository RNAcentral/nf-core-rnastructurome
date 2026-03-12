#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable


DEFAULT_CONTAINER_IMAGE = "docker.io/rnastructurome/rnaframework:2.9.6-r1"
DEFAULT_SINGULARITY_IMAGE = "/hps/nobackup/agb/rnacentral/chemprob/nf-core-rnastructurome/work/singularity/img/depot.galaxyproject.org-singularity-samtools-1.22.1--h96c455f_0.img"
DEFAULT_QUANTILES = (0.2, 0.5, 0.8)
SAM_FLAG_UNMAPPED = 0x4
SAM_FLAG_SECONDARY = 0x100
SAM_FLAG_DUPLICATE = 0x400
SAM_FLAG_SUPPLEMENTARY = 0x800


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select low/medium/high coverage transcripts from transcript-aligned BAMs "
            "and subset matching FASTQ reads for each sample in an existing samplesheet."
        )
    )
    parser.add_argument("--samplesheet", required=True, help="Existing pipeline samplesheet CSV.")
    parser.add_argument(
        "--bam",
        action="append",
        required=True,
        metavar="SAMPLE=/path/to/sample.bam",
        help="Map a samplesheet sample or sample_id value to a transcript-aligned BAM. Repeat per sample.",
    )
    parser.add_argument("--output-dir", required=True, help="Directory for subset FASTQs and reports.")
    parser.add_argument(
        "--min-alignments",
        type=int,
        default=10,
        help="Minimum primary alignment records per transcript to consider during auto-selection.",
    )
    parser.add_argument(
        "--transcript",
        action="append",
        default=[],
        help="Override transcript selection. Repeat exactly three times to skip auto-selection.",
    )
    parser.add_argument(
        "--quantiles",
        default="0.2,0.5,0.8",
        help="Ascending quantiles used to choose low, medium, and high coverage transcripts.",
    )
    parser.add_argument(
        "--container-engine",
        default="auto",
        choices=("auto", "host", "singularity", "docker"),
        help="How to run samtools: use host PATH, singularity, docker, or auto-detect. Default: auto",
    )
    parser.add_argument(
        "--container-image",
        default=DEFAULT_CONTAINER_IMAGE,
        help=f"Docker image used to run samtools. Default: {DEFAULT_CONTAINER_IMAGE}",
    )
    parser.add_argument(
        "--singularity-image",
        default=DEFAULT_SINGULARITY_IMAGE,
        help=f"Singularity/Apptainer image used to run samtools. Default: {DEFAULT_SINGULARITY_IMAGE}",
    )
    parser.add_argument(
        "--container-platform",
        default="linux/amd64",
        help="Container platform passed to docker run. Default: linux/amd64",
    )
    parser.add_argument(
        "--sample-id-column",
        default="sample_id",
        help="Preferred samplesheet column for matching --bam mappings. Falls back to sample.",
    )
    return parser.parse_args()


class SamtoolsRunner:
    def __init__(self, engine: str, image: str, singularity_image: str, platform: str) -> None:
        self.engine = engine
        self.image = image
        self.singularity_image = singularity_image
        self.platform = platform
        self.uid = os.getuid()
        self.gid = os.getgid()

    def _resolve_command(self, bam_path: Path) -> list[str]:
        bam_path = bam_path.resolve()
        if self.engine == "host":
            return ["samtools", "view", str(bam_path)]

        if self.engine == "singularity":
            image_path = Path(self.singularity_image).expanduser().resolve()
            if not image_path.exists():
                raise FileNotFoundError(f"Singularity image not found: {image_path}")
            singularity_bin = shutil.which("singularity") or shutil.which("apptainer")
            if singularity_bin is None:
                raise RuntimeError("Neither singularity nor apptainer is available on PATH.")
            return [singularity_bin, "exec", str(image_path), "samtools", "view", str(bam_path)]

        parent = bam_path.parent
        container_bam = f"/input/{bam_path.name}"
        return [
            "docker",
            "run",
            "--rm",
            "--platform",
            self.platform,
            "-u",
            f"{self.uid}:{self.gid}",
            "-v",
            f"{parent}:/input:ro",
            self.image,
            "samtools",
            "view",
            container_bam,
        ]

    def view_lines(self, bam_path: Path) -> Iterable[str]:
        command = self._resolve_command(bam_path)
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert process.stdout is not None
        try:
            for line in process.stdout:
                if line.strip():
                    yield line
        finally:
            stderr = process.stderr.read() if process.stderr is not None else ""
            return_code = process.wait()
            if return_code != 0:
                raise RuntimeError(f"samtools view failed for {bam_path}: {stderr.strip()}")


def resolve_container_engine(requested: str) -> str:
    if requested != "auto":
        return requested
    if shutil.which("samtools"):
        return "host"
    if shutil.which("singularity") or shutil.which("apptainer"):
        return "singularity"
    if shutil.which("docker"):
        return "docker"
    raise RuntimeError("No usable samtools backend found. Install samtools, singularity/apptainer, or docker.")


def parse_bam_mappings(values: list[str]) -> dict[str, Path]:
    mappings: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Invalid --bam value '{value}'. Expected SAMPLE=/path/to/sample.bam")
        sample_key, bam_raw = value.split("=", 1)
        sample_key = sample_key.strip()
        bam_path = Path(bam_raw.strip()).expanduser().resolve()
        if not sample_key:
            raise ValueError(f"Invalid --bam value '{value}'. Sample key is empty.")
        if sample_key in mappings:
            raise ValueError(f"Duplicate --bam mapping for sample '{sample_key}'.")
        if not bam_path.exists():
            raise FileNotFoundError(f"BAM path does not exist for sample '{sample_key}': {bam_path}")
        mappings[sample_key] = bam_path
    return mappings


def parse_quantiles(raw: str) -> tuple[float, float, float]:
    parts = [part.strip() for part in raw.split(",") if part.strip()]
    if len(parts) != 3:
        raise ValueError("--quantiles must contain exactly three comma-separated values.")
    quantiles = tuple(float(part) for part in parts)
    if any(value < 0.0 or value > 1.0 for value in quantiles):
        raise ValueError("--quantiles values must be between 0 and 1.")
    if list(quantiles) != sorted(quantiles):
        raise ValueError("--quantiles must be in ascending order.")
    return quantiles  # type: ignore[return-value]


def load_samplesheet(path: Path, sample_id_column: str) -> tuple[list[dict[str, str]], str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Samplesheet is missing a header row: {path}")
        rows = list(reader)

    if not rows:
        raise ValueError(f"Samplesheet contains no samples: {path}")

    if sample_id_column in reader.fieldnames:
        return rows, sample_id_column
    if "sample" in reader.fieldnames:
        return rows, "sample"
    raise ValueError(
        f"Samplesheet must contain either '{sample_id_column}' or 'sample' for BAM mapping: {path}"
    )


def fastq_reader(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8")


def fastq_writer(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "wt", encoding="utf-8")


def normalize_fastq_name(header: str) -> str:
    token = header.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def should_skip_alignment(flag: int, rname: str) -> bool:
    if rname == "*":
        return True
    return bool(flag & (SAM_FLAG_UNMAPPED | SAM_FLAG_SECONDARY | SAM_FLAG_DUPLICATE | SAM_FLAG_SUPPLEMENTARY))


def scan_alignment_counts(samtools: SamtoolsRunner, bam_by_sample: dict[str, Path]) -> dict[str, dict[str, int]]:
    counts_by_sample: dict[str, dict[str, int]] = {}
    for sample_key, bam_path in bam_by_sample.items():
        sample_counts: dict[str, int] = defaultdict(int)
        for line in samtools.view_lines(bam_path):
            fields = line.split("\t", 4)
            if len(fields) < 4:
                continue
            flag = int(fields[1])
            rname = fields[2]
            if should_skip_alignment(flag, rname):
                continue
            sample_counts[rname] += 1
        counts_by_sample[sample_key] = dict(sample_counts)
    return counts_by_sample


def aggregate_counts(counts_by_sample: dict[str, dict[str, int]]) -> dict[str, int]:
    aggregated: dict[str, int] = defaultdict(int)
    for sample_counts in counts_by_sample.values():
        for transcript, count in sample_counts.items():
            aggregated[transcript] += count
    return dict(aggregated)


def choose_transcripts(
    aggregated_counts: dict[str, int],
    min_alignments: int,
    quantiles: tuple[float, float, float],
) -> list[tuple[str, int, str]]:
    eligible = sorted(
        ((transcript, count) for transcript, count in aggregated_counts.items() if count >= min_alignments),
        key=lambda item: (item[1], item[0]),
    )
    if len(eligible) < 3:
        raise ValueError(
            f"Need at least 3 transcripts with >= {min_alignments} alignments, found {len(eligible)}."
        )

    labels = ("low", "medium", "high")
    selected: list[tuple[str, int, str]] = []
    used: set[str] = set()
    last_index = len(eligible) - 1
    for label, quantile in zip(labels, quantiles):
        target_index = min(max(round(last_index * quantile), 0), last_index)
        candidate_indices = sorted(range(len(eligible)), key=lambda idx: (abs(idx - target_index), idx))
        chosen_index = next(idx for idx in candidate_indices if eligible[idx][0] not in used)
        transcript, count = eligible[chosen_index]
        selected.append((transcript, count, label))
        used.add(transcript)
    return selected


def validate_requested_transcripts(requested: list[str], aggregated_counts: dict[str, int]) -> list[tuple[str, int, str]]:
    unique = []
    seen: set[str] = set()
    for transcript in requested:
        if transcript in seen:
            raise ValueError(f"Duplicate --transcript value '{transcript}'.")
        seen.add(transcript)
        unique.append(transcript)
    if len(unique) != 3:
        raise ValueError("Provide exactly three --transcript values when overriding auto-selection.")
    labels = ("low", "medium", "high")
    selected: list[tuple[str, int, str]] = []
    for label, transcript in zip(labels, unique):
        if transcript not in aggregated_counts:
            raise ValueError(f"Requested transcript '{transcript}' was not observed in the BAM inputs.")
        selected.append((transcript, aggregated_counts[transcript], label))
    return selected


def collect_read_names(
    samtools: SamtoolsRunner,
    bam_by_sample: dict[str, Path],
    transcripts: set[str],
) -> tuple[dict[str, set[str]], dict[str, dict[str, int]]]:
    read_names_by_sample: dict[str, set[str]] = {}
    transcript_counts_by_sample: dict[str, dict[str, int]] = {}
    for sample_key, bam_path in bam_by_sample.items():
        sample_names: set[str] = set()
        sample_counts: dict[str, int] = defaultdict(int)
        for line in samtools.view_lines(bam_path):
            fields = line.split("\t", 4)
            if len(fields) < 4:
                continue
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            if should_skip_alignment(flag, rname) or rname not in transcripts:
                continue
            sample_names.add(qname)
            sample_counts[rname] += 1
        read_names_by_sample[sample_key] = sample_names
        transcript_counts_by_sample[sample_key] = dict(sample_counts)
    return read_names_by_sample, transcript_counts_by_sample


def subset_fastq(source: Path, destination: Path, selected_names: set[str]) -> int:
    written = 0
    with fastq_reader(source) as reader, fastq_writer(destination) as writer:
        while True:
            header = reader.readline()
            if not header:
                break
            sequence = reader.readline()
            plus = reader.readline()
            quality = reader.readline()
            if not quality:
                raise ValueError(f"Malformed FASTQ record in {source}")
            if normalize_fastq_name(header) in selected_names:
                writer.write(header)
                writer.write(sequence)
                writer.write(plus)
                writer.write(quality)
                written += 1
    return written


def build_output_fastq_path(output_dir: Path, sample_key: str, original_path: Path, suffix: str) -> Path:
    base_name = original_path.name
    if base_name.endswith(".fastq.gz"):
        stem = base_name[:-9]
        extension = ".fastq.gz"
    elif base_name.endswith(".fq.gz"):
        stem = base_name[:-6]
        extension = ".fq.gz"
    else:
        stem = original_path.stem
        extension = "".join(original_path.suffixes) or ".fastq"
    return output_dir / f"{sample_key}.{suffix}.{stem}.subset{extension}"


def write_summary(
    output_dir: Path,
    selected_transcripts: list[tuple[str, int, str]],
    aggregated_counts: dict[str, int],
    transcript_counts_by_sample: dict[str, dict[str, int]],
    subset_counts: dict[str, dict[str, int]],
) -> None:
    transcript_path = output_dir / "selected_transcripts.tsv"
    with transcript_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["label", "transcript", "aggregate_primary_alignments"])
        for transcript, count, label in selected_transcripts:
            writer.writerow([label, transcript, aggregated_counts[transcript]])

    per_sample_path = output_dir / "subset_summary.tsv"
    with per_sample_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "selected_read_names",
                "fastq_1_records_written",
                "fastq_2_records_written",
                "transcript_alignment_breakdown",
            ]
        )
        for sample_key in sorted(subset_counts):
            transcript_summary = ",".join(
                f"{transcript}:{transcript_counts_by_sample.get(sample_key, {}).get(transcript, 0)}"
                for transcript, _count, _label in selected_transcripts
            )
            writer.writerow(
                [
                    sample_key,
                    subset_counts[sample_key]["selected_read_names"],
                    subset_counts[sample_key]["fastq_1_records_written"],
                    subset_counts[sample_key]["fastq_2_records_written"],
                    transcript_summary,
                ]
            )


def write_subset_samplesheet(
    output_dir: Path,
    rows: list[dict[str, str]],
    match_column: str,
    output_paths: dict[str, dict[str, Path]],
    source_name: str,
) -> None:
    destination = output_dir / f"{Path(source_name).stem}.subset.csv"
    fieldnames = list(rows[0].keys())
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            sample_key = row[match_column]
            updated = dict(row)
            updated["fastq_1"] = str(output_paths[sample_key]["fastq_1"])
            if "fastq_2" in updated:
                updated["fastq_2"] = str(output_paths[sample_key].get("fastq_2", ""))
            writer.writerow(updated)


def main() -> int:
    args = parse_args()
    samplesheet_path = Path(args.samplesheet).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    bam_by_sample = parse_bam_mappings(args.bam)
    rows, match_column = load_samplesheet(samplesheet_path, args.sample_id_column)
    quantiles = parse_quantiles(args.quantiles)

    sample_keys = {row[match_column] for row in rows}
    missing_bams = sorted(sample_keys - set(bam_by_sample))
    unknown_bams = sorted(set(bam_by_sample) - sample_keys)
    if missing_bams:
        raise ValueError(f"Missing --bam mappings for samples: {', '.join(missing_bams)}")
    if unknown_bams:
        raise ValueError(f"--bam mappings not present in samplesheet: {', '.join(unknown_bams)}")

    engine = resolve_container_engine(args.container_engine)
    samtools = SamtoolsRunner(engine, args.container_image, args.singularity_image, args.container_platform)
    counts_by_sample = scan_alignment_counts(samtools, bam_by_sample)
    aggregated_counts = aggregate_counts(counts_by_sample)
    if args.transcript:
        selected_transcripts = validate_requested_transcripts(args.transcript, aggregated_counts)
    else:
        selected_transcripts = choose_transcripts(aggregated_counts, args.min_alignments, quantiles)

    selected_names_by_sample, transcript_counts_by_sample = collect_read_names(
        samtools,
        bam_by_sample,
        {transcript for transcript, _count, _label in selected_transcripts},
    )

    output_paths: dict[str, dict[str, Path]] = {}
    subset_counts: dict[str, dict[str, int]] = {}
    for row in rows:
        sample_key = row[match_column]
        fastq_1 = Path(row["fastq_1"]).expanduser().resolve()
        if not fastq_1.exists():
            raise FileNotFoundError(f"FASTQ path does not exist for sample '{sample_key}': {fastq_1}")
        selected_names = selected_names_by_sample[sample_key]
        fastq_1_output = build_output_fastq_path(output_dir, sample_key, fastq_1, "R1")
        fastq_1_written = subset_fastq(fastq_1, fastq_1_output, selected_names)
        sample_output_paths = {"fastq_1": fastq_1_output}

        fastq_2_written = 0
        fastq_2_value = row.get("fastq_2", "").strip()
        if fastq_2_value:
            fastq_2 = Path(fastq_2_value).expanduser().resolve()
            if not fastq_2.exists():
                raise FileNotFoundError(f"FASTQ path does not exist for sample '{sample_key}': {fastq_2}")
            fastq_2_output = build_output_fastq_path(output_dir, sample_key, fastq_2, "R2")
            fastq_2_written = subset_fastq(fastq_2, fastq_2_output, selected_names)
            sample_output_paths["fastq_2"] = fastq_2_output

        output_paths[sample_key] = sample_output_paths
        subset_counts[sample_key] = {
            "selected_read_names": len(selected_names),
            "fastq_1_records_written": fastq_1_written,
            "fastq_2_records_written": fastq_2_written,
        }

    write_summary(output_dir, selected_transcripts, aggregated_counts, transcript_counts_by_sample, subset_counts)
    write_subset_samplesheet(output_dir, rows, match_column, output_paths, samplesheet_path.name)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
