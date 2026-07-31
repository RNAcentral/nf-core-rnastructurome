#!/usr/bin/env python3
"""Sort FASTA records alphabetically by transcript ID, normalising yeast isoform names.

Usage: fasta_sort.py <input_fasta> <output_prefix> <organism>
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

YEAST_ISOFORM_PATTERN = re.compile(
    r'^(Y[A-P][LR][0-9]{3}[CW])-([A-Z])(_(?:mRNA|ncRNA|snRNA|snoRNA|rRNA|tRNA))$'
)
# Keep in sync with sanitize_gtf_ids.py: strip regex/shell-unsafe chars (e.g. yeast tRNA tK(UUU)K).
UNSAFE_ID_CHARS = re.compile(r'[^A-Za-z0-9._-]')


def open_fasta(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def normalize_header(header: str, organism: str) -> str:
    token, sep, remainder = header.partition(" ")
    if organism == "saccharomyces_cerevisiae":
        token = YEAST_ISOFORM_PATTERN.sub(r'\1_\2\3', token)
    token = UNSAFE_ID_CHARS.sub("_", token)
    return f"{token}{sep}{remainder}" if sep else token


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_fasta> <output_prefix> <organism>", file=sys.stderr)
        sys.exit(1)

    input_path = Path(sys.argv[1])
    output_path = Path(f"{sys.argv[2]}.sorted.fa")
    organism = sys.argv[3].strip().lower().replace(" ", "_")

    records = []
    current_header = None
    current_sequence: list[str] = []

    with open_fasta(input_path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_sequence)))
                current_header = normalize_header(line[1:], organism)
                current_sequence = []
            else:
                current_sequence.append(line)

    if current_header is not None:
        records.append((current_header, "".join(current_sequence)))

    records.sort(key=lambda record: record[0].split()[0])

    with output_path.open("wt", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for idx in range(0, len(sequence), 80):
                handle.write(sequence[idx:idx + 80] + "\n")


if __name__ == "__main__":
    main()
