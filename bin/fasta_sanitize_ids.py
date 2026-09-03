#!/usr/bin/env python3
"""Normalize bespoke organism-specific FASTA header quirks before FASTA_SORT.

Bundles two workarounds needed only for organisms with problematic transcript IDs:
- yeast systematic ORF isoform names (e.g. ``YAL016C-A_mRNA``) have their isoform
  hyphen normalized to an underscore.
- any character outside a conservative safe set (e.g. parentheses in yeast tRNA
  names like ``tK(UUU)K``) is replaced with ``_`` so downstream RNAframework
  tooling that shells out on transcript IDs does not choke.

Keep in sync with sanitize_gtf_ids.py, which applies the same unsafe-character
rule to GTF transcript_id/gene_id attributes.

Usage: fasta_sanitize_ids.py <input_fasta[.gz]> <output_fasta> <organism>
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

YEAST_ISOFORM_PATTERN = re.compile(
    r'^(Y[A-P][LR][0-9]{3}[CW])-([A-Z])(_(?:mRNA|ncRNA|snRNA|snoRNA|rRNA|tRNA))$'
)
UNSAFE_ID_CHARS = re.compile(r'[^A-Za-z0-9._-]')


def open_fasta(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def sanitize_header(header: str, organism: str) -> str:
    token, sep, remainder = header.partition(" ")
    if organism == "saccharomyces_cerevisiae":
        token = YEAST_ISOFORM_PATTERN.sub(r'\1_\2\3', token)
    token = UNSAFE_ID_CHARS.sub("_", token)
    return f"{token}{sep}{remainder}" if sep else token


def main() -> int:
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_fasta> <output_fasta> <organism>", file=sys.stderr)
        return 1

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    organism = sys.argv[3].strip().lower().replace(" ", "_")

    with open_fasta(input_path) as handle, output_path.open("wt", encoding="utf-8") as out_handle:
        for line in handle:
            if line.startswith(">"):
                header = sanitize_header(line[1:].rstrip("\n"), organism)
                out_handle.write(f">{header}\n")
            else:
                out_handle.write(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
