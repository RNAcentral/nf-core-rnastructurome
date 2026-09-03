#!/usr/bin/env python3
"""Sanitize GTF transcript/gene IDs so they are safe for RNAframework and shell use.

RNAframework's XML parser interpolates attribute values straight into a regex
(lib/Data/IO/XML.pm), and rf-fold builds shell commands from transcript IDs.  IDs
containing regex/shell metacharacters — e.g. yeast tRNA names like ``tK(UUU)K`` —
break both: the parser spins forever on the un-escaped parentheses.  Replace any
character outside a conservative safe set with ``_`` in the ``transcript_id`` and
``gene_id`` attribute values.  ``.`` and ``-`` are preserved so version suffixes
and the yeast ORF-isoform normalisation (see fasta_sanitize_ids.py) keep working.

Usage: sanitize_gtf_ids.py <input_gtf[.gz]> <output_gtf>
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

# Characters kept verbatim; everything else in an ID becomes "_".
UNSAFE = re.compile(r"[^A-Za-z0-9._-]")
# Match an attribute (transcript_id / gene_id) and capture its quoted value.
ATTR = re.compile(r'((?:transcript_id|gene_id)\s+")([^"]*)(")')


def sanitize_id(value: str) -> str:
    return UNSAFE.sub("_", value)


def _open(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def sanitize_line(line: str) -> str:
    return ATTR.sub(lambda m: f"{m.group(1)}{sanitize_id(m.group(2))}{m.group(3)}", line)


def main() -> int:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_gtf[.gz]> <output_gtf>", file=sys.stderr)
        return 1

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    changed = 0
    with _open(input_path, "rt") as fin, output_path.open("wt", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            new_line = sanitize_line(line)
            if new_line != line:
                changed += 1
            fout.write(new_line)

    print(f"[sanitize_gtf_ids] Rewrote IDs on {changed} line(s) of {input_path.name}.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
