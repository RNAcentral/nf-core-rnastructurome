#!/usr/bin/env python3
"""Concatenate per-transcript rf-wiggle WIG files into a single track for a sample group.

WIG files are expected in the current working directory (staged flat by Nextflow).
The output keeps one leading "track type=wiggle_0" header and drops the per-file headers.

Usage: merge_wig.py <prefix>
"""
from __future__ import annotations

import sys
from pathlib import Path


def main(prefix: str) -> int:
    wig_files = sorted(Path(".").glob("*.wig"))
    if not wig_files:
        print("No WIG files were provided for merging", file=sys.stderr)
        return 1

    merged_wig = Path(f"{prefix}.merged.wig")
    with merged_wig.open("wt", encoding="utf-8") as out_handle:
        out_handle.write("track type=wiggle_0\n")
        for wig_file in wig_files:
            for line in wig_file.read_text(encoding="utf-8").splitlines(keepends=True):
                if not line.startswith("track "):
                    out_handle.write(line)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: merge_wig.py <prefix>", file=sys.stderr)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
