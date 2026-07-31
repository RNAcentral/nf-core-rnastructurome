#!/usr/bin/env python3
"""Merge per-transcript Shannon entropy WIG files into a single WIG file.

Usage: merge_shannon_wig.py <prefix>
"""

import sys
from pathlib import Path


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <prefix>", file=sys.stderr)
        sys.exit(1)

    prefix = sys.argv[1]

    wig_files = sorted(Path(".").glob("*.wig"))
    if not wig_files:
        print("No Shannon entropy WIG files were provided for merging", file=sys.stderr)
        sys.exit(1)

    merged_wig = Path(f"{prefix}.merged.wig")
    with merged_wig.open("wt", encoding="utf-8") as out_handle:
        out_handle.write("track type=wiggle_0\n")
        for wig_file in wig_files:
            for line in wig_file.read_text(encoding="utf-8").splitlines(keepends=True):
                if not line.startswith("track "):
                    out_handle.write(line)


if __name__ == "__main__":
    main()
