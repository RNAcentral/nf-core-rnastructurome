#!/usr/bin/env python3
"""Extract chromosome sizes from a fixedStep WIG file.

Counts the number of data lines per fixedStep block to determine each
chromosome (transcript) length, then writes a tab-separated chrom.sizes
file suitable for wigToBigWig.
"""
from __future__ import annotations

import argparse
import re
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wig", help="Input fixedStep WIG file")
    parser.add_argument("-o", "--output", default="-", help="Output chrom.sizes file (default: stdout)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    chrom: str | None = None
    count = 0
    sizes: dict[str, int] = {}

    with open(args.wig) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            m = re.match(r"fixedStep\s+chrom=(\S+)", line)
            if m:
                if chrom is not None:
                    sizes[chrom] = count
                chrom = m.group(1)
                count = 0
            else:
                count += 1

    if chrom is not None:
        sizes[chrom] = count

    out = open(args.output, "w") if args.output != "-" else sys.stdout
    try:
        for name in sorted(sizes):
            out.write(f"{name}\t{sizes[name]}\n")
    finally:
        if args.output != "-":
            out.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
