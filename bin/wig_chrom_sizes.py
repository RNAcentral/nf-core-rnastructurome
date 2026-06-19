#!/usr/bin/env python3
"""Extract chromosome sizes from a WIG file (fixedStep or variableStep).

For fixedStep blocks, computes end = start + (count-1)*step + span - 1.
For variableStep blocks, uses max(position) + span - 1.
Writes a tab-separated chrom.sizes file suitable for wigToBigWig.
"""
from __future__ import annotations

import argparse
import re
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wig", help="Input WIG file")
    parser.add_argument("-o", "--output", default="-", help="Output chrom.sizes file (default: stdout)")
    return parser.parse_args()


def _int_header(line: str, key: str, default: int) -> int:
    m = re.search(rf"{key}=(\d+)", line)
    return int(m.group(1)) if m else default


def main() -> int:
    args = parse_args()

    sizes: dict[str, int] = {}

    # State for the current block
    chrom: str | None = None
    mode: str | None = None   # "fixed" or "variable"
    start = 1
    step = 1
    span = 1
    count = 0
    max_pos = 0

    def _flush() -> None:
        if chrom is None:
            return
        if mode == "fixed":
            end = start + (count - 1) * step + span - 1
        else:
            end = max_pos + span - 1
        sizes[chrom] = max(sizes.get(chrom, 0), end)

    with open(args.wig) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue

            if line.startswith("fixedStep"):
                _flush()
                chrom = re.search(r"chrom=(\S+)", line).group(1)
                mode  = "fixed"
                start = _int_header(line, "start", 1)
                step  = _int_header(line, "step",  1)
                span  = _int_header(line, "span",  1)
                count = 0
                max_pos = 0

            elif line.startswith("variableStep"):
                _flush()
                chrom = re.search(r"chrom=(\S+)", line).group(1)
                mode  = "variable"
                span  = _int_header(line, "span", 1)
                count = 0
                max_pos = 0

            else:
                if mode == "fixed":
                    count += 1
                elif mode == "variable":
                    pos = int(line.split()[0])
                    if pos > max_pos:
                        max_pos = pos

    _flush()

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
