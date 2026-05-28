#!/usr/bin/env python3
"""Average per-replicate merged WIG files; pass through unchanged when only one replicate.

WIG files are expected in an inputs/ subdirectory (Nextflow stageAs "inputs/*.wig" pattern).

Usage: average_wig.py <prefix>
"""
from __future__ import annotations

import sys
from pathlib import Path


def parse_wig(filepath: Path) -> dict[str, dict[int, float]]:
    data: dict[str, dict[int, float]] = {}
    chrom = None
    is_fixed = False
    fixed_start = 1
    fixed_step = 1
    pos_counter = 0
    for line in filepath.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("track"):
            continue
        if line.startswith("variableStep"):
            parts = dict(p.split("=") for p in line.split() if "=" in p)
            chrom = parts.get("chrom")
            is_fixed = False
            data.setdefault(chrom, {})
        elif line.startswith("fixedStep"):
            parts = dict(p.split("=") for p in line.split() if "=" in p)
            chrom = parts.get("chrom")
            fixed_start = int(parts.get("start", 1))
            fixed_step = int(parts.get("step", 1))
            is_fixed = True
            pos_counter = 0
            data.setdefault(chrom, {})
        else:
            if is_fixed:
                pos = fixed_start + pos_counter * fixed_step
                pos_counter += 1
                data[chrom][pos] = float(line)
            else:
                p, v = line.split()
                data[chrom][int(p)] = float(v)
    return data


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <prefix>", file=sys.stderr)
        sys.exit(1)

    prefix = sys.argv[1]
    wig_files = sorted(Path("inputs").glob("*.wig"))
    output_path = Path(f"{prefix}.merged.wig")

    if len(wig_files) == 1:
        output_path.write_bytes(wig_files[0].read_bytes())
        return

    all_data = [parse_wig(f) for f in wig_files]
    all_chroms = set().union(*(d.keys() for d in all_data))

    with output_path.open("w") as fh:
        fh.write("track type=wiggle_0\n")
        for chrom in sorted(all_chroms):
            all_pos = set().union(*(d[chrom].keys() for d in all_data if chrom in d))
            fh.write(f"variableStep chrom={chrom}\n")
            for pos in sorted(all_pos):
                vals = [d[chrom][pos] for d in all_data if chrom in d and pos in d[chrom]]
                fh.write(f"{pos} {sum(vals)/len(vals):.6g}\n")


if __name__ == "__main__":
    main()
