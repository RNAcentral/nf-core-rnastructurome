#!/usr/bin/env python3
"""Merge per-transcript Shannon entropy WIG files and build a transcript chrom.sizes from XML metadata.

XML files are expected in xml*/ subdirectories (Nextflow stageAs "xml*/*" pattern).
Duplicates across replicates are resolved by keeping the first occurrence per filename.

Usage: merge_shannon_wig.py <prefix>
"""

import re
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <prefix>", file=sys.stderr)
        sys.exit(1)

    prefix = sys.argv[1]

    seen_xml = {}
    for f in sorted(Path(".").glob("xml*/*.xml")):
        seen_xml.setdefault(f.name, f)

    chrom_sizes = {}
    for xml_file in sorted(seen_xml.values(), key=lambda f: f.name):
        content = xml_file.read_text(encoding="utf-8")
        match = re.search(r'<transcript\b[^>]+\bid="([^"]+)"[^>]*\blength="([0-9]+)"', content)
        if not match:
            match = re.search(r'<transcript\b[^>]+\blength="([0-9]+)"[^>]*\bid="([^"]+)"', content)
            if match:
                length_str, transcript_id = match.group(1), match.group(2)
            else:
                continue
        else:
            transcript_id, length_str = match.group(1), match.group(2)
        chrom_sizes[transcript_id] = int(length_str)

    if not chrom_sizes:
        print("No transcript lengths found in XML files", file=sys.stderr)
        sys.exit(1)

    sizes_path = Path(f"{prefix}_chrom.sizes")
    with sizes_path.open("wt", encoding="utf-8") as handle:
        for transcript_id, length in sorted(chrom_sizes.items()):
            handle.write(f"{transcript_id}\t{length}\n")

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
