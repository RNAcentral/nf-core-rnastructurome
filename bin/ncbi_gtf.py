#!/usr/bin/env python3
"""Generate a synthetic transcript-level GTF from an NCBI FASTA download.

Each sequence becomes a single-exon transcript spanning its full length, with
transcript_id and gene_id set to the accession ID.  This gives
rnaframework_dotplot2bp a valid GTF lookup for viral genomes where RNAframework
names dotplot files after the FASTA sequence IDs (e.g. NC_002023.1).
"""
from __future__ import annotations

import argparse
import gzip
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a synthetic single-exon GTF from an NCBI FASTA file."
    )
    parser.add_argument("--fasta", required=True, help="Input gzip-compressed FASTA")
    parser.add_argument("--output", required=True, help="Output GTF file path")
    return parser.parse_args()


def iter_fasta_lengths(fasta_path: str):
    """Yield (accession, length) for each sequence in a (gzip-compressed) FASTA."""
    opener = gzip.open if fasta_path.endswith(".gz") else open
    accession = None
    length = 0
    with opener(fasta_path, "rt", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if accession is not None:
                    yield accession, length
                accession = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if accession is not None:
        yield accession, length


def main() -> int:
    args = parse_args()
    records = list(iter_fasta_lengths(args.fasta))

    if not records:
        print(f"[NCBI_GTF] Error: no sequences found in {args.fasta!r}", file=sys.stderr)
        return 1

    with open(args.output, "w", encoding="utf-8") as out:
        for accession, length in records:
            attrs = f'gene_id "{accession}"; transcript_id "{accession}";'
            out.write(f'{accession}\tncbi\ttranscript\t1\t{length}\t.\t+\t.\t{attrs}\n')
            out.write(f'{accession}\tncbi\texon\t1\t{length}\t.\t+\t.\t{attrs}\n')

    print(
        f"[NCBI_GTF] Written {len(records)} transcript record(s) to {args.output!r}.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
