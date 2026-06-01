#!/usr/bin/env python3
"""
Extract FASTA sequences whose IDs appear in a dotbracket fold directory.

Usage
-----
    r2dt_extract_sequences.py --fold-dir <dir> --fasta <ref.fa> --out <out.fa>
"""

import argparse
import sys
from pathlib import Path


def main():
    """Extract FASTA sequences for R2DT input, filtering by length and dotbracket IDs."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fold-dir', required=True, type=Path,
                    help='Directory containing dotbracket/*.db files')
    ap.add_argument('--fasta',    required=True, type=Path,
                    help='Reference FASTA to extract sequences from')
    ap.add_argument('--out',      required=True, type=Path,
                    help='Output FASTA path')
    args = ap.parse_args()

    ids = {p.stem for p in (args.fold_dir / 'dotbracket').glob('*.db')}
    if not ids:
        print(f'[R2DT] No .db files found in {args.fold_dir}/dotbracket', file=sys.stderr)
        args.out.write_text('')
        return

    # Longest R2DT model is the human LSU (HS_LSU_3D, ~3305 nt); skip anything longer.
    # If R2DT is updated and adds longer templates, raise this value accordingly.
    _max_len = 3500

    extracted = []
    skipped   = 0
    hdr = None
    seq = []

    def flush(h, s):
        nonlocal skipped
        if h and h in ids:
            joined = ''.join(s)
            if len(joined) <= _max_len:
                extracted.append((h, joined))
            else:
                skipped += 1

    with open(args.fasta, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                flush(hdr, seq)
                hdr = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    flush(hdr, seq)

    with open(args.out, 'w', encoding='utf-8') as fh:
        for h, s in extracted:
            fh.write(f'>{h}\n{s}\n')

    print(
        f'[R2DT] Extracted {len(extracted)}/{len(ids)} sequences for template search'
        + (f' ({skipped} skipped: >{_max_len} nt)' if skipped else ''),
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
