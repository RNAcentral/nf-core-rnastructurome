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

    extracted = []
    hdr = None
    seq = []

    def flush(h, s):
        if h and h in ids:
            extracted.append((h, ''.join(s)))

    with open(args.fasta) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                flush(hdr, seq)
                hdr = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    flush(hdr, seq)

    with open(args.out, 'w') as fh:
        for h, s in extracted:
            fh.write(f'>{h}\n{s}\n')

    print(f'[R2DT] Extracted {len(extracted)}/{len(ids)} sequences for template search',
          file=sys.stderr)


if __name__ == '__main__':
    main()
