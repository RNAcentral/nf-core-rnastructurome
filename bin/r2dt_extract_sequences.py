#!/usr/bin/env python3
"""
Extract FASTA sequences from rf-fold dotbracket output for R2DT template search.

Usage
-----
    r2dt_extract_sequences.py --fold-dir <dir> --fasta <ref.fa> --out <out.fa>

--fasta is accepted for backward compatibility but is not used.
"""

import argparse
import sys
from pathlib import Path


def main():
    """Extract FASTA sequences for R2DT input from dotbracket .db files."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fold-dir', required=True, type=Path,
                    help='Directory containing dotbracket/*.db files')
    ap.add_argument('--fasta',    required=False, type=Path,
                    help='Accepted for backward compatibility; not used')
    ap.add_argument('--out',      required=True, type=Path,
                    help='Output FASTA path')
    args = ap.parse_args()

    db_files = sorted((args.fold_dir / 'dotbracket').glob('*.db'))
    if not db_files:
        print(f'[R2DT] No .db files found in {args.fold_dir}/dotbracket', file=sys.stderr)
        args.out.write_text('')
        return

    # Longest R2DT model is the human LSU (HS_LSU_3D, ~3305 nt); skip anything longer.
    _max_len = 3500

    extracted = []
    skipped   = 0

    for db_path in db_files:
        lines = db_path.read_text(encoding='utf-8').splitlines()
        if len(lines) < 2 or not lines[0].startswith('>'):
            continue
        tid = lines[0][1:].split()[0]
        seq = lines[1].strip()
        if not seq:
            continue
        if len(seq) <= _max_len:
            extracted.append((tid, seq))
        else:
            skipped += 1

    with open(args.out, 'w', encoding='utf-8') as fh:
        for h, s in extracted:
            fh.write(f'>{h}\n{s}\n')

    print(
        f'[R2DT] Extracted {len(extracted)}/{len(db_files)} sequences for template search'
        + (f' ({skipped} skipped: >{_max_len} nt)' if skipped else ''),
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
