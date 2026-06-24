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

    # Build FASTA index keyed by both the full ID and the version-stripped base ID.
    # This lets genome-route fold IDs (from GTF, e.g. ENST00000389680) match FASTA
    # entries that carry a version suffix (e.g. ENST00000389680.2), and vice versa.
    import re as _re
    _ver_re = _re.compile(r'^(.+)\.\d+$')

    fasta_index: dict = {}  # id -> (display_id, sequence)
    hdr = None
    seq: list = []

    def _index(h, s):
        if not h:
            return
        joined = ''.join(s)
        fasta_index[h] = (h, joined)
        m = _ver_re.match(h)
        if m:
            base = m.group(1)
            fasta_index.setdefault(base, (h, joined))

    with open(args.fasta, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                _index(hdr, seq)
                hdr = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    _index(hdr, seq)

    extracted = []
    skipped   = 0
    for fold_id in ids:
        entry = fasta_index.get(fold_id)
        if entry is None:
            continue
        display_id, joined = entry
        if len(joined) <= _max_len:
            extracted.append((display_id, joined))
        else:
            skipped += 1

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
