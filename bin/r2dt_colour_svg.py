#!/usr/bin/env python3
"""
Overlay SHAPE / chemical-probing reactivities onto R2DT SVG diagrams.

Reads RNAframework rf-norm XML files and R2DT SVG outputs, then colours
each nucleotide according to its normalised reactivity value using the
standard SHAPE colouring scheme:

  NaN / missing / negative  →  grey     (#B1B3B6)
  0.00 – 0.40               →  black    (#000000)
  0.40 – 0.70               →  yellow   (#FFCD2F)
  > 0.70                    →  dark red (#9A2322)

When multiple XML files contain data for the same transcript (e.g. biological
replicates), per-position reactivities are averaged, ignoring NaN values.

Usage
-----
    r2dt_colour_svg.py \\
        --svg-dir       /path/to/r2dt/results/svg \\
        --xml-search-dir /path/to/work/dir \\
        --out-dir       /path/to/coloured_svg
"""

import argparse
import math
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

# SVG namespace used by R2DT / Traveler
_SVG_NS = 'http://www.w3.org/2000/svg'
_TITLE_RE = re.compile(r'^(\d+)\s')   # leading integer = 1-based position


# ── SHAPE colour scheme ──────────────────────────────────────────────────────

def _shape_colour(r) -> str:
    if r is None or (isinstance(r, float) and math.isnan(r)) or r < 0:
        return '#B1B3B6'   # grey     – no data
    if r <= 0.40:
        return '#000000'   # black    – low reactivity
    if r < 0.70:
        return '#FFCD2F'   # yellow   – medium reactivity
    return '#9A2322'       # dark red – high reactivity


# ── rf-norm XML parser ───────────────────────────────────────────────────────

def _parse_rfnorm_xml(path: Path) -> dict:
    """Return {transcript_id: [float|NaN, ...]} from one rf-norm XML file."""
    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        print(f'[WARN] Could not parse {path}: {exc}', file=sys.stderr)
        return {}

    result = {}
    for transcript in tree.getroot().iter('transcript'):
        tid = transcript.get('id') or transcript.findtext('id')
        if not tid:
            continue
        raw = transcript.findtext('reactivity') or transcript.findtext('values') or ''
        values = []
        for token in raw.strip().split(','):
            token = token.strip()
            if token.lower() in ('nan', 'na', '', 'null', 'none'):
                values.append(float('nan'))
            else:
                try:
                    values.append(float(token))
                except ValueError:
                    values.append(float('nan'))
        if values:
            result[tid] = values
    return result


def load_reactivities(xml_search_dir: Path, tids_needed: set) -> dict:
    """
    Load rf-norm XML reactivities only for transcripts in tids_needed.
    Builds a filename index first (fast), then parses only the relevant files.
    Multiple XML files for the same transcript are averaged position-by-position.
    """
    # Index all XML paths by transcript ID (filename stem = transcript ID)
    xml_index: dict = defaultdict(list)
    for xml_path in xml_search_dir.rglob('xml_input*/*.xml'):
        xml_index[xml_path.stem].append(xml_path)

    if not xml_index:
        print('[WARN] No reactivity data found in xml_input* directories.', file=sys.stderr)
        return {}

    raw: dict = defaultdict(list)   # tid -> list of reactivity vectors
    for tid in tids_needed:
        for xml_path in xml_index.get(tid, []):
            for _, values in _parse_rfnorm_xml(xml_path).items():
                raw[tid].append(values)

    averaged = {}
    for tid, reps in raw.items():
        if len(reps) == 1:
            averaged[tid] = reps[0]
        else:
            length = max(len(r) for r in reps)
            avg = []
            for i in range(length):
                vals = [
                    r[i] for r in reps
                    if i < len(r) and r[i] is not None and not math.isnan(r[i])
                ]
                avg.append(sum(vals) / len(vals) if vals else float('nan'))
            averaged[tid] = avg

    return averaged


# ── SVG colouring ────────────────────────────────────────────────────────────

def colour_svg(svg_path: Path, reactivities: list, out_path: Path) -> int:
    """
    Colour nucleotide text elements in an R2DT SVG by SHAPE reactivity.
    Nucleotides are identified by the leading integer in their <title> text,
    which gives the 1-based sequence position.
    Returns the number of nucleotides coloured.
    """
    ET.register_namespace('', _SVG_NS)
    ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')

    try:
        tree = ET.parse(svg_path)
    except ET.ParseError as exc:
        print(f'[WARN] Could not parse SVG {svg_path}: {exc}', file=sys.stderr)
        return 0

    coloured = 0
    for g in tree.getroot().iter(f'{{{_SVG_NS}}}g'):
        title_el = g.find(f'{{{_SVG_NS}}}title')
        if title_el is None or not title_el.text:
            continue
        m = _TITLE_RE.match(title_el.text.strip())
        if not m:
            continue
        pos = int(m.group(1)) - 1          # 0-based index
        if pos < 0 or pos >= len(reactivities):
            continue
        colour = _shape_colour(reactivities[pos])

        text_el = g.find(f'{{{_SVG_NS}}}text')
        if text_el is None:
            continue
        # Inline style takes precedence over class-based CSS fill
        existing = text_el.get('style', '')
        existing = re.sub(r'fill\s*:[^;]+;?\s*', '', existing).strip()
        new_style = f'fill: {colour};'
        if existing:
            new_style = f'{new_style} {existing}'
        text_el.set('style', new_style)
        coloured += 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(str(out_path), xml_declaration=True, encoding='unicode')
    return coloured


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    """Overlay SHAPE reactivities onto R2DT SVGs and report colouring statistics."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--svg-dir',        required=True, type=Path,
                    help='Directory containing R2DT SVG output files')
    ap.add_argument('--xml-search-dir', required=True, type=Path,
                    help='Root directory to search recursively for rf-norm XML files '
                         '(looks for xml_input*/*.xml pattern)')
    ap.add_argument('--out-dir',        required=True, type=Path,
                    help='Output directory for reactivity-coloured SVGs')
    args = ap.parse_args()

    tids_needed = {
        svg_path.stem.split('-')[0]
        for svg_path in args.svg_dir.glob('*.colored.svg')
    }
    if not tids_needed:
        print('[R2DT colour] 0 SVGs coloured, 0 skipped (no SVGs found)', file=sys.stderr)
        return

    reactivities = load_reactivities(args.xml_search_dir, tids_needed)
    if not reactivities:
        print('[ERROR] No reactivity data loaded — cannot colour SVGs.', file=sys.stderr)
        sys.exit(1)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    n_coloured = 0
    n_skipped  = 0

    for svg_path in sorted(args.svg_dir.glob('*.colored.svg')):
        # R2DT names SVGs as {URS_ID}-{TEMPLATE}.colored.svg; extract the URS ID
        tid = svg_path.stem.split('-')[0]
        if tid not in reactivities:
            n_skipped += 1
            continue
        n = colour_svg(svg_path, reactivities[tid], args.out_dir / (tid + '.svg'))
        if n:
            n_coloured += 1
        else:
            n_skipped += 1

    print(
        f'[R2DT colour] {n_coloured} SVGs coloured, {n_skipped} skipped '
        f'(no reactivity data or empty SVG)',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
