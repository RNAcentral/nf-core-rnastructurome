#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import statistics
import sys
from pathlib import Path
from xml.etree import ElementTree



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert RNAframework XML reactivity files + rf-fold .db structure files to RDAT format."
    )
    parser.add_argument("--xml-dir", required=True, help="Directory containing rf-norm XML files")
    parser.add_argument("--structures-dir", required=True, help="Directory containing rf-fold .db files")
    parser.add_argument("--prefix", required=True, help="Output prefix; files written to <prefix>_rdat/")
    return parser.parse_args()


def parse_xml(xml_path: Path) -> tuple[str, str, list[float | None]]:
    """Return (transcript_id, dna_sequence, reactivities).

    Reactivities are per-base floats; None indicates a missing/NaN value.
    """
    tree = ElementTree.parse(xml_path)
    root = tree.getroot()

    transcript_el = root.find("transcript")
    if transcript_el is None:
        raise ValueError(f"No <transcript> element found in {xml_path}")

    transcript_id = transcript_el.get("id", xml_path.stem)

    seq_el = transcript_el.find("sequence")
    dna_sequence = "".join((seq_el.text or "").split()) if seq_el is not None else ""

    react_el = transcript_el.find("reactivity")
    reactivities: list[float | None] = []
    if react_el is not None and react_el.text:
        for token in react_el.text.split(","):
            token = token.strip()
            if not token or token.lower() == "nan":
                reactivities.append(None)
            else:
                try:
                    reactivities.append(float(token))
                except ValueError:
                    reactivities.append(None)

    return transcript_id, dna_sequence, reactivities


def parse_db(db_path: Path) -> tuple[str, str]:
    """Return (rna_sequence, dot_bracket) from an rf-fold .db file.

    .db format:
        Line 1: >transcript_id
        Line 2: RNA sequence (may use T or U)
        Line 3: dot-bracket string, optionally followed by whitespace and MFE score
    """
    lines = [line.rstrip("\n") for line in db_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 3:
        raise ValueError(f"Malformed .db file (fewer than 3 non-empty lines): {db_path}")

    rna_sequence = lines[1].strip().upper().replace("T", "U")

    # Strip optional trailing MFE score like " (-770.27)"
    dot_bracket_field = lines[2].strip()
    dot_bracket = re.split(r"\s", dot_bracket_field)[0]

    return rna_sequence, dot_bracket


def resolve_db_path(structures_dir: Path, transcript_id: str, xml_stems: list[str]) -> Path | None:
    candidates: list[str] = [transcript_id, transcript_id.rsplit(".", 1)[0]]
    for stem in xml_stems:
        candidates.extend([stem, stem.rsplit(".", 1)[0]])

    seen: set[str] = set()
    for candidate in candidates:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        db_path = structures_dir / f"{candidate}.db"
        if db_path.exists():
            return db_path

    return None


def aggregate_reactivities(
    reactivity_sets: list[list[float | None]],
    length: int,
) -> tuple[list[float | None], list[float | None] | None]:
    aggregated: list[float | None] = []
    errors: list[float | None] | None = [] if len(reactivity_sets) > 1 else None

    for idx in range(length):
        values = [
            reactivities[idx]
            for reactivities in reactivity_sets
            if idx < len(reactivities) and reactivities[idx] is not None
        ]
        aggregated.append(statistics.fmean(values) if values else None)

        if errors is not None:
            if len(values) >= 2:
                errors.append(statistics.stdev(values) / math.sqrt(len(values)))
            else:
                errors.append(None)

    return aggregated, errors


def write_rdat(
    out_path: Path,
    transcript_id: str,
    rna_sequence: str,
    dot_bracket: str,
    reactivities: list[float | None],
    reactivity_errors: list[float | None] | None = None,
) -> None:
    length = len(rna_sequence)

    # REACTIVITY: one value per position; NaN for missing
    reactivity_values = []
    for i in range(length):
        val = reactivities[i] if i < len(reactivities) else None
        reactivity_values.append("NaN" if val is None else f"{val:.6g}")

    error_values = []
    if reactivity_errors is not None:
        for i in range(length):
            val = reactivity_errors[i] if i < len(reactivity_errors) else None
            error_values.append("NaN" if val is None else f"{val:.6g}")

    with out_path.open("wt", encoding="utf-8") as fh:
        fh.write(f"NAME\t{transcript_id}\n")
        fh.write(f"SEQUENCE\t{rna_sequence}\n")
        fh.write(f"STRUCTURE\t{dot_bracket}\n")
        fh.write(f"ANNOTATION_DATA:1\tmodifier:DMS\n")
        fh.write(f"REACTIVITY\t{' '.join(reactivity_values)}\n")
        if reactivity_errors is not None:
            fh.write(f"REACTIVITY_ERROR\t{' '.join(error_values)}\n")


def main() -> int:
    args = parse_args()
    xml_dir = Path(args.xml_dir)
    structures_dir = Path(args.structures_dir)
    out_dir = Path(f"{args.prefix}_rdat")
    out_dir.mkdir(parents=True, exist_ok=True)

    xml_files = sorted(path for path in xml_dir.rglob("*.xml") if path.is_file())
    if not xml_files:
        print(f"No XML files found in {xml_dir}", file=sys.stderr)
        return 1

    written = 0
    skipped = 0
    transcript_reactivities: dict[str, list[list[float | None]]] = {}
    transcript_xml_stems: dict[str, list[str]] = {}

    for xml_path in xml_files:
        try:
            transcript_id, _dna_sequence, reactivities = parse_xml(xml_path)
        except Exception as exc:
            print(f"WARNING: Could not parse {xml_path.name}: {exc}", file=sys.stderr)
            skipped += 1
            continue

        transcript_reactivities.setdefault(transcript_id, []).append(reactivities)
        transcript_xml_stems.setdefault(transcript_id, []).append(xml_path.stem)

    for transcript_id in sorted(transcript_reactivities):
        xml_stems = transcript_xml_stems.get(transcript_id, [])
        db_path = resolve_db_path(structures_dir, transcript_id, xml_stems)
        if db_path is None:
            print(
                f"WARNING: No .db file found for {transcript_id} in {structures_dir}; skipping.",
                file=sys.stderr,
            )
            skipped += 1
            continue

        try:
            rna_sequence, dot_bracket = parse_db(db_path)
        except Exception as exc:
            print(f"WARNING: Could not parse {db_path.name}: {exc}", file=sys.stderr)
            skipped += 1
            continue

        aggregated_reactivities, reactivity_errors = aggregate_reactivities(
            transcript_reactivities[transcript_id],
            len(rna_sequence),
        )
        out_path = out_dir / f"{transcript_id}.rdat"
        try:
            write_rdat(
                out_path,
                transcript_id,
                rna_sequence,
                dot_bracket,
                aggregated_reactivities,
                reactivity_errors,
            )
            written += 1
        except Exception as exc:
            print(f"WARNING: Could not write RDAT for {transcript_id}: {exc}", file=sys.stderr)
            skipped += 1

    print(f"Wrote {written} RDAT file(s), skipped {skipped}.", file=sys.stderr)

    if written == 0 and xml_files:
        print("ERROR: No RDAT files were produced.", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
