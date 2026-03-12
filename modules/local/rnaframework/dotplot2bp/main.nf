process RNAFRAMEWORK_DOTPLOT2BP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12.11'

    input:
    tuple val(meta), path(fold_dir), path(gtf)

    output:
    tuple val(meta), path("${meta.id}_bp/dotplot/*.bp"), optional: true, emit: bp
    tuple val(meta), path("${meta.id}_bp/conversion_warnings.log"), optional: true, emit: warnings
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python - <<'PY'
import gzip
import re
from pathlib import Path

organism = "${meta.organism ?: meta.id}".strip().lower().replace(" ", "_")
prefix = "${meta.id}"
fold_dir = Path("${fold_dir}")
gtf_path = Path("${gtf}")
out_dir = Path(f"{prefix}_bp")
dotplot_out_dir = out_dir / "dotplot"
warnings_path = out_dir / "conversion_warnings.log"
dotplot_in_dir = fold_dir / "dotplot"

dotplot_out_dir.mkdir(parents=True, exist_ok=True)
warnings = []

yeast_isoform_pattern = re.compile(r'^(Y[A-P][LR][0-9]{3}[CW])-([A-Z])(_(?:mRNA|ncRNA|snRNA|snoRNA|rRNA|tRNA))\$')
transcript_id_pattern = re.compile(r'transcript_id "([^"]+)"')

def normalize_transcript_id(transcript_id: str) -> str:
    if organism == "saccharomyces_cerevisiae":
        return yeast_isoform_pattern.sub(r'\\1_\\2\\3', transcript_id)
    return transcript_id

def transcript_id_candidates(transcript_id: str):
    normalized = normalize_transcript_id(transcript_id)
    candidates = [normalized]
    if "." in normalized:
        candidates.append(normalized.split(".", 1)[0])
    return candidates

def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")

transcripts = {}
with open_text(gtf_path) as handle:
    for raw_line in handle:
        if not raw_line or raw_line.startswith("#"):
            continue
        fields = raw_line.rstrip("\\n").split("\\t")
        if len(fields) < 9 or fields[2] != "exon":
            continue
        match = transcript_id_pattern.search(fields[8])
        if not match:
            continue
        transcript_id = normalize_transcript_id(match.group(1))
        seqname = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        entry = transcripts.setdefault(transcript_id, {"seqname": seqname, "strand": strand, "exons": []})
        entry["exons"].append((start, end))

for transcript_id, entry in transcripts.items():
    exons = entry["exons"]
    if entry["strand"] == "-":
        exons.sort(key=lambda exon: exon[0], reverse=True)
    else:
        exons.sort(key=lambda exon: exon[0])

def map_transcript_pos(entry, position: int) -> int:
    cursor = 1
    for start, end in entry["exons"]:
        exon_len = end - start + 1
        if cursor <= position < cursor + exon_len:
            offset = position - cursor
            if entry["strand"] == "-":
                return end - offset
            return start + offset
        cursor += exon_len
    raise ValueError(f"Transcript position {position} exceeds transcript length")

dotplot_paths = sorted(dotplot_in_dir.glob("*.dp")) if dotplot_in_dir.is_dir() else []
bp_count = 0

for dotplot_path in dotplot_paths:
    transcript_id = normalize_transcript_id(dotplot_path.stem)
    entry = None
    matched_transcript_id = None
    for candidate in transcript_id_candidates(transcript_id):
        entry = transcripts.get(candidate)
        if entry:
            matched_transcript_id = candidate
            break
    if not entry:
        warnings.append(f"Missing transcript_id '{transcript_id}' in annotation for {dotplot_path.name}; skipping.")
        continue
    if matched_transcript_id != transcript_id:
        warnings.append(
            f"Matched dotplot transcript_id '{transcript_id}' to annotation transcript_id '{matched_transcript_id}' for {dotplot_path.name}."
        )

    output_path = dotplot_out_dir / f"{transcript_id}.bp"
    with dotplot_path.open("rt", encoding="utf-8") as reader, output_path.open("wt", encoding="utf-8") as writer:
        first_line = reader.readline()
        header_line = reader.readline()
        if not first_line or not header_line:
            warnings.append(f"Malformed dotplot file {dotplot_path.name}; skipping.")
            output_path.unlink(missing_ok=True)
            continue
        writer.write("color\\t31\\t119\\t180\\tRNAFramework dotplot\\n")
        converted_any = False
        for raw_line in reader:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\\t")
            if len(fields) < 3:
                continue
            left_pos = int(fields[0])
            right_pos = int(fields[1])
            left_genome = map_transcript_pos(entry, left_pos)
            right_genome = map_transcript_pos(entry, right_pos)
            start, end = sorted((left_genome, right_genome))
            writer.write(f"{entry['seqname']}\\t{start}\\t{start}\\t{end}\\t{end}\\t0\\n")
            converted_any = True

    if converted_any:
        bp_count += 1
    else:
        output_path.unlink(missing_ok=True)
        warnings.append(f"No convertible base-pair records found in {dotplot_path.name}; skipping.")

if warnings:
    warnings_path.write_text("\\n".join(warnings) + "\\n", encoding="utf-8")
elif warnings_path.exists():
    warnings_path.unlink()

if dotplot_paths and bp_count == 0:
    raise SystemExit("No .bp files were generated from the available .dp files.")
PY

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """

    stub:
    """
    mkdir -p ${meta.id}_bp/dotplot
    touch ${meta.id}_bp/dotplot/stub.bp

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """
}
