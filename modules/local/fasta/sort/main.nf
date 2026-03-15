process FASTA_SORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12.11'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.sorted.fa"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python - <<'PY'
import gzip
import re
from pathlib import Path

input_path = Path("${fasta}")
output_path = Path("${prefix}.sorted.fa")
organism = "${meta.organism ?: meta.id}".strip().lower().replace(" ", "_")

yeast_isoform_pattern = re.compile(r'^(Y[A-P][LR][0-9]{3}[CW])-([A-Z])(_(?:mRNA|ncRNA|snRNA|snoRNA|rRNA|tRNA))\$')

def open_fasta(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")

def normalize_header(header: str) -> str:
    token, sep, remainder = header.partition(" ")
    if organism == "saccharomyces_cerevisiae":
        token = yeast_isoform_pattern.sub(r'\\1_\\2\\3', token)
    return f"{token}{sep}{remainder}" if sep else token

records = []
current_header = None
current_sequence = []

with open_fasta(input_path) as handle:
    for raw_line in handle:
        line = raw_line.rstrip("\\n")
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(current_sequence)))
            current_header = normalize_header(line[1:])
            current_sequence = []
        else:
            current_sequence.append(line)

if current_header is not None:
    records.append((current_header, "".join(current_sequence)))

records.sort(key=lambda record: record[0].split()[0])

with output_path.open("wt", encoding="utf-8") as handle:
    for header, sequence in records:
        handle.write(f">{header}\\n")
        for idx in range(0, len(sequence), 80):
            handle.write(sequence[idx:idx + 80] + "\\n")
PY

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.fa

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """
}
