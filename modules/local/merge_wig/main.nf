process MERGE_WIG {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), path(wig), path(xml, stageAs: "xml/*")

    output:
    tuple val(meta), path("*.merged.wig"),  emit: merged_wig
    tuple val(meta), path("*_chrom.sizes"), emit: chrom_sizes
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python - <<'PY'
import re
import sys
from pathlib import Path

chrom_sizes = {}
for xml_file in sorted(Path("xml").glob("*.xml")):
    content = xml_file.read_text(encoding="utf-8")
    match = re.search(r'<transcript\\b[^>]+\\bid="([^"]+)"[^>]*\\blength="([0-9]+)"', content)
    if not match:
        match = re.search(r'<transcript\\b[^>]+\\blength="([0-9]+)"[^>]*\\bid="([^"]+)"', content)
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

sizes_path = Path("${prefix}_chrom.sizes")
with sizes_path.open("wt", encoding="utf-8") as handle:
    for transcript_id, length in sorted(chrom_sizes.items()):
        handle.write(f"{transcript_id}\\t{length}\\n")

wig_files = sorted(Path(".").glob("*.wig"))
if not wig_files:
    print("No WIG files were provided for merging", file=sys.stderr)
    sys.exit(1)

merged_wig = Path("${prefix}.merged.wig")
with merged_wig.open("wt", encoding="utf-8") as out_handle:
    out_handle.write("track type=wiggle_0\\n")
    for wig_file in wig_files:
        for line in wig_file.read_text(encoding="utf-8").splitlines(keepends=True):
            if not line.startswith("track "):
                out_handle.write(line)
PY

    printf '"%s":\\n    python: %s\\n' \\
        "${task.process}" \\
        "\$(python --version 2>&1 | sed 's/^Python //')" \\
        > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.wig
    touch ${prefix}_chrom.sizes

    printf '"%s":\\n    python: %s\\n' \\
        "${task.process}" \\
        "\$(python --version 2>&1 | sed 's/^Python //')" \\
        > versions.yml
    """
}
