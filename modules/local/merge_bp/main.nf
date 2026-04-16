process MERGE_BP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), path(bp_files, stageAs: "inputs/*.bp")

    output:
    tuple val(meta), path("${meta.id}_merged.bp"), optional: true, emit: bp
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python - <<'PY'
import sys
from pathlib import Path

bp_files = sorted(Path("inputs").glob("*.bp"))
output_path = Path("${prefix}_merged.bp")

EXPECTED_DATA_FIELDS = 6

def is_valid_data_line(fields):
    # Return True only if all six fields are non-empty, non-None strings
    # and the five numeric fields (cols 1-5) are integers.
    if len(fields) != EXPECTED_DATA_FIELDS:
        return False
    for f in fields:
        if not f or f.lower() == "none":
            return False
    # fields: chr  start  start  end  end  color_index
    try:
        int(fields[1])
        int(fields[2])
        int(fields[3])
        int(fields[4])
        int(fields[5])
    except ValueError:
        return False
    return True

color_lines = []
data_lines = []
color_written = False
skipped = 0

for bp_file in bp_files:
    for line in bp_file.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("color"):
            if not color_written:
                color_lines.append(stripped)
        else:
            fields = stripped.split("\\t")
            if is_valid_data_line(fields):
                data_lines.append(stripped)
            else:
                skipped += 1
                print(f"WARNING: skipping malformed/null bp line in {bp_file.name}: {stripped!r}", file=sys.stderr)
    color_written = True

if skipped:
    print(f"WARNING: {skipped} line(s) were dropped due to missing or non-numeric fields.", file=sys.stderr)

if data_lines:
    with output_path.open("w") as fh:
        for line in color_lines:
            fh.write(line + "\\n")
        for line in data_lines:
            fh.write(line + "\\n")
PY

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged.bp

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """
}
