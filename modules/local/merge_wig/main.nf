process MERGE_WIG {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), path(wig)

    output:
    tuple val(meta), path("*.merged.wig"), emit: merged_wig
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python - <<'PY'
from pathlib import Path

wig_files = sorted(Path(".").glob("*.wig"))
if not wig_files:
    import sys
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

    printf '"%s":\\n    python: %s\\n' \\
        "${task.process}" \\
        "\$(python --version 2>&1 | sed 's/^Python //')" \\
        > versions.yml
    """
}
