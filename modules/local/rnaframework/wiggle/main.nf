process RNAFRAMEWORK_RFWIGGLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../fold/environment.yml"
    container 'ghcr.io/vicbeg/rnaframework:2.9.6-r5-runtime'

    input:
    tuple val(meta), path(xml, stageAs: "xml/*")

    output:
    tuple val(meta), path("${prefix}_wiggle/*.wig"), emit: wig
    tuple val(meta), path("${prefix}.merged.wig"),   emit: merged_wig
    tuple val(meta), path("${prefix}_chrom.sizes"),  emit: chrom_sizes
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    rf-wiggle \\
        -p ${task.cpus} \\
        -o ${prefix}_wiggle \\
        -ow \\
        ${args} \\
        xml/*.xml

    # Build chrom.sizes from XML length attributes
    python3 - <<'PY'
import re
import sys
from pathlib import Path

xml_files = sorted(Path("xml").glob("*.xml"))
chrom_sizes = {}
for xml_file in xml_files:
    content = xml_file.read_text(encoding="utf-8")
    m = re.search(r'<transcript\\b[^>]+\\bid="([^"]+)"[^>]*\\blength="([0-9]+)"', content)
    if not m:
        m = re.search(r'<transcript\\b[^>]+\\blength="([0-9]+)"[^>]*\\bid="([^"]+)"', content)
        if m:
            length_str, transcript_id = m.group(1), m.group(2)
        else:
            continue
    else:
        transcript_id, length_str = m.group(1), m.group(2)
    chrom_sizes[transcript_id] = int(length_str)

if not chrom_sizes:
    print("No transcript lengths found in XML files", file=sys.stderr)
    sys.exit(1)

sizes_path = Path("${prefix}_chrom.sizes")
with sizes_path.open("wt") as fh:
    for tid, length in sorted(chrom_sizes.items()):
        fh.write(f"{tid}\\t{length}\\n")

# Merge per-transcript WIG files into one, stripping per-file track headers
wig_files = sorted(Path("${prefix}_wiggle").glob("*.wig"))
merged_wig = Path("${prefix}.merged.wig")
with merged_wig.open("wt") as out_fh:
    out_fh.write("track type=wiggle_0\\n")
    for wig_file in wig_files:
        for line in wig_file.read_text(encoding="utf-8").splitlines(keepends=True):
            if not line.startswith("track "):
                out_fh.write(line)
PY

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-wiggle 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_wiggle
    touch ${prefix}_wiggle/stub.wig
    touch ${prefix}.merged.wig
    touch ${prefix}_chrom.sizes

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-wiggle 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
