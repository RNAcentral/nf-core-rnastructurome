process AVERAGE_WIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), path(wigs, stageAs: "inputs/*.wig"), path(chrom_sizes, stageAs: "sizes/*.sizes")

    output:
    tuple val(meta), path("${prefix}.merged.wig"),  emit: merged_wig
    tuple val(meta), path("${prefix}_chrom.sizes"), emit: chrom_sizes
    path "versions.yml",                            emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python - <<'PY'
import sys
from pathlib import Path

wig_files = sorted(Path("inputs").glob("*.wig"))
output_path = Path("${prefix}.merged.wig")
sizes_out   = Path("${prefix}_chrom.sizes")

sizes_map = {}
for sf in sorted(Path("sizes").glob("*.sizes")):
    for line in sf.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        tid, length = line.split("\\t", 1)
        length = int(length)
        if tid not in sizes_map or length > sizes_map[tid]:
            sizes_map[tid] = length
with sizes_out.open("w") as fh:
    for tid, length in sorted(sizes_map.items()):
        fh.write(f"{tid}\\t{length}\\n")

if len(wig_files) == 1:
    output_path.write_bytes(wig_files[0].read_bytes())
else:
    def parse_wig(filepath):
        data = {}
        chrom = None
        is_fixed = False
        fixed_start = 1
        fixed_step  = 1
        pos_counter  = 0
        for line in filepath.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("track"):
                continue
            if line.startswith("variableStep"):
                parts = dict(p.split("=") for p in line.split() if "=" in p)
                chrom = parts.get("chrom")
                is_fixed = False
                data.setdefault(chrom, {})
            elif line.startswith("fixedStep"):
                parts = dict(p.split("=") for p in line.split() if "=" in p)
                chrom       = parts.get("chrom")
                fixed_start = int(parts.get("start", 1))
                fixed_step  = int(parts.get("step",  1))
                is_fixed    = True
                pos_counter  = 0
                data.setdefault(chrom, {})
            else:
                if is_fixed:
                    pos = fixed_start + pos_counter * fixed_step
                    pos_counter += 1
                    data[chrom][pos] = float(line)
                else:
                    p, v = line.split()
                    data[chrom][int(p)] = float(v)
        return data

    all_data = [parse_wig(f) for f in wig_files]
    all_chroms = set().union(*(d.keys() for d in all_data))

    with output_path.open("w") as fh:
        fh.write("track type=wiggle_0\\n")
        for chrom in sorted(all_chroms):
            all_pos = set().union(*(d[chrom].keys() for d in all_data if chrom in d))
            fh.write(f"variableStep chrom={chrom}\\n")
            for pos in sorted(all_pos):
                vals = [d[chrom][pos] for d in all_data if chrom in d and pos in d[chrom]]
                fh.write(f"{pos} {sum(vals)/len(vals):.6g}\\n")
PY

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.wig
    touch ${prefix}_chrom.sizes

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """
}
