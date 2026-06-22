process RNAFRAMEWORK_RFRCTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(treated_rc), path(untreated_rc), path(rci)

    output:
    tuple val(meta), path("chunks/treated/${prefix}_chunk_*.rc"), emit: treated_chunks
    tuple val(meta), path("chunks/untreated/${prefix}_chunk_*.rc"), optional: true, emit: untreated_chunks
    path "versions.yml",                                            emit: versions

    script:
    def chunk_size    = task.ext.chunk_size ?: 5000
    prefix            = task.ext.prefix ?: "${meta.id}"
    """
    export TERM="\${TERM:-xterm}"

    # The .rci sidecar is staged alongside treated_rc by the channel (path(rci) input).
    # rf-rctools stats auto-discovers it; no explicit index call needed here.
    rf-rctools stats ${treated_rc} > rc_stats.txt

    python3 << 'PYEOF'
import os, sys

chunk_size = ${chunk_size}
entries = []  # list of (tx_id, length)

with open('rc_stats.txt') as f:
    for line in f:
        line = line.rstrip('\\n')
        if not line or line.startswith('#'):
            continue
        cols = line.split('\\t') if '\\t' in line else line.split()
        if len(cols) < 2:
            continue
        tx_id = cols[0]
        try:
            length = int(cols[1])
            if length > 0:
                entries.append((tx_id, length))
        except ValueError:
            continue  # skip header rows with non-numeric length column

if not entries:
    sys.exit("ERROR: could not parse transcript IDs and lengths from rf-rctools stats output. "
             "Check rc_stats.txt in the work directory.")

os.makedirs('chunk_beds', exist_ok=True)
for chunk_idx in range(0, len(entries), chunk_size):
    chunk = entries[chunk_idx:chunk_idx + chunk_size]
    idx   = chunk_idx // chunk_size
    with open(f'chunk_beds/chunk_{idx:04d}.bed', 'w') as f:
        for tx_id, length in chunk:
            f.write(f'{tx_id}\\t0\\t{length}\\n')
PYEOF

    mkdir -p chunks/treated chunks/untreated
    for bed in \$(ls chunk_beds/chunk_*.bed | sort); do
        chunk_name=\$(basename "\${bed}" .bed)
        rf-rctools extract \\
            -a "\${bed}" \\
            -o "chunks/treated/${prefix}_\${chunk_name}.rc" \\
            -ow \\
            ${treated_rc}
        rf-rctools index "chunks/treated/${prefix}_\${chunk_name}.rc"

        if [[ -n "${untreated_rc}" && -f "${untreated_rc}" ]]; then
            rf-rctools extract \\
                -a "\${bed}" \\
                -o "chunks/untreated/${prefix}_\${chunk_name}_untreated.rc" \\
                -ow \\
                ${untreated_rc}
            rf-rctools index "chunks/untreated/${prefix}_\${chunk_name}_untreated.rc"
        fi
    done

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-rctools 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    def chunk_size = task.ext.chunk_size ?: 5000
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p chunks/treated chunks/untreated
    touch chunks/treated/${prefix}_chunk_0000.rc

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-rctools 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
