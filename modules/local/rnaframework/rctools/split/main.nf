process RNAFRAMEWORK_RFRCTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(treated_rc), path(untreated_rc), path(gtf)

    output:
    tuple val(meta), path("chunks/treated/${prefix}_chunk_*.rc"), emit: treated_chunks
    tuple val(meta), path("chunks/untreated/${prefix}_chunk_*.rc"), optional: true, emit: untreated_chunks
    path "versions.yml",                                            emit: versions

    script:
    def chunk_size    = task.ext.chunk_size ?: 5000
    def gtf_feature   = task.ext.gtf_feature ?: 'exon'
    def gtf_attribute = task.ext.gtf_attribute ?: 'transcript_id'
    def min_coverage  = task.ext.min_coverage != null ? task.ext.min_coverage : 1
    prefix            = task.ext.prefix ?: "${meta.id}"
    """
    export TERM="\${TERM:-xterm}"

    # Build the RCI index locally. rf-rctools index writes <file>.rci into the (writable)
    # work dir even when the RC is a staged symlink, and rf-rctools view/extract need it.
    set -e
    rf-rctools index ${treated_rc}
    if [[ -n "${untreated_rc}" && -f "${untreated_rc}" ]]; then
        rf-rctools index ${untreated_rc}
    fi
    set +e

    # Pre-filter: keep only transcripts with real coverage in the treated RC. In genome mode
    # the RC contains EVERY GTF transcript, the vast majority with zero coverage, so chunking
    # and norming them all is wasted work. The view 'coverage' track (4th line per transcript)
    # is the ground truth (rf-rctools stats does not report usable per-transcript coverage).
    # A transcript is kept if any base has coverage >= min_coverage. min_coverage=0 keeps all.
    rf-rctools view ${treated_rc} 2>/dev/null | awk -v mc=${min_coverage} '
        NR % 4 == 1 { id = \$0; next }
        NR % 4 == 0 {
            n = split(\$0, cov, ",")
            for (i = 1; i <= n; i++) if (cov[i] + 0 >= mc) { print id; break }
        }
    ' > covered_ids.txt

    if [[ ! -s covered_ids.txt ]]; then
        echo "ERROR: no covered transcripts in ${treated_rc} (min_coverage=${min_coverage})." >&2
        exit 1
    fi

    # Derive per-transcript lengths from the GTF (spliced length = sum of exon lengths per
    # transcript). The RC was built by rf-rctools extract from this same GTF, so these lengths
    # match the RC exactly. rf-rctools stats does NOT report transcript lengths, so the GTF is
    # the authoritative source. Restrict to the covered transcripts identified above.
    python3 << 'PYEOF'
import os, re, sys

chunk_size = ${chunk_size}
feature    = "${gtf_feature}"
attribute  = "${gtf_attribute}"

attr_re = re.compile(r'${gtf_attribute}\\s+"([^"]+)"')

with open("covered_ids.txt") as f:
    covered = {line.strip() for line in f if line.strip()}

# transcript_id -> spliced length; dict preserves first-seen order (py3.7+).
lengths = {}
with open("${gtf}") as f:
    for line in f:
        if not line or line.startswith('#'):
            continue
        cols = line.rstrip('\\n').split('\\t')
        if len(cols) < 9 or cols[2] != feature:
            continue
        m = attr_re.search(cols[8])
        if not m:
            continue
        tx_id = m.group(1)
        if tx_id not in covered:
            continue
        try:
            start = int(cols[3])
            end   = int(cols[4])
        except ValueError:
            continue
        lengths[tx_id] = lengths.get(tx_id, 0) + (end - start + 1)

entries = [(tx_id, length) for tx_id, length in lengths.items() if length > 0]
if not entries:
    sys.exit("ERROR: no covered transcript lengths parsed from GTF '${gtf}' "
             "(feature='%s', attribute='%s'). Covered IDs may not match GTF %s values."
             % (feature, attribute, attribute))

# 4-column BED: chrom=transcript_id, 0, length, name=transcript_id.
# The name column makes rf-rctools extract preserve the clean transcript ID instead of
# renaming the region to '<id>_0-<end>', which would break treated/untreated pairing.
os.makedirs('chunk_beds', exist_ok=True)
for chunk_idx in range(0, len(entries), chunk_size):
    chunk = entries[chunk_idx:chunk_idx + chunk_size]
    idx   = chunk_idx // chunk_size
    with open(f'chunk_beds/chunk_{idx:04d}.bed', 'w') as bed:
        for tx_id, length in chunk:
            bed.write(f'{tx_id}\\t0\\t{length}\\t{tx_id}\\n')
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
