process R2DT {
    tag "$meta.id"
    label 'process_single'
    label 'process_long'

    container params.r2dt_container

    input:
    tuple val(meta), path(fold_dir), path(xml_files, stageAs: "xml_input*/*"), path(fasta), path(gtf)

    output:
    tuple val(meta), path("${prefix}_r2dt/"), optional: true, emit: diagrams
    tuple val(meta), path("r2dt_drawn_ids.txt"),              emit: drawn_ids
    tuple val(meta), path("${prefix}_r2dt.log"),                        emit: log
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    // Restrict template search to transcripts R2DT actually models (by GTF biotype); without a GTF
    // (e.g. viral synthetic references) fall back to drawing everything.
    def biotypes = task.ext.templatable_biotypes ?: ''
    def gtf_arg  = gtf ? "--gtf ${gtf} --allowed-biotypes \"${biotypes}\"" : ''
    """
    # ── 1. Extract sequences for transcripts present in fold dotbracket output ──
    r2dt_extract_sequences.py \\
        --fold-dir ${fold_dir} \\
        --fasta    ${fasta} \\
        ${gtf_arg} \\
        --out      r2dt_input.fa \\
        2>&1 | tee ${prefix}_r2dt.log

    if [[ ! -s r2dt_input.fa ]]; then
        echo "[R2DT] No sequences extracted — skipping." | tee -a ${prefix}_r2dt.log
        touch r2dt_drawn_ids.txt
        cat <<END_VERSIONS > versions.yml
"${task.process}":
    r2dt: \$(r2dt.py version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "unknown")
END_VERSIONS
        exit 0
    fi

    # ── 2. Run R2DT template-based layout ──────────────────────────────────────
    # Keep R2DT's FULL output + exit status: r2dt_status below is the authoritative
    # success/failure signal (see step 3) — the log itself is for diagnostics only.
    mkdir -p r2dt_raw
    set +e
    r2dt.py draw \\
        --skip_ribovore_filters \\
        $args \\
        r2dt_input.fa \\
        r2dt_raw \\
        > r2dt_draw.out 2>&1
    r2dt_status=\$?
    set -e
    grep -E '^(Analysing|Elapsed time|Visualising|Traveler crashed|Failed cmalign|Traceback|OSError|Errno|[Ee]rror|usage:)' \\
        r2dt_draw.out | tee -a ${prefix}_r2dt.log || true

    # ── 3. Overlay reactivities onto SVGs ──────────────────────────────────────
    # R2DT logs Tracebacks/"error" text to stderr for individual sequences that briefly
    # mismatch a template (e.g. depaired Infernal mapping, see r2dt-bio/R2DT#93) — routine
    # noise from a successful run. More importantly, r2dt.py returns a NON-ZERO exit when a
    # SINGLE sequence fails (e.g. Failed esl-sfetch, traveler json2svg NoneType) even after
    # drawing hundreds of others, so a non-zero status is NOT proof that nothing usable was
    # produced. Salvage whatever R2DT drew rather than discarding it all.
    if [[ \${r2dt_status} -ne 0 ]]; then
        echo "[R2DT] r2dt.py draw exited \${r2dt_status} — keeping any structures it did draw; error context:" | tee -a ${prefix}_r2dt.log
        grep -inE 'traceback|error|errno|exception|no such|read-only|permission|denied|cannot|traveler' r2dt_draw.out | tail -n 40 | tee -a ${prefix}_r2dt.log >&2
    fi

    # R2DT assembles results/svg only on a clean exit; a fatal crash on a single sequence (e.g. the
    # esl-sfetch failure above) aborts before assembly, stranding the per-template *.colored.svg it
    # already drew. Backfill results/svg from those so a late crash doesn't discard good diagrams.
    # No-op on a clean run: cp -n keeps R2DT's own assembled results/svg.
    mkdir -p r2dt_raw/results/svg
    find r2dt_raw -path r2dt_raw/results -prune -o -name '*.colored.svg' -print0 \\
        | xargs -0 -r -I{} cp -n {} r2dt_raw/results/svg/ || true

    # Always run the colour step over whatever SVGs R2DT produced. Don't pre-check
    # `find results/svg` in bash: that was racy — results/svg can hold hundreds of real SVGs
    # (confirmed on this cluster) yet read as empty to a `find` immediately afterwards, even
    # with a retry loop. r2dt_colour_svg.py globs --svg-dir itself (reliable) and no-ops
    # cleanly with "0 SVGs coloured" if nothing is there.
    mkdir -p ${prefix}_r2dt
    r2dt_colour_svg.py \\
        --svg-dir       r2dt_raw/results/svg \\
        --xml-search-dir . \\
        --out-dir       ${prefix}_r2dt \\
        2>&1 | tee -a ${prefix}_r2dt.log

    # Write the list of transcript IDs coloured, and record whether we produced any.
    # Remove the output directory if empty so optional: true suppresses publishing.
    : > r2dt_drawn_ids.txt
    drew_any=0
    for _f in ${prefix}_r2dt/*.svg; do
        [[ -f "\$_f" ]] || continue
        basename "\$_f" .svg >> r2dt_drawn_ids.txt
        drew_any=1
    done
    [[ \${drew_any} -eq 1 ]] || rm -rf ${prefix}_r2dt

    # R2DT is a best-effort template renderer and is NEVER fatal to the pipeline: ViennaRNA
    # renders every transcript R2DT does not cover. r2dt.py returns non-zero if it crashes on a
    # single sequence (e.g. traveler failing on one tRNA), and on such a crash it may leave
    # results/svg empty (nothing to salvage). Warn loudly but let the run finish via ViennaRNA.
    if [[ \${r2dt_status} -ne 0 && \${drew_any} -eq 0 ]]; then
        echo "[R2DT] r2dt.py draw FAILED (exit \${r2dt_status}) and produced no usable structures — falling back to ViennaRNA for this group." | tee -a ${prefix}_r2dt.log >&2
    fi

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    r2dt: \$(r2dt.py version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "unknown")
    python: \$(python3 --version | cut -d' ' -f2)
END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_r2dt
    touch ${prefix}_r2dt/stub_URS000035F234.svg
    echo "stub_URS000035F234" > r2dt_drawn_ids.txt
    echo "[R2DT] Extracted 1/1 sequences for template search" > ${prefix}_r2dt.log
    echo "[R2DT colour] 1 SVGs coloured, 0 skipped (no reactivity data or empty SVG)" >> ${prefix}_r2dt.log
    cat <<END_VERSIONS > versions.yml
"${task.process}":
    r2dt: 2.2.0
    python: 3.11.0
END_VERSIONS
    """
}
