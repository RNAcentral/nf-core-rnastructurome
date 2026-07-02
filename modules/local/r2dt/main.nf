process R2DT {
    tag "$meta.id"
    label 'process_single'
    label 'process_long'

    container params.r2dt_container

    input:
    tuple val(meta), path(fold_dir), path(xml_files, stageAs: "xml_input*/*"), path(fasta)
    path colour_script
    path extract_script

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
    """
    # ── 1. Extract sequences for transcripts present in fold dotbracket output ──
    python3 ${extract_script} \\
        --fold-dir ${fold_dir} \\
        --fasta    ${fasta} \\
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
    # mismatch a template (e.g. depaired Infernal mapping, see r2dt-bio/R2DT#93) — this is
    # routine noise from a successful run, NOT proof of failure. r2dt_status (the process's
    # own exit code) is the only reliable fatal signal; grepping for error keywords in the
    # full log previously caused runs with hundreds of good SVGs to be marked as failed.
    # The `find` below also retries briefly: on this filesystem, results/svg can take a
    # moment to become visible right after r2dt.py exits (same lag `scratch = false` was
    # meant to fix — it resurfaced here).
    svg_found=""
    for _attempt in 1 2 3 4 5; do
        if find r2dt_raw/results/svg -maxdepth 1 -name '*.svg' 2>/dev/null | grep -q .; then  # find avoids ARG_MAX with many per-transcript SVGs
            svg_found=1
            break
        fi
        sleep 2
    done

    if [[ -n "\${svg_found}" ]]; then
        mkdir -p ${prefix}_r2dt
        python3 ${colour_script} \\
            --svg-dir       r2dt_raw/results/svg \\
            --xml-search-dir . \\
            --out-dir       ${prefix}_r2dt \\
            2>&1 | tee -a ${prefix}_r2dt.log
    elif [[ \${r2dt_status} -ne 0 ]]; then
        echo "[R2DT] r2dt.py draw FAILED (exit \${r2dt_status}) — error context:" | tee -a ${prefix}_r2dt.log
        grep -inE 'traceback|error|errno|exception|no such|read-only|permission|denied|cannot|traveler' r2dt_draw.out | tail -n 40 | tee -a ${prefix}_r2dt.log >&2
        exit 1
    else
        echo "[R2DT] No template matches for any transcript — skipping." | tee -a ${prefix}_r2dt.log
    fi

    # Write list of transcript IDs that R2DT successfully drew
    # Remove output directory if empty so optional: true suppresses publishing
    if [[ -d ${prefix}_r2dt ]]; then
        for _f in ${prefix}_r2dt/*.svg; do
            [[ -f "\$_f" ]] || continue
            basename "\$_f" .svg
        done > r2dt_drawn_ids.txt
        find ${prefix}_r2dt -maxdepth 1 -name '*.svg' | grep -q . || rm -rf ${prefix}_r2dt
    else
        touch r2dt_drawn_ids.txt
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
