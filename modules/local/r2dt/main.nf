process R2DT {
    tag "$meta.id"
    label 'process_single'
    label 'process_long'

    container 'docker.io/rnacentral/r2dt:latest'

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
    mkdir -p r2dt_raw
    r2dt.py draw \\
        --skip_ribovore_filters \\
        $args \\
        r2dt_input.fa \\
        r2dt_raw \\
        2>&1 | grep -E '^(Analysing|Elapsed time|Traveler crashed|Failed cmalign|[Ee]rror|usage:)' | tee -a ${prefix}_r2dt.log || true

    # ── 3. Overlay reactivities onto SVGs ──────────────────────────────────────
    if find r2dt_raw/results/svg -maxdepth 1 -name '*.svg' 2>/dev/null | grep -q .; then  # find avoids ARG_MAX with many per-transcript SVGs
        mkdir -p ${prefix}_r2dt
        python3 ${colour_script} \\
            --svg-dir       r2dt_raw/results/svg \\
            --xml-search-dir . \\
            --out-dir       ${prefix}_r2dt \\
            2>&1 | tee -a ${prefix}_r2dt.log
    else
        echo "[R2DT] No template matches — skipping." | tee -a ${prefix}_r2dt.log
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
