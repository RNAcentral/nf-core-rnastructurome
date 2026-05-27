process R2DT {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/rnacentral/r2dt:2.2'

    input:
    tuple val(meta), path(fold_dir), path(xml_files, stageAs: "xml_input*/*"), path(fasta)
    path colour_script
    path extract_script

    output:
    tuple val(meta), path("${prefix}_r2dt/"), optional: true, emit: diagrams
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
        --out      r2dt_input.fa

    if [[ ! -s r2dt_input.fa ]]; then
        echo "[R2DT] No sequences extracted — skipping." >&2
        mkdir -p ${prefix}_r2dt
        cat <<END_VERSIONS > versions.yml
"${task.process}":
    r2dt: \$(r2dt.py version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "unknown")
END_VERSIONS
        exit 0
    fi

    # ── 2. Run R2DT template-based layout ──────────────────────────────────────
    mkdir -p r2dt_raw
    r2dt.py draw \\
        --skip_ribovore_filters \\
        $args \\
        r2dt_input.fa \\
        r2dt_raw

    # ── 3. Overlay reactivities onto SVGs ──────────────────────────────────────
    mkdir -p ${prefix}_r2dt/svg

    if ls r2dt_raw/results/svg/*.svg 1>/dev/null 2>&1; then
        python3 ${colour_script} \\
            --svg-dir       r2dt_raw/results/svg \\
            --xml-search-dir . \\
            --out-dir       ${prefix}_r2dt/svg
    else
        echo "[R2DT] No template matches — output will be empty." >&2
    fi

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    r2dt: \$(r2dt.py version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "unknown")
    python: \$(python3 --version | cut -d' ' -f2)
END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_r2dt/svg
    touch ${prefix}_r2dt/svg/stub_URS000035F234.svg
    cat <<END_VERSIONS > versions.yml
"${task.process}":
    r2dt: 2.2.0
    python: 3.11.0
END_VERSIONS
    """
}
