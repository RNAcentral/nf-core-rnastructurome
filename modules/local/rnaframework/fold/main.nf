process RNAFRAMEWORK_RFFOLD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'dincarnato/rnaframework:2.9.6'

    input:
    tuple val(meta), path(xml)

    output:
    tuple val(meta), path("${prefix}_fold/"), emit: structures
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def xml_list = xml instanceof List ? xml.join(' ') : "${xml}"
    def auto_window_min_len = params.rffold_auto_window_min_len != null ? params.rffold_auto_window_min_len as Integer : 10000
    """
    auto_window_arg=""
    if [[ ! " ${args} " =~ [[:space:]]-w([[:space:]]|$) ]]; then
        max_xml_len=\$(grep -hoE 'length="[0-9]+"' ${xml_list} | sed -E 's/.*"([0-9]+)".*/\\1/' | sort -nr | head -1 || true)
        max_xml_len=\${max_xml_len:-0}
        if [[ "\${max_xml_len}" -ge ${auto_window_min_len} ]]; then
            echo "[RNAFRAMEWORK_RFFOLD] Auto-enabling windowed folding (-w): max transcript length \${max_xml_len} >= ${auto_window_min_len}" >&2
            auto_window_arg="-w"
        fi
    fi

    rf-fold \\
        -p ${task.cpus} \\
        -o ${prefix}_fold \\
        -ow \\
        ${args} \\
        \${auto_window_arg} \\
        ${xml_list}

    # rf-fold can return exit 0 even when all folds fail and details are written to error.out.
    # Treat this as a hard failure so the pipeline does not continue with empty fold outputs.
    if [[ -s ${prefix}_fold/error.out ]]; then
        echo "[RNAFRAMEWORK_RFFOLD] rf-fold reported errors:" >&2
        cat ${prefix}_fold/error.out >&2
        exit 1
    fi

    if [[ ! -d ${prefix}_fold/structures ]] || ! find ${prefix}_fold/structures -type f -print -quit | grep -q .; then
        echo "[RNAFRAMEWORK_RFFOLD] No structure files were produced in ${prefix}_fold/structures." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-fold 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_fold

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-fold 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
