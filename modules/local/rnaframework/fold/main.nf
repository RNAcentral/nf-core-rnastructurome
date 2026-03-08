process RNAFRAMEWORK_RFFOLD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'rnastructurome/rnaframework:2.9.6-r1'

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
    """
    rf-fold \\
        -p ${task.cpus} \\
        -o ${prefix}_fold \\
        -ow \\
        ${args} \\
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
