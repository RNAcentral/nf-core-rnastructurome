process VIENNARNA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(fold_dir), path(xml, stageAs: "xml_input*/*"), path(drawn_ids)
    path xml_script
    path colour_script

    output:
    tuple val(meta), path("${prefix}_2D_structures/*.svg"), optional: true, emit: plots
    path "versions.yml",                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def rnaplot = task.ext.rnaplot ?: 'RNAplot'
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_2D_structures

    sort "${drawn_ids}" > .r2dt_drawn.txt

    for _db in ${fold_dir}/dotbracket/*.db; do
        [[ -f "\$_db" ]] || continue
        _id=\$(basename "\$_db" .db)

        if grep -qxF "\$_id" .r2dt_drawn.txt 2>/dev/null; then
            continue
        fi

        python3 "${xml_script}" "\${_id}" xml_input*/*.xml > "\${_id}.shape" || true

        "${rnaplot}" --output-format=svg < "\$_db" || true

        if [[ -f "\${_id}_ss.svg" ]]; then
            mv "\${_id}_ss.svg" "${prefix}_2D_structures/\${_id}.svg"

            if [[ -s "\${_id}.shape" ]]; then
                python3 "${colour_script}" "\${_id}.shape" "${prefix}_2D_structures/\${_id}.svg" || true
            fi
        fi
    done

    printf '"%s":\\n    viennarna: %s\\n' \\
        "${task.process}" \\
        "\$("${rnaplot}" --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo 'unknown')" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_2D_structures
    touch ${prefix}_2D_structures/stub_ENST00000000001.svg
    printf '"%s":\\n    viennarna: 2.6.4\\n' "${task.process}" > versions.yml
    """
}
