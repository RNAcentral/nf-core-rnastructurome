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
    tuple val(meta), path("${prefix}_structures/*.svg"), optional: true, emit: plots
    path "versions.yml",                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def rnaplot = task.ext.rnaplot ?: 'RNAplot'
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_structures

    sort "${drawn_ids}" > .r2dt_drawn.txt

    _process_db() {
        local _db="\$1"
        local _id=\$(basename "\$_db" .db)
        grep -qxF "\$_id" .r2dt_drawn.txt 2>/dev/null && return
        python3 "${xml_script}" "\${_id}" xml_input*/"\${_id}".xml > "\${_id}.shape" || true
        "${rnaplot}" --output-format=svg < "\$_db" || true
        if [[ -f "\${_id}_ss.svg" ]]; then
            mv "\${_id}_ss.svg" "${prefix}_structures/\${_id}.svg"
            python3 "${colour_script}" "\${_id}.shape" "${prefix}_structures/\${_id}.svg" || true
        fi
    }

    _n_jobs=0
    for _db in ${fold_dir}/dotbracket/*.db; do
        [[ -f "\$_db" ]] || continue
        _process_db "\$_db" &
        (( ++_n_jobs ))
        if (( _n_jobs >= ${task.cpus} )); then
            wait -n
            (( --_n_jobs ))
        fi
    done
    wait

    printf '"%s":\\n    viennarna: %s\\n' \\
        "${task.process}" \\
        "\$("${rnaplot}" --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo 'unknown')" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_structures
    touch ${prefix}_structures/stub_ENST00000000001.svg
    printf '"%s":\\n    viennarna: 2.6.4\\n' "${task.process}" > versions.yml
    """
}
