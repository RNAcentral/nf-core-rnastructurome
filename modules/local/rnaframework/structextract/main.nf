process RNAFRAMEWORK_RFSTRUCTEXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(fold_dir), path(xmls, stageAs: 'xml_input/*')

    output:
    tuple val(meta), path("${prefix}_structextract/"),       optional: true, emit: motifs
    tuple val(meta), path("${prefix}.rfstructextract.log"),  optional: true, emit: log
    path "versions.yml",                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export TERM="\${TERM:-xterm}"

    rf-structextract \\
        -p ${task.cpus} \\
        -ro ${fold_dir} \\
        -xf xml_input \\
        -o ${prefix}_structextract \\
        -ow \\
        ${args} 2>&1 | tee "${prefix}.rfstructextract.log"

    rnaframework_version=\$(rf-structextract -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_structextract
    touch ${prefix}_structextract/${prefix}.db
    touch ${prefix}.rfstructextract.log

    rnaframework_version=\$(rf-structextract -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """
}
