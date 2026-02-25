process RNAFRAMEWORK_RFFOLD {
    tag "$meta.id"
    label 'process_high'

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
    """
    rf-fold \\
        -p ${task.cpus} \\
        -o ${prefix}_fold \\
        -ow \\
        ${args} \\
        ${xml_list}

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
