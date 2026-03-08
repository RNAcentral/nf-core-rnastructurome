process RNAFRAMEWORK_RFNORM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'rnastructurome/rnaframework:2.9.6-r1'

    input:
    tuple val(meta), path(treated), path(untreated), path(denatured), path(rci_files)

    output:
    tuple val(meta), path("${prefix}_norm/*.xml"), emit: xml
    tuple val(meta), path("${prefix}_norm/plots/*.pdf"), optional: true, emit: plots
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    prefix            = task.ext.prefix ?: "${meta.id}"
    def untreated_arg = untreated ? "-u ${untreated}" : ''
    def denatured_arg = denatured ? "-d ${denatured}" : ''
    def treated_list  = treated instanceof List ? treated.join(' ') : "${treated}"
    """
    rf-norm \\
        -p ${task.cpus} \\
        -o ${prefix}_norm \\
        -ow \\
        ${args} \\
        -t ${treated_list} \\
        ${untreated_arg} \\
        ${denatured_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-norm 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_norm
    touch ${prefix}_norm/stub.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-norm 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
