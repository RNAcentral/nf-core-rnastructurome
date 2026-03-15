process RNAFRAMEWORK_RFNORM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/vicbeg/rnaframework:2.9.6-r2-runtime'

    input:
    tuple val(meta), path(treated), path(untreated), path(denatured), path(rci_files)

    output:
    tuple val(meta), path("${prefix}_norm/*.xml"), emit: xml
    tuple val(meta), path("${prefix}_norm/rfnorm.log"), optional: true, emit: log
    tuple val(meta), path("${prefix}_norm/plots/*.pdf"), optional: true, emit: plots
    path "versions.yml"                           , emit: versions

    script:
    def args          = task.ext.args ?: ''
    prefix            = task.ext.prefix ?: "${meta.id}"
    def untreated_arg = untreated ? "-u ${untreated}" : ''
    def denatured_arg = denatured ? "-d ${denatured}" : ''
    def treated_list  = treated instanceof List ? treated.join(' ') : "${treated}"
    """
    export TERM="\${TERM:-xterm}"
    rfnorm_log_tmp="${prefix}.rfnorm.log"

    rf-norm \\
        -p ${task.cpus} \\
        -o ${prefix}_norm \\
        -ow \\
        ${args} \\
        -t ${treated_list} \\
        ${untreated_arg} \\
        ${denatured_arg} 2>&1 | tee "\${rfnorm_log_tmp}"

    if [[ -d ${prefix}_norm ]]; then
        mv "\${rfnorm_log_tmp}" ${prefix}_norm/rfnorm.log
    fi

    printf '"%s":\n    rnaframework: %s\n' \\
        "${task.process}" \\
        "\$(rf-norm 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_norm
    touch ${prefix}_norm/stub.xml
    touch ${prefix}_norm/rfnorm.log

    printf '"%s":\n    rnaframework: %s\n' \\
        "${task.process}" \\
        "\$(rf-norm 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
