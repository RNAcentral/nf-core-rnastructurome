process RNAFRAMEWORK_RFNORMFACTOR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(treated), path(untreated), path(denatured), path(rci_files)

    output:
    tuple val(meta), path("${prefix}.norm_factors.txt"), emit: factors
    tuple val(meta), path("${prefix}.rfnormfactor.log"), optional: true, emit: log
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    prefix            = task.ext.prefix ?: "${meta.id}"
    // rf-normfactor expects comma-separated file lists for -t/-u/-d (unlike rf-norm).
    def untreated_arg = untreated ? "-u ${untreated instanceof List ? untreated.join(',') : untreated}" : ''
    def denatured_arg = denatured ? "-d ${denatured instanceof List ? denatured.join(',') : denatured}" : ''
    def treated_list  = treated instanceof List ? treated.join(',') : "${treated}"
    """
    export TERM="\${TERM:-xterm}"
    rfnormfactor_log_tmp="${prefix}.rfnormfactor.log"

    rf-normfactor \\
        -p ${task.cpus} \\
        -o ${prefix}.norm_factors.txt \\
        -ow \\
        ${args} \\
        -t ${treated_list} \\
        ${untreated_arg} \\
        ${denatured_arg} 2>&1 | tee "\${rfnormfactor_log_tmp}"

    mv "\${rfnormfactor_log_tmp}" ${prefix}.rfnormfactor.log

    # rf-normfactor prints its version banner only under -h and exits non-zero; capture with '|| true'
    # so 'set -o pipefail' does not abort, and fall back to 'unknown' if no version is found.
    rnaframework_version=\$(rf-normfactor -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.norm_factors.txt
    touch ${prefix}.rfnormfactor.log

    rnaframework_version=\$(rf-normfactor -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """
}
