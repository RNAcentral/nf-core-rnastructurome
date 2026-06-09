process RNAFRAMEWORK_RFJACKKNIFE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(xml, stageAs: "input*/*")
    path reference

    output:
    tuple val(meta), path("${prefix}_jackknife/*.csv"),          emit: csv
    tuple val(meta), path("${prefix}_jackknife/*.pdf"), optional: true, emit: heatmap
    tuple val(meta), path("${prefix}_jackknife/rfjackknife.log"), optional: true, emit: log
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export TERM="\${TERM:-xterm}"
    log_tmp="\$(mktemp "${prefix}_jackknife.XXXXXX.log")"

    rf-jackknife \\
        -p ${task.cpus} \\
        -o ${prefix}_jackknife \\
        -ow \\
        -r ${reference} \\
        ${args} \\
        input*/ 2>&1 | tee "\${log_tmp}"

    mkdir -p ${prefix}_jackknife
    mv "\${log_tmp}" ${prefix}_jackknife/rfjackknife.log

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-jackknife 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_jackknife

    cat <<-END_CSV > ${prefix}_jackknife/${prefix}_fmi.csv
    slope,intercept,FMI
    1.0,-1.5,0.72
    1.0,-1.0,0.75
    1.2,-1.0,0.78
    END_CSV

    touch ${prefix}_jackknife/rfjackknife.log

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-jackknife 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
