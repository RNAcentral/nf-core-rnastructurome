process RNAFRAMEWORK_RFRCTOOLS_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(rc, stageAs: "input/*"), path(rci, stageAs: "input/*")
    tuple val(meta_ref), path(gtf)

    output:
    tuple val(meta), path("${prefix}_rctools_extract/${prefix}.rc"),                       emit: rc
    tuple val(meta), path("${prefix}_rctools_extract/${prefix}.rc.rci"), optional: true,   emit: rci
    path "versions.yml",                                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rctools_extract"
    """
    export TERM="\${TERM:-xterm}"
    mkdir -p ${outdir}

    rf-rctools extract \\
        -a ${gtf} \\
        -o ${outdir}/${prefix}.rc \\
        -ow \\
        ${args} \\
        input/*.rc

    rf-rctools index ${outdir}/${prefix}.rc

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-rctools 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rctools_extract"
    """
    mkdir -p ${outdir}
    touch ${outdir}/${prefix}.rc
    touch ${outdir}/${prefix}.rc.rci

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-rctools 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
