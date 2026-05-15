process RNAFRAMEWORK_RFWIGGLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../fold/environment.yml"
    container 'ghcr.io/vicbeg/rnaframework:2.9.6-r5-runtime'

    input:
    tuple val(meta), path(xml, stageAs: "xml/*")

    output:
    tuple val(meta), path("${prefix}_wiggle/*.wig"), emit: wig
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    rf-wiggle \\
        -p ${task.cpus} \\
        -o ${prefix}_wiggle \\
        -ow \\
        ${args} \\
        xml/

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-wiggle 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_wiggle
    touch ${prefix}_wiggle/stub.wig

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-wiggle 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
