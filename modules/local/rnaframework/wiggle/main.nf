process RNAFRAMEWORK_RFWIGGLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/../fold/environment.yml"
    container params.rnaframework_container

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

    # rf-wiggle names its single output WIG after the input directory (xml.wig); rename it to the
    # group id so the published track is self-describing.
    produced_wig=\$(find ${prefix}_wiggle -maxdepth 1 -type f -name '*.wig' | head -1)
    if [[ -n "\${produced_wig}" && "\${produced_wig}" != "${prefix}_wiggle/${prefix}.wig" ]]; then
        mv "\${produced_wig}" "${prefix}_wiggle/${prefix}.wig"
    fi

    rnaframework_version=\$(rf-wiggle -h 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_wiggle
    touch ${prefix}_wiggle/${prefix}.wig

    rnaframework_version=\$(rf-wiggle -h 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """
}
