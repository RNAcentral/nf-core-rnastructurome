process RNAFRAMEWORK_DOTPLOT2BP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12.11'

    input:
    tuple val(meta), path(fold_dir), path(gtf)
    path dotplot2bp_script

    output:
    tuple val(meta), path("${meta.id}_bp/dotplot/*.bp"), optional: true, emit: bp
    tuple val(meta), path("${meta.id}_bp/conversion_warnings.log"), optional: true, emit: warnings
    path "versions.yml", emit: versions

    script:
    """
    python "${dotplot2bp_script}" \
        --organism "${meta.organism ?: meta.id}" \
        --prefix "${meta.id}" \
        --fold-dir "${fold_dir}" \
        --gtf "${gtf}"

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """

    stub:
    """
    mkdir -p ${meta.id}_bp/dotplot
    touch ${meta.id}_bp/dotplot/stub.bp

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """
}
