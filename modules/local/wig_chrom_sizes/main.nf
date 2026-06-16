process WIG_CHROM_SIZES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), path(wig)
    path chrom_sizes_script

    output:
    tuple val(meta), path("${prefix}.chrom.sizes"), emit: sizes
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python "${chrom_sizes_script}" \\
        "${wig}" \\
        -o "${prefix}.chrom.sizes"

    printf '"%s":\\n    python: %s\\n' \\
        "${task.process}" \\
        "\$(python --version 2>&1 | sed 's/^Python //')" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.chrom.sizes"
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
