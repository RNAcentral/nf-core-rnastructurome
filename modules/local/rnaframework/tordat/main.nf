process RNAFRAMEWORK_TORDAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), path(xml, stageAs: "xml/*"), path(fold_dir, stageAs: "fold_dir")
    path tordat_script

    output:
    tuple val(meta), path("${prefix}_rdat/*.rdat"), optional: true, emit: rdat
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python "${tordat_script}" \\
        --xml-dir xml \\
        --structures-dir fold_dir/structures \\
        --prefix "${prefix}"

    printf '"%s":\n    python: %s\n' \\
        "${task.process}" \\
        "\$(python --version 2>&1 | sed 's/^Python //')" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_rdat
    touch ${prefix}_rdat/stub.rdat

    printf '"%s":\n    python: %s\n' \\
        "${task.process}" \\
        "\$(python --version 2>&1 | sed 's/^Python //')" \\
        > versions.yml
    """
}
