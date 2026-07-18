process AVERAGE_WIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(wigs, stageAs: "inputs/*.wig")
    path average_script

    output:
    tuple val(meta), path("${prefix}.merged.wig"), emit: merged_wig
    path "versions.yml",                           emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 "${average_script}" "${prefix}"

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python3 --version | cut -d' ' -f2)" \
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.wig

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python3 --version | cut -d' ' -f2)" \
        > versions.yml
    """
}
