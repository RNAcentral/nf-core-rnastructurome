process MERGE_BP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(bp_files, stageAs: "inputs/*.bp")
    path merge_script

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.bp"), optional: true, emit: bp
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 "${merge_script}" "${prefix}"

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python3 --version | cut -d' ' -f2)" \
        > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bp

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python3 --version | cut -d' ' -f2)" \
        > versions.yml
    """
}
