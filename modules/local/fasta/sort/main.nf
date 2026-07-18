process FASTA_SORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(fasta)
    path sort_script

    output:
    tuple val(meta), path("*.sorted.fa"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 "${sort_script}" "${fasta}" "${prefix}" "${meta.organism ?: meta.id}"

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python3 --version | cut -d' ' -f2)" \
        > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.fa

    printf '"%s":\n    python: %s\n' \
        "${task.process}" \
        "\$(python3 --version | cut -d' ' -f2)" \
        > versions.yml
    """
}
