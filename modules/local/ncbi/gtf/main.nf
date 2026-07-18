process NCBI_GTF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(fasta)
    path ncbi_gtf_script

    output:
    tuple val(meta), path("${meta.id}.annotation.gtf"), emit: gtf
    path "versions.yml",                                emit: versions

    script:
    """
    python "${ncbi_gtf_script}" \
        --fasta "${fasta}" \
        --output "${meta.id}.annotation.gtf"

    printf '"%s":\\n    python: %s\\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """

    stub:
    """
    touch ${meta.id}.annotation.gtf
    printf '"%s":\\n    python: %s\\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """
}
