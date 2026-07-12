process MERGE_SOURCE_URLS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(url_files, stageAs: "inputs/url_??.txt")

    output:
    tuple val(meta), path("ensembl_source_url.txt"), emit: urls
    path "versions.yml",                               emit: versions

    script:
    """
    cat inputs/* > ensembl_source_url.txt

    printf '"%s":\\n    python: %s\\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """

    stub:
    """
    touch ensembl_source_url.txt
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
