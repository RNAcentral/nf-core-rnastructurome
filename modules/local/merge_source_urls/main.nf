process MERGE_SOURCE_URLS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(url_files, stageAs: "inputs/url_??.txt")

    output:
    tuple val(meta), path("ensembl_source_url.txt"), emit: urls
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    """
    cat inputs/* > ensembl_source_url.txt
    """

    stub:
    """
    touch ensembl_source_url.txt
    """
}
