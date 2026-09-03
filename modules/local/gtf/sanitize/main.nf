process GTF_SANITIZE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${meta.id}.sanitized.gtf"), emit: gtf
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    """
    sanitize_gtf_ids.py "${gtf}" "${meta.id}.sanitized.gtf"
    """

    stub:
    """
    touch ${meta.id}.sanitized.gtf
    """
}
