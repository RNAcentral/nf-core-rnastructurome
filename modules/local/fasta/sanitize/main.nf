process FASTA_SANITIZE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.sanitized.fa"), emit: fasta
    tuple val("${task.process}"), val('python'), eval("python3 --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fasta_sanitize_ids.py "${fasta}" "${prefix}.sanitized.fa" "${meta.organism ?: meta.id}"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sanitized.fa
    """
}
