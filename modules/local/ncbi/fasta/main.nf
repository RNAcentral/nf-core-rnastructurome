process NCBI_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), val(accessions)
    path ncbi_fasta_script

    output:
    tuple val(meta), path("${meta.id}.transcripts.fa.gz"), emit: fasta
    tuple val(meta), path("ncbi_source_accessions.txt"),   emit: source_accessions
    path "versions.yml",                                   emit: versions

    script:
    // accessions may be null for unseen organisms — the script falls back to esearch
    def acc_arg = accessions ? "--accessions \"${accessions}\"" : ""
    def args = task.ext.args ?: ''
    """
    python "${ncbi_fasta_script}" \
        ${acc_arg} \
        --organism "${meta.original_organism}" \
        --output "${meta.id}.transcripts.fa.gz" \
        --source-accessions "ncbi_source_accessions.txt" \
        ${args}

    printf '"%s":\\n    python: %s\\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """

    stub:
    """
    touch ${meta.id}.transcripts.fa.gz
    printf '%s\\n' "stub://${meta.id}" > ncbi_source_accessions.txt
    printf '"%s":\\n    python: %s\\n' \
        "${task.process}" \
        "\$(python --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """
}
