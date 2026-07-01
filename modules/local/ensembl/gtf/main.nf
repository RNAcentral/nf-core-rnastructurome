process ENSEMBL_GTF {
    tag "${meta.id}:${ensembl_species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), val(ensembl_species)
    val ensembl_config_input
    path ensembl_gtf_script

    output:
    tuple val(meta), path("${meta.id}.annotation.gtf"),   optional: true, emit: gtf
    tuple val(meta), path("${meta.id}.annotation.gtf.gz"), optional: true, emit: gtf_gz
    tuple val(meta), path("ensembl_source_url.txt"),       optional: true, emit: source_urls
    tuple val(meta), path("${meta.id}.not_found"),          optional: true, emit: not_found
    path "versions.yml",                                                     emit: versions

    script:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    python "${ensembl_gtf_script}" \
        --species "${ensembl_species}" \
        --release "${ensembl_config.ensembl_release}" \
        --base-url "${ensembl_config.ensembl_base_url}" \
        --output "${meta.id}.annotation.gtf.gz" \
        --source-urls "ensembl_source_url.txt" \
        --not-found-file "${meta.id}.not_found"

    [[ -f "${meta.id}.annotation.gtf.gz" ]] && gzip -dc "${meta.id}.annotation.gtf.gz" > "${meta.id}.annotation.gtf" || true

    printf '%s\n' \
        '"${task.process}":' \
        '    ensembl_release: "${ensembl_config.ensembl_release}"' \
        > versions.yml
    """

    stub:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    touch ${meta.id}.annotation.gtf.gz
    touch ${meta.id}.annotation.gtf
    printf '%s\n' \
        "stub://${meta.id}.annotation.gtf.gz" \
        > ensembl_source_url.txt
    printf '%s\n' \
        '"${task.process}":' \
        '    ensembl_release: "${ensembl_config.ensembl_release}"' \
        > versions.yml
    """
}

def defaultEnsemblConfig() {
    [
        ensembl_release : 'current',
        ensembl_base_url: 'https://ftp.ensembl.org/pub'
    ]
}

process GTF_SANITIZE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), path(gtf)
    path sanitize_script

    output:
    tuple val(meta), path("${meta.id}.sanitized.gtf"), emit: gtf
    path "versions.yml",                               emit: versions

    script:
    """
    python3 "${sanitize_script}" "${gtf}" "${meta.id}.sanitized.gtf"

    printf '"%s":\\n    python: %s\\n' \
        "${task.process}" \
        "\$(python3 --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """

    stub:
    """
    touch ${meta.id}.sanitized.gtf
    printf '"%s":\\n    python: %s\\n' \
        "${task.process}" \
        "\$(python3 --version 2>&1 | sed 's/^Python //')" \
        > versions.yml
    """
}
