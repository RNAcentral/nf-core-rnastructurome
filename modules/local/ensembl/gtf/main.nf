process ENSEMBL_GTF {
    tag "${meta.id}:${ensembl_species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12.11'

    input:
    tuple val(meta), val(ensembl_species)
    val ensembl_config_input

    output:
    tuple val(meta), path("${meta.id}.annotation.gtf.gz"), emit: gtf
    path "ensembl_source_urls.txt", emit: source_urls
    path "versions.yml", emit: versions

    script:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    python "${projectDir}/bin/ensembl_gtf.py" \
        --species "${ensembl_species}" \
        --release "${ensembl_config.ensembl_release}" \
        --base-url "${ensembl_config.ensembl_base_url}" \
        --output "${meta.id}.annotation.gtf.gz" \
        --source-urls "ensembl_source_urls.txt"

    printf '%s\n' \
        '"${task.process}":' \
        '    ensembl_release: "${ensembl_config.ensembl_release}"' \
        > versions.yml
    """

    stub:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    touch ${meta.id}.annotation.gtf.gz
    printf '%s\n' \
        "stub://${meta.id}.annotation.gtf.gz" \
        > ensembl_source_urls.txt
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
