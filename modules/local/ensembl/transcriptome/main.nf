process ENSEMBL_TRANSCRIPTOME {
    tag "${meta.id}:${ensembl_species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12.11'

    input:
    tuple val(meta), val(ensembl_species)
    val ensembl_config_input

    output:
    tuple val(meta), path("${meta.id}.transcripts.fa.gz"), emit: fasta
    path "ensembl_source_urls.txt", emit: source_urls
    path "ensembl_warnings.log", optional: true, emit: warnings
    path "versions.yml", emit: versions

    script:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    python "${projectDir}/bin/ensembl_transcriptome.py" \
        --species "${ensembl_species}" \
        --release "${ensembl_config.ensembl_release}" \
        --base-url "${ensembl_config.ensembl_base_url}" \
        --output "${meta.id}.transcripts.fa.gz" \
        --source-urls "ensembl_source_urls.txt" \
        --warnings-log "ensembl_warnings.log"

    printf '%s\n' \
        '"${task.process}":' \
        '    ensembl_release: "${ensembl_config.ensembl_release}"' \
        > versions.yml
    """

    stub:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    touch ${meta.id}.transcripts.fa.gz
    printf '%s\n' \
        "stub://${meta.id}.transcripts.fa.gz" \
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
