process ENSEMBL_TRANSCRIPTOME {
    tag "${meta.id}:${ensembl_species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), val(ensembl_species)
    val ensembl_config_input
    path ensembl_transcriptome_script

    output:
    tuple val(meta), path("${meta.id}.transcripts.fa.gz"), optional: true, emit: fasta
    tuple val(meta), path("ensembl_source_url.txt"),      optional: true, emit: source_urls
    tuple val(meta), path("${meta.id}.not_found"),         optional: true, emit: not_found
    path "ensembl_warnings.log",                           optional: true, emit: warnings
    path "versions.yml",                                                    emit: versions

    script:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    python "${ensembl_transcriptome_script}" \
        --species "${ensembl_species}" \
        --release "${ensembl_config.ensembl_release}" \
        --base-url "${ensembl_config.ensembl_base_url}" \
        --output "${meta.id}.transcripts.fa.gz" \
        --source-urls "ensembl_source_url.txt" \
        --warnings-log "ensembl_warnings.log" \
        --not-found-file "${meta.id}.not_found"

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
