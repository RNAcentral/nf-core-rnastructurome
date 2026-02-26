process RNAFRAMEWORK_RFCOUNT {
    tag "$meta.id"
    label 'process_medium'

    container 'dincarnato/rnaframework:2.9.6'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta)

    output:
    tuple val(meta), path("*.rc"), emit: rc
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rf-count \\
        -p ${task.cpus} \\
        -f ${fasta} \\
        -o ./ \\
        -ow \\
        ${args} \\
        "${prefix}:${bam}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-count 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-count 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
