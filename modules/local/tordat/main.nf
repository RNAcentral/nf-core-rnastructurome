process RNAFRAMEWORK_TORDAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/python:3.12.11'

    input:
    tuple val(meta), path(xml, stageAs: "xml_inputs/rep??/*"), path(fold_dir, stageAs: "fold_dir")
    path tordat_script

    output:
    tuple val(meta), path("${prefix}_rdat/*.rdat"), optional: true, emit: rdat
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_str  = meta.fasta_name ?: 'unknown'
    def gtf_str    = meta.gtf_name   ?: ''
    def principle  = meta.principle  ?: ''
    def scoring_sm = meta.rfnorm_scoring_method != null ? meta.rfnorm_scoring_method.toString() : ''
    def norm_nm    = meta.rfnorm_norm_method    != null ? meta.rfnorm_norm_method.toString()    : ''
    def extra_args = task.ext.args ?: ''
    def gtf_arg    = gtf_str ? "--gtf \"${gtf_str}\"" : ''
    """
    python "${tordat_script}" \\
        --xml-dir xml_inputs \\
        --structures-dir fold_dir/dotbracket \\
        --prefix "${prefix}" \\
        --fasta "${fasta_str}" \\
        ${gtf_arg} \\
        --principle "${principle}" \\
        --rfnorm-scoring-method "${scoring_sm}" \\
        --rfnorm-norm-method "${norm_nm}" \\
        ${extra_args}

    printf '"%s":\n    python: %s\n' \\
        "${task.process}" \\
        "\$(python --version 2>&1 | sed 's/^Python //')" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_rdat
    touch ${prefix}_rdat/stub.rdat

    printf '"%s":\n    python: %s\n' \\
        "${task.process}" \\
        "\$(python --version 2>&1 | sed 's/^Python //')" \\
        > versions.yml
    """
}
