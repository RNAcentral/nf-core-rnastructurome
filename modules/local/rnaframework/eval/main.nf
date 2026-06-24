process RNAFRAMEWORK_RFEVAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(xml, stageAs: "xml_input*/*")
    path structures

    output:
    tuple val(meta), path("${prefix}_rfeval/*.csv"),          optional: true, emit: csv
    tuple val(meta), path("${prefix}_rfeval/plots/*.pdf"),    optional: true, emit: plots
    tuple val(meta), path("${prefix}_rfeval/rfeval.log"),     optional: true, emit: log
    path "versions.yml",                                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export TERM="\${TERM:-xterm}"
    log_tmp="\$(mktemp "${prefix}_rfeval.XXXXXX.log")"

    rf-eval \\
        -p ${task.cpus} \\
        -o ${prefix}_rfeval \\
        -ow \\
        -s ${structures} \\
        -r xml_input*/ \\
        ${args} 2>&1 | tee "\${log_tmp}"

    # Strip ANSI/CR progress noise; drop [+] status and | progress-bar lines
    perl -pe 's/\\r/\\n/g; s/\\e\\[[0-9;]*[A-Za-z]//g' "\${log_tmp}" \\
        | grep -vE '^\\[+\\]|^\\|' \\
        | cat -s \\
        > "\${log_tmp}.clean"
    mv "\${log_tmp}.clean" "\${log_tmp}"

    mkdir -p ${prefix}_rfeval
    mv "\${log_tmp}" ${prefix}_rfeval/rfeval.log

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-eval 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_rfeval/plots

    cat <<-END_CSV > ${prefix}_rfeval/${prefix}_eval.csv
    transcript,unpaired_coeff,DSCI,AUROC
    ENST00000000001,0.82,0.79,0.85
    ENST00000000002,0.74,0.71,0.78
    END_CSV

    touch ${prefix}_rfeval/rfeval.log

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-eval 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
