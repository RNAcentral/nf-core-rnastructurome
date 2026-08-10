// Version probe must not abort the task or leak a shell error into versions.yml,
// so failures fall back to 'unknown' the same way the rf-eval probe does.
def pythonVersionCmd() {
    '''python_version=$(python3 --version 2>/dev/null | sed 's/^Python //') || true
    printf '    python: %s\\n' "${python_version:-unknown}" >> versions.yml'''
}

process RNAFRAMEWORK_RFEVAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:affb2f7a4bac9a7a' :
        'community.wave.seqera.io/library/rnaframework_findutils:3db7cd7277dc8f08' }"

    input:
    tuple val(meta), path(xml, stageAs: "xml_input*/*")
    path structures
    path windows
    path window_script

    output:
    tuple val(meta), path("${prefix}_rfeval/*.csv"),          optional: true, emit: csv
    tuple val(meta), path("${prefix}_rfeval/plots/*.pdf"),    optional: true, emit: plots
    tuple val(meta), path("${prefix}_rfeval/rfeval.log"),     optional: true, emit: log
    tuple val(meta), path("${prefix}_rfeval_windows/*.xml"),  optional: true, emit: windows
    path "versions.yml",                                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    // With a windows manifest, slice reactivities to each reference region first so a
    // sub-region structure is scored against a matching XML, not the full transcript.
    def window_cmd = windows
        ? "python3 \"${window_script}\" --windows \"${windows}\" --xml-glob 'xml_input*/*.xml' --outdir ${prefix}_rfeval_windows"
        : ''
    def reactivity_dir = windows ? "${prefix}_rfeval_windows/" : 'xml_input*/'
    def python_version = windows ? pythonVersionCmd() : ''
    """
    export TERM="\${TERM:-xterm}"
    log_tmp="\$(mktemp "${prefix}_rfeval.XXXXXX.log")"

    ${window_cmd}

    rf-eval \\
        -p ${task.cpus} \\
        -o ${prefix}_rfeval \\
        -ow \\
        -s ${structures} \\
        -r ${reactivity_dir} \\
        ${args} 2>&1 | tee "\${log_tmp}"

    # Strip ANSI/CR progress noise; drop [+] status and | progress-bar lines
    perl -pe 's/\\r/\\n/g; s/\\e\\[[0-9;]*[A-Za-z]//g' "\${log_tmp}" \\
        | grep -vE '^\\[+\\]|^\\|' \\
        | cat -s \\
        > "\${log_tmp}.clean"
    mv "\${log_tmp}.clean" "\${log_tmp}"

    mkdir -p ${prefix}_rfeval
    mv "\${log_tmp}" ${prefix}_rfeval/rfeval.log

    rnaframework_version=\$(rf-eval -h 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    ${python_version}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def stub_windows = windows
        ? "mkdir -p ${prefix}_rfeval_windows && touch ${prefix}_rfeval_windows/stub_element.xml"
        : ''
    def stub_python = windows ? pythonVersionCmd() : ''
    """
    mkdir -p ${prefix}_rfeval/plots
    ${stub_windows}

    cat <<-END_CSV > ${prefix}_rfeval/${prefix}_eval.csv
    transcript,unpaired_coeff,DSCI,AUROC
    ENST00000000001,0.82,0.79,0.85
    ENST00000000002,0.74,0.71,0.78
    END_CSV

    touch ${prefix}_rfeval/rfeval.log

    rnaframework_version=\$(rf-eval -h 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    ${stub_python}
    """
}
