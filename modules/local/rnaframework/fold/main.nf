process RNAFRAMEWORK_RFFOLD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/vicbeg/rnaframework:2.9.6-r2-runtime'

    input:
    tuple val(meta), path(xml, stageAs: "input*/*")

    output:
    tuple val(meta), path("${prefix}_fold/"),                                        emit: structures
    tuple val(meta), path("${prefix}_fold_publish/dotbracket/*"),  optional: true,   emit: dotbracket
    tuple val(meta), path("${prefix}_fold_publish/2D-structures/*"), optional: true, emit: structure_plots
    tuple val(meta), path("${prefix}_fold_publish/summaries/*"),   optional: true,   emit: summaries
    tuple val(meta), path("${prefix}_fold_publish/dotplot/*"),     optional: true,   emit: dotplot
    tuple val(meta), path("${prefix}_fold_publish/shannon/*.wig"), optional: true,   emit: shannon_wig
    tuple val(meta), path("${prefix}_fold_publish/rffold.log"),    optional: true,   emit: log
    tuple val(meta), path("${prefix}_fold_publish/missing_transcripts.txt"), optional: true, emit: missing_transcripts
    tuple val(meta), path("${prefix}_fold_publish/partial_fold_warning.log"), optional: true, emit: partial_warning
    path "versions.yml",                                                             emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def isArm64 = ((System.properties['os.arch'] ?: '').toLowerCase() in ['aarch64', 'arm64'])
    // Only clear conflicting Perl env vars in container mode — in conda mode these vars
    // point to the conda-installed ViennaRNA Perl bindings (RNA.pm) and must be preserved.
    def inContainer = workflow.containerEngine && workflow.containerEngine != 'none'
    def perlEnvCleanup = (isArm64 && inContainer) ? 'unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT' : ''
    """
    export TERM="\${TERM:-xterm}"
    ${perlEnvCleanup}

    rffold_dedup_xml.sh

    log_tmp=\$(mktemp "${prefix}_fold.XXXXXX.log")

    rf-fold \\
        -p ${task.cpus} \\
        -o ${prefix}_fold \\
        -ow \\
        ${args} \\
        unique_xml/ 2>&1 | tee "\${log_tmp}"

    mkdir -p ${prefix}_fold
    mv "\${log_tmp}" ${prefix}_fold/rffold.log

    if [[ -s ${prefix}_fold/error.out ]]; then
        if grep -qv "Unable to open RNAplot" ${prefix}_fold/error.out; then
            echo "[RNAFRAMEWORK_RFFOLD] rf-fold reported errors:" >&2
            cat ${prefix}_fold/error.out >&2
            exit 1
        else
            echo "[RNAFRAMEWORK_RFFOLD] rf-fold reported RNAplot warnings (non-fatal):" >&2
            cat ${prefix}_fold/error.out >&2
        fi
    fi

    if [[ ! -d ${prefix}_fold/structures ]] || ! find ${prefix}_fold/structures -type f -print -quit | grep -q .; then
        echo "[RNAFRAMEWORK_RFFOLD] No structure files were produced in ${prefix}_fold/structures." >&2
        exit 1
    fi

    rffold_check_missing.sh ${prefix}_fold unique_xml

    mv ${prefix}_fold/structures ${prefix}_fold/dotbracket
    [[ -d ${prefix}_fold/plots/structures ]] && mv ${prefix}_fold/plots/structures ${prefix}_fold/2D-structures
    [[ -d ${prefix}_fold/plots/summaries  ]] && mv ${prefix}_fold/plots/summaries  ${prefix}_fold/summaries
    rmdir ${prefix}_fold/plots 2>/dev/null || true

    rffold_publish.sh ${prefix}_fold ${prefix}_fold_publish

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-fold 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_fold/dotbracket ${prefix}_fold/2D-structures ${prefix}_fold/summaries ${prefix}_fold/shannon
    touch ${prefix}_fold/rffold.log ${prefix}_fold/dotbracket/example.db ${prefix}_fold/shannon/example.wig

    rffold_publish.sh ${prefix}_fold ${prefix}_fold_publish

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-fold 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
