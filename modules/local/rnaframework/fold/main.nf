process RNAFRAMEWORK_RFFOLD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/vicbeg/rnaframework:2.9.6-r2-runtime'

    input:
    tuple val(meta), path(xml)

    output:
    tuple val(meta), path("${prefix}_fold/"),             emit: structures
    tuple val(meta), path("${prefix}_fold/shannon/*.wig"), optional: true, emit: shannon_wig
    tuple val(meta), path("${prefix}_fold/rffold.log"),   optional: true, emit: log
    path "versions.yml"                                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def xml_list = xml instanceof List ? xml.join(' ') : "${xml}"
    def isArm64 = ((System.properties['os.arch'] ?: '').toLowerCase() in ['aarch64', 'arm64'])
    // Only clear conflicting Perl env vars in container mode — in conda mode these vars
    // point to the conda-installed ViennaRNA Perl bindings (RNA.pm) and must be preserved.
    def inContainer = workflow.containerEngine && workflow.containerEngine != 'none'
    def perlEnvCleanup = (isArm64 && inContainer) ? 'unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT' : ''
    """
    export TERM="\${TERM:-xterm}"
    ${perlEnvCleanup}

    log_tmp=\$(mktemp "${prefix}_fold.XXXXXX.log")

    rf-fold \\
        -p ${task.cpus} \\
        -o ${prefix}_fold \\
        -ow \\
        ${args} \\
        ${xml_list} 2>&1 | tee "\${log_tmp}"

    mkdir -p ${prefix}_fold
    mv "\${log_tmp}" ${prefix}_fold/rffold.log

    # rf-fold can return exit 0 even when all folds fail and details are written to error.out.
    # Treat this as a hard failure so the pipeline does not continue with empty fold outputs.
    if [[ -s ${prefix}_fold/error.out ]]; then
        echo "[RNAFRAMEWORK_RFFOLD] rf-fold reported errors:" >&2
        cat ${prefix}_fold/error.out >&2
        exit 1
    fi

    if [[ ! -d ${prefix}_fold/structures ]] || ! find ${prefix}_fold/structures -type f -print -quit | grep -q .; then
        echo "[RNAFRAMEWORK_RFFOLD] No structure files were produced in ${prefix}_fold/structures." >&2
        exit 1
    fi

    missing_list="${prefix}_fold/missing_transcripts.txt"
    warning_log="${prefix}_fold/partial_fold_warning.log"

    expected_list=\$(mktemp)
    folded_list=\$(mktemp)

    printf '%s\n' ${xml_list} | sed 's#.*/##; s#\\.xml\$##' | sort -u >| "\${expected_list}"
    find ${prefix}_fold/structures -maxdepth 1 -type f -name '*.db' -print \\
        | sed 's#.*/##; s#\\.db\$##' | sort -u >| "\${folded_list}"
    comm -23 "\${expected_list}" "\${folded_list}" >| "\${missing_list}"

    expected_count=\$(wc -l < "\${expected_list}" | tr -d ' ')
    folded_count=\$(wc -l < "\${folded_list}" | tr -d ' ')
    missing_count=\$(wc -l < "\${missing_list}" | tr -d ' ')

    if [[ "\${missing_count}" -gt 0 ]]; then
        {
            printf '[RNAFRAMEWORK_RFFOLD] Partial fold output detected.\\n'
            printf '[RNAFRAMEWORK_RFFOLD] Expected %s transcript(s), folded %s, missing %s.\\n' "\${expected_count}" "\${folded_count}" "\${missing_count}"
            printf '[RNAFRAMEWORK_RFFOLD] Missing transcript IDs:\\n'
            cat "\${missing_list}"
        } | tee "\${warning_log}" >&2
    else
        rm -f "\${warning_log}" "\${missing_list}"
    fi

    rm -f "\${expected_list}" "\${folded_list}"

    mv ${prefix}_fold/structures ${prefix}_fold/dotbracket
    if [[ -d ${prefix}_fold/plots/structures ]]; then
        mv ${prefix}_fold/plots/structures ${prefix}_fold/2D-structures
    fi
    if [[ -d ${prefix}_fold/plots/summaries ]]; then
        mv ${prefix}_fold/plots/summaries ${prefix}_fold/summaries
    fi
    rmdir ${prefix}_fold/plots 2>/dev/null || true

    printf '"%s":\n    rnaframework: %s\n' \\
        "${task.process}" \\
        "\$(rf-fold 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_fold/dotbracket
    mkdir -p ${prefix}_fold/2D-structures
    mkdir -p ${prefix}_fold/summaries
    mkdir -p ${prefix}_fold/shannon
    touch ${prefix}_fold/rffold.log
    touch ${prefix}_fold/dotbracket/example.db
    touch ${prefix}_fold/shannon/example.wig

    printf '"%s":\n    rnaframework: %s\n' \\
        "${task.process}" \\
        "\$(rf-fold 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
