process RNAFRAMEWORK_RFFOLD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(xml, stageAs: "input*/*")

    output:
    tuple val(meta), path("${prefix}_fold/"),                                        emit: structures
    tuple val(meta), path("${prefix}_fold_publish/dotbracket/*"),  optional: true,   emit: dotbracket
    tuple val(meta), path("${prefix}_fold_publish/structures/*"), optional: true, emit: structure_plots
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
    def sl_m          = (args =~ /-sl\s+(\S+)/)
    def in_m          = (args =~ /-in\s+(\S+)/)
    def slope_log     = sl_m ? sl_m[0][1] : 'not set'
    def intercept_log = in_m ? in_m[0][1] : 'not set'
    """
    export TERM="\${TERM:-xterm}"
    ${perlEnvCleanup}

    # Rebuild one experiment directory per replicate so rf-fold folds replicates together via
    # majority voting (rf-fold experiment1/ experiment2/ ...). The XMLs were staged one-per-dir as
    # input*/ in replicate-grouped order; fold_replicate_sizes says how many consecutive input dirs
    # belong to each replicate. Without this, a single merged dir makes rf-fold see only 1 sample
    # and the -oc (only-common) consensus parameter is rejected.
    sizes=(${meta.fold_replicate_sizes.join(' ')})
    mapfile -t _input_dirs < <(ls -d input*/ 2>/dev/null | sort -V)
    if [[ \${#_input_dirs[@]} -ne ${meta.fold_replicate_sizes.sum()} ]]; then
        echo "[RNAFRAMEWORK_RFFOLD] staged input dir count (\${#_input_dirs[@]}) != expected (${meta.fold_replicate_sizes.sum()})." >&2
        exit 1
    fi
    idx=0
    for k in "\${!sizes[@]}"; do
        expdir="experiment\$((k + 1))"
        mkdir -p "\${expdir}"
        for ((j = 0; j < sizes[k]; j++)); do
            for x in "\${_input_dirs[idx]}"*.xml; do
                [[ -e "\${x}" ]] && ln -sf "\$(readlink -f "\${x}")" "\${expdir}/\$(basename "\${x}")"
            done
            idx=\$((idx + 1))
        done
    done

    # Expected transcripts for the missing-check = those rf-fold will attempt = present in ALL
    # replicate experiments (matches default -oc = replicate count). Building the intersection here
    # stops the check from flagging transcripts that -oc intentionally skips (present in only some
    # replicates). For a single-replicate group this is just that experiment's transcript set.
    mkdir -p expected_xml
    _intersect=\$(mktemp)
    _this=\$(mktemp)
    _first=1
    for d in experiment*/; do
        ls "\${d}"*.xml 2>/dev/null | sed 's#.*/##' | sort -u >| "\${_this}"
        if [[ \${_first} -eq 1 ]]; then
            cp "\${_this}" "\${_intersect}"
            _first=0
        else
            comm -12 "\${_intersect}" "\${_this}" >| "\${_intersect}.new" && mv "\${_intersect}.new" "\${_intersect}"
        fi
    done
    while read -r _b; do
        [[ -n "\${_b}" ]] && ln -sf "\$(readlink -f "experiment1/\${_b}")" "expected_xml/\${_b}"
    done < "\${_intersect}"
    rm -f "\${_intersect}" "\${_this}"

    log_tmp=\$(mktemp "${prefix}_fold.XXXXXX.log")
    printf 'slope=%s intercept=%s\n' "${slope_log}" "${intercept_log}" >> "\${log_tmp}"

    rf-fold \\
        -p ${task.cpus} \\
        -o ${prefix}_fold \\
        -ow \\
        ${args} \\
        experiment*/ 2>&1 | tee -a "\${log_tmp}"

    mkdir -p ${prefix}_fold
    mv "\${log_tmp}" ${prefix}_fold/rffold.log

    if [[ -s ${prefix}_fold/error.out ]]; then
        exception_count=\$(grep -Fc "[!] Exception" ${prefix}_fold/error.out || true)
        rnaplot_count=\$(grep -Fc "Unable to open RNAplot" ${prefix}_fold/error.out || true)
        if [[ \$exception_count -gt 0 && \$exception_count -ne \$rnaplot_count ]]; then
            echo "[RNAFRAMEWORK_RFFOLD] rf-fold reported errors:" >&2
            cat ${prefix}_fold/error.out >&2
            exit 1
        elif [[ \$exception_count -gt 0 ]]; then
            echo "[RNAFRAMEWORK_RFFOLD] rf-fold reported RNAplot warnings (non-fatal):" >&2
            cat ${prefix}_fold/error.out >&2
        fi
    fi

    if [[ ! -d ${prefix}_fold/structures ]] || ! find ${prefix}_fold/structures -type f -print -quit | grep -q .; then
        echo "[RNAFRAMEWORK_RFFOLD] No structure files were produced in ${prefix}_fold/structures." >&2
        exit 1
    fi

    rffold_check_missing.sh ${prefix}_fold expected_xml

    mv ${prefix}_fold/structures ${prefix}_fold/dotbracket
    [[ -d ${prefix}_fold/plots/structures ]] && mv ${prefix}_fold/plots/structures ${prefix}_fold/structures
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
    mkdir -p ${prefix}_fold/dotbracket ${prefix}_fold/structures ${prefix}_fold/summaries ${prefix}_fold/shannon
    touch ${prefix}_fold/rffold.log ${prefix}_fold/dotbracket/example.db ${prefix}_fold/shannon/example.wig

    rffold_publish.sh ${prefix}_fold ${prefix}_fold_publish

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-fold 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
