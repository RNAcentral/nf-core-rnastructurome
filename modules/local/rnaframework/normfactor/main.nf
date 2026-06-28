process RNAFRAMEWORK_RFNORMFACTOR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(treated), path(untreated), path(denatured), path(rci_files)

    output:
    tuple val(meta), path("${prefix}.norm_factors.txt"), optional: true, emit: factors
    tuple val(meta), path("${prefix}.rfnormfactor.log"), optional: true, emit: log
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    prefix             = task.ext.prefix ?: "${meta.id}"
    def treatedNames   = (treated instanceof List ? treated : [treated]).collect { it.name }
    def untreatedNames = untreated ? (untreated instanceof List ? untreated : [untreated]).collect { it.name } : []
    def denaturedNames = denatured ? (denatured instanceof List ? denatured : [denatured]).collect { it.name } : []
    """
    export TERM="\${TERM:-xterm}"
    rfnormfactor_log_tmp="${prefix}.rfnormfactor.log"

    treated_files=(${treatedNames.collect { "\"${it}\"" }.join(' ')})
    untreated_files=(${untreatedNames ? untreatedNames.collect { "\"${it}\"" }.join(' ') : ''})
    denatured_files=(${denaturedNames ? denaturedNames.collect { "\"${it}\"" }.join(' ') : ''})

    # rf-normfactor requires every input RC file to contain the SAME transcript set (it errors with
    # "Provided RC files have unequal sizes" otherwise). On the genome route, the per-sample coverage
    # pre-filter in rf-rctools extract leaves treated/untreated covering different transcripts. Align
    # all inputs to their common transcript set by re-extracting them (see bin/rfnormfactor_align_rc.sh);
    # transcript order is identical across outputs, so the positional -t/-u/-d pairing is preserved.
    all_files=( "\${treated_files[@]}" )
    [[ \${#untreated_files[@]} -gt 0 ]] && all_files+=( "\${untreated_files[@]}" )
    [[ \${#denatured_files[@]} -gt 0 ]] && all_files+=( "\${denatured_files[@]}" )

    rfnormfactor_align_rc.sh "${prefix}" aligned "\${all_files[@]}"

    aligned_list() { local out=""; for f in "\$@"; do out+="aligned/\${f},"; done; echo "\${out%,}"; }
    treated_arg="-t \$(aligned_list "\${treated_files[@]}")"
    untreated_arg=""
    [[ \${#untreated_files[@]} -gt 0 ]] && untreated_arg="-u \$(aligned_list "\${untreated_files[@]}")"
    denatured_arg=""
    [[ \${#denatured_files[@]} -gt 0 ]] && denatured_arg="-d \$(aligned_list "\${denatured_files[@]}")"

    # Best-effort: rf-normfactor can legitimately fail to produce factors (e.g. "No bases covered
    # across all samples" when no base meets -mc across every sample). Cross-experiment normalisation
    # is auto-enabled, so a shortfall must NOT break the run — capture the status, warn, and emit no
    # factor file. The subworkflow then leaves rf-norm to normalise this reference per-sample.
    set +e
    rf-normfactor \\
        -p ${task.cpus} \\
        -o ${prefix}.norm_factors.txt \\
        -ow \\
        ${args} \\
        \${treated_arg} \\
        \${untreated_arg} \\
        \${denatured_arg} 2>&1 | tee "\${rfnormfactor_log_tmp}"
    nf_status="\${PIPESTATUS[0]}"
    set -e

    mv "\${rfnormfactor_log_tmp}" ${prefix}.rfnormfactor.log || true

    if [[ "\${nf_status}" -ne 0 || ! -s "${prefix}.norm_factors.txt" ]]; then
        echo "[RNAFRAMEWORK_RFNORMFACTOR] rf-normfactor did not produce normalisation factors for '${prefix}' (exit \${nf_status}); rf-norm will fall back to per-sample normalisation for this reference. If coverage is the cause, lower --rfnorm_normfactor_min_coverage." >&2
        rm -f "${prefix}.norm_factors.txt"
    fi

    # rf-normfactor prints its version banner only under -h and exits non-zero; capture with '|| true'
    # so 'set -o pipefail' does not abort, and fall back to 'unknown' if no version is found.
    rnaframework_version=\$(rf-normfactor -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.norm_factors.txt
    touch ${prefix}.rfnormfactor.log

    rnaframework_version=\$(rf-normfactor -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    printf '"%s":\\n    rnaframework: %s\\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """
}
