process RNAFRAMEWORK_RFNORM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(treated), path(untreated), path(denatured), path(rci_files), path(norm_factor)

    output:
    tuple val(meta), path("${prefix}_norm/*.xml"), emit: xml
    tuple val(meta), path("${prefix}_norm/rfnorm.log"), optional: true, emit: log
    tuple val(meta), path("${prefix}_norm/plots/*.pdf"), optional: true, emit: plots
    path "versions.yml"                           , emit: versions

    script:
    def args          = task.ext.args ?: ''
    prefix            = task.ext.prefix ?: "${meta.id}"
    def untreated_arg = untreated ? "-u ${untreated}" : ''
    def denatured_arg = denatured ? "-d ${denatured}" : ''
    // Cross-experiment normalisation factor from rf-normfactor (opt-in via --rfnorm_use_normfactor).
    // IMPORTANT: rf-norm's -nf takes a NUMERIC factor value, NOT a path. rf-normfactor writes a
    // per-experiment table (col 1 = treated RC basename, col 2 = factor); we look up THIS group's
    // factor by its treated RC name below and pass the value. Only attempted for single-treated
    // groups (cross-experiment normalisation emits one treated per group); a missing entry or a
    // multi-treated group falls back to rf-norm's internal per-sample normalisation.
    def treatedNames  = (treated instanceof List ? treated : [treated]).collect { f -> f.name }
    def nfKeyName     = (norm_factor && treatedNames.size() == 1) ? treatedNames[0] : ''
    def treated_list  = treated instanceof List ? treated.join(' ') : "${treated}"
    """
    export TERM="\${TERM:-xterm}"
    rfnorm_log_tmp="${prefix}.rfnorm.log"

    # Resolve the numeric -nf value for this group from the rf-normfactor table (see note above).
    norm_factor_arg=""
    nf_key_name="${nfKeyName}"
    if [[ -n "\${nf_key_name}" && -s "${norm_factor}" ]]; then
        nf_key="\${nf_key_name%.rc}"
        nf_value="\$(awk -v k="\${nf_key}" '\$1 == k { print \$2; exit }' "${norm_factor}")"
        if [[ -n "\${nf_value}" ]]; then
            norm_factor_arg="-nf \${nf_value}"
        else
            echo "[RNAFRAMEWORK_RFNORM] no rf-normfactor entry for '\${nf_key}' in ${norm_factor}; using rf-norm internal per-sample normalisation" >&2
        fi
    fi

    rf-norm \\
        -p ${task.cpus} \\
        -o ${prefix}_norm \\
        -ow \\
        ${args} \\
        \${norm_factor_arg} \\
        -t ${treated_list} \\
        ${untreated_arg} \\
        ${denatured_arg} 2>&1 | tee "\${rfnorm_log_tmp}"

    if [[ -d ${prefix}_norm ]]; then
        mv "\${rfnorm_log_tmp}" ${prefix}_norm/rfnorm.log
    fi

    if [[ -d "${prefix}_norm/plots" ]]; then
        find "${prefix}_norm/plots" -maxdepth 1 -type f -name '*.pdf' | while IFS= read -r plot; do
            transcript_id="\$(basename "\${plot}" .pdf)"
            if [[ ! -f "${prefix}_norm/\${transcript_id}.xml" ]]; then
                rm -f "\${plot}"
            fi
        done
    fi

    rnaframework_version=\$(rf-norm -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    printf '"%s":\n    rnaframework: %s\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_norm
    mkdir -p ${prefix}_norm/plots
    touch ${prefix}_norm/${prefix}.xml
    touch ${prefix}_norm/plots/${prefix}.pdf
    touch ${prefix}_norm/rfnorm.log

    rnaframework_version=\$(rf-norm -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    printf '"%s":\n    rnaframework: %s\n' "${task.process}" "\${rnaframework_version:-unknown}" > versions.yml
    """
}
