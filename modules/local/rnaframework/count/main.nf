process RNAFRAMEWORK_RFCOUNT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:affb2f7a4bac9a7a' :
        'community.wave.seqera.io/library/rnaframework_findutils:3db7cd7277dc8f08' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta)

    output:
    tuple val(meta), path("*_rfcount/*.rc"), emit: rc
    tuple val(meta), path("*_rfcount/*.rc.rci"), optional: true, emit: rci
    tuple val(meta), path("*_rfcount/index.rci"), optional: true, emit: index_rci
    tuple val(meta), path("*_rfcount/error.out"), optional: true, emit: error_log
    tuple val(meta), path("*_rfcount/samtools.log"), optional: true, emit: samtools_log
    tuple val(meta), path("*_rfcount/*.rfcount_summary.tsv"), optional: true, emit: summary
    tuple val(meta), path("*_rfcount/*.rfcount.log"),         optional: true, emit: log
    tuple val(meta), path("*_rfcount/plots/*.pdf"),           optional: true, emit: plots
    path "versions.yml"          , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rfcount"
    // Mirrors the -m (mutation mode) decision in conf/modules.config, so the summary parser knows
    // which rf-count table layout to expect (MaP has a Mutated-alignments column; RT-stop does not).
    def is_map = ((meta.principle ?: '').toLowerCase() == 'map') ? '1' : '0'
    """
    FASTA_PATH="${fasta}"

    export TERM="\${TERM:-xterm}"

    mkdir -p ${outdir}
    rfcount_outdir="${outdir}"
    rfcount_log_tmp="${prefix}.rfcount.log"

    set -o pipefail
    set +e
    rf-count \\
        -p ${task.cpus} \\
        -f "\${FASTA_PATH}" \\
        -o ${outdir} \\
        -ow \\
        ${args} \\
        "${prefix}:${bam}" 2>&1 | tee "\${rfcount_log_tmp}"
    pipeline_statuses=( "\${PIPESTATUS[@]}" )
    rfcount_status="\${pipeline_statuses[0]}"
    tee_status="\${pipeline_statuses[1]}"
    set -e

    if [[ "\${tee_status}" -ne 0 ]]; then
        echo "[RNAFRAMEWORK_RFCOUNT] tee failed while writing rf-count log for sample '${prefix}'." >&2
        exit "\${tee_status}"
    fi

    cleaned_log="${prefix}.rfcount.clean.log"
    sed -E 's/\\x1b\\[[0-9;]*[A-Za-z]//g' "\${rfcount_log_tmp}" | tr '\\r' '\\n' > "\${cleaned_log}"

    rfcount_completed_with_nonzero=0
    if [[ "\${rfcount_status}" -ne 0 ]]; then
        if grep -Fq '[+] All done.' "\${cleaned_log}"; then
            rfcount_completed_with_nonzero=1
        else
            echo "[RNAFRAMEWORK_RFCOUNT] rf-count failed with exit status \${rfcount_status} for sample '${prefix}'." >&2
            exit "\${rfcount_status}"
        fi
    fi

    summary_tsv="${outdir}/${prefix}.rfcount_summary.tsv"
    rfcount_parse_summary.awk \\
        -v sample="${prefix}" -v is_map="${is_map}" -v match_mode="exact" \\
        "\${cleaned_log}" > "\${summary_tsv}"

    covered=\$(awk -F'\\t' 'NR == 2 {print \$2}' "\${summary_tsv}")
    case "\${covered}" in
        ''|*[!0-9]*)
            ;;
        0)
            echo "[RNAFRAMEWORK_RFCOUNT] rf-count reported zero covered transcripts for sample '${prefix}'." >&2
            echo "[RNAFRAMEWORK_RFCOUNT] No usable signal was available for rf-norm/rf-fold. Use deeper or structure-probing compatible input." >&2
            exit 1
            ;;
    esac

    rc_count=\$(find "\${rfcount_outdir}" -type f -name '*.rc' 2>/dev/null | wc -l)
    if [[ "\${rc_count}" -eq 0 ]]; then
        echo "[RNAFRAMEWORK_RFCOUNT] rf-count produced no RC files for sample '${prefix}'." >&2
        exit 1
    fi
    if [[ "\${rfcount_completed_with_nonzero}" -eq 1 ]]; then
        echo "[RNAFRAMEWORK_RFCOUNT] rf-count exited with status \${rfcount_status} after reporting completion; continuing because RC files were produced." >&2
    fi

    # Publish only the final statistics section (from "[+] Statistics:" to the end); the verbose
    # per-transcript progress above it is noise. Fall back to the full log if the marker is absent,
    # e.g. an early failure, so nothing useful is lost.
    if grep -Fq '[+] Statistics:' "\${cleaned_log}"; then
        awk '/^\\[\\+\\] Statistics:/{p=1} p' "\${cleaned_log}" > "${outdir}/${prefix}.rfcount.log"
    else
        cp "\${cleaned_log}" "${outdir}/${prefix}.rfcount.log"
    fi
    rm -f "\${cleaned_log}" "\${rfcount_log_tmp}"

    rnaframework_version=\$(rf-count -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \${rnaframework_version:-unknown}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rfcount"
    """
    mkdir -p ${outdir}
    touch ${outdir}/${prefix}.rc
    touch ${outdir}/${prefix}.rc.rci
    touch ${outdir}/index.rci
    touch ${outdir}/error.out
    touch ${outdir}/samtools.log
    cat <<-END_SUMMARY > ${outdir}/${prefix}.rfcount_summary.tsv
    sample	covered	pct_mutated	pct_a_muts	pct_c_muts	pct_g_muts	pct_u_muts
    ${prefix}	1	25.0	25.0	25.0	25.0
    END_SUMMARY

    rnaframework_version=\$(rf-count -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \${rnaframework_version:-unknown}
    END_VERSIONS
    """
}
