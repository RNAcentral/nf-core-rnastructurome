process RNAFRAMEWORK_RFCOUNT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

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
    tuple val(meta), path("*_rfcount/plots/*.pdf"), optional: true, emit: plots
    path "versions.yml"          , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rfcount"
    """
    FASTA_PATH="${fasta}"

    export TERM="\${TERM:-xterm}"

    mkdir -p ${outdir}
    rfcount_log_tmp="${prefix}.rfcount.log"

    set -o pipefail
    rf-count \\
        -p ${task.cpus} \\
        -f "\${FASTA_PATH}" \\
        -o ${outdir} \\
        -ow \\
        ${args} \\
        "${prefix}:${bam}" 2>&1 | tee "\${rfcount_log_tmp}"

    cleaned_log="${prefix}.rfcount.clean.log"
    sed -E 's/\\x1b\\[[0-9;]*[A-Za-z]//g' "\${rfcount_log_tmp}" | tr '\\r' '\\n' > "\${cleaned_log}"

    summary_tsv="${outdir}/${prefix}.rfcount_summary.tsv"
    {
        printf 'sample\\tcovered\\tpct_a_stops\\tpct_c_stops\\tpct_g_stops\\tpct_u_stops\\n'
        awk -v sample="${prefix}" '
            \$1 == sample {
                if (\$3 ~ /\\// && \$4 ~ /^\(/) {
                    pct_a = \$4; gsub(/[()%]/, "", pct_a)
                    print \$1 "\\t" \$2 "\\t" pct_a "\\t" \$5 "\\t" \$6 "\\t" \$7
                } else {
                    print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6
                }
            }
        ' "\${cleaned_log}" | tail -n 1
    } > "\${summary_tsv}"

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

    rm -f "\${rfcount_log_tmp}" "\${cleaned_log}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-count 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1 || echo "unknown")
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
    sample	covered	pct_a_stops	pct_c_stops	pct_g_stops	pct_u_stops
    ${prefix}	1	25.0	25.0	25.0	25.0
    END_SUMMARY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-count 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
