process RNAFRAMEWORK_RFCOUNT_GENOME {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta)

    output:
    tuple val(meta), path("*_rfcount_genome/*.rc"),             optional: true, emit: rc
    tuple val(meta), path("*_rfcount_genome/*.rc.rci"),         optional: true, emit: rci
    tuple val(meta), path("*_rfcount_genome/index.rci"),        optional: true, emit: index_rci
    tuple val(meta), path("*_rfcount_genome/error.out"),        optional: true, emit: error_log
    tuple val(meta), path("*_rfcount_genome/samtools.log"),     optional: true, emit: samtools_log
    tuple val(meta), path("*_rfcount_genome/*.rfcount_genome_summary.tsv"), optional: true, emit: summary
    tuple val(meta), path("*_rfcount_genome/*.rfcount_genome.log"),         optional: true, emit: log
    tuple val(meta), path("*_rfcount_genome/plots/*.pdf"),                  optional: true, emit: plots
    path "versions.yml",                                                        emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rfcount_genome"
    """
    FASTA_PATH="${fasta}"

    export TERM="\${TERM:-xterm}"

    mkdir -p ${outdir}
    rfcount_outdir="${outdir}"
    rfcount_log_tmp="${prefix}.rfcount_genome.log"

    set -o pipefail
    set +e
    rf-count-genome \\
        -p ${task.cpus} \\
        -f "\${FASTA_PATH}" \\
        -o ${outdir} \\
        -ow \\
        ${args} \\
        "${bam}" 2>&1 | tee "\${rfcount_log_tmp}"
    pipeline_statuses=( "\${PIPESTATUS[@]}" )
    rfcount_status="\${pipeline_statuses[0]}"
    tee_status="\${pipeline_statuses[1]}"
    set -e

    if [[ "\${tee_status}" -ne 0 ]]; then
        echo "[RNAFRAMEWORK_RFCOUNT_GENOME] tee failed while writing rf-count-genome log for sample '${prefix}'." >&2
        exit "\${tee_status}"
    fi

    cleaned_log="${prefix}.rfcount_genome.clean.log"
    sed -E 's/\\x1b\\[[0-9;]*[A-Za-z]//g' "\${rfcount_log_tmp}" | tr '\\r' '\\n' > "\${cleaned_log}"

    rfcount_completed_with_nonzero=0
    if [[ "\${rfcount_status}" -ne 0 ]]; then
        if grep -Fq '[+] All done.' "\${cleaned_log}"; then
            rfcount_completed_with_nonzero=1
        else
            echo "[RNAFRAMEWORK_RFCOUNT_GENOME] rf-count-genome failed with exit status \${rfcount_status} for sample '${prefix}'." >&2
            exit "\${rfcount_status}"
        fi
    fi

    summary_tsv="${outdir}/${prefix}.rfcount_genome_summary.tsv"
    {
        printf 'sample\\tcovered\\tmutated_alignments\\tpct_mutated\\tpct_a_muts\\tpct_c_muts\\tpct_g_muts\\tpct_u_muts\\n'
        # Sample column carries the staged BAM's filename suffix (e.g. "${prefix}.sorted"), so match by prefix.
        awk -v sample="${prefix}" '
            index(\$1, sample) == 1 {
                if (index(\$3, "/") > 0 && substr(\$4, 1, 1) == "(") {
                    pct_mut = \$4; gsub("[()%]", "", pct_mut)
                    print sample "\\t" \$2 "\\t" \$3 "\\t" pct_mut "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8
                } else {
                    print sample "\\t" \$2 "\\t" "" "\\t" "" "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6
                }
            }
        ' "\${cleaned_log}" | tail -n 1
    } > "\${summary_tsv}"

    # rf-count-genome reports "Covered: 0" when run without a -a annotation file
    # (genome-wide mode); that is expected here — rf-rctools extract handles transcript
    # extraction using the GTF in the next step.  Fail only if no RC files were produced.
    rc_count=\$(find "\${rfcount_outdir}" -type f -name '*.rc' 2>/dev/null | wc -l)
    if [[ "\${rc_count}" -eq 0 ]]; then
        echo "[RNAFRAMEWORK_RFCOUNT_GENOME] rf-count-genome produced no RC files for sample '${prefix}'." >&2
        exit 1
    fi
    if [[ "\${rfcount_completed_with_nonzero}" -eq 1 ]]; then
        echo "[RNAFRAMEWORK_RFCOUNT_GENOME] rf-count-genome exited with status \${rfcount_status} after reporting completion; continuing because RC files were produced." >&2
    fi

    mv "\${cleaned_log}" "${outdir}/${prefix}.rfcount_genome.log"
    rm -f "\${rfcount_log_tmp}"

    rnaframework_version=\$(rf-count-genome -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \${rnaframework_version:-unknown}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rfcount_genome"
    """
    mkdir -p ${outdir}
    touch ${outdir}/${prefix}.plus.rc
    touch ${outdir}/${prefix}.plus.rc.rci
    touch ${outdir}/index.rci
    touch ${outdir}/error.out
    touch ${outdir}/samtools.log
    cat <<-END_SUMMARY > ${outdir}/${prefix}.rfcount_genome_summary.tsv
    sample	covered	pct_mutated	pct_a_muts	pct_c_muts	pct_g_muts	pct_u_muts
    ${prefix}	1	25.0	25.0	25.0	25.0
    END_SUMMARY

    rnaframework_version=\$(rf-count-genome -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1) || true
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \${rnaframework_version:-unknown}
    END_VERSIONS
    """
}
