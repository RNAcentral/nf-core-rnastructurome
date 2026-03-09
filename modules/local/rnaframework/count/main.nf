process RNAFRAMEWORK_RFCOUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'rnastructurome/rnaframework:2.9.6-r1'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta)

    output:
    tuple val(meta), path("*_rfcount/*.rc"), emit: rc
    tuple val(meta), path("*_rfcount/*.rc.rci"), optional: true, emit: rci
    tuple val(meta), path("*_rfcount/index.rci"), optional: true, emit: index_rci
    tuple val(meta), path("*_rfcount/error.out"), optional: true, emit: error_log
    tuple val(meta), path("*_rfcount/samtools.log"), optional: true, emit: samtools_log
    tuple val(meta), path("*_rfcount/plots/*.pdf"), optional: true, emit: plots
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fallback_fasta = params.fasta ?: (
        meta_ref?.organism && params.genomes?.containsKey(meta_ref.organism)
            ? (params.genomes[meta_ref.organism]?.transcript_fasta ?: params.genomes[meta_ref.organism]?.transcriptome ?: params.genomes[meta_ref.organism]?.cdna ?: '')
            : ''
    )
    def outdir = "${prefix}_rfcount"
    """
    FASTA_PATH="${fasta}"
    if [[ ! -f "\${FASTA_PATH}" && -n "${fallback_fasta}" ]]; then
        FASTA_PATH="${fallback_fasta}"
    fi

    export TERM="\${TERM:-xterm}"

    set -o pipefail
    rf-count \\
        -p ${task.cpus} \\
        -f "\${FASTA_PATH}" \\
        -o ${outdir} \\
        -ow \\
        ${args} \\
        "${prefix}:${bam}" 2>&1 | tee rfcount.log

    covered=\$(sed -E 's/\\x1b\\[[0-9;]*[A-Za-z]//g' rfcount.log | awk -v sample="${prefix}" '\$1 == sample {print \$2}' | tail -n 1)
    case "\${covered}" in
        ''|*[!0-9]*)
            ;;
        0)
            echo "[RNAFRAMEWORK_RFCOUNT] rf-count reported zero covered transcripts for sample '${prefix}'." >&2
            echo "[RNAFRAMEWORK_RFCOUNT] No usable signal was available for rf-norm/rf-fold. Use deeper or structure-probing compatible input." >&2
            exit 1
            ;;
    esac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-count 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaframework: \$(rf-count 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
