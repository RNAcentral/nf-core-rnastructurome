process RNAFRAMEWORK_RFRCTOOLS_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(rc, stageAs: "input/*"), path(rci, stageAs: "input/*")
    tuple val(meta_ref), path(gtf)

    output:
    tuple val(meta), path("${prefix}_rctools_extract/${prefix}.rc"),                       emit: rc
    tuple val(meta), path("${prefix}_rctools_extract/${prefix}.rc.rci"), optional: true,   emit: rci
    path "versions.yml",                                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def outdir    = "${prefix}_rctools_extract"
    """
    export TERM="\${TERM:-xterm}"
    mkdir -p ${outdir}

    # Build per-file RCI indexes so rf-rctools can do strand-aware extraction.
    # set -e here so index failures surface rather than silently producing a
    # successful cached task with no .rci output.
    set -e
    for f in input/*.rc; do
        rf-rctools index "\${f}"
    done
    set +e

    # rf-rctools extract requires the BASENAME (no extension) to activate strand-aware
    # extraction: it auto-discovers Sample.plus.rc + Sample.minus.rc and uses the GTF
    # strand field to select the correct file per transcript.  Passing individual *.rc
    # paths disables this logic and produces all-zero counts.
    _rc_first=\$(ls input/*.rc 2>/dev/null | head -1)
    _rc_filename=\$(basename "\${_rc_first}")
    _rc_basename="\${_rc_filename%.plus.rc}"
    _rc_basename="\${_rc_basename%.minus.rc}"
    rc_basename="\${_rc_basename%.rc}"

    rf-rctools extract \\
        -a ${gtf} \\
        -o ${outdir}/${prefix}.rc \\
        -ow \\
        ${args} \\
        input/\${rc_basename}

    set -e
    rf-rctools index ${outdir}/${prefix}.rc
    set +e

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-rctools 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rctools_extract"
    """
    mkdir -p ${outdir}
    touch ${outdir}/${prefix}.rc
    touch ${outdir}/${prefix}.rc.rci

    printf '"%s":\\n    rnaframework: %s\\n' \\
        "${task.process}" \\
        "\$(rf-rctools 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 || echo "unknown")" \\
        > versions.yml
    """
}
