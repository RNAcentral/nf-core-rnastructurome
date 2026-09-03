process MERGE_SHANNON_WIG {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(wig)

    output:
    tuple val(meta), path("*.merged.wig"), emit: merged_wig
    tuple val("${task.process}"), val('bash'), eval("bash --version | head -n1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+'"), topic: versions, emit: versions_bash

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    shopt -s nullglob
    wig_files=(*.wig)
    if [[ \${#wig_files[@]} -eq 0 ]]; then
        echo "No Shannon entropy WIG files were provided for merging" >&2
        exit 1
    fi

    {
        echo "track type=wiggle_0"
        for f in "\${wig_files[@]}"; do
            while IFS= read -r line || [[ -n "\$line" ]]; do
                [[ "\$line" == "track "* ]] && continue
                printf '%s\\n' "\$line"
            done < "\$f"
        done
    } > "${prefix}.merged.wig"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.wig
    """
}
