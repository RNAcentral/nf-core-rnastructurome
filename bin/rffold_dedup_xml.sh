#!/usr/bin/env bash
set -euo pipefail

# Deduplicate XML files staged from multiple replicates into unique_xml/.
# Nextflow stages each file under input<N>/ with the original filename preserved.
# Keeps the first occurrence of each transcript XML when replicates overlap.

mkdir -p unique_xml
declare -A _seen

for _f in input*/*.xml; do
    [[ -f "$_f" ]] || continue
    _base=$(basename "$_f")
    if [[ -z "${_seen[$_base]:-}" ]]; then
        _seen[$_base]=1
        ln -sf "$(readlink -f "$_f")" "unique_xml/$_base"
    fi
done
