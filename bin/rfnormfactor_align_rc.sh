#!/usr/bin/env bash
# rfnormfactor_align_rc.sh — align a set of RNA Framework RC files to their common transcript set.
#
# rf-normfactor requires every input RC file to contain the SAME transcript set (it errors with
# "Provided RC files have unequal sizes" otherwise). On the genome route, the per-sample coverage
# pre-filter in rf-rctools extract leaves treated/untreated covering different transcripts. This
# script aligns all inputs to their common transcript set (intersection) by re-extracting with
# rf-rctools, writing aligned copies into <aligned_dir>/<basename>. Transcript order is identical
# across outputs (one shared BED), so the caller's positional -t/-u/-d pairing is preserved.
#
# Usage: rfnormfactor_align_rc.sh <label> <aligned_dir> <rc_file>...
#   <label>        identifier used only in the "no common transcripts" error message
#   <aligned_dir>  directory to write the aligned RC files into (created if missing)
#   <rc_file>...   the RC files to align (treated, untreated, denatured — order-preserving)
#
# Exits 1 if the RC files share no transcripts (cross-experiment normalisation impossible).
set -euo pipefail

label="$1"; shift
aligned_dir="$1"; shift

# Per-transcript id list, and id<TAB>length, parsed from rf-rctools view's 4-line-per-transcript output.
idlist() { rf-rctools view "$1" | awk 'NF==0{l=0;next}{l++} l==1{print $0} l==4{l=0}' | sort -u; }
idlen()  { rf-rctools view "$1" | awk 'NF==0{l=0;next}{l++} l==1{id=$0} l==2{len=length($0)} l==4{print id"\t"len; l=0}' | sort -u; }

first=1
for rc in "$@"; do
    rf-rctools index "$rc" >/dev/null 2>&1 || true
    if [[ ${first} -eq 1 ]]; then
        idlist "$rc" > common_ids.txt
        idlen  "$rc" > lengths.txt
        first=0
    else
        comm -12 common_ids.txt <(idlist "$rc") > common_ids.new
        mv common_ids.new common_ids.txt
    fi
done

if [[ ! -s common_ids.txt ]]; then
    echo "[RNAFRAMEWORK_RFNORMFACTOR] no transcripts common to all RC files for '${label}'; cannot compute cross-experiment normalisation factors." >&2
    exit 1
fi

# BED of the common transcripts (full length), used to re-extract every RC to the same transcript set.
awk 'NR==FNR{keep[$1]=1; next} ($1 in keep){print $1"\t0\t"$2"\t"$1}' common_ids.txt lengths.txt | sort -k1,1 > common.bed

mkdir -p "${aligned_dir}"
for rc in "$@"; do
    rf-rctools extract -a common.bed -o "${aligned_dir}/$(basename "${rc}")" -ow "${rc}"
    rf-rctools index "${aligned_dir}/$(basename "${rc}")"
done
