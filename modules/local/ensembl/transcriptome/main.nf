process ENSEMBL_TRANSCRIPTOME {
    tag "${meta.id}:${ensembl_species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12.11'

    input:
    tuple val(meta), val(ensembl_species)

    output:
    tuple val(meta), path("${meta.id}.transcripts.fa.gz"), emit: fasta
    path "ensembl_source_urls.txt", emit: source_urls
    path "ensembl_warnings.log", optional: true, emit: warnings
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python - <<'PY'
import gzip
import os
import re
import shutil
import sys
import urllib.request

species = "${ensembl_species}".strip().lower().replace(" ", "_")
release = "${(params.ensembl_release ?: 'current').toString()}".strip()
base_url = "${(params.ensembl_base_url ?: 'https://ftp.ensembl.org/pub').toString()}".rstrip("/")
out_gz = "${meta.id}.transcripts.fa.gz"
warnings_log = "ensembl_warnings.log"

if release in ("current", "latest"):
    release_path = "current_fasta"
elif release.startswith("release-"):
    release_path = f"{release}/fasta"
else:
    release_path = f"release-{release}/fasta"

species_root = f"{base_url}/{release_path}/{species}"
cdna_dir = f"{species_root}/cdna/"
ncrna_dir = f"{species_root}/ncrna/"

def fetch_text(url: str) -> str:
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            return response.read().decode("utf-8", errors="ignore")
    except Exception:
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        with urllib.request.urlopen(fallback_url, timeout=60) as response:
            return response.read().decode("utf-8", errors="ignore")

def find_ensembl_file(listing_url: str, pattern: str) -> str:
    listing = fetch_text(listing_url)
    matches = re.findall(r'href="([^"]+)"', listing)
    filtered = [m for m in matches if re.search(pattern, m)]
    if not filtered:
        raise RuntimeError(f"No file matching /{pattern}/ found at {listing_url}")
    return filtered[0]

def find_optional_ensembl_file(listing_url: str, pattern: str):
    try:
        listing = fetch_text(listing_url)
    except Exception:
        return None
    matches = re.findall(r'href="([^"]+)"', listing)
    filtered = [m for m in matches if re.search(pattern, m)]
    return filtered[0] if filtered else None

cdna_name = find_ensembl_file(cdna_dir, r'\\.cdna\\.all\\.fa\\.gz')
ncrna_name = find_optional_ensembl_file(ncrna_dir, r'\\.ncrna\\.fa\\.gz')
cdna_url = f"{cdna_dir}{cdna_name}"
ncrna_url = f"{ncrna_dir}{ncrna_name}" if ncrna_name else None

cdna_local = "cdna.fa.gz"
ncrna_local = "ncrna.fa.gz"
urllib.request.urlretrieve(cdna_url, cdna_local)
if ncrna_url:
    urllib.request.urlretrieve(ncrna_url, ncrna_local)
else:
    warning = f"[ENSEMBL_TRANSCRIPTOME] Warning: no ncrna FASTA found for species '{species}' at {ncrna_dir}. Continuing with cdna only."
    print(warning, file=sys.stderr)
    with open(warnings_log, "w", encoding="utf-8") as warning_handle:
        warning_handle.write(f"{warning}\\n")

with gzip.open(out_gz, "wb") as out_handle:
    with gzip.open(cdna_local, "rb") as in_handle:
        shutil.copyfileobj(in_handle, out_handle)
    if ncrna_url:
        with gzip.open(ncrna_local, "rb") as in_handle:
            shutil.copyfileobj(in_handle, out_handle)

with open("ensembl_source_urls.txt", "w", encoding="utf-8") as handle:
    handle.write(f"{cdna_url}\\n")
    if ncrna_url:
        handle.write(f"{ncrna_url}\\n")
PY

    printf '%s\n' \
        '"${task.process}":' \
        '    ensembl_release: "${(params.ensembl_release ?: 'current').toString()}"' \
        > versions.yml
    """

    stub:
    """
    touch ${meta.id}.transcripts.fa.gz
    printf '%s\n' \
        "stub://${meta.id}.transcripts.fa.gz" \
        > ensembl_source_urls.txt
    printf '%s\n' \
        '"${task.process}":' \
        '    ensembl_release: "${(params.ensembl_release ?: 'current').toString()}"' \
        > versions.yml
    """
}
