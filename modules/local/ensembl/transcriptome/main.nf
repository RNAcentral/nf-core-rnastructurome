process ENSEMBL_TRANSCRIPTOME {
    tag "${meta.id}:${ensembl_species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12-slim'

    input:
    tuple val(meta), val(ensembl_species)

    output:
    tuple val(meta), path("${meta.id}.transcripts.fa.gz"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def release = (params.ensembl_release ?: 'current').toString()
    def baseUrl = (params.ensembl_base_url ?: 'https://ftp.ensembl.org/pub').toString()
    """
    python - <<'PY'
    import gzip
    import os
    import re
    import shutil
    import urllib.request

    species = "${ensembl_species}".strip().lower()
    release = "${release}".strip()
    base_url = "${baseUrl}".rstrip("/")
    out_gz = "${meta.id}.transcripts.fa.gz"

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
        with urllib.request.urlopen(url, timeout=60) as response:
            return response.read().decode("utf-8", errors="ignore")

    def find_ensembl_file(listing_url: str, pattern: str) -> str:
        listing = fetch_text(listing_url)
        matches = re.findall(r'href="([^"]+)"', listing)
        filtered = [m for m in matches if re.search(pattern, m)]
        if not filtered:
            raise RuntimeError(f"No file matching /{pattern}/ found at {listing_url}")
        return filtered[0]

    cdna_name = find_ensembl_file(cdna_dir, r'\\.cdna\\.all\\.fa\\.gz')
    ncrna_name = find_ensembl_file(ncrna_dir, r'\\.ncrna\\.fa\\.gz')
    cdna_url = f"{cdna_dir}{cdna_name}"
    ncrna_url = f"{ncrna_dir}{ncrna_name}"

    cdna_local = "cdna.fa.gz"
    ncrna_local = "ncrna.fa.gz"
    urllib.request.urlretrieve(cdna_url, cdna_local)
    urllib.request.urlretrieve(ncrna_url, ncrna_local)

    with gzip.open(out_gz, "wb") as out_handle:
        with gzip.open(cdna_local, "rb") as in_handle:
            shutil.copyfileobj(in_handle, out_handle)
        with gzip.open(ncrna_local, "rb") as in_handle:
            shutil.copyfileobj(in_handle, out_handle)

    with open("ensembl_source_urls.txt", "w", encoding="utf-8") as handle:
        handle.write(f"{cdna_url}\\n{ncrna_url}\\n")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl_release: "${release}"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.transcripts.fa.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl_release: "${release}"
    END_VERSIONS
    """
}
