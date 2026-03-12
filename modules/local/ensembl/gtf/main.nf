process ENSEMBL_GTF {
    tag "${meta.id}:${ensembl_species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12.11'

    input:
    tuple val(meta), val(ensembl_species)

    output:
    tuple val(meta), path("${meta.id}.annotation.gtf.gz"), emit: gtf
    path "ensembl_source_urls.txt", emit: source_urls
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python - <<'PY'
import re
import urllib.request

species = "${ensembl_species}".strip().lower().replace(" ", "_")
release = "${(params.ensembl_release ?: 'current').toString()}".strip()
base_url = "${(params.ensembl_base_url ?: 'https://ftp.ensembl.org/pub').toString()}".rstrip("/")
out_gz = "${meta.id}.annotation.gtf.gz"

if release in ("current", "latest"):
    species_root = f"{base_url}/current_gtf/{species}/"
elif release.startswith("release-"):
    species_root = f"{base_url}/{release}/gtf/{species}/"
else:
    species_root = f"{base_url}/release-{release}/gtf/{species}/"

def fetch_text(url: str) -> str:
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            return response.read().decode("utf-8", errors="ignore")
    except Exception:
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        with urllib.request.urlopen(fallback_url, timeout=60) as response:
            return response.read().decode("utf-8", errors="ignore")

listing = fetch_text(species_root)
matches = [match for match in re.findall(r'href="([^"]+)"', listing) if re.search(r'\\.gtf\\.gz\$', match)]
preferred_matches = [match for match in matches if 'abinitio' not in match.lower()]
gtf_name = (preferred_matches or matches or [None])[0]
if not gtf_name:
    raise RuntimeError(f"No .gtf.gz file found at {species_root}")

gtf_url = f"{species_root}{gtf_name}"
urllib.request.urlretrieve(gtf_url, out_gz)

with open("ensembl_source_urls.txt", "w", encoding="utf-8") as handle:
    handle.write(f"{gtf_url}\\n")
PY

    printf '%s\n' \
        '"${task.process}":' \
        '    ensembl_release: "${(params.ensembl_release ?: 'current').toString()}"' \
        > versions.yml
    """

    stub:
    """
    touch ${meta.id}.annotation.gtf.gz
    printf '%s\n' \
        "stub://${meta.id}.annotation.gtf.gz" \
        > ensembl_source_urls.txt
    printf '%s\n' \
        '"${task.process}":' \
        '    ensembl_release: "${params.ensembl_release ?: 'current'}"' \
        > versions.yml
    """
}
