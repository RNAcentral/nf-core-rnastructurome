#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.error
import urllib.request
import re


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download and merge Ensembl transcript FASTA files for one species."
    )
    parser.add_argument("--species", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--base-url", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-urls", required=True)
    parser.add_argument("--warnings-log", required=True)
    return parser.parse_args()


def fetch_text(url: str) -> str:
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            return response.read().decode("utf-8", errors="ignore")
    except urllib.error.URLError:
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        with urllib.request.urlopen(fallback_url, timeout=60) as response:
            return response.read().decode("utf-8", errors="ignore")


def find_ensembl_file(listing_url: str, pattern: str) -> str:
    listing = fetch_text(listing_url)
    matches = re.findall(r'href="([^"]+)"', listing)
    filtered = [match for match in matches if re.search(pattern, match)]
    if not filtered:
        raise RuntimeError(f"No file matching /{pattern}/ found at {listing_url}")
    return filtered[0]


def find_optional_ensembl_file(listing_url: str, pattern: str) -> str | None:
    try:
        listing = fetch_text(listing_url)
    except urllib.error.URLError:
        return None
    matches = re.findall(r'href="([^"]+)"', listing)
    filtered = [match for match in matches if re.search(pattern, match)]
    return filtered[0] if filtered else None


def release_path_for_value(release: str) -> str:
    if release in ("current", "latest"):
        return "current_fasta"
    if release.startswith("release-"):
        return f"{release}/fasta"
    return f"release-{release}/fasta"


def main() -> int:
    args = parse_args()
    species = args.species.strip().lower().replace(" ", "_")
    release = args.release.strip()
    base_url = args.base_url.rstrip("/")
    warnings_log = args.warnings_log

    species_root = f"{base_url}/{release_path_for_value(release)}/{species}"
    cdna_dir = f"{species_root}/cdna/"
    ncrna_dir = f"{species_root}/ncrna/"

    cdna_name = find_ensembl_file(cdna_dir, r"\.cdna\.all\.fa\.gz")
    ncrna_name = find_optional_ensembl_file(ncrna_dir, r"\.ncrna\.fa\.gz")
    cdna_url = f"{cdna_dir}{cdna_name}"
    ncrna_url = f"{ncrna_dir}{ncrna_name}" if ncrna_name else None

    cdna_local = "cdna.fa.gz"
    ncrna_local = "ncrna.fa.gz"
    urllib.request.urlretrieve(cdna_url, cdna_local)
    if ncrna_url:
        urllib.request.urlretrieve(ncrna_url, ncrna_local)
    else:
        warning = (
            f"[ENSEMBL_TRANSCRIPTOME] Warning: no ncrna FASTA found for species "
            f"'{species}' at {ncrna_dir}. Continuing with cdna only."
        )
        print(warning, file=sys.stderr)
        with open(warnings_log, "w", encoding="utf-8") as warning_handle:
            warning_handle.write(f"{warning}\n")

    with gzip.open(args.output, "wb") as out_handle:
        with gzip.open(cdna_local, "rb") as in_handle:
            shutil.copyfileobj(in_handle, out_handle)
        if ncrna_url:
            with gzip.open(ncrna_local, "rb") as in_handle:
                shutil.copyfileobj(in_handle, out_handle)

    with open(args.source_urls, "w", encoding="utf-8") as handle:
        handle.write(f"{cdna_url}\n")
        if ncrna_url:
            handle.write(f"{ncrna_url}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
