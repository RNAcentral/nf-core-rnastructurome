#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
import urllib.error
import urllib.request


class EnsemblSpeciesNotFound(Exception):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download one Ensembl GTF annotation file for one species."
    )
    parser.add_argument("--species", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--base-url", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-urls", required=True)
    parser.add_argument(
        "--not-found-file",
        required=True,
        help="Path to write (empty) when the species is absent from Ensembl FTP; "
             "process exits 0 and GTF is silently skipped for this reference.",
    )
    return parser.parse_args()


def fetch_text(url: str) -> str:
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            return response.read().decode("utf-8", errors="ignore")
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise EnsemblSpeciesNotFound(f"HTTP 404 at {url}")
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        try:
            with urllib.request.urlopen(fallback_url, timeout=60) as response:
                return response.read().decode("utf-8", errors="ignore")
        except urllib.error.HTTPError as fallback_exc:
            if fallback_exc.code == 404:
                raise EnsemblSpeciesNotFound(f"HTTP 404 at {fallback_url}")
            raise
    except urllib.error.URLError:
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        with urllib.request.urlopen(fallback_url, timeout=60) as response:
            return response.read().decode("utf-8", errors="ignore")


def species_root_for_release(base_url: str, release: str, species: str) -> str:
    if release in ("current", "latest"):
        return f"{base_url}/current_gtf/{species}/"
    if release.startswith("release-"):
        return f"{base_url}/{release}/gtf/{species}/"
    return f"{base_url}/release-{release}/gtf/{species}/"


def main() -> int:
    args = parse_args()
    species = args.species.strip().lower().replace(" ", "_")
    release = args.release.strip()
    base_url = args.base_url.rstrip("/")

    species_root = species_root_for_release(base_url, release, species)

    try:
        listing = fetch_text(species_root)
    except EnsemblSpeciesNotFound as exc:
        print(
            f"[ENSEMBL_GTF] Species '{species}' not found on Ensembl FTP "
            f"— dotplot-to-bp conversion will be skipped: {exc}",
            file=sys.stderr,
        )
        open(args.not_found_file, "w").close()
        return 0

    matches = [match for match in re.findall(r'href="([^"]+)"', listing) if re.search(r"\.gtf\.gz$", match)]
    preferred_matches = [match for match in matches if "abinitio" not in match.lower()]
    gtf_name = (preferred_matches or matches or [None])[0]
    if not gtf_name:
        raise RuntimeError(f"No .gtf.gz file found at {species_root}")

    gtf_url = f"{species_root}{gtf_name}"
    urllib.request.urlretrieve(gtf_url, args.output)

    with open(args.source_urls, "w", encoding="utf-8") as handle:
        handle.write(f"{gtf_url}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
