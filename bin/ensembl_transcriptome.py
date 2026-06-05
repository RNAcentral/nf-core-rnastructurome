#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.error
import urllib.parse
import urllib.request
import re
from ftplib import FTP


class EnsemblSpeciesNotFound(Exception):
    pass


# EnsemblGenomes divisions with a flat species-directory structure (same as main Ensembl).
_EG_FLAT_DIVISIONS = [
    ("metazoa",  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa"),
    ("fungi",    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi"),
    ("plants",   "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants"),
    ("protists", "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists"),
]
# Bacteria: HTTPS to ftp.ensemblgenomes.ebi.ac.uk is blocked at EBI; use FTP instead.
_BACTERIA_FTP_HOST = "bacteria.ensembl.org"


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
    parser.add_argument(
        "--not-found-file",
        required=True,
        help="Path to write (empty) when the species is absent from all Ensembl FTPs; "
             "process exits 0 and the NCBI fallback is triggered.",
    )
    return parser.parse_args()


def fetch_text(url: str, timeout: int = 60) -> str:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="ignore")
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise EnsemblSpeciesNotFound(f"HTTP 404 at {url}") from exc
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        try:
            with urllib.request.urlopen(fallback_url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="ignore")
        except urllib.error.HTTPError as fallback_exc:
            if fallback_exc.code == 404:
                raise EnsemblSpeciesNotFound(f"HTTP 404 at {fallback_url}") from fallback_exc
            raise
    except urllib.error.URLError:
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        with urllib.request.urlopen(fallback_url, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="ignore")


def find_ensembl_file(listing_url: str, pattern: str, timeout: int = 60) -> str:
    listing = fetch_text(listing_url, timeout=timeout)
    matches = re.findall(r'href="([^"]+)"', listing)
    filtered = [match for match in matches if re.search(pattern, match)]
    if not filtered:
        raise RuntimeError(f"No file matching /{pattern}/ found at {listing_url}")
    return filtered[0]


def find_optional_ensembl_file(listing_url: str, pattern: str) -> str | None:
    if listing_url.startswith("ftp://"):
        parsed = urllib.parse.urlparse(listing_url)
        try:
            names = _list_ftp_names(parsed.netloc, parsed.path.rstrip("/"))
        except Exception:
            return None
        filtered = [n for n in names if re.search(pattern, n)]
        return filtered[0] if filtered else None
    try:
        listing = fetch_text(listing_url)
    except (urllib.error.URLError, EnsemblSpeciesNotFound):
        return None
    matches = re.findall(r'href="([^"]+)"', listing)
    filtered = [match for match in matches if re.search(pattern, match)]
    return filtered[0] if filtered else None


def _list_ftp_names(host: str, path: str) -> list[str]:
    """Return bare filenames in an FTP directory."""
    with FTP(host) as ftp:
        ftp.login()
        return [name.split("/")[-1] for name in ftp.nlst(path)]


def _find_bacteria_ftp(release: str, species: str) -> tuple[str, str, str]:
    """Locate bacteria FASTA files via FTP on bacteria.ensembl.org.

    Returns (cdna_dir_ftp_url, cdna_filename, ncrna_dir_ftp_url).
    Raises EnsemblSpeciesNotFound if the species is absent.
    """
    release_dir = release_path_for_value(release)  # e.g. "current_fasta"
    base_path = f"/pub/bacteria/{release_dir}"
    prefix = f"{species}_"

    try:
        top_names = _list_ftp_names(_BACTERIA_FTP_HOST, base_path)
    except Exception as exc:
        raise EnsemblSpeciesNotFound(
            f"FTP listing of {_BACTERIA_FTP_HOST}{base_path} failed: {exc}"
        ) from exc

    ftp_root = f"ftp://{_BACTERIA_FTP_HOST}"

    def _try_species(species_path: str) -> tuple[str, str, str] | None:
        try:
            cdna_names = _list_ftp_names(_BACTERIA_FTP_HOST, f"{species_path}/cdna")
            cdna_files = [n for n in cdna_names if re.search(r"\.cdna\.all\.fa\.gz$", n)]
            if cdna_files:
                return (
                    f"{ftp_root}{species_path}/cdna/",
                    cdna_files[0],
                    f"{ftp_root}{species_path}/ncrna/",
                )
        except Exception:
            pass
        return None

    # Direct match (no collection)
    for name in sorted(top_names):
        if name == species or name.startswith(prefix):
            result = _try_species(f"{base_path}/{name}")
            if result:
                return result

    # Collection-based layout (bacteria_N_collection/)
    collections = sorted(n for n in top_names if re.match(r"[a-z]+_\d+_collection$", n))
    for collection in collections:
        try:
            coll_names = _list_ftp_names(_BACTERIA_FTP_HOST, f"{base_path}/{collection}")
        except Exception:
            continue
        for name in sorted(coll_names):
            if name == species or name.startswith(prefix):
                result = _try_species(f"{base_path}/{collection}/{name}")
                if result:
                    print(
                        f"[ENSEMBL_TRANSCRIPTOME] '{species}' matched '{name}' "
                        f"in {collection} on EnsemblBacteria FTP.",
                        file=sys.stderr,
                    )
                    return result

    raise EnsemblSpeciesNotFound(
        f"Species '{species}' not found on {_BACTERIA_FTP_HOST} FTP"
    )


def release_path_for_value(release: str) -> str:
    if release in ("current", "latest"):
        return "current_fasta"
    if release.startswith("release-"):
        return f"{release}/fasta"
    return f"release-{release}/fasta"


def _try_flat_division(base_url: str, release: str, species: str) -> tuple[str, str, str]:
    """Return (cdna_dir, cdna_name, ncrna_dir) for a flat Ensembl/EnsemblGenomes division.

    Raises EnsemblSpeciesNotFound if the species directory does not exist.
    """
    species_root = f"{base_url}/{release_path_for_value(release)}/{species}"
    cdna_dir = f"{species_root}/cdna/"
    ncrna_dir = f"{species_root}/ncrna/"
    cdna_name = find_ensembl_file(cdna_dir, r"\.cdna\.all\.fa\.gz")
    return cdna_dir, cdna_name, ncrna_dir


def _prefix_candidates(base_url: str, release: str, prefix: str) -> list[tuple[str | None, str]]:
    """Return (subcollection_or_None, species) pairs under base_url matching prefix.

    Handles both flat FTP layouts (species directly under the release directory)
    and collection-based layouts where species are nested inside
    name_N_collection/ subdirectories (as used by EnsemblBacteria).
    prefix should be '{species}_' so that both exact names and
    strain/subspecies-suffixed variants are matched.
    """
    root = f"{base_url}/{release_path_for_value(release)}/"
    try:
        listing = fetch_text(root, timeout=30)
    except Exception:
        return []
    all_dirs = re.findall(r'href="([a-z][a-z0-9_]+)/"', listing)
    exact = prefix.rstrip("_")
    results: list[tuple[str | None, str]] = [
        (None, d) for d in sorted(all_dirs) if d == exact or d.startswith(prefix)
    ]
    collections = [d for d in all_dirs if re.match(r'[a-z]+_\d+_collection$', d)]
    for collection in sorted(collections):
        try:
            coll_listing = fetch_text(f"{root}{collection}/", timeout=30)
        except Exception:
            continue
        coll_dirs = re.findall(r'href="([a-z][a-z0-9_]+)/"', coll_listing)
        for d in sorted(d for d in coll_dirs if d == exact or d.startswith(prefix)):
            results.append((collection, d))
    return results


def _find_species(base_url: str, release: str, species: str) -> tuple[str, str, str, str]:
    """Try every Ensembl source in priority order.

    Returns (cdna_dir, cdna_name, ncrna_dir, division_label).
    Raises EnsemblSpeciesNotFound if the species is absent from all sources.
    """
    # 1. Primary URL (main Ensembl — eukaryotes)
    try:
        cdna_dir, cdna_name, ncrna_dir = _try_flat_division(base_url, release, species)
        return cdna_dir, cdna_name, ncrna_dir, "Ensembl"
    except EnsemblSpeciesNotFound:
        pass

    # 2. EnsemblGenomes flat divisions (metazoa, fungi, plants, protists)
    for division_name, division_base in _EG_FLAT_DIVISIONS:
        try:
            cdna_dir, cdna_name, ncrna_dir = _try_flat_division(division_base, release, species)
            return cdna_dir, cdna_name, ncrna_dir, f"EnsemblGenomes/{division_name}"
        except EnsemblSpeciesNotFound:
            continue

    # 3. Prefix fallback across HTTPS sources — handles strain/subspecies suffixes
    #    (e.g. 'escherichia_coli' matching 'escherichia_coli_k_12').
    prefix = f"{species}_"
    https_sources = [("Ensembl", base_url)] + [
        (f"EnsemblGenomes/{n}", b) for n, b in _EG_FLAT_DIVISIONS
    ]
    for label, src_url in https_sources:
        for subcollection, candidate in _prefix_candidates(src_url, release, prefix):
            release_path = release_path_for_value(release)
            species_root = (
                f"{src_url}/{release_path}/{subcollection}/{candidate}"
                if subcollection
                else f"{src_url}/{release_path}/{candidate}"
            )
            cdna_dir = f"{species_root}/cdna/"
            try:
                cdna_name = find_ensembl_file(cdna_dir, r"\.cdna\.all\.fa\.gz", timeout=10)
                ncrna_dir = f"{species_root}/ncrna/"
                print(
                    f"[ENSEMBL_TRANSCRIPTOME] '{species}' matched '{candidate}' on {label}.",
                    file=sys.stderr,
                )
                return cdna_dir, cdna_name, ncrna_dir, label
            except (EnsemblSpeciesNotFound, RuntimeError):
                continue

    # 4. EnsemblBacteria via FTP — HTTPS to ftp.ensemblgenomes.ebi.ac.uk is blocked at EBI.
    try:
        cdna_dir, cdna_name, ncrna_dir = _find_bacteria_ftp(release, species)
        return cdna_dir, cdna_name, ncrna_dir, "EnsemblBacteria"
    except EnsemblSpeciesNotFound:
        pass

    raise EnsemblSpeciesNotFound(
        f"Species '{species}' not found on Ensembl, EnsemblGenomes, or EnsemblBacteria FTP"
    )


def main() -> int:
    args = parse_args()
    species = args.species.strip().lower().replace(" ", "_")
    release = args.release.strip()
    base_url = args.base_url.rstrip("/")
    warnings_log = args.warnings_log

    try:
        cdna_dir, cdna_name, ncrna_dir, division = _find_species(base_url, release, species)
    except EnsemblSpeciesNotFound as exc:
        print(
            f"[ENSEMBL_TRANSCRIPTOME] Species '{species}' not found on any Ensembl FTP "
            f"— will fall back to NCBI: {exc}",
            file=sys.stderr,
        )
        with open(args.not_found_file, "w", encoding="utf-8"):
            pass
        return 0

    print(
        f"[ENSEMBL_TRANSCRIPTOME] Resolved '{species}' via {division}.",
        file=sys.stderr,
    )

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
