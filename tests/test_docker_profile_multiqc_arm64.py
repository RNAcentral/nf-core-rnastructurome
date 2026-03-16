from __future__ import annotations

from pathlib import Path


def test_docker_profile_falls_back_to_conda_for_multiqc_on_arm64() -> None:
    repo = Path(__file__).resolve().parents[1]
    config_text = (repo / "nextflow.config").read_text(encoding="utf-8")

    assert "conda.enabled           = (System.properties['os.arch'] ?: '').toLowerCase() in ['aarch64', 'arm64']" in config_text
    assert "MultiQC's current Wave image crashes under amd64 emulation on arm64 Docker hosts." in config_text
    assert "withName: '(^|.*:)MULTIQC$'" in config_text
    assert "container = ((System.properties['os.arch'] ?: '').toLowerCase() in ['aarch64', 'arm64']) ? null : 'community.wave.seqera.io/library/multiqc:1.32--d58f60e4deb769bf'" in config_text
    assert "containerOptions = null" in config_text
