from __future__ import annotations

from pathlib import Path


def test_utils_nfschema_plugin_has_local_plugin_config() -> None:
    repo = Path(__file__).resolve().parents[1]
    config_path = repo / "subworkflows" / "nf-core" / "utils_nfschema_plugin" / "nextflow.config"

    text = config_path.read_text(encoding="utf-8")

    assert 'id "nf-schema@2.6.1"' in text


def test_nf_core_parent_scope_has_plugin_config() -> None:
    repo = Path(__file__).resolve().parents[1]
    config_path = repo / "subworkflows" / "nf-core" / "nextflow.config"

    text = config_path.read_text(encoding="utf-8")

    assert 'id "nf-schema@2.6.1"' in text
