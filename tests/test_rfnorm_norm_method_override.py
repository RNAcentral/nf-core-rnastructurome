from __future__ import annotations

from pathlib import Path


def test_rfnorm_norm_method_override_is_declared_and_validated() -> None:
    repo = Path(__file__).resolve().parents[1]
    workflow_text = (repo / "workflows" / "rnastructurome.nf").read_text(encoding="utf-8")
    config_text = (repo / "nextflow.config").read_text(encoding="utf-8")
    docs_text = (repo / "docs" / "usage.md").read_text(encoding="utf-8")

    assert "rfnorm_norm_method                : null" in workflow_text
    assert "resolveRfNormNormMethod(pipeline_config, scoringMethod)" in workflow_text
    assert "Unsupported rf-norm normalization method" in workflow_text
    assert "rfnorm_norm_method          = null" in config_text
    assert "--rfnorm_norm_method 2|3" in docs_text
