from __future__ import annotations

import json
from pathlib import Path


def test_rffold_window_is_enabled_by_default() -> None:
    repo = Path(__file__).resolve().parents[1]
    config_text = (repo / "nextflow.config").read_text(encoding="utf-8")
    schema = json.loads((repo / "nextflow_schema.json").read_text(encoding="utf-8"))
    docs_text = (repo / "docs" / "usage.md").read_text(encoding="utf-8")

    assert "rffold_window              = true" in config_text
    assert (
        schema["$defs"]["rnaframework_options"]["properties"]["rffold_window"]["default"]
        is True
    )
    assert "Windowed folding is enabled by default (`-w`); disable with `--rffold_window false`." in docs_text
