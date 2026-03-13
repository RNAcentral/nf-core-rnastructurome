from __future__ import annotations

import json
from pathlib import Path


def test_workspace_settings_exclude_venv_from_nextflow_scan() -> None:
    repo = Path(__file__).resolve().parents[1]
    settings_path = repo / ".vscode" / "settings.json"

    settings = json.loads(settings_path.read_text(encoding="utf-8"))

    assert ".venv" in settings["nextflow.files.exclude"]
    assert settings["files.exclude"]["**/.venv"] is True
    assert settings["search.exclude"]["**/.venv"] is True
    assert settings["files.watcherExclude"]["**/.venv/**"] is True
