from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_nf_metro_animation_injection(tmp_path: Path) -> None:
    repo = Path(__file__).resolve().parents[1]
    source_svg = repo / "nf-metro" / "rendered" / "dag.metro.2.multiqc_last.svg"
    script = repo / "nf-metro" / "animate_svg.py"
    output_svg = tmp_path / "animated.svg"

    subprocess.run(
        [sys.executable, str(script), str(source_svg), str(output_svg)],
        check=True,
    )

    text = output_svg.read_text(encoding="utf-8")
    assert "animateMotion" in text
    assert "motion-path-rtstop-" in text
    assert "motion-path-map-" in text
