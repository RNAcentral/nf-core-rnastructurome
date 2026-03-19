from __future__ import annotations

from pathlib import Path


def test_dotplot2bp_emits_probability_binned_bp_colors() -> None:
    repo = Path(__file__).resolve().parents[1]
    module_text = (
        repo / "modules" / "local" / "rnaframework" / "dotplot2bp" / "main.nf"
    ).read_text(encoding="utf-8")

    assert '(189, 189, 189, "5-10% probability")' in module_text
    assert '(242, 204, 84, "10-40% probability")' in module_text
    assert '(120, 182, 220, "40-70% probability")' in module_text
    assert '(126, 198, 143, "70-100% probability")' in module_text
    assert "probability = math.pow(10.0, -score)" in module_text
    assert "writer.write(f\"{entry['seqname']}\\\\t{start}\\\\t{start}\\\\t{end}\\\\t{end}\\\\t{color_index}\\\\n\")" in module_text
