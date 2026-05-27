import math
import subprocess
import sys
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "bin" / "r2dt_colour_svg.py"

_SVG_NS = "http://www.w3.org/2000/svg"

# ── helpers ──────────────────────────────────────────────────────────────────

sys.path.insert(0, str(REPO_ROOT / "bin"))
from r2dt_colour_svg import _shape_colour, _parse_rfnorm_xml, load_reactivities, colour_svg


def _make_svg(tmp_path: Path, name: str, positions: list[int]) -> Path:
    """Create a minimal R2DT-style SVG with nucleotide <g> elements."""
    groups = "\n".join(
        f"  <g>\n    <title>{pos} (position.label in template: {pos}.A)</title>\n    <text x=\"0\" y=\"0\">A</text>\n  </g>"
        for pos in positions
    )
    svg = (
        f'<?xml version="1.0" encoding="utf-8"?>\n'
        f'<svg xmlns="{_SVG_NS}" xmlns:xlink="http://www.w3.org/1999/xlink">\n'
        f'{groups}\n'
        f'</svg>\n'
    )
    p = tmp_path / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(svg)
    return p


def _make_xml(tmp_path: Path, tid: str, values: list) -> Path:
    """Create a minimal rf-norm XML file."""
    reactivity_str = ",".join(
        "nan" if (v is None or (isinstance(v, float) and math.isnan(v))) else str(v)
        for v in values
    )
    xml = textwrap.dedent(f"""\
        <?xml version="1.0" encoding="utf-8"?>
        <data>
            <transcript id="{tid}">
                <length>{len(values)}</length>
                <reactivity>{reactivity_str}</reactivity>
            </transcript>
        </data>
    """)
    p = tmp_path / f"{tid}.xml"
    p.write_text(xml)
    return p


# ── _shape_colour ─────────────────────────────────────────────────────────────

def test_shape_colour_nan():
    assert _shape_colour(float("nan")) == "#808080"


def test_shape_colour_none():
    assert _shape_colour(None) == "#808080"


def test_shape_colour_negative():
    assert _shape_colour(-0.1) == "#808080"


def test_shape_colour_low():
    assert _shape_colour(0.0) == "#000000"
    assert _shape_colour(0.39) == "#000000"


def test_shape_colour_medium():
    assert _shape_colour(0.40) == "#E07B54"
    assert _shape_colour(0.84) == "#E07B54"


def test_shape_colour_high():
    assert _shape_colour(0.85) == "#C0392B"
    assert _shape_colour(1.0) == "#C0392B"


# ── _parse_rfnorm_xml ─────────────────────────────────────────────────────────

def test_parse_rfnorm_xml_basic(tmp_path):
    xml_path = _make_xml(tmp_path, "URS001", [0.1, float("nan"), 0.9])
    result = _parse_rfnorm_xml(xml_path)
    assert "URS001" in result
    vals = result["URS001"]
    assert vals[0] == pytest.approx(0.1)
    assert math.isnan(vals[1])
    assert vals[2] == pytest.approx(0.9)


def test_parse_rfnorm_xml_invalid_file(tmp_path):
    bad = tmp_path / "bad.xml"
    bad.write_text("not valid xml <<<")
    result = _parse_rfnorm_xml(bad)
    assert result == {}


def test_parse_rfnorm_xml_multiple_transcripts(tmp_path):
    xml_content = textwrap.dedent("""\
        <?xml version="1.0" encoding="utf-8"?>
        <data>
            <transcript id="A"><length>2</length><reactivity>0.1,0.2</reactivity></transcript>
            <transcript id="B"><length>2</length><reactivity>nan,0.5</reactivity></transcript>
        </data>
    """)
    p = tmp_path / "multi.xml"
    p.write_text(xml_content)
    result = _parse_rfnorm_xml(p)
    assert set(result.keys()) == {"A", "B"}
    assert result["A"] == pytest.approx([0.1, 0.2])
    assert math.isnan(result["B"][0])
    assert result["B"][1] == pytest.approx(0.5)


# ── load_reactivities ─────────────────────────────────────────────────────────

def test_load_reactivities_single_replicate(tmp_path):
    xml_dir = tmp_path / "xml_input_rep1"
    xml_dir.mkdir()
    _make_xml(xml_dir, "URS001", [0.1, 0.5, 0.9])
    result = load_reactivities(tmp_path)
    assert "URS001" in result
    assert result["URS001"] == pytest.approx([0.1, 0.5, 0.9])


def test_load_reactivities_averages_replicates(tmp_path):
    for i, vals in enumerate(([0.0, 0.6], [0.2, 0.8])):
        d = tmp_path / f"xml_input_rep{i}"
        d.mkdir()
        _make_xml(d, "URS001", vals)
    result = load_reactivities(tmp_path)
    assert result["URS001"] == pytest.approx([0.1, 0.7])


def test_load_reactivities_nan_ignored_in_average(tmp_path):
    # rep1 has NaN at position 0; rep2 has 0.4 — average should be 0.4, not NaN
    d1 = tmp_path / "xml_input_rep1"
    d1.mkdir()
    _make_xml(d1, "URS001", [float("nan"), 0.2])
    d2 = tmp_path / "xml_input_rep2"
    d2.mkdir()
    _make_xml(d2, "URS001", [0.4, 0.6])
    result = load_reactivities(tmp_path)
    assert result["URS001"][0] == pytest.approx(0.4)
    assert result["URS001"][1] == pytest.approx(0.4)


def test_load_reactivities_empty_dir(tmp_path):
    result = load_reactivities(tmp_path)
    assert result == {}


def test_load_reactivities_ignores_non_xml_input_dirs(tmp_path):
    # XML in a plain subdir (not matching xml_input*) should be ignored
    d = tmp_path / "some_other_dir"
    d.mkdir()
    _make_xml(d, "URS001", [0.5])
    result = load_reactivities(tmp_path)
    assert result == {}


# ── colour_svg ────────────────────────────────────────────────────────────────

def test_colour_svg_colours_nucleotides(tmp_path):
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", [1, 2, 3])
    reactivities = [0.2, 0.6, 0.9]  # black, orange, red
    out_path = tmp_path / "out" / "URS001.svg"

    n = colour_svg(svg_path, reactivities, out_path)

    assert n == 3
    tree = ET.parse(out_path)
    texts = list(tree.getroot().iter(f"{{{_SVG_NS}}}text"))
    assert len(texts) == 3
    styles = [t.get("style", "") for t in texts]
    assert any("#000000" in s for s in styles)
    assert any("#E07B54" in s for s in styles)
    assert any("#C0392B" in s for s in styles)


def test_colour_svg_nan_reactivity_gives_grey(tmp_path):
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", [1])
    out_path = tmp_path / "out" / "URS001.svg"
    colour_svg(svg_path, [float("nan")], out_path)
    tree = ET.parse(out_path)
    text = next(tree.getroot().iter(f"{{{_SVG_NS}}}text"))
    assert "#808080" in text.get("style", "")


def test_colour_svg_out_of_range_position_skipped(tmp_path):
    # SVG has position 5 but reactivity list only has 3 entries
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", [5])
    out_path = tmp_path / "out" / "URS001.svg"
    n = colour_svg(svg_path, [0.1, 0.2, 0.3], out_path)
    assert n == 0


def test_colour_svg_invalid_svg(tmp_path):
    bad = tmp_path / "bad.svg"
    bad.write_text("not xml <<<")
    n = colour_svg(bad, [0.5], tmp_path / "out.svg")
    assert n == 0


def test_colour_svg_creates_output_dir(tmp_path):
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", [1])
    nested_out = tmp_path / "deep" / "nested" / "URS001.svg"
    colour_svg(svg_path, [0.5], nested_out)
    assert nested_out.exists()


# ── CLI integration ───────────────────────────────────────────────────────────

def test_cli_end_to_end(tmp_path):
    svg_dir = tmp_path / "svg"
    svg_dir.mkdir()
    _make_svg(svg_dir, "URS001.svg", [1, 2, 3])

    xml_input_dir = tmp_path / "xml_input_rep1"
    xml_input_dir.mkdir()
    _make_xml(xml_input_dir, "URS001", [0.1, 0.5, 0.9])

    out_dir = tmp_path / "coloured"

    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--svg-dir", str(svg_dir),
         "--xml-search-dir", str(tmp_path),
         "--out-dir", str(out_dir)],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    assert (out_dir / "URS001.svg").exists()


def test_cli_exits_nonzero_when_no_reactivities(tmp_path):
    svg_dir = tmp_path / "svg"
    svg_dir.mkdir()
    _make_svg(svg_dir, "URS001.svg", [1])
    out_dir = tmp_path / "coloured"

    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--svg-dir", str(svg_dir),
         "--xml-search-dir", str(tmp_path),
         "--out-dir", str(out_dir)],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0


def test_cli_skips_svg_with_no_matching_transcript(tmp_path):
    svg_dir = tmp_path / "svg"
    svg_dir.mkdir()
    _make_svg(svg_dir, "UNKNOWN.svg", [1])

    xml_input_dir = tmp_path / "xml_input_rep1"
    xml_input_dir.mkdir()
    _make_xml(xml_input_dir, "URS001", [0.5])

    out_dir = tmp_path / "coloured"

    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--svg-dir", str(svg_dir),
         "--xml-search-dir", str(tmp_path),
         "--out-dir", str(out_dir)],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert not (out_dir / "UNKNOWN.svg").exists()
