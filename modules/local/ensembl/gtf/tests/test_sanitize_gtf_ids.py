"""Guard tests for bin/sanitize_gtf_ids.py (GTF transcript/gene ID sanitization)."""
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[5] / "bin" / "sanitize_gtf_ids.py"


def _run(tmp_path, gtf_text):
    src = tmp_path / "in.gtf"
    out = tmp_path / "out.gtf"
    src.write_text(gtf_text)
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(src), str(out)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    return out.read_text()


def test_parenthesised_trna_ids_are_sanitized(tmp_path):
    line = 'chrI\tSGD\texon\t1\t72\t.\t+\t.\tgene_id "tK(UUU)K"; transcript_id "tK(UUU)K_tRNA"\n'
    out = _run(tmp_path, line)
    assert 'transcript_id "tK_UUU_K_tRNA"' in out
    assert 'gene_id "tK_UUU_K"' in out
    assert "(" not in out and ")" not in out


def test_dash_and_dot_and_other_attributes_are_preserved(tmp_path):
    # Dash (yeast isoform) and dot (version suffix) are kept; non-target attributes are untouched.
    attrs = 'gene_id "YNL042W-B"; transcript_id "YNL042W-B.1_mRNA"; note "keep(this)"'
    out = _run(tmp_path, f"chrI\tSGD\texon\t1\t9\t.\t+\t.\t{attrs}\n")
    assert 'transcript_id "YNL042W-B.1_mRNA"' in out
    assert 'note "keep(this)"' in out


def test_comment_lines_pass_through(tmp_path):
    out = _run(tmp_path, '#!genome-build R64\nchrI\tSGD\texon\t1\t9\t.\t+\t.\ttranscript_id "a b"\n')
    assert out.startswith("#!genome-build R64\n")
    assert 'transcript_id "a_b"' in out
