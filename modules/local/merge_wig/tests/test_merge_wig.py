import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[4] / "bin" / "merge_wig.py"
FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_merge_wig_concatenates_wig_files_and_drops_per_file_track_headers(tmp_path):
    shutil.copy(FIXTURES / "TX1.wig", tmp_path / "TX1.wig")
    shutil.copy(FIXTURES / "TX2.wig", tmp_path / "TX2.wig")

    result = subprocess.run(
        [sys.executable, str(SCRIPT), "group1"],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "group1.merged.wig").read_text() == "\n".join(
        [
            "track type=wiggle_0",
            "fixedStep chrom=TX1 start=1 step=1",
            "1",
            "2",
            "variableStep chrom=TX2",
            "1 5",
            "",
        ]
    )


def test_merge_wig_fails_when_no_wig_files_are_present(tmp_path):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "empty"],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 1
    assert "No WIG files were provided for merging" in result.stderr
    assert not (tmp_path / "empty.merged.wig").exists()
