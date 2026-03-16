from __future__ import annotations

from pathlib import Path


def test_test_profile_uses_human_subset_and_raw_rfnorm() -> None:
    repo = Path(__file__).resolve().parents[1]
    test_config = (repo / "conf" / "test.config").read_text(encoding="utf-8")
    default_nf_test = (repo / "tests" / "default.nf.test").read_text(encoding="utf-8")

    assert 'input        = "${projectDir}/assets/testdata/human/rnastruct00001_samplesheet.first2.trim_subset_top3.csv"' in test_config
    assert 'fasta        = "${projectDir}/assets/testdata/human/Homo_sapiens.GRCh38.selected_top3_subset.fa"' in test_config
    assert "rfnorm_raw   = true" in test_config
    assert 'fasta = "${projectDir}/tests/data/ref.fa"' in default_nf_test
    assert "rfnorm_raw = false" in default_nf_test
