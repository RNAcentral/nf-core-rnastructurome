from __future__ import annotations

from pathlib import Path


def test_utils_nextflow_pipeline_uses_explicit_entry_params() -> None:
    repo = Path(__file__).resolve().parents[1]
    util_text = (
        repo / "subworkflows" / "nf-core" / "utils_nextflow_pipeline" / "main.nf"
    ).read_text(encoding="utf-8")
    caller_text = (
        repo / "subworkflows" / "local" / "utils_nfcore_rnastructurome_pipeline" / "main.nf"
    ).read_text(encoding="utf-8")
    entry_text = (repo / "main.nf").read_text(encoding="utf-8")

    assert "entry_params         //    map: params captured at the entry workflow" in util_text
    assert "dumpParametersToJSON(outdir, entry_params ?: [:])" in util_text
    assert "JsonOutput.toJson(params)" not in util_text
    assert "entry_params_input // map: raw params captured at the entry workflow" in caller_text
    assert "entry_params" in entry_text
    assert "buildEntryParamsMap(params)" in entry_text


def test_completion_email_uses_explicit_max_multiqc_email_size() -> None:
    repo = Path(__file__).resolve().parents[1]
    util_text = (
        repo / "subworkflows" / "nf-core" / "utils_nfcore_pipeline" / "main.nf"
    ).read_text(encoding="utf-8")
    caller_text = (
        repo / "subworkflows" / "local" / "utils_nfcore_rnastructurome_pipeline" / "main.nf"
    ).read_text(encoding="utf-8")
    entry_text = (repo / "main.nf").read_text(encoding="utf-8")

    assert "max_multiqc_email_size_input=0" in util_text
    assert "params.containsKey('max_multiqc_email_size')" not in util_text
    assert "max_multiqc_email_size // string|memoryunit: maximum MultiQC size to attach to completion emails" in caller_text
    assert "params.max_multiqc_email_size" in entry_text
