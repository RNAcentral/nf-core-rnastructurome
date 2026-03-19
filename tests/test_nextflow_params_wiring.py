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


def test_genome_param_and_igenomes_config_are_removed() -> None:
    repo = Path(__file__).resolve().parents[1]
    main_text = (repo / "main.nf").read_text(encoding="utf-8")
    config_text = (repo / "nextflow.config").read_text(encoding="utf-8")
    schema_text = (repo / "nextflow_schema.json").read_text(encoding="utf-8")
    docs_text = (repo / "docs" / "usage.md").read_text(encoding="utf-8")
    test_config_text = (repo / "conf" / "test.config").read_text(encoding="utf-8")
    test_full_config_text = (repo / "conf" / "test_full.config").read_text(encoding="utf-8")

    assert "genome                           : all_params.genome" not in main_text
    assert "igenomes_base" not in config_text
    assert "igenomes_ignore" not in config_text
    assert "params.genome" not in config_text
    assert "conf/igenomes.config" not in config_text
    assert '"genome": {' not in schema_text
    assert '"igenomes_ignore": {' not in schema_text
    assert '"igenomes_base": {' not in schema_text
    assert "--genome" not in docs_text
    assert "genome" not in test_config_text
    assert "genome" not in test_full_config_text
    assert not (repo / "conf" / "igenomes.config").exists()
    assert not (repo / "conf" / "igenomes_ignored.config").exists()
