from __future__ import annotations

from pathlib import Path


def test_main_workflow_test_targets_nfcore_rnastructurome() -> None:
    repo = Path(__file__).resolve().parents[1]
    text = (repo / "tests" / "main.nf.test").read_text(encoding="utf-8")

    assert 'workflow "NFCORE_RNASTRUCTUROME"' in text


def test_colocated_workflow_and_subworkflow_tests_exist() -> None:
    repo = Path(__file__).resolve().parents[1]

    workflow_test = (repo / "workflows" / "tests" / "rnastructurome.nf.test").read_text(encoding="utf-8")
    initialisation_test = (
        repo
        / "subworkflows"
        / "local"
        / "utils_nfcore_rnastructurome_pipeline"
        / "tests"
        / "pipeline_initialisation.workflow.nf.test"
    ).read_text(encoding="utf-8")
    completion_test = (
        repo
        / "subworkflows"
        / "local"
        / "utils_nfcore_rnastructurome_pipeline"
        / "tests"
        / "pipeline_completion.workflow.nf.test"
    ).read_text(encoding="utf-8")

    assert 'workflow "RNASTRUCTUROME"' in workflow_test
    assert 'workflow "PIPELINE_INITIALISATION"' in initialisation_test
    assert 'workflow "PIPELINE_COMPLETION"' in completion_test
