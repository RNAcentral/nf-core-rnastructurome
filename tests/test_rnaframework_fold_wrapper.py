from __future__ import annotations

from pathlib import Path


def test_rnaframework_fold_logs_outside_output_directory_until_completion() -> None:
    repo = Path(__file__).resolve().parents[1]
    text = (repo / "modules" / "local" / "rnaframework" / "fold" / "main.nf").read_text(
        encoding="utf-8"
    )

    assert 'log_tmp=\\$(mktemp "${prefix}_fold.XXXXXX.log")' in text
    assert '${xml_list} 2>&1 | tee "\\${log_tmp}"' in text
    assert 'tee ${prefix}_fold/rffold.log' not in text
    assert 'mv "\\${log_tmp}" ${prefix}_fold/rffold.log' in text


def test_rnaframework_docker_container_options_are_arch_aware() -> None:
    repo = Path(__file__).resolve().parents[1]
    text = (repo / "conf" / "modules.config").read_text(encoding="utf-8")

    assert text.count("def isArm64 = ((System.properties['os.arch'] ?: '').toLowerCase() in ['aarch64', 'arm64'])") >= 3
    assert text.count("return isArm64 ? '-u 1000:1000' : '--platform linux/amd64 -u 1000:1000'") >= 3
