import os
from pathlib import Path
from tactic_pipeline import version


def test_version_from_git_tag_env_var(monkeypatch):
    monkeypatch.setenv("GIT_TAG", "v1.0.0-rc.2")
    assert version.get_version() == "v1.0.0-rc.2"


def test_version_from_git_tag_without_v_prefix(monkeypatch):
    monkeypatch.setenv("GIT_TAG", "1.0.0-rc.2")
    assert version.get_version() == "v1.0.0-rc.2"


def test_version_from_legacy_tag_git(monkeypatch):
    monkeypatch.delenv("GIT_TAG", raising=False)
    monkeypatch.setenv("TAG_GIT", "v1.0.0-rc.1")
    assert version.get_version() == "v1.0.0-rc.1"


def test_version_from_tactic_version(monkeypatch):
    monkeypatch.delenv("GIT_TAG", raising=False)
    monkeypatch.delenv("TAG_GIT", raising=False)
    monkeypatch.setenv("TACTIC_VERSION", "v1.2.3")
    assert version.get_version() == "v1.2.3"


def test_version_ignores_na_and_empty(monkeypatch):
    monkeypatch.setenv("GIT_TAG", "NA")
    monkeypatch.delenv("TAG_GIT", raising=False)
    monkeypatch.delenv("TACTIC_VERSION", raising=False)
    # Should ignore 'NA' and fall back to git or default, not return 'vNA'
    ver = version.get_version()
    assert ver != "vNA"
    assert ver.startswith("v")


def test_version_from_version_file(tmp_path, monkeypatch):
    monkeypatch.delenv("GIT_TAG", raising=False)
    monkeypatch.delenv("TAG_GIT", raising=False)
    monkeypatch.delenv("TACTIC_VERSION", raising=False)

    version_file = Path(version.__file__).resolve().parent / "VERSION"
    try:
        version_file.write_text("v1.5.0-beta.1\n", encoding="utf-8")
        assert version.get_version() == "v1.5.0-beta.1"
    finally:
        if version_file.exists():
            version_file.unlink()

