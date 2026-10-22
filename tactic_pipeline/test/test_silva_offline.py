import logging
import shutil
import pytest
import requests

from tactic_pipeline import processing_helper, silva_remote
from tactic_pipeline.processing_helper import (
    find_local_silva_database,
    prepare_silva_database,
    DatabaseUnavailableException,
    SILVA_OFFLINE_PREFIX,
    SUCCESS_MARKER,
)


def _create_release(base_dir, release_name, *, has_arb=True, marker=SUCCESS_MARKER):
    release_dir = base_dir / "SILVA" / release_name
    release_dir.mkdir(parents=True, exist_ok=True)
    if has_arb:
        (release_dir / "SILVA.arb").write_bytes(b"mock arb content")
    if marker is not None:
        (release_dir / "SILVA.sidx.success").write_text(
            f"Verified release {release_name}\n{marker}\n",
            encoding="utf-8",
        )
    return release_dir


class FakeResponse:
    def __init__(self, text="", status_code=200):
        self.text = text
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(str(self.status_code))


def test_offline_fallback_picks_newest_verified_release(tmp_path, monkeypatch, caplog):
    _create_release(tmp_path, "138.1")
    _create_release(tmp_path, "138.2")
    # 139 is newer but lacks success marker
    _create_release(tmp_path, "139", marker=None)

    monkeypatch.setenv("SILVA_SOURCE", "ftp")
    monkeypatch.setattr(
        requests, "get",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            requests.ConnectionError("Connection refused by host")
        ),
    )
    monkeypatch.setattr(
        requests, "head",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            requests.ConnectionError("Connection refused by host")
        ),
    )

    logger = logging.getLogger("test_offline_fallback")
    with caplog.at_level(logging.WARNING, logger=logger.name):
        dir_path, version = prepare_silva_database(
            dest_dir=str(tmp_path),
            logger_obj=logger,
        )

    assert dir_path == str(tmp_path / "SILVA" / "138.2")
    assert version == "138.2"

    warn_msgs = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warn_msgs) >= 2
    assert warn_msgs[-2].startswith("Could not reach the SILVA server: ")
    assert warn_msgs[-1] == f"{SILVA_OFFLINE_PREFIX}138.2"


def test_offline_fallback_with_worker_source(tmp_path, monkeypatch, caplog):
    _create_release(tmp_path, "138.1")
    _create_release(tmp_path, "138.2")
    _create_release(tmp_path, "139", marker=None)

    monkeypatch.setenv("SILVA_SOURCE", "worker")

    def _worker_fails(**kwargs):
        raise silva_remote.SilvaRemoteUnavailable("Worker endpoint unreachable")

    monkeypatch.setattr(
        processing_helper,
        "prepare_prebuilt_silva_database",
        _worker_fails,
    )
    monkeypatch.setattr(
        requests, "get",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            requests.ConnectionError("Connection refused by host")
        ),
    )
    monkeypatch.setattr(
        requests, "head",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            requests.ConnectionError("Connection refused by host")
        ),
    )

    logger = logging.getLogger("test_worker_offline_fallback")
    with caplog.at_level(logging.WARNING, logger=logger.name):
        dir_path, version = prepare_silva_database(
            dest_dir=str(tmp_path),
            logger_obj=logger,
        )

    assert dir_path == str(tmp_path / "SILVA" / "138.2")
    assert version == "138.2"

    warn_msgs = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warn_msgs) >= 2
    assert warn_msgs[-1] == f"{SILVA_OFFLINE_PREFIX}138.2"


def test_offline_marker_lacking_success_marker_is_ignored(tmp_path, monkeypatch):
    # 140 has marker file, but lacks SUCCESS_MARKER
    _create_release(tmp_path, "140", marker="INCOMPLETE_DOWNLOAD")
    _create_release(tmp_path, "138.2")

    monkeypatch.setenv("SILVA_SOURCE", "ftp")
    monkeypatch.setattr(
        requests, "get",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            requests.ConnectionError("Network down")
        ),
    )
    monkeypatch.setattr(
        requests, "head",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            requests.ConnectionError("Network down")
        ),
    )

    dir_path, version = prepare_silva_database(dest_dir=str(tmp_path))
    assert dir_path == str(tmp_path / "SILVA" / "138.2")
    assert version == "138.2"

    # If 140 is the only release, it must be ignored and raise DatabaseUnavailableException
    shutil.rmtree(tmp_path / "SILVA" / "138.2")
    with pytest.raises(DatabaseUnavailableException):
        prepare_silva_database(dest_dir=str(tmp_path))


def test_offline_no_local_copy_raises_database_unavailable(tmp_path, monkeypatch):
    monkeypatch.setenv("SILVA_SOURCE", "ftp")
    conn_err = requests.ConnectionError("Network unreachable")
    monkeypatch.setattr(
        requests, "get",
        lambda *args, **kwargs: (_ for _ in ()).throw(conn_err),
    )
    monkeypatch.setattr(
        requests, "head",
        lambda *args, **kwargs: (_ for _ in ()).throw(conn_err),
    )

    with pytest.raises(DatabaseUnavailableException) as exc_info:
        prepare_silva_database(dest_dir=str(tmp_path))

    assert str(tmp_path) in str(exc_info.value)
    assert "internet connection" in str(exc_info.value).lower()
    # Ensure {exc} is included in the exception message
    assert str(conn_err) in str(exc_info.value)
    assert exc_info.value.__cause__ is conn_err


def test_offline_fallback_on_index_error_when_arb_requests_fail(tmp_path, monkeypatch, caplog):
    _create_release(tmp_path, "138.2")
    monkeypatch.setenv("SILVA_SOURCE", "ftp")

    # Base URL succeeds with release directory, but ARB_files fetch fails
    def mock_get(url, **kwargs):
        if "ARB_files" in url:
            raise requests.ConnectionError("ARB_files index unreachable")
        return FakeResponse('<a href="release_138.2/">release_138.2/</a>')

    monkeypatch.setattr(requests, "get", mock_get)

    logger = logging.getLogger("test_index_error")
    with caplog.at_level(logging.WARNING, logger=logger.name):
        dir_path, version = prepare_silva_database(
            dest_dir=str(tmp_path),
            logger_obj=logger,
        )

    assert dir_path == str(tmp_path / "SILVA" / "138.2")
    assert version == "138.2"
    warn_msgs = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert warn_msgs[-1] == f"{SILVA_OFFLINE_PREFIX}138.2"

    # When no local copy exists, raises DatabaseUnavailableException with {exc}
    shutil.rmtree(tmp_path / "SILVA" / "138.2")
    with pytest.raises(DatabaseUnavailableException) as exc_info:
        prepare_silva_database(dest_dir=str(tmp_path))
    assert "No accessible SILVA releases found" in str(exc_info.value)


def test_offline_fallback_on_value_error_captive_portal(tmp_path, monkeypatch, caplog):
    _create_release(tmp_path, "138.2")
    monkeypatch.setenv("SILVA_SOURCE", "ftp")

    # Captive portal returns HTML login page with no release_ links
    captive_html = "<html><body>Please log in to Hotel Guest Wi-Fi</body></html>"
    monkeypatch.setattr(requests, "get", lambda *a, **kw: FakeResponse(captive_html))

    logger = logging.getLogger("test_captive_portal")
    with caplog.at_level(logging.WARNING, logger=logger.name):
        dir_path, version = prepare_silva_database(
            dest_dir=str(tmp_path),
            logger_obj=logger,
        )

    assert dir_path == str(tmp_path / "SILVA" / "138.2")
    assert version == "138.2"
    warn_msgs = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert "No release_ directories found" in warn_msgs[-2]
    assert warn_msgs[-1] == f"{SILVA_OFFLINE_PREFIX}138.2"

    # When no local copy exists, raises DatabaseUnavailableException with {exc}
    shutil.rmtree(tmp_path / "SILVA" / "138.2")
    with pytest.raises(DatabaseUnavailableException) as exc_info:
        prepare_silva_database(dest_dir=str(tmp_path))
    assert "No release_ directories found" in str(exc_info.value)


def test_online_download_success_no_offline_warning(tmp_path, monkeypatch, caplog):
    _create_release(tmp_path, "138.2")

    monkeypatch.setenv("SILVA_SOURCE", "ftp")
    downloaded_arb = tmp_path / "SILVA" / "144" / "SILVA.arb"

    def mock_download(**kwargs):
        downloaded_arb.parent.mkdir(parents=True, exist_ok=True)
        downloaded_arb.write_bytes(b"fresh arb")
        return str(downloaded_arb), "144"

    def mock_index(**kwargs):
        pass

    monkeypatch.setattr(processing_helper, "download_silva_databases", mock_download)
    monkeypatch.setattr(processing_helper, "index_silva_database", mock_index)

    logger = logging.getLogger("test_online_success")
    with caplog.at_level(logging.WARNING, logger=logger.name):
        dir_path, version = prepare_silva_database(
            dest_dir=str(tmp_path),
            logger_obj=logger,
        )

    assert dir_path == str(downloaded_arb.parent)
    assert version == "144"
    assert SILVA_OFFLINE_PREFIX not in caplog.text


def test_version_ordering_and_underscore_handling(tmp_path):
    # None when no directory
    assert find_local_silva_database(str(tmp_path / "nonexistent")) is None
    assert find_local_silva_database("") is None

    _create_release(tmp_path, "138.1")
    _create_release(tmp_path, "138.2")
    _create_release(tmp_path, "139")

    # 139 > 138.2 > 138.1
    best = find_local_silva_database(str(tmp_path))
    assert best == (str(tmp_path / "SILVA" / "139"), "139")

    shutil.rmtree(tmp_path / "SILVA" / "139")
    best = find_local_silva_database(str(tmp_path))
    assert best == (str(tmp_path / "SILVA" / "138.2"), "138.2")

    shutil.rmtree(tmp_path / "SILVA" / "138.2")
    best = find_local_silva_database(str(tmp_path))
    assert best == (str(tmp_path / "SILVA" / "138.1"), "138.1")

    # Underscore version handling: 138_2 > 138.1
    _create_release(tmp_path, "138_2")
    best = find_local_silva_database(str(tmp_path))
    assert best == (str(tmp_path / "SILVA" / "138_2"), "138_2")

    # 139 > 138_2
    _create_release(tmp_path, "139")
    best = find_local_silva_database(str(tmp_path))
    assert best == (str(tmp_path / "SILVA" / "139"), "139")


def test_find_local_silva_database_handles_unreadable_file_safely(tmp_path, monkeypatch):
    release_dir = _create_release(tmp_path, "138.2")
    success_file = release_dir / "SILVA.sidx.success"

    # Make reading fail with OSError
    original_open = open

    def failing_open(file, *args, **kwargs):
        if str(file) == str(success_file):
            raise OSError("Unreadable disk sector")
        return original_open(file, *args, **kwargs)

    monkeypatch.setattr("builtins.open", failing_open)
    # Must not raise, just return None
    assert find_local_silva_database(str(tmp_path)) is None


def test_run_tactic_offline_no_db_exits_172(tmp_path, monkeypatch):
    import run_tactic

    def failing_prepare(**kwargs):
        raise DatabaseUnavailableException("Could not reach SILVA server and no local DB")

    monkeypatch.setattr(
        run_tactic.proc_helper,
        "prepare_silva_database",
        failing_prepare,
    )

    with pytest.raises(SystemExit) as exc_info:
        run_tactic.prepare_silva_db(str(tmp_path / "dbs"))

    assert exc_info.value.code == 172


def test_run_tactic_prepare_silva_db_success(tmp_path, monkeypatch):
    import run_tactic

    monkeypatch.setattr(
        run_tactic.proc_helper,
        "prepare_silva_database",
        lambda **kwargs: ("/mock/db/dir", "138.2"),
    )

    db_dir, version = run_tactic.prepare_silva_db(str(tmp_path / "dbs"))
    assert db_dir == "/mock/db/dir"
    assert version == "138.2"
