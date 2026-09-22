import hashlib
import json
import logging
import sys

import pytest
import zstandard

from tactic_pipeline import processing_helper, silva_remote


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self.payload = payload
        self.status_code = status_code

    def json(self):
        return self.payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise silva_remote.requests.HTTPError(str(self.status_code), response=self)

    def iter_content(self, size):
        data = self.payload
        for offset in range(0, len(data), size):
            yield data[offset:offset + size]

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _sha(data):
    return hashlib.sha256(data).hexdigest()


@pytest.fixture(autouse=True)
def application_token(monkeypatch):
    monkeypatch.setenv("SILVA_APP_TOKEN", "test-app-token")


def _fake_r2(monkeypatch, binary, *, corrupt_index=False, incompatible_latest=False,
             compressed=False, wrong_raw_hash=False):
    arb = b"prebuilt SILVA ARB"
    index = b"prebuilt SINA index"
    payloads = {"SILVA.arb": arb, "SILVA.sidx": index}
    if compressed:
        compressor = zstandard.ZstdCompressor()
        payloads["SILVA.arb.zst"] = compressor.compress(arb)
        payloads["SILVA.sidx.zst"] = compressor.compress(index)
    downloads = []
    api_calls = []

    def manifest(release):
        data = {
            "platform": "linux-arm64" if incompatible_latest and release == "144"
                        else silva_remote._native_platform(),
            "byte_order": sys.byteorder,
            "sina_version": "1.6.0",
            "sina_binary_sha256": _sha(binary.read_bytes()),
            "arb_sha256": _sha(arb),
            "index_sha256": _sha(b"wrong" if wrong_raw_hash else index),
        }
        if compressed:
            data["arb_zst_sha256"] = _sha(payloads["SILVA.arb.zst"])
            data["index_zst_sha256"] = _sha(payloads["SILVA.sidx.zst"])
        return json.dumps(data).encode()

    def fake_get(url, **kwargs):
        if url.startswith(silva_remote.DEFAULT_WORKER_URL):
            assert kwargs.get("headers") == {"Authorization": "Bearer test-app-token"}
        else:
            assert not kwargs.get("headers"), "R2 downloads must not receive application tokens"
        if url.endswith("/v1/releases"):
            api_calls.append("releases")
            return FakeResponse({"releases": ["138.2", "144"]})
        if url.endswith("/v1/download-url"):
            release = kwargs["params"]["release"]
            filename = kwargs["params"]["file"]
            api_calls.append((release, filename))
            data = manifest(release) if filename == "manifest.json" else payloads[filename]
            return FakeResponse({
                "release": release,
                "file": filename,
                "url": f"https://test.r2.cloudflarestorage.com/{release}/{filename}",
                "size_bytes": len(data),
            })
        if url.startswith("https://test.r2.cloudflarestorage.com/"):
            release, filename = url.split("/")[-2:]
            data = manifest(release) if filename == "manifest.json" else payloads[filename]
            if filename in ("SILVA.sidx", "SILVA.sidx.zst") and corrupt_index:
                data = b"X" + data[1:]
            downloads.append((release, filename))
            return FakeResponse(data)
        raise AssertionError(f"Unexpected request: {url}")

    monkeypatch.setattr(silva_remote.requests, "get", fake_get)
    return downloads, api_calls


@pytest.mark.parametrize("compressed", [False, True])
def test_prebuilt_install_verifies_files_and_skips_local_index(tmp_path, monkeypatch, compressed):
    binary = tmp_path / "sina"
    binary.write_bytes(b"known sina binary")
    downloads, api_calls = _fake_r2(monkeypatch, binary, compressed=compressed)
    directory, release = silva_remote.prepare_prebuilt_silva_database(
        tmp_path, "latest", binary, silva_remote.DEFAULT_WORKER_URL,
        logging.getLogger("silva.test"),
    )
    assert release == "144"
    assert (tmp_path / "SILVA" / "144" / "SILVA.arb").read_bytes() == b"prebuilt SILVA ARB"
    assert (tmp_path / "SILVA" / "144" / "SILVA.sidx").read_bytes() == b"prebuilt SINA index"
    assert (tmp_path / "SILVA" / "144" / "SILVA.sidx.success").exists()
    assert ("144", "SILVA.arb") in api_calls
    assert ("144", "SILVA.sidx") in api_calls
    assert not any(filename.endswith(".zst") for _, filename in downloads)

    def unexpected_local_index(*args, **kwargs):
        raise AssertionError("SINA must not rebuild the verified index")

    monkeypatch.setattr(processing_helper, "system_sub", unexpected_local_index)
    assert processing_helper.index_silva_database(
        str(tmp_path / "SILVA" / "144" / "SILVA.arb"), sina_bin=str(binary),
    ) == str(tmp_path / "SILVA" / "144" / "SILVA.sidx")

    downloads.clear()
    assert silva_remote.prepare_prebuilt_silva_database(
        tmp_path, "latest", binary, silva_remote.DEFAULT_WORKER_URL,
    ) == (directory, release)
    assert downloads == [("144", "manifest.json")]


@pytest.mark.parametrize("compressed", [False, True])
def test_corrupt_download_never_creates_success_marker(tmp_path, monkeypatch, compressed):
    binary = tmp_path / "sina"
    binary.write_bytes(b"known sina binary")
    downloads, _ = _fake_r2(monkeypatch, binary, corrupt_index=True, compressed=compressed)
    with pytest.raises(ValueError, match="SHA-256"):
        silva_remote.prepare_prebuilt_silva_database(
            tmp_path, "latest", binary, silva_remote.DEFAULT_WORKER_URL,
        )
    directory = tmp_path / "SILVA" / "144"
    assert not (directory / "SILVA.sidx").exists()
    assert not (directory / "SILVA.sidx.success").exists()
    assert not list(directory.glob("*.part"))
    assert not list(directory.glob(".silva-download-*"))
    assert not any(filename.endswith(".zst") for _, filename in downloads)


def test_latest_skips_incompatible_hosted_release(tmp_path, monkeypatch):
    binary = tmp_path / "sina"
    binary.write_bytes(b"known sina binary")
    downloads, _ = _fake_r2(monkeypatch, binary, incompatible_latest=True)
    _, release = silva_remote.prepare_prebuilt_silva_database(
        tmp_path, "latest", binary, silva_remote.DEFAULT_WORKER_URL,
    )
    assert release == "138.2"
    assert ("144", "SILVA.arb") not in downloads


def test_pipeline_prepares_database_with_application_token(tmp_path, monkeypatch):
    monkeypatch.delenv("SILVA_ACCESS_TOKEN", raising=False)
    monkeypatch.delenv("SILVA_SOURCE", raising=False)
    monkeypatch.delenv("SILVA_WORKER_URL", raising=False)
    binary = tmp_path / "sina"
    binary.write_bytes(b"known sina binary")
    _fake_r2(monkeypatch, binary, compressed=True)
    assert processing_helper.prepare_silva_database(
        dest_dir=str(tmp_path), sina_binary=str(binary),
    ) == (str(tmp_path / "SILVA" / "144"), "144")


@pytest.mark.parametrize("status", [302, 401, 403])
def test_protected_worker_explains_required_dashboard_change(monkeypatch, status):
    monkeypatch.setattr(
        silva_remote.requests, "get",
        lambda *args, **kwargs: FakeResponse({}, status_code=status),
    )
    client = silva_remote.SilvaWorkerClient(silva_remote.DEFAULT_WORKER_URL)
    with pytest.raises(silva_remote.SilvaAuthenticationError, match="Cloudflare Access"):
        client.releases()


def test_application_token_can_be_loaded_from_image_file(tmp_path, monkeypatch):
    monkeypatch.delenv("SILVA_APP_TOKEN")
    token_file = tmp_path / "app-token"
    token_file.write_text("file-token\n", encoding="ascii")
    monkeypatch.setattr(silva_remote, "DEFAULT_TOKEN_FILE", token_file)
    seen_headers = []

    def fake_get(url, **kwargs):
        seen_headers.append(kwargs["headers"])
        return FakeResponse({"releases": []})

    monkeypatch.setattr(silva_remote.requests, "get", fake_get)
    client = silva_remote.SilvaWorkerClient(silva_remote.DEFAULT_WORKER_URL)
    assert client.releases() == []
    assert seen_headers == [{"Authorization": "Bearer file-token"}]


def test_missing_application_token_is_not_silent_ftp_fallback(tmp_path, monkeypatch):
    monkeypatch.delenv("SILVA_APP_TOKEN")
    monkeypatch.setattr(silva_remote, "DEFAULT_TOKEN_FILE", tmp_path / "missing")

    def unexpected_ftp(**kwargs):
        raise AssertionError("Missing application token must not trigger local indexing")

    monkeypatch.setattr(processing_helper, "download_silva_databases", unexpected_ftp)
    with pytest.raises(ValueError, match="application token is missing"):
        processing_helper.prepare_silva_database(dest_dir=str(tmp_path))


def test_rejected_application_token_is_not_silent_ftp_fallback(tmp_path, monkeypatch):
    binary = tmp_path / "sina"
    binary.write_bytes(b"known sina binary")
    monkeypatch.setattr(
        silva_remote.requests, "get",
        lambda *args, **kwargs: FakeResponse({}, status_code=401),
    )

    def unexpected_ftp(**kwargs):
        raise AssertionError("Rejected application token must not trigger local indexing")

    monkeypatch.setattr(processing_helper, "download_silva_databases", unexpected_ftp)
    with pytest.raises(silva_remote.SilvaAuthenticationError, match="application token"):
        processing_helper.prepare_silva_database(
            dest_dir=str(tmp_path), sina_binary=str(binary),
        )


def test_ftp_fallback_is_explicit(tmp_path, monkeypatch):
    monkeypatch.setenv("SILVA_SOURCE", "ftp")
    monkeypatch.delenv("SILVA_ACCESS_TOKEN", raising=False)
    arb_path = tmp_path / "SILVA" / "144" / "SILVA.arb"
    calls = []

    def fake_download(**kwargs):
        calls.append("download")
        return str(arb_path), "144"

    def fake_index(**kwargs):
        calls.append("index")

    monkeypatch.setattr(processing_helper, "download_silva_databases", fake_download)
    monkeypatch.setattr(processing_helper, "index_silva_database", fake_index)
    assert processing_helper.prepare_silva_database(
        dest_dir=str(tmp_path), sina_binary="sina-test",
    ) == (str(arb_path.parent), "144")
    assert calls == ["download", "index"]


def test_extracted_checksum_is_verified_before_index_install(tmp_path, monkeypatch):
    binary = tmp_path / "sina"
    binary.write_bytes(b"known sina binary")
    downloads, api_calls = _fake_r2(
        monkeypatch, binary, compressed=True, wrong_raw_hash=True,
    )
    original_get = silva_remote.requests.get

    def missing_raw_index(url, **kwargs):
        if (url.endswith("/v1/download-url") and
                kwargs["params"]["file"] == "SILVA.sidx"):
            return FakeResponse({}, status_code=404)
        return original_get(url, **kwargs)

    monkeypatch.setattr(silva_remote.requests, "get", missing_raw_index)
    with pytest.raises(ValueError, match="extracted SHA-256"):
        silva_remote.prepare_prebuilt_silva_database(
            tmp_path, "latest", binary, silva_remote.DEFAULT_WORKER_URL,
        )
    directory = tmp_path / "SILVA" / "144"
    assert not (directory / "SILVA.sidx").exists()
    assert not (directory / "SILVA.sidx.success").exists()
    assert not list(directory.glob(".silva-download-*"))
    assert ("144", "SILVA.sidx.zst") in api_calls
    assert ("144", "SILVA.sidx.zst") in downloads


def test_compressed_fallback_when_uncompressed_objects_are_missing(tmp_path, monkeypatch):
    binary = tmp_path / "sina"
    binary.write_bytes(b"known sina binary")
    downloads, api_calls = _fake_r2(monkeypatch, binary, compressed=True)
    original_get = silva_remote.requests.get
    missing_requests = []

    def missing_uncompressed(url, **kwargs):
        if (url.endswith("/v1/download-url") and
                kwargs["params"]["file"] in ("SILVA.arb", "SILVA.sidx")):
            missing_requests.append(kwargs["params"]["file"])
            return FakeResponse({}, status_code=404)
        return original_get(url, **kwargs)

    monkeypatch.setattr(silva_remote.requests, "get", missing_uncompressed)
    silva_remote.prepare_prebuilt_silva_database(
        tmp_path, "latest", binary, silva_remote.DEFAULT_WORKER_URL,
    )
    assert missing_requests == ["SILVA.arb", "SILVA.sidx"]
    assert ("144", "SILVA.arb.zst") in api_calls
    assert ("144", "SILVA.sidx.zst") in api_calls
    assert ("144", "SILVA.arb.zst") in downloads
    assert ("144", "SILVA.sidx.zst") in downloads
    assert ("144", "SILVA.arb") not in downloads
    assert ("144", "SILVA.sidx") not in downloads
    directory = tmp_path / "SILVA" / "144"
    assert (directory / "SILVA.arb").read_bytes() == b"prebuilt SILVA ARB"
    assert (directory / "SILVA.sidx").read_bytes() == b"prebuilt SINA index"
    assert (directory / "SILVA.sidx.success").exists()


@pytest.mark.parametrize("failure", [
    "worker_timeout", "worker_503", "manifest_missing", "index_missing",
    "r2_interrupted", "no_releases", "incompatible",
])
def test_pipeline_automatically_falls_back_to_ftp(tmp_path, monkeypatch, caplog, failure):
    monkeypatch.delenv("SILVA_SOURCE", raising=False)
    binary = tmp_path / "sina"
    binary.write_bytes(b"known sina binary")
    _fake_r2(monkeypatch, binary, compressed=True)
    normal_get = silva_remote.requests.get
    interrupted_downloads = []

    class InterruptedResponse(FakeResponse):
        def iter_content(self, size):
            yield b"partial index"
            raise silva_remote.requests.ConnectionError("broken stream X-Amz-Signature=secret")

    def unavailable(url, **kwargs):
        if url.endswith("/v1/releases"):
            if failure == "worker_timeout":
                raise silva_remote.requests.Timeout("Worker did not respond")
            if failure == "worker_503":
                return FakeResponse({}, status_code=503)
            if failure == "no_releases":
                return FakeResponse({"releases": []})
        if url.endswith("/v1/download-url"):
            filename = kwargs["params"]["file"]
            if failure == "manifest_missing" and filename == "manifest.json":
                return FakeResponse({}, status_code=404)
            if failure == "index_missing" and filename.startswith("SILVA.sidx"):
                return FakeResponse({}, status_code=404)
        if failure == "r2_interrupted" and url.endswith("/SILVA.sidx"):
            interrupted_downloads.append(url)
            return InterruptedResponse(b"")
        return normal_get(url, **kwargs)

    monkeypatch.setattr(silva_remote.requests, "get", unavailable)
    if failure == "incompatible":
        monkeypatch.setattr(silva_remote, "SINA_VERSION", "unsupported-version")

    ftp_path = tmp_path / "SILVA" / "138.2" / "SILVA.arb.gz"
    calls = []

    def ftp_download(**kwargs):
        calls.append(("ftp", kwargs))
        return str(ftp_path), "138.2"

    def local_index(**kwargs):
        calls.append(("index", kwargs))

    monkeypatch.setattr(processing_helper, "download_silva_databases", ftp_download)
    monkeypatch.setattr(processing_helper, "index_silva_database", local_index)
    logger = logging.getLogger("silva.fallback.test")
    with caplog.at_level(logging.WARNING, logger=logger.name):
        result = processing_helper.prepare_silva_database(
            version="138_2", md5_check=False, dest_dir=str(tmp_path),
            sina_binary=str(binary), sina_index_threads=2, logger_obj=logger,
        )
    assert result == (str(ftp_path.parent), "138.2")
    assert [name for name, _ in calls] == ["ftp", "index"]
    assert calls[0][1] == {
        "download_dir": str(tmp_path), "version": "138_2",
        "md5_check": False, "expected_file_basename": "SILVA", "logger_obj": logger,
    }
    assert calls[1][1] == {
        "silva_arb_file": str(ftp_path), "logger_obj": logger,
        "sina_bin": str(binary), "threads": 2,
    }
    assert "Falling back to SILVA FTP" in caplog.text
    assert "substantial RAM" in caplog.text
    assert "X-Amz-Signature" not in caplog.text
    assert not list(tmp_path.rglob("*.part"))
    assert not list(tmp_path.rglob(".silva-download-*"))
    if failure == "r2_interrupted":
        assert len(interrupted_downloads) == 2


def test_pipeline_does_not_fallback_on_corrupt_index(tmp_path, monkeypatch):
    monkeypatch.delenv("SILVA_SOURCE", raising=False)
    binary = tmp_path / "sina"
    binary.write_bytes(b"known sina binary")
    _fake_r2(monkeypatch, binary, compressed=True, corrupt_index=True)

    def unexpected_ftp(**kwargs):
        raise AssertionError("An integrity failure must not trigger FTP fallback")

    monkeypatch.setattr(processing_helper, "download_silva_databases", unexpected_ftp)
    with pytest.raises(ValueError, match="SHA-256"):
        processing_helper.prepare_silva_database(
            dest_dir=str(tmp_path), sina_binary=str(binary),
        )


def test_pipeline_does_not_fallback_on_local_io_failure(tmp_path, monkeypatch):
    monkeypatch.delenv("SILVA_SOURCE", raising=False)

    def disk_full(**kwargs):
        raise OSError("No space left on device")

    def unexpected_ftp(**kwargs):
        raise AssertionError("A local filesystem failure must not trigger FTP fallback")

    monkeypatch.setattr(processing_helper, "prepare_prebuilt_silva_database", disk_full)
    monkeypatch.setattr(processing_helper, "download_silva_databases", unexpected_ftp)
    with pytest.raises(OSError, match="No space left"):
        processing_helper.prepare_silva_database(dest_dir=str(tmp_path))
