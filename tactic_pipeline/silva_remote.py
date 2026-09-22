"""Install verified, prebuilt SILVA/SINA files through the download Worker."""

import hashlib
import json
import logging
import os
import platform
import re
import sys
import tempfile
from pathlib import Path
from urllib.parse import urlparse

import requests
import zstandard


DEFAULT_WORKER_URL = "https://silva-downloads-worker.mohsenpm50.workers.dev"
SINA_VERSION = "1.6.0"
SUCCESS_MARKER = "d3fa056627656a96abf4006b9f5d928b"
CHUNK_SIZE = 1024 * 1024
MAX_MANIFEST_BYTES = 65536
DEFAULT_TOKEN_FILE = Path("/base/.silva-app-token")


class SilvaRemoteUnavailable(RuntimeError):
    """Remote service or compatible prebuilt files are unavailable."""

    def __init__(self, message, status_code=None):
        super().__init__(message)
        self.status_code = status_code


class SilvaAuthenticationError(RuntimeError):
    """The Worker rejected the application token or requires Cloudflare Access."""


def _application_token():
    token = os.environ.get("SILVA_APP_TOKEN")
    if token is None:
        try:
            token = DEFAULT_TOKEN_FILE.read_text(encoding="ascii")
        except FileNotFoundError as exc:
            raise ValueError(
                "SILVA application token is missing; set SILVA_APP_TOKEN or use an "
                "image built with the token"
            ) from exc
    token = token.strip()
    if not token:
        raise ValueError("SILVA application token is empty")
    return token


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _expected_hash(manifest, field):
    value = manifest.get(field)
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", value):
        raise ValueError(f"SILVA manifest has no valid {field}")
    return value.lower()


def _native_platform():
    machine = platform.machine().lower()
    architecture = {"x86_64": "amd64", "amd64": "amd64",
                    "aarch64": "arm64", "arm64": "arm64"}.get(machine, machine)
    return f"{platform.system().lower()}-{architecture}"


class SilvaWorkerClient:
    def __init__(self, worker_url):
        parsed = urlparse(worker_url)
        if parsed.scheme != "https" or not parsed.netloc or parsed.path not in ("", "/"):
            raise ValueError("SILVA_WORKER_URL must be an HTTPS Worker origin")
        self.worker_url = worker_url.rstrip("/")
        self._application_token = _application_token()

    def _api(self, route, params=None):
        try:
            response = requests.get(
                self.worker_url + route,
                params=params,
                headers={"Authorization": f"Bearer {self._application_token}"},
                timeout=(10, 30),
                allow_redirects=False,
            )
            if response.status_code in (301, 302, 303, 307, 308, 401, 403):
                raise SilvaAuthenticationError(
                    "SILVA Worker rejected the application token or redirected to "
                    "Cloudflare Access. Check the image token and Worker secret."
                )
            response.raise_for_status()
        except requests.RequestException as exc:
            status = exc.response.status_code if exc.response is not None else None
            detail = f" (HTTP {status})" if status is not None else ""
            raise SilvaRemoteUnavailable(f"SILVA Worker request failed{detail}", status) from None
        try:
            return response.json()
        except ValueError as exc:
            raise ValueError("SILVA Worker returned invalid JSON") from exc

    def releases(self):
        data = self._api("/v1/releases")
        releases = data.get("releases") if isinstance(data, dict) else None
        if not isinstance(releases, list) or any(
            not isinstance(item, str) or not re.fullmatch(r"[0-9]+(?:\.[0-9]+)*", item)
            for item in releases
        ):
            raise ValueError("SILVA Worker returned an invalid release list")
        return sorted(set(releases), key=lambda item: tuple(map(int, item.split("."))), reverse=True)

    def download_info(self, release, filename):
        data = self._api("/v1/download-url", {"release": release, "file": filename})
        if not isinstance(data, dict) or data.get("release") != release or data.get("file") != filename:
            raise ValueError("SILVA Worker returned inconsistent download metadata")
        size = data.get("size_bytes")
        url = data.get("url")
        parsed = urlparse(url) if isinstance(url, str) else None
        if (not isinstance(size, int) or isinstance(size, bool) or size <= 0 or
                parsed is None or parsed.scheme != "https" or
                not parsed.hostname or not parsed.hostname.endswith(".r2.cloudflarestorage.com")):
            raise ValueError("SILVA Worker returned an invalid R2 download URL")
        return url, size

    def manifest(self, release):
        url, size = self.download_info(release, "manifest.json")
        if size > MAX_MANIFEST_BYTES:
            raise ValueError("SILVA manifest is unexpectedly large")
        try:
            with requests.get(url, stream=True, timeout=(10, 120), allow_redirects=False) as response:
                response.raise_for_status()
                content = bytearray()
                for chunk in response.iter_content(CHUNK_SIZE):
                    content.extend(chunk)
                    if len(content) > MAX_MANIFEST_BYTES:
                        raise ValueError("SILVA manifest exceeds size limit")
        except requests.RequestException:
            raise SilvaRemoteUnavailable("SILVA manifest download failed") from None
        if len(content) != size:
            raise ValueError("SILVA manifest size does not match R2 metadata")
        try:
            manifest = json.loads(content)
        except (ValueError, UnicodeDecodeError) as exc:
            raise ValueError("SILVA manifest is invalid JSON") from exc
        if not isinstance(manifest, dict):
            raise ValueError("SILVA manifest must be a JSON object")
        return manifest

    def artifact(self, release, filename, destination, expected_sha256, logger):
        if destination.is_file() and _sha256_file(destination) == expected_sha256:
            logger.info("Verified existing SILVA file %s", destination)
            return
        for attempt in range(2):
            url, size = self.download_info(release, filename)
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(
                    mode="wb", dir=str(destination.parent),
                    prefix=f".{filename}.", suffix=".part", delete=False,
                ) as output:
                    temporary = Path(output.name)
                    digest = hashlib.sha256()
                    bytes_read = 0
                    with requests.get(url, stream=True, timeout=(10, 120), allow_redirects=False) as response:
                        response.raise_for_status()
                        for chunk in response.iter_content(CHUNK_SIZE):
                            if chunk:
                                output.write(chunk)
                                digest.update(chunk)
                                bytes_read += len(chunk)
                    if bytes_read != size or digest.hexdigest() != expected_sha256:
                        raise ValueError(f"SILVA {filename} failed size or SHA-256 verification")
                os.replace(temporary, destination)
                logger.info("Installed verified SILVA file %s", destination)
                return
            except requests.RequestException as exc:
                if attempt == 1:
                    status = exc.response.status_code if exc.response is not None else None
                    raise SilvaRemoteUnavailable(f"SILVA {filename} download failed", status) from None
                logger.warning("SILVA download interrupted; requesting a fresh URL")
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)

    def install_artifact(self, release, filename, destination, manifest, hash_field, logger):
        expected = _expected_hash(manifest, hash_field)
        if destination.is_file() and _sha256_file(destination) == expected:
            logger.info("Verified existing SILVA file %s", destination)
            return
        compressed_field = hash_field.replace("_sha256", "_zst_sha256")
        try:
            self.artifact(release, filename, destination, expected, logger)
            return
        except SilvaRemoteUnavailable as exc:
            # Only a missing raw object calls for the compressed variant.
            # A service outage or interrupted transfer should use the FTP fallback.
            if exc.status_code != 404 or compressed_field not in manifest:
                raise

        compressed_hash = _expected_hash(manifest, compressed_field)
        with tempfile.TemporaryDirectory(dir=str(destination.parent), prefix=".silva-download-") as staging:
            compressed = Path(staging) / (filename + ".zst")
            logger.info("SILVA %s is not available uncompressed; trying .zst", filename)
            self.artifact(release, filename + ".zst", compressed, compressed_hash, logger)

            extracted = Path(staging) / filename
            digest = hashlib.sha256()
            with open(compressed, "rb") as source, open(extracted, "wb") as output:
                with zstandard.ZstdDecompressor().stream_reader(source, read_across_frames=True) as reader:
                    for chunk in iter(lambda: reader.read(CHUNK_SIZE), b""):
                        output.write(chunk)
                        digest.update(chunk)
            if digest.hexdigest() != expected:
                raise ValueError(f"SILVA {filename} failed extracted SHA-256 verification")
            os.replace(extracted, destination)
            logger.info("Installed verified SILVA file %s", destination)


def prepare_prebuilt_silva_database(
    download_dir, version, sina_binary, worker_url,
    logger_obj=logging.getLogger(__name__),
):
    """Select a compatible hosted release and install its ARB and SIDX atomically."""
    client = SilvaWorkerClient(worker_url)
    binary = Path(sina_binary)
    if not binary.is_file():
        raise FileNotFoundError(f"SINA binary is missing: {binary}")
    binary_hash = _sha256_file(binary)
    releases = client.releases()
    requested = version.replace("_", ".") if version != "latest" else version
    if requested != "latest":
        releases = [item for item in releases if item == requested]
    if not releases:
        raise SilvaRemoteUnavailable(f"No prebuilt SILVA release available for {version}")

    selected = None
    for release in releases:
        manifest = client.manifest(release)
        compatible = (
            manifest.get("platform") == _native_platform() and
            manifest.get("byte_order") == sys.byteorder and
            manifest.get("sina_version") == SINA_VERSION and
            isinstance(manifest.get("sina_binary_sha256"), str) and
            manifest["sina_binary_sha256"].lower() == binary_hash
        )
        if compatible:
            selected = (release, manifest)
            break
        logger_obj.info("Skipping incompatible prebuilt SILVA release %s", release)
    if selected is None:
        raise SilvaRemoteUnavailable("No prebuilt SILVA index matches this machine and SINA binary")

    release, manifest = selected
    _expected_hash(manifest, "arb_sha256")
    _expected_hash(manifest, "index_sha256")
    directory = Path(download_dir).absolute() / "SILVA" / release
    directory.mkdir(parents=True, exist_ok=True)
    arb_path = directory / "SILVA.arb"
    index_path = directory / "SILVA.sidx"
    marker_path = directory / "SILVA.sidx.success"
    marker_path.unlink(missing_ok=True)
    client.install_artifact(release, "SILVA.arb", arb_path, manifest, "arb_sha256", logger_obj)
    client.install_artifact(release, "SILVA.sidx", index_path, manifest, "index_sha256", logger_obj)

    # SINA rejects an index older than the ARB, regardless of its checksum.
    if index_path.stat().st_mtime_ns < arb_path.stat().st_mtime_ns:
        new_time = arb_path.stat().st_mtime_ns + 1
        os.utime(index_path, ns=(new_time, new_time))
    temporary_marker = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=str(directory),
            prefix=".SILVA.sidx.success.", suffix=".tmp", delete=False,
        ) as output:
            temporary_marker = Path(output.name)
            output.write(f"Verified prebuilt SILVA {release} index\n{SUCCESS_MARKER}\n")
        os.replace(temporary_marker, marker_path)
    finally:
        if temporary_marker is not None:
            temporary_marker.unlink(missing_ok=True)
    return str(directory), release
