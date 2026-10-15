import os
import re
import subprocess
from pathlib import Path

DEFAULT_VERSION = "v1.0.0"


def _format_v_tag(tag: str) -> str:
    """
    Ensures the tag conforms to 'vX.X.X' format while preserving
    any prerelease suffixes (e.g. -alpha-1, -beta.2).
    """
    tag = str(tag).strip()
    if not tag or tag in ("NA", "None"):
        return DEFAULT_VERSION
    return tag if tag.startswith("v") else f"v{tag}"


def get_version() -> str:
    """
    Resolves the pipeline version as 'vX.X.X'.
    Priority:
      1. os.environ['GIT_TAG'], os.environ['TAG_GIT'], or os.environ['TACTIC_VERSION'] (injected by GitHub Actions/Docker)
      2. tactic_pipeline/VERSION file (if generated during build)
      3. git describe --tags (when running in a local cloned git repo)
      4. DEFAULT_VERSION fallback ('v1.0.0')
    """
    # 1. Environment variable (e.g. in Docker container or CI)
    for var in ("GIT_TAG", "TAG_GIT", "TACTIC_VERSION"):
        env_tag = os.environ.get(var)
        if env_tag and env_tag.strip() not in ("NA", "None", ""):
            return _format_v_tag(env_tag)

    # 2. Check for a VERSION file in the package directory
    version_file = Path(__file__).resolve().parent / "VERSION"
    if version_file.is_file():
        try:
            content = version_file.read_text(encoding="utf-8").strip()
            if content and content not in ("NA", "None", ""):
                return _format_v_tag(content)
        except Exception:
            pass

    # 3. Check git repository if .git directory exists
    repo_root = Path(__file__).resolve().parent.parent
    if (repo_root / ".git").exists():
        try:
            res = subprocess.run(
                ["git", "describe", "--tags", "--match", "v[0-9]*", "--abbrev=0"],
                cwd=str(repo_root),
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                check=False,
            )
            if res.returncode == 0 and res.stdout.strip():
                return _format_v_tag(res.stdout)
        except Exception:
            pass

    # 4. Fallback default
    return _format_v_tag(DEFAULT_VERSION)


__version__ = get_version()

