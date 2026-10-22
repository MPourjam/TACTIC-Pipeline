"""Test the portable commands from a relocated copy of the small collection."""
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--full-datasets", nargs="*", choices=["silva", "ltp", "gg2"], default=[])
    args = parser.parse_args()
    workdir = args.workdir.resolve()
    if workdir.exists():
        raise RuntimeError("Use a new --workdir for an isolated reproduction test")
    relocated = workdir / "collection"
    relocated.mkdir(parents=True)
    for name in ("code", "config", "tools", "scripts", "environment", "metadata"):
        shutil.copytree(ROOT / name, relocated / name, ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(ROOT / "results/archive_run", relocated / "results/archive_run")
    shutil.copy2(ROOT / "README.md", relocated / "README.md")
    # The runtime sees only the selected Conda environment's programs and OS utilities.
    environment = dict(os.environ)
    environment["PATH"] = str(Path(sys.executable).parent) + ":/usr/bin:/bin"
    environment.pop("PYTHONPATH", None)
    environment["PYTHONNOUSERSITE"] = "1"
    report = {"relocated_collection": str(relocated), "commands": [], "passed": True}
    jobs = [
        ("check", ["check", "--outdir", str(workdir / "check")]),
        ("plots", ["plots", "--outdir", str(workdir / "plots")]),
        ("smoke", ["reproduce", "--data-dir", str(args.data_dir.resolve()),
                   "--outdir", str(workdir / "smoke"), "--smoke-records", "64"]),
    ]
    if args.full_datasets:
        jobs.append(("full", ["reproduce", "--data-dir", str(args.data_dir.resolve()),
                              "--outdir", str(workdir / "full"), "--datasets", *args.full_datasets,
                              "--threads", "19", "--workers", "3"]))
    for name, arguments in jobs:
        command = [sys.executable, str(relocated / "scripts/run_stepwise.py"), *arguments]
        print("Running " + name, flush=True)
        with (workdir / (name + ".log")).open("w") as log:
            result = subprocess.run(command, cwd=relocated, env=environment, stdout=log, stderr=subprocess.STDOUT)
        report["commands"].append({"name": name, "command": command, "returncode": result.returncode})
        report["passed"] = report["passed"] and result.returncode == 0
        (workdir / "test_report.json").write_text(json.dumps(report, indent=2) + "\n")
        if result.returncode:
            print(f"Test failed: {name}; see {workdir / (name + '.log')}")
            return 1
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
