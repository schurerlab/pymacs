#!/usr/bin/env python3
"""Safely resume explicitly opted-in legacy PyMACS MD runs after LSF run limits.

This program is intentionally conservative.  It never scans arbitrary MD
folders or restarts a job merely because a log is old.  A project must contain
``.portal-auto-resume.json`` with ``enabled: true`` and must have both a
GROMACS checkpoint and a confirmed LSF timeout signal.  Every restart is
recorded locally in that marker before the command returns.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
from pathlib import Path


ACTIVE = {"PEND", "RUN", "PSUSP", "USUSP", "SSUSP"}
JOB_ID = re.compile(r"<(\d+)>")
TIMEOUT_MARKER = "User defined signal 2"
LOCK_NAME = ".portal-autoresume.lock"
MARKER_NAME = ".portal-auto-resume.json"


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def scheduler_jobs() -> list[tuple[str, str, str]]:
    result = subprocess.run(
        ["bjobs", "-a", "-noheader", "-o", "jobid stat job_name"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=30,
    )
    rows: list[tuple[str, str, str]] = []
    for line in result.stdout.splitlines():
        fields = line.split(maxsplit=2)
        if len(fields) >= 2:
            rows.append((fields[0], fields[1], fields[2] if len(fields) == 3 else ""))
    return rows


def active_md_job(project: str, jobs: list[tuple[str, str, str]]) -> bool:
    prefixes = (f"{project}_md", f"{project}_resume")
    return any(state in ACTIVE and job_name.startswith(prefixes) for _, state, job_name in jobs)


def timeout_evidence(project: Path) -> bool:
    # Consider only the most recent production-MD stderr.  An older run-limit
    # message must not mask a newer scientific/setup failure, and an analysis
    # or report stderr must not control production recovery.
    errors = sorted(project.glob(f"{project.name}_md.*.err"), key=lambda path: path.stat().st_mtime, reverse=True)
    if not errors:
        return False
    try:
        return TIMEOUT_MARKER in errors[0].read_text(errors="replace")[-20000:]
    except OSError:
        return False


def parse_job_id(output: str) -> str:
    match = JOB_ID.search(output)
    if not match:
        raise RuntimeError(f"bsub did not return a job id: {output.strip()}")
    return match.group(1)


def submit_from_lsf(path: Path, extra: list[str]) -> str:
    with path.open("rb") as handle:
        result = subprocess.run(["bsub", *extra], stdin=handle, text=False, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, check=False, timeout=60)
    if result.returncode:
        raise RuntimeError(result.stderr.decode(errors="replace").strip() or "bsub failed")
    return parse_job_id(result.stdout.decode(errors="replace"))


def submit_report(project: Path, name: str, analysis_id: str) -> str:
    command = " && ".join((
        "source /etc/bashrc",
        "set -euo pipefail",
        f"cd {project}",
        "module purge",
        "module load gcc-runtime/13.4.0-none-none-gh6rpjx",
        "module load gromacs/2025.1-gcc-13.4.0-6vq7xfo",
        "source /projectnb/triton/home/jxs794/miniforge3/etc/profile.d/conda.sh",
        "conda activate mdanalysis",
        "python /projectnb/triton/home/jxs794/Pymacs/5_PyMACS_Result_Dashboard.py --root . --out Analysis_Results/PyMACS_Result_Dashboard.html",
    ))
    result = subprocess.run(
        ["bsub", "-w", f"done({analysis_id})", "-J", f"{name}_report", "-q", "normal", "-P", "brd",
         "-W", "2:00", "-n", "4", "-R", "span[hosts=1]", "-oo", str(project / f"{name}_report.%J.out"),
         "-eo", str(project / f"{name}_report.%J.err"), "bash", "-lc", command],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, timeout=60,
    )
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "report bsub failed")
    return parse_job_id(result.stdout)


def process(root: Path, execute: bool) -> list[dict[str, object]]:
    jobs = scheduler_jobs()
    outcomes: list[dict[str, object]] = []
    for project in sorted(path for path in root.iterdir() if path.is_dir()):
        marker_path = project / MARKER_NAME
        try:
            marker = json.loads(marker_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(marker, dict) or marker.get("enabled") is not True:
            continue
        name = project.name
        state: dict[str, object] = {"project": name}
        if (project / "md_0_1.log").is_file() and "Finished mdrun" in (project / "md_0_1.log").read_text(errors="replace"):
            state["status"] = "complete"
        elif active_md_job(name, jobs):
            state["status"] = "active"
        elif not ((project / "md_0_1.tpr").is_file() and (project / "md_0_1.cpt").is_file()):
            state["status"] = "missing-checkpoint"
        elif not timeout_evidence(project):
            state["status"] = "not-a-confirmed-timeout"
        else:
            attempts = int(marker.get("restart_count", 0))
            max_attempts = int(marker.get("max_restarts", 8))
            if attempts >= max_attempts:
                state["status"] = "restart-limit-reached"
            elif not execute:
                state["status"] = "would-resume"
            else:
                lock = project / LOCK_NAME
                try:
                    os.mkdir(lock)
                except FileExistsError:
                    state["status"] = "locked"
                else:
                    try:
                        md = submit_from_lsf(project / f"run_{name}_md.lsf", ["-J", f"{name}_resume{attempts + 1}"])
                        analysis = submit_from_lsf(project / f"run_{name}_analysis.lsf", ["-w", f"done({md})", "-J", f"{name}_analysis{attempts + 1}"])
                        report = submit_report(project, name, analysis)
                        marker["restart_count"] = attempts + 1
                        marker.setdefault("history", []).append({"at": now(), "md": md, "analysis": analysis, "report": report})
                        marker_path.write_text(json.dumps(marker, indent=2) + "\n")
                        state.update({"status": "resubmitted", "md": md, "analysis": analysis, "report": report})
                    except Exception as error:
                        state.update({"status": "error", "error": str(error)})
                    finally:
                        lock.rmdir()
        outcomes.append(state)
    return outcomes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--execute", action="store_true", help="Submit eligible resume jobs; default is dry-run.")
    args = parser.parse_args()
    print(json.dumps({"checked_at": now(), "execute": args.execute, "projects": process(args.root, args.execute)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
