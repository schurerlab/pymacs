#!/usr/bin/env python3
"""Publish completed ligand dossier HTML files from Triton into Kyle's cache.

This is intentionally transport-only: it never invokes PyMACS, modifies a
simulation, or replaces a cache file until the incoming HTML validates.
"""

from __future__ import annotations

import json
import shlex
import subprocess
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PLAN = Path("/home/jxs794/ligand-backfill-20261029-plan.json")
CACHE_ROOT = Path("/home/jxs794/lab-job-portal/jobs/published-md-results")
REPORT = Path("/home/jxs794/ligand-backfill-20261029-cache-report.json")
RETRIES = 4


def triton_stream(command: str, *, stdout=None) -> subprocess.CompletedProcess:
    """Use the verified Kyle -> token -> Triton fallback route."""
    return subprocess.run(
        ["ssh", "token", f"ssh triton {shlex.quote(command)}"],
        stdout=stdout,
        stderr=subprocess.PIPE,
        timeout=180,
    )


def publish_one(item: dict) -> list[dict]:
    source_html = item["source"] + "/Analysis_Results/PyMACS_Result_Dashboard.html"
    outcomes: list[dict] = []
    for relative_cache_dir in item.get("cache_directories", []):
        destination = CACHE_ROOT / relative_cache_dir / "Analysis_Results"
        destination.mkdir(parents=True, exist_ok=True)
        partial = destination / "PyMACS_Result_Dashboard.html.backfill-partial"
        final = destination / "PyMACS_Result_Dashboard.html"
        result = None
        for attempt in range(1, RETRIES + 1):
            try:
                with partial.open("wb") as output:
                    result = triton_stream(f"cat {shlex.quote(source_html)}", stdout=output)
            except (subprocess.TimeoutExpired, OSError) as exc:
                detail = str(exc)
            else:
                detail = result.stderr.decode(errors="replace")[-400:]
                if result.returncode == 0:
                    body = partial.read_bytes()
                    if len(body) >= 100_000 and b"PyMACS" in body and b"Strongest" in body:
                        partial.replace(final)
                        outcomes.append({
                            "job_id": item["job_id"], "source": item["source"],
                            "cache": str(destination), "result": "cached", "size": len(body),
                        })
                        break
                    detail = f"received invalid dashboard ({len(body)} bytes)"
            partial.unlink(missing_ok=True)
            if attempt < RETRIES:
                time.sleep(attempt * 3)
        else:
            outcomes.append({
                "job_id": item["job_id"], "source": item["source"],
                "cache": str(destination), "result": "unavailable", "detail": detail,
            })
    return outcomes


def main() -> None:
    items = json.loads(PLAN.read_text())["submitted"]
    report: list[dict] = []
    # Keep fallback traffic conservative; retries make transient Acorn failures safe.
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(publish_one, item) for item in items]
        for index, future in enumerate(as_completed(futures), start=1):
            report.extend(future.result())
            print(f"published check {index}/{len(items)}", flush=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(Counter(record["result"] for record in report), flush=True)


if __name__ == "__main__":
    main()
