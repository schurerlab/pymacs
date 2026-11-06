#!/usr/bin/env python3
"""Regenerate every published protein--protein PyMACS dossier on a portal host.

Run this on Kyle after deploying ``5_PyMACS_Result_Dashboard.py``.  The portal
index is the source of truth for published directories and recorded MD rate, so
historical systems receive the same scorecard provenance as newly analysed jobs.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, required=True, help="Portal md_results_index.json")
    parser.add_argument("--published-root", type=Path, required=True, help="Portal published-md-results directory")
    parser.add_argument("--dashboard", type=Path, required=True, help="5_PyMACS_Result_Dashboard.py")
    parser.add_argument("--skip", action="append", default=[], help="System identifier to skip; repeatable")
    parser.add_argument("--only", action="append", default=[], help="System identifier to render; repeatable")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def is_ppi(row: dict) -> bool:
    return row.get("kind") == "md" and row.get("molecule_type") == "protein_protein" and row.get("analysis_ready") is True


def finite_rate(row: dict) -> float | None:
    value = row.get("rate_ns_day")
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) and value > 0 else None


def main() -> int:
    args = parse_args()
    skip, only = set(args.skip), set(args.only)
    systems = json.loads(args.index.read_text())["systems"]
    selected = [row for row in systems if is_ppi(row)]
    if only:
        selected = [row for row in selected if row["system"] in only]
    selected = [row for row in selected if row["system"] not in skip]
    selected.sort(key=lambda row: row["system"].lower())
    print(f"Selected {len(selected)} published protein--protein dossier(s).", flush=True)

    failures: list[str] = []
    for number, row in enumerate(selected, start=1):
        system = row["system"]
        root = args.published_root / row["directory"]
        output = root / "Analysis_Results" / "PyMACS_Result_Dashboard.html"
        command = [sys.executable, str(args.dashboard), "--root", "."]
        rate = finite_rate(row)
        if rate is not None:
            command.extend(["--rate-ns-day", f"{rate:.6f}"])
        print(f"[{number}/{len(selected)}] {system}: {' '.join(command)}", flush=True)
        if args.dry_run:
            continue
        if not root.is_dir():
            failures.append(f"{system}: missing published directory {root}")
            continue
        completed = subprocess.run(command, cwd=root, text=True)
        if completed.returncode != 0 or not output.is_file() or output.stat().st_size < 100_000:
            failures.append(f"{system}: renderer failed or output is incomplete")
        else:
            print(f"[{number}/{len(selected)}] {system}: complete ({output.stat().st_size:,} bytes)", flush=True)

    if failures:
        print("Failures:\n" + "\n".join(failures), file=sys.stderr)
        return 1
    print("All selected PPI dossiers completed successfully.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
