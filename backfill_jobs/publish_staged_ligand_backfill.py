#!/usr/bin/env python3
"""Publish staged, read-only-source ligand backfills after their LSF jobs finish."""

from __future__ import annotations

import json
import shlex
import subprocess
import time
from pathlib import Path

PLAN = Path("/home/jxs794/ligand-backfill-20261029-plan.json")
CACHE_ROOT = Path("/home/jxs794/lab-job-portal/jobs/published-md-results")
REPORT = Path("/home/jxs794/ligand-backfill-20261029-staged-cache-report.json")
STAGE_ROOT = "/projectnb/triton/home/jxs794/Pymacs/backfill_staging/ligands_20261029"
STAGE_BY_SOURCE = {
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/LEDGF_9g94/rxn513b____EMOL299971444____EMOL719228____EMOL978659": "1d8ae69e3f",
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/LEDGF_9g94/Grubbs____45907639____96488367": "2362494eeb",
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/LEDGF_9g94/rxn110b____EMOL45477804____EMOL24630141": "240a57bb6d",
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/LEDGF_9g94/WXVL022____HU0022____AG0572": "27454898d5",
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/LEDGF_9g94/m_207bcb____26285856____10860938": "7e1874c0ea",
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/Nahuel/MOUSE_UMF": "a29dc1228b",
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/LEDGF_9g94/rxn512a____EMOL56771003____EMOL37531308____EMOL48572214": "a8ef95d0d2",
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/LEDGF_9g94/rxn512b____EMOL498026____EMOL93404083____EMOL888772": "aa2d6f6666",
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/LEDGF_9g94/Reamb_c009be____Amb16233958____Amb35863514____Amb33542963": "ac9f9f6f63",
    "/projectnb/cheminfo/home/jxs794/MD-Backups/KYLE_MD_JOBS/LEDGF_9g94/rxn512b____EMOL56771003____EMOL828220____EMOL888772": "e2033bce3b",
}


def remote_cat(source: str, target: Path) -> bool:
    for attempt in range(4):
        try:
            with target.open("wb") as output:
                completed = subprocess.run(
                    ["ssh", "token", f"ssh triton {shlex.quote('cat ' + shlex.quote(source))}"],
                    stdout=output,
                    stderr=subprocess.PIPE,
                    timeout=180,
                )
            body = target.read_bytes()
            if completed.returncode == 0 and len(body) >= 100_000 and b"PyMACS" in body and b"Strongest" in body:
                return True
        except (OSError, subprocess.TimeoutExpired):
            pass
        target.unlink(missing_ok=True)
        time.sleep(3 * (attempt + 1))
    return False


def main() -> None:
    cache_by_source = {x["source"]: x.get("cache_directories", []) for x in json.loads(PLAN.read_text())["submitted"]}
    pending = {
        source: (f"{STAGE_ROOT}/{slug}", cache_by_source[source])
        for source, slug in STAGE_BY_SOURCE.items()
    }
    report = []
    # The jobs are short, but a 2-hour window survives normal queue delays.
    for _ in range(60):
        for source, (stage, destinations) in list(pending.items()):
            remote_html = stage + "/Analysis_Results/PyMACS_Result_Dashboard.html"
            remaining = []
            for relative in destinations:
                directory = CACHE_ROOT / relative / "Analysis_Results"
                directory.mkdir(parents=True, exist_ok=True)
                partial = directory / "PyMACS_Result_Dashboard.html.backfill-partial"
                final = directory / "PyMACS_Result_Dashboard.html"
                if remote_cat(remote_html, partial):
                    partial.replace(final)
                    report.append({"source": source, "stage": stage, "cache": str(directory), "result": "cached"})
                else:
                    remaining.append(relative)
            if not remaining:
                del pending[source]
        if not pending:
            break
        time.sleep(120)
    for source, (stage, destinations) in pending.items():
        report.append({"source": source, "stage": stage, "cache_directories": destinations, "result": "unavailable"})
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"cached={sum(x['result'] == 'cached' for x in report)} unavailable={sum(x['result'] == 'unavailable' for x in report)}", flush=True)


if __name__ == "__main__":
    main()
