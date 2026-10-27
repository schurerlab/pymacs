"""Durable, job-scoped progress reporting for PyMACS analysis pipelines.

The portal reads this small JSON document during its scheduler poll.  It is
written atomically so an interrupted write can never make a live job look
failed or complete.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Iterable


STANDARD_ANALYSIS_STAGES = (
    "A_B", "C", "D1", "D2", "D3", "D4", "D4B", "D5", "D6", "D6B",
    "D6C", "E_INTERFACE", "E1", "E1Q", "E1QP", "E1A", "E1B", "E2",
    "E3", "E4", "E4B", "E5", "E6", "E7", "E7B", "E8", "E9", "E10",
    "E11", "D6D", "FINALIZE",
)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def standard_stage_id(label: str) -> str | None:
    """Normalize the established ``STEP D1`` / ``STEP E4B`` labels."""
    upper = str(label).upper()
    if "STEP A+B" in upper:
        return "A_B"
    if "STEP E-INTERFACE" in upper:
        return "E_INTERFACE"
    match = re.search(r"STEP\s+(D6D|D6C|D6B|D4B|E1QP|E1Q|E1A|E1B|E4B|E7B|D[1-6]|E(?:[1-9]|1[01])|C)\b", upper)
    return match.group(1) if match else None


class ProgressReporter:
    """Maintain a single analysis-attempt manifest below ``Analysis_Results``."""

    def __init__(self, outdir: str | Path, *, pipeline: str, stages: Iterable[str]):
        self.path = Path(outdir) / "portal_progress.json"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.job_id = str(os.environ.get("LSB_JOBID") or os.environ.get("PYMACS_ANALYSIS_JOB_ID") or "local")
        self.data = {
            "schema_version": 1,
            "job_id": self.job_id,
            "attempt_id": f"{self.job_id}-{utc_now().replace(':', '').replace('+00:00', 'Z')}",
            "pipeline": pipeline,
            "state": "running",
            "updated_at": utc_now(),
            "current_stage": None,
            "stages": {stage: {"status": "pending"} for stage in dict.fromkeys(stages)},
        }
        self._write()

    def _write(self) -> None:
        stages = self.data["stages"]
        resolved = sum(1 for item in stages.values() if item.get("status") in {"complete", "skipped", "reused"})
        total = len(stages)
        self.data.update({
            "updated_at": utc_now(), "completed": resolved, "total": total,
            "percent": round((resolved / total * 100) if total else 0.0, 1),
        })
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=self.path.parent,
                                         prefix=".portal_progress.", suffix=".tmp", delete=False) as handle:
            json.dump(self.data, handle, sort_keys=True, indent=2)
            handle.write("\n")
            temporary = Path(handle.name)
        os.replace(temporary, self.path)

    def start(self, stage: str | None, label: str | None = None) -> None:
        if not stage:
            return
        stages = self.data["stages"]
        if stage not in stages:
            stages[stage] = {"status": "pending"}
        previous = self.data.get("current_stage")
        if previous and previous != stage and stages.get(previous, {}).get("status") == "running":
            stages[previous].update({"status": "complete", "completed_at": utc_now(), "inferred_from_next_stage": True})
        stages[stage].update({"status": "running", "started_at": utc_now(), "label": label or stage})
        self.data["current_stage"] = stage
        self._write()

    def complete(self, stage: str | None) -> None:
        if not stage:
            return
        stages = self.data["stages"]
        stages.setdefault(stage, {})
        stages[stage].update({"status": "complete", "completed_at": utc_now()})
        if self.data.get("current_stage") == stage:
            self.data["current_stage"] = None
        self._write()

    def frame_progress(self, stage: str, done: int, total: int) -> None:
        self.start(stage)
        item = self.data["stages"][stage]
        item.update({"frames_done": int(done), "frames_total": int(total)})
        self._write()

    def fail(self, error: str) -> None:
        self.data["state"] = "failed"
        self.data["error"] = str(error)[-1000:]
        self._write()

    def finish(self) -> None:
        for stage, item in self.data["stages"].items():
            if item.get("status") == "pending":
                item.update({"status": "skipped", "completed_at": utc_now()})
        current = self.data.get("current_stage")
        if current and self.data["stages"].get(current, {}).get("status") == "running":
            self.complete(current)
        self.data["state"] = "complete"
        self.data["current_stage"] = None
        self._write()
