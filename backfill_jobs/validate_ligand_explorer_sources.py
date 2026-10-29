#!/usr/bin/env python3
"""Validate analysis-only 3B backfill inputs without modifying source runs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--eligible-output", required=True, type=Path)
    parser.add_argument("--blocked-output", required=True, type=Path)
    args = parser.parse_args()
    eligible, blocked = [], []
    for entry in json.loads(args.plan.read_text(encoding="utf-8")):
        source = Path(entry["source"])
        events = source / "Analysis_Results" / "InteractionTypes_Framewise.csv"
        chemistry = sorted(source.glob("*.cgenff.mol2"))
        if source.is_dir() and events.is_file() and events.stat().st_size and chemistry:
            eligible.append(entry)
        else:
            blocked.append({
                "cache": entry["cache"],
                "source": str(source),
                "source_directory": source.is_dir(),
                "interaction_events": events.is_file() and bool(events.stat().st_size),
                "cgenff_mol2_count": len(chemistry),
            })
    args.eligible_output.write_text(json.dumps(eligible, indent=2) + "\n", encoding="utf-8")
    args.blocked_output.write_text(json.dumps(blocked, indent=2) + "\n", encoding="utf-8")
    print(f"eligible={len(eligible)} blocked={len(blocked)}")


if __name__ == "__main__":
    main()
