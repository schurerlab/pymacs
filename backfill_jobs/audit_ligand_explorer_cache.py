#!/usr/bin/env python3
"""List ligand dashboards that have typed events but lack a NETWORX explorer."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    for dashboard in args.cache_root.glob("**/Analysis_Results/PyMACS_Result_Dashboard.html"):
        try:
            html = dashboard.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        analysis = dashboard.parent
        ligand_mode = "Ligand / PROTAC:" in html or "Ligand interaction type over time" in html
        events = analysis / "InteractionTypes_Framewise.csv"
        explorer = list((analysis / "NETWORX").glob("*_interactive.html"))
        if ligand_mode and events.is_file() and events.stat().st_size and not explorer:
            rows.append({
                "directory": str(analysis.parent.relative_to(args.cache_root)),
                "events_bytes": events.stat().st_size,
            })
    args.output.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    print(f"{len(rows)} eligible ligand dashboards without a NETWORX explorer")


if __name__ == "__main__":
    main()
