#!/usr/bin/env python3
"""Embed already-rendered NETWORX explorers into cached ligand dossiers.

This recovery helper is deliberately narrow: it only touches ligand reports
that already have both typed interaction events and a generated explorer.
It leaves data-incomplete reports unchanged and saves the prior HTML once.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path


def load_dashboard_builder(path: Path):
    spec = importlib.util.spec_from_file_location("pymacs_dashboard", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load dashboard builder: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-root", required=True, type=Path)
    parser.add_argument("--dashboard-script", required=True, type=Path)
    args = parser.parse_args()
    builder = load_dashboard_builder(args.dashboard_script)
    updated: list[Path] = []
    for report in args.cache_root.glob("**/Analysis_Results/PyMACS_Result_Dashboard.html"):
        body = report.read_text(encoding="utf-8", errors="ignore")
        ligand_mode = "Ligand / PROTAC:" in body or "Ligand interaction type over time" in body
        if not ligand_mode or "Interactive ligand interaction explorer" in body:
            continue
        analysis = report.parent
        if not (analysis / "InteractionTypes_Framewise.csv").is_file():
            continue
        if not list((analysis / "NETWORX").glob("*_interactive.html")):
            continue
        explorer = builder.interactive_networx_html(analysis, ligand_expected=True)
        if "Interactive ligand interaction explorer" not in explorer:
            raise RuntimeError(f"could not build explorer for {report}")
        candidate = body.replace("</main>", explorer + "</main>", 1)
        if candidate == body or candidate.count("Interactive ligand interaction explorer") != 1:
            raise RuntimeError(f"unsafe explorer insertion for {report}")
        backup = report.with_name(report.name + ".pre-explorer-20261029")
        if not backup.exists():
            backup.write_text(body, encoding="utf-8")
        temporary = report.with_name(report.name + ".explorer-partial")
        temporary.write_text(candidate, encoding="utf-8")
        if temporary.stat().st_size < 100_000:
            raise RuntimeError(f"invalid dashboard candidate for {report}")
        temporary.replace(report)
        updated.append(report)
    for report in updated:
        print(report.relative_to(args.cache_root))
    print(f"embedded={len(updated)}")


if __name__ == "__main__":
    main()
