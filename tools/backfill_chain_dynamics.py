#!/usr/bin/env python3
"""Recover chain-resolved Cα RMSD/RMSF and protein Rg from a finished GROMACS run.

This is deliberately narrow: it never changes the production trajectory and
does not rerun contact chemistry.  It is suitable for historical dossiers that
already have ``md_0_1.tpr``, ``md_0_1.xtc`` and PyMACS ``atomIndex.txt``.
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


def components(path: Path):
    found = []
    for line in path.read_text(encoding="utf-8").splitlines():
        body = line.split("#", 1)[0].strip()
        match = re.match(r"(\S+)\s+(\d+)\s*-\s*(\d+)$", body)
        if match:
            found.append((match.group(1), int(match.group(2)), int(match.group(3))))
    if not found:
        raise ValueError(f"No PyMACS component ranges in {path}")
    return found


def numeric_xvg(path: Path):
    rows = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.lstrip().startswith(("#", "@")):
            continue
        fields = raw.split()
        if len(fields) >= 2:
            try:
                rows.append((float(fields[0]), float(fields[1])))
            except ValueError:
                pass
    if len(rows) < 2:
        raise ValueError(f"No numeric trace in {path}")
    return rows


def run(command, stdin=""):
    subprocess.run(command, input=stdin, text=True, check=True)


def write_csv(path: Path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--gmx", default="gmx_mpi")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    tpr, xtc, mapping = (run_dir / "md_0_1.tpr", run_dir / "md_0_1.xtc", run_dir / "atomIndex.txt")
    for source in (tpr, xtc, mapping):
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(source)
    outdir = run_dir / "Analysis_Results"
    outdir.mkdir(exist_ok=True)
    for label, start, end in components(mapping):
        stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
        index, rmsd_xvg, rmsf_xvg, rg_xvg = (
            outdir / f".chain_{stem}.ndx", outdir / f".chain_{stem}_rmsd.xvg",
            outdir / f".chain_{stem}_rmsf.xvg", outdir / f".chain_{stem}_rg.xvg",
        )
        # atomIndex ranges are GROMACS' one-based inclusive atom numbering.
        run([args.gmx, "select", "-s", str(tpr), "-select", f"atomnr {start} to {end} and name CA", "-on", str(index)])
        run([args.gmx, "rms", "-s", str(tpr), "-f", str(xtc), "-n", str(index), "-o", str(rmsd_xvg), "-tu", "ns"], "0\n0\n")
        run([args.gmx, "rmsf", "-s", str(tpr), "-f", str(xtc), "-n", str(index), "-o", str(rmsf_xvg), "-res", "-fit"], "0\n")
        # All-atom component Rg complements the Cα RMSD/RMSF measures above.
        # The atomIndex component ranges are contiguous production atom ranges.
        run([args.gmx, "select", "-s", str(tpr), "-select", f"atomnr {start} to {end}", "-on", str(index)])
        run([args.gmx, "gyrate", "-s", str(tpr), "-f", str(xtc), "-n", str(index), "-o", str(rg_xvg), "-tu", "ns"], "0\n")
        write_csv(outdir / f"RMSD_{stem}.csv", ("Time_ns", "RMSD_Å", "Chain_Component", "Atom_Selection", "Reference_Frame"),
                  ((time, value * 10.0, label, "CA (one per residue)", 0) for time, value in numeric_xvg(rmsd_xvg)))
        write_csv(outdir / f"RMSF_{stem}_per_residue.csv", ("Residue", "Original_Residue", "RMSF_Å", "Chain_Component", "Anchor_Atom"),
                  ((int(residue), int(residue), value * 10.0, label, "CA") for residue, value in numeric_xvg(rmsf_xvg)))
        write_csv(outdir / f"Radius_of_Gyration_{stem}.csv", ("Time_ns", "Rg_Component_Å", "Chain_Component", "Atom_Selection"),
                  ((time, value * 10.0, label, "all component atoms") for time, value in numeric_xvg(rg_xvg)))
        for temporary in (index, rmsd_xvg, rmsf_xvg, rg_xvg):
            temporary.unlink(missing_ok=True)
        print(f"Recovered Cα dynamics for {label} ({start}-{end})")
    # Existing index.ndx has the production Protein selection.  Rg is not a
    # backbone metric: it describes compactness of all protein atoms.
    index = run_dir / "index.ndx"
    if index.is_file():
        rg_xvg = outdir / ".protein_rg.xvg"
        run([args.gmx, "gyrate", "-s", str(tpr), "-f", str(xtc), "-n", str(index), "-o", str(rg_xvg), "-tu", "ns"], "Protein\n")
        write_csv(outdir / "Radius_of_Gyration_Protein.csv", ("Time_ns", "Rg_Protein_Å"),
                  ((time, value * 10.0) for time, value in numeric_xvg(rg_xvg)))
        rg_xvg.unlink(missing_ok=True)
        print("Recovered protein radius of gyration")


if __name__ == "__main__":
    main()
