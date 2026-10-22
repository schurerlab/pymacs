#!/usr/bin/env python3
"""Write reference-defined ligand--protein native-contact Q for a finished MD run.

This is intentionally analysis-only: it reads an existing topology/trajectory and
adds Q CSV/PNG/summary files beneath the selected output directory.  The
reference is frame 0 of the supplied trajectory, so Q describes retention of
contacts present in that submitted/production-starting model; it is not a
binding-affinity estimate.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd

from pymacs_native_contacts import contact_fraction, reference_contact_pairs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default="Final_Trajectory.pdb")
    parser.add_argument("--trajectory", default="Final_Trajectory.xtc")
    parser.add_argument("--ligand", required=True, help="Ligand residue name, e.g. LIG or PTC")
    parser.add_argument("--outdir", default="Analysis_Results")
    parser.add_argument("--cutoff", type=float, default=4.5, help="Heavy-atom contact cutoff in angstroms")
    parser.add_argument("--stride", type=int, default=1, help="Analyze every Nth frame; default is every frame")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.cutoff <= 0 or args.stride <= 0:
        raise SystemExit("--cutoff and --stride must be positive")
    topology = Path(args.topology)
    trajectory = Path(args.trajectory)
    if not topology.is_file() or not trajectory.is_file():
        raise SystemExit("Finished topology and trajectory are both required")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    universe = mda.Universe(str(topology), str(trajectory))
    ligand = universe.select_atoms(f"resname {args.ligand} and not name H*")
    protein = universe.select_atoms("protein and not name H*")
    if not ligand.n_atoms:
        raise SystemExit(f"No heavy atoms found for ligand residue {args.ligand!r}")
    if not protein.n_atoms:
        raise SystemExit("No protein heavy atoms found")

    universe.trajectory[0]
    pairs = reference_contact_pairs(protein.positions, ligand.positions, args.cutoff)
    if not len(pairs):
        raise SystemExit(f"No frame-0 protein--{args.ligand} heavy-atom contacts within {args.cutoff:.2f} A")

    rows: list[dict[str, float | int]] = []
    frame_indices = range(0, len(universe.trajectory), args.stride)
    total = len(range(0, len(universe.trajectory), args.stride))
    for ordinal, frame_index in enumerate(frame_indices, start=1):
        ts = universe.trajectory[frame_index]
        rows.append({
            "Frame": int(frame_index),
            "Time_ns": float(ts.time) / 1000.0,
            "Q_native_contact_fraction": contact_fraction(
                protein.positions, ligand.positions, pairs, args.cutoff
            ),
        })
        if ordinal % 1000 == 0 or ordinal == total:
            print(f"Q progress: {ordinal}/{total}", flush=True)

    frame_data = pd.DataFrame(rows)
    csv_path = outdir / "native_contact_q.csv"
    frame_data.to_csv(csv_path, index=False)
    mean_q = float(frame_data["Q_native_contact_fraction"].mean())
    final_q = float(frame_data["Q_native_contact_fraction"].iloc[-1])
    summary_path = outdir / "native_contact_q_summary.txt"
    summary_path.write_text(
        "Reference-defined native-contact retention Q\n"
        f"Ligand residue: {args.ligand}\n"
        f"Reference: frame 0 of {trajectory.name}\n"
        f"Cutoff (A): {args.cutoff:.3f}\n"
        f"Reference pairs: {len(pairs)}\n"
        f"Frames analyzed: {len(frame_data)}\n"
        f"Mean Q: {mean_q:.4f}\n"
        f"Final Q: {final_q:.4f}\n"
    )
    plt.figure(figsize=(9, 4.5))
    plt.plot(frame_data["Time_ns"], frame_data["Q_native_contact_fraction"], color="#18ff8b", linewidth=1.25)
    plt.ylim(-0.02, 1.02)
    plt.xlabel("Time (ns)")
    plt.ylabel("Native-contact Q")
    plt.title(f"{args.ligand} native-contact retention relative to frame 0")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(outdir / "native_contact_q.png", dpi=180)
    print(f"Q complete: mean={mean_q:.4f}, final={final_q:.4f}, frames={len(frame_data)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
