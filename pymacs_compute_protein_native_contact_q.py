#!/usr/bin/env python3
"""Write frame-0 Cα native-contact retention for an apo protein trajectory.

This is intentionally distinct from PyMACS' ligand--protein Q(x).  It reports
the fraction of nonlocal Cα contacts present in the first trajectory frame
that remain within the documented cutoff, and is appropriate when no ligand
is present in the archived topology.
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default="Final_Trajectory.pdb")
    parser.add_argument("--trajectory", default="Final_Trajectory.xtc")
    parser.add_argument("--outdir", default="Analysis_Results")
    parser.add_argument("--cutoff", type=float, default=8.0, help="Cα contact cutoff in angstroms")
    parser.add_argument("--sequence-separation", type=int, default=3,
                        help="Exclude residue pairs within this sequence separation")
    parser.add_argument("--stride", type=int, default=1, help="Analyze every Nth frame")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.cutoff <= 0 or args.stride <= 0 or args.sequence_separation < 0:
        raise SystemExit("--cutoff and --stride must be positive; --sequence-separation must be non-negative")
    topology, trajectory = Path(args.topology), Path(args.trajectory)
    if not topology.is_file() or not trajectory.is_file():
        raise SystemExit("Finished topology and trajectory are both required")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    universe = mda.Universe(str(topology), str(trajectory))
    alpha_carbons = universe.select_atoms("protein and name CA")
    if alpha_carbons.n_atoms < 2:
        raise SystemExit("At least two protein Cα atoms are required")
    sequence_gap = np.abs(alpha_carbons.resids[:, None] - alpha_carbons.resids[None, :])
    upper_nonlocal = np.triu(sequence_gap > args.sequence_separation, k=1)
    universe.trajectory[0]
    reference_distances = np.linalg.norm(
        alpha_carbons.positions[:, None, :] - alpha_carbons.positions[None, :, :], axis=-1
    )
    native_pairs = np.argwhere(upper_nonlocal & (reference_distances <= args.cutoff))
    if not len(native_pairs):
        raise SystemExit("No nonlocal frame-0 Cα contacts found at the selected cutoff")

    rows: list[dict[str, float | int]] = []
    frame_indices = range(0, len(universe.trajectory), args.stride)
    total = len(range(0, len(universe.trajectory), args.stride))
    for ordinal, frame_index in enumerate(frame_indices, start=1):
        ts = universe.trajectory[frame_index]
        distances = np.linalg.norm(
            alpha_carbons.positions[native_pairs[:, 0]] - alpha_carbons.positions[native_pairs[:, 1]], axis=1
        )
        rows.append({"Frame": int(frame_index), "Time_ns": float(ts.time) / 1000.0,
                     "Q_protein_native_contact_fraction": float(np.mean(distances <= args.cutoff))})
        if ordinal % 1000 == 0 or ordinal == total:
            print(f"Protein native-contact Q progress: {ordinal}/{total}", flush=True)

    frame_data = pd.DataFrame(rows)
    frame_data.to_csv(outdir / "protein_native_contact_q.csv", index=False)
    mean_q = float(frame_data["Q_protein_native_contact_fraction"].mean())
    final_q = float(frame_data["Q_protein_native_contact_fraction"].iloc[-1])
    (outdir / "protein_native_contact_q_summary.txt").write_text(
        "Protein Cα native-contact retention Q\n"
        f"Reference: frame 0 of {trajectory.name}\n"
        f"Contact cutoff (A): {args.cutoff:.3f}\n"
        f"Minimum sequence separation: {args.sequence_separation}\n"
        f"Reference Cα contact pairs: {len(native_pairs)}\n"
        f"Frames analyzed: {len(frame_data)}\n"
        f"Mean Q_protein: {mean_q:.4f}\n"
        f"Final Q_protein: {final_q:.4f}\n",
        encoding="utf-8",
    )
    plt.figure(figsize=(9, 4.5))
    plt.plot(frame_data["Time_ns"], frame_data["Q_protein_native_contact_fraction"], color="#18ff8b", linewidth=1.25)
    plt.ylim(-0.02, 1.02)
    plt.xlabel("Time (ns)")
    plt.ylabel("Protein native-contact Q")
    plt.title("Protein Cα native-contact retention relative to frame 0")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(outdir / "protein_native_contact_q.png", dpi=180)
    print(f"Protein native-contact Q complete: mean={mean_q:.4f}, final={final_q:.4f}, frames={len(frame_data)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
