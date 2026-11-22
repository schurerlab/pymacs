#!/usr/bin/env python3
"""Publish the standard PyMACS native-contact Q suite for a finished run.

The suite always measures protein-fold retention.  When a ligand or peptide
selection is supplied it also measures partner--protein retention; for
multi-chain protein systems ``--interface-q`` adds inter-chain retention.
All measurements use frame zero as an explicitly declared reference.
"""
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd

from pymacs_native_contacts import (contact_fraction, intragroup_contact_fraction,
                                    reference_contact_pairs, reference_protein_ca_pairs)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default="Final_Trajectory.pdb")
    parser.add_argument("--trajectory", default="Final_Trajectory.xtc")
    parser.add_argument("--outdir", default="Analysis_Results")
    parser.add_argument("--partner-selection", default="", help="MDAnalysis selection for a ligand or peptide")
    parser.add_argument("--partner-kind", choices=("ligand", "peptide"), default="ligand")
    parser.add_argument("--partner-label", default="")
    parser.add_argument("--interface-q", action="store_true", help="Measure protein--protein inter-chain Q")
    parser.add_argument("--contact-cutoff", type=float, default=4.5)
    parser.add_argument("--protein-contact-cutoff", type=float, default=8.0)
    parser.add_argument("--protein-sequence-separation", type=int, default=3)
    parser.add_argument("--stride", type=int, default=1)
    return parser.parse_args()


def frame_indices(universe: mda.Universe, stride: int):
    return range(0, len(universe.trajectory), stride)


def write_metric(path: Path, rows: list[dict], value_column: str) -> None:
    frame_data = pd.DataFrame(rows)
    frame_data.to_csv(path, index=False)
    if not rows:
        return
    values = pd.to_numeric(frame_data[value_column], errors="coerce").dropna()
    path.with_name(path.stem + "_summary.txt").write_text(
        f"Reference-defined native-contact retention Q\nMean Q: {values.mean():.4f}\nFinal Q: {values.iloc[-1]:.4f}\n",
        encoding="utf-8",
    )
    plt.figure(figsize=(9, 4.5))
    plt.plot(frame_data["Time_ns"], frame_data[value_column], color="#18ff8b", linewidth=1.25)
    plt.ylim(-0.02, 1.02)
    plt.xlabel("Time (ns)")
    plt.ylabel("Native-contact Q")
    plt.title(path.stem.replace("_", " "))
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=180)
    plt.close()


def append_normalized(rows: list[dict], metric: str, component_a: str, component_b: str,
                      series: list[dict], value_column: str, cutoff: float, pairs: int) -> None:
    for row in series:
        rows.append({"Metric": metric, "ComponentA": component_a, "ComponentB": component_b,
                     "Frame": row["Frame"], "Time_ns": row["Time_ns"], "Q": row[value_column],
                     "ReferenceFrame": 0, "ReferenceType": "first_trajectory_frame",
                     "Cutoff_A": cutoff, "ReferencePairs": pairs})


def main() -> int:
    args = parse_args()
    if args.stride <= 0 or args.contact_cutoff <= 0 or args.protein_contact_cutoff <= 0:
        raise SystemExit("Q cutoffs and --stride must be positive.")
    topology, trajectory, outdir = Path(args.topology), Path(args.trajectory), Path(args.outdir)
    if not topology.is_file() or not trajectory.is_file():
        raise SystemExit("Finished topology and trajectory are both required for Q analysis.")
    outdir.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe(str(topology), str(trajectory))
    protein_ca = universe.select_atoms("protein and name CA")
    protein_heavy = universe.select_atoms("protein and not name H*")
    if not protein_ca.n_atoms or not protein_heavy.n_atoms:
        raise SystemExit("Q analysis requires protein atoms and C-alpha atoms.")
    universe.trajectory[0]
    normalized: list[dict] = []

    # Protein fold Q is always present, including apo systems.
    fold_pairs = reference_protein_ca_pairs(protein_ca.positions.copy(), protein_ca.resids,
                                            args.protein_contact_cutoff, args.protein_sequence_separation)
    if len(fold_pairs):
        fold_rows = []
        for index in frame_indices(universe, args.stride):
            ts = universe.trajectory[index]
            fold_rows.append({"Frame": int(ts.frame), "Time_ns": float(ts.time) / 1000.0,
                              "Q_protein_native_contact_fraction": intragroup_contact_fraction(
                                  protein_ca.positions, fold_pairs, args.protein_contact_cutoff),
                              "ReferenceFrame": 0, "ReferenceType": "first_trajectory_frame",
                              "ContactCutoff_A": args.protein_contact_cutoff,
                              "MinimumSequenceSeparation": args.protein_sequence_separation,
                              "NativeProteinCAPairs": int(len(fold_pairs))})
        write_metric(outdir / "protein_native_contact_q.csv", fold_rows, "Q_protein_native_contact_fraction")
        append_normalized(normalized, "protein_fold", "protein", "protein", fold_rows,
                          "Q_protein_native_contact_fraction", args.protein_contact_cutoff, len(fold_pairs))

    # Ligand and peptide engagement Q share a calculation, but retain explicit labels.
    if args.partner_selection:
        partner = universe.select_atoms(args.partner_selection + " and not name H*")
        if not partner.n_atoms:
            raise SystemExit(f"No partner atoms matched: {args.partner_selection!r}")
        universe.trajectory[0]
        pairs = reference_contact_pairs(protein_heavy.positions, partner.positions, args.contact_cutoff)
        if len(pairs):
            partner_rows = []
            for index in frame_indices(universe, args.stride):
                ts = universe.trajectory[index]
                partner_rows.append({"Frame": int(ts.frame), "Time_ns": float(ts.time) / 1000.0,
                                     "Q_native_contact_fraction": contact_fraction(
                                         protein_heavy.positions, partner.positions, pairs, args.contact_cutoff),
                                     "ReferenceFrame": 0, "ReferenceType": "first_trajectory_frame",
                                     "ContactCutoff_A": args.contact_cutoff, "NativeHeavyAtomContactPairs": int(len(pairs))})
            write_metric(outdir / "native_contact_q.csv", partner_rows, "Q_native_contact_fraction")
            if args.partner_kind == "peptide":
                peptide_rows = [dict(row, Q_peptide_native_contact_fraction=row["Q_native_contact_fraction"])
                                for row in partner_rows]
                write_metric(outdir / "peptide_native_contact_q.csv", peptide_rows, "Q_peptide_native_contact_fraction")
            append_normalized(normalized, f"{args.partner_kind}_protein", args.partner_label or args.partner_kind,
                              "protein", partner_rows, "Q_native_contact_fraction", args.contact_cutoff, len(pairs))

    if args.interface_q:
        chain_ids = sorted({str(atom.chainID).strip() for atom in protein_ca if str(atom.chainID).strip()})
        interface_rows = []
        universe.trajectory[0]
        for left, right in combinations(chain_ids, 2):
            left_ca, right_ca = universe.select_atoms(f"protein and name CA and chain {left}"), universe.select_atoms(f"protein and name CA and chain {right}")
            pairs = reference_contact_pairs(left_ca.positions, right_ca.positions, args.protein_contact_cutoff)
            if not len(pairs):
                continue
            for index in frame_indices(universe, args.stride):
                ts = universe.trajectory[index]
                interface_rows.append({"Frame": int(ts.frame), "Time_ns": float(ts.time) / 1000.0,
                                       "ChainA": left, "ChainB": right,
                                       "Q_interface_native_contact_fraction": contact_fraction(
                                           left_ca.positions, right_ca.positions, pairs, args.protein_contact_cutoff),
                                       "ReferenceFrame": 0, "ReferenceType": "first_trajectory_frame",
                                       "ContactCutoff_A": args.protein_contact_cutoff, "ReferenceCAPairs": int(len(pairs))})
            append_normalized(normalized, "protein_interface", left, right, interface_rows[-len(frame_indices(universe, args.stride)):],
                              "Q_interface_native_contact_fraction", args.protein_contact_cutoff, len(pairs))
        if interface_rows:
            write_metric(outdir / "protein_interface_native_contact_q.csv", interface_rows, "Q_interface_native_contact_fraction")

    if normalized:
        pd.DataFrame(normalized).to_csv(outdir / "native_contact_q_metrics.csv", index=False)
    print(f"Q suite complete: {len(normalized)} metric frames written to {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
