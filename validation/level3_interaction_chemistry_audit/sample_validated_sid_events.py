#!/usr/bin/env python3
"""Generate a fast, trajectory-backed validated-event sample for SID review.

This is a display-validation job, not a replacement for the full Level-3
analysis: it samples evenly spaced frames and limits chemistry calculation to
residues observed in an existing contact table.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd

from pymacs_interaction_chemistry import (
    EVENT_COLUMNS,
    classify_ligand_residue_events,
    events_to_legacy_framewise,
    legacy_summary,
)


def parse_residue(label: str):
    text = str(label).strip().upper()
    letters = "".join(char for char in text if char.isalpha())[:3]
    digits = "".join(char for char in text if char.isdigit() or char == "-")
    return letters, int(digits)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", default="binding_pocket_only.pdb")
    parser.add_argument("--xtc", default="binding_pocket_only.xtc")
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--contacts", default="Analysis_Results/FilteredContacts_Framewise.csv")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--frames", type=int, default=500)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    contacts = pd.read_csv(args.contacts, usecols=["Residue"])
    requested = {parse_residue(value) for value in contacts["Residue"].dropna().unique()}
    universe = mda.Universe(args.pdb, args.xtc)
    ligand = universe.select_atoms(f"resname {args.ligand}")
    residues = [res for res in universe.residues if (str(res.resname).upper(), int(res.resid)) in requested]
    if ligand.n_atoms == 0 or not residues:
        raise SystemExit("No ligand or historical contact residues found in pocket trajectory.")
    frame_indices = np.linspace(0, len(universe.trajectory) - 1, min(args.frames, len(universe.trajectory)), dtype=int)
    rows = []
    for ordinal, frame in enumerate(np.unique(frame_indices), 1):
        ts = universe.trajectory[int(frame)]
        for residue in residues:
            rows.extend(classify_ligand_residue_events(ligand, residue, int(frame), float(ts.time) / 1000.0))
        if ordinal % 25 == 0 or ordinal == len(frame_indices):
            print(f"validated SID sample: {ordinal}/{len(frame_indices)} frames", flush=True)
    events = pd.DataFrame(rows, columns=EVENT_COLUMNS)
    events.to_csv(outdir / "ValidatedInteractionEvents.csv", index=False)
    legacy = events_to_legacy_framewise(events)
    legacy.to_csv(outdir / "InteractionTypes_Framewise.csv", index=False)
    legacy_summary(events).to_csv(outdir / "InteractionTypes_Summary.csv", index=False)
    legacy.to_csv(outdir / "FilteredContacts_Framewise.csv", index=False)
    (outdir / "README.txt").write_text(
        f"Trajectory-backed SID validation sample: {len(np.unique(frame_indices))} evenly spaced frames; "
        f"{len(residues)} historically contacting residues. Not a full-trajectory occupancy report.\n"
    )
    print(f"wrote {len(events)} events", flush=True)


if __name__ == "__main__":
    main()
