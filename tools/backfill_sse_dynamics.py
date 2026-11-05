#!/usr/bin/env python3
"""Recover DSSP/SSE time series for a completed protein-complex trajectory.

The regular chunked 3A workflow already emits these files.  This small,
read-only recovery utility brings historical complex dossiers to the same
standard without rerunning MD or changing the production trajectory.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd


def components(path: Path) -> list[tuple[str, int, int]]:
    found = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"(\S+)\s+(\d+)\s*-\s*(\d+)$", line.split("#", 1)[0].strip())
        if match:
            found.append((match.group(1), int(match.group(2)), int(match.group(3))))
    if not found:
        raise ValueError(f"No PyMACS component ranges in {path}")
    return found


def encode(dssp: np.ndarray) -> np.ndarray:
    result = np.zeros(dssp.shape, dtype=np.int8)
    result[np.isin(dssp, ["E", "B"])] = 1
    result[np.isin(dssp, ["H", "G", "I"])] = 2
    return result


def write_fraction(path: Path, values: np.ndarray, time_ns: np.ndarray) -> None:
    pd.DataFrame({
        "Time_ns": time_ns,
        "Helix": (values == 2).mean(axis=1),
        "Sheet": (values == 1).mean(axis=1),
        "Coil": (values == 0).mean(axis=1),
    }).to_csv(path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, default=None)
    parser.add_argument("--topology", type=Path, default=None)
    parser.add_argument("--chunk-size", type=int, default=500)
    args = parser.parse_args()
    root = args.run_dir.resolve()
    trajectory = (args.trajectory or root / "Final_Trajectory.xtc").resolve()
    topology = (args.topology or root / "Final_Trajectory.pdb").resolve()
    mapping = root / "atomIndex.txt"
    for source in (trajectory, topology, mapping):
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(source)

    outdir = root / "Analysis_Results"
    outdir.mkdir(exist_ok=True)
    topology_model = md.load_frame(str(trajectory), 0, top=str(topology)).topology
    protein_residues = [residue for residue in topology_model.residues if residue.is_protein]
    if not protein_residues:
        raise RuntimeError("No protein residues were found for DSSP")
    labels = components(mapping)
    component_for_residue: list[str] = []
    for residue in protein_residues:
        atom_numbers = [atom.index + 1 for atom in residue.atoms]
        label = next((name for name, start, end in labels if atom_numbers and start <= atom_numbers[0] <= end), None)
        component_for_residue.append(label or f"chain_{residue.chain.index + 1}")

    matrices, times = [], []
    for chunk in md.iterload(str(trajectory), top=str(topology), chunk=args.chunk_size):
        matrices.append(encode(md.compute_dssp(chunk)))
        times.append(np.asarray(chunk.time, dtype=float) / 1000.0)
    matrix = np.vstack(matrices)
    time_ns = np.concatenate(times)
    # mdtraj's DSSP output describes protein residues; reject a mismatched
    # topology rather than assigning residue rows to the wrong chain.
    if matrix.shape[1] != len(protein_residues):
        raise RuntimeError(f"DSSP has {matrix.shape[1]} residues but topology has {len(protein_residues)} protein residues")

    pd.DataFrame(matrix).to_csv(outdir / "DSSP_encoded_numeric.csv", index=False)
    write_fraction(outdir / "SSE_Fractions.csv", matrix, time_ns)
    persistence = pd.DataFrame({
        "Residue_Index": np.arange(matrix.shape[1]),
        "Helix_Persistence": (matrix == 2).mean(axis=0),
        "Sheet_Persistence": (matrix == 1).mean(axis=0),
        "Coil_Persistence": (matrix == 0).mean(axis=0),
    })
    persistence.to_csv(outdir / "Residue_SSE_Persistence.csv", index=False)
    for label in dict.fromkeys(component_for_residue):
        indices = np.array([i for i, value in enumerate(component_for_residue) if value == label], dtype=int)
        if not len(indices):
            continue
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
        component_matrix = matrix[:, indices]
        pd.DataFrame(component_matrix).to_csv(outdir / f"DSSP_chain_{safe}.csv", index=False)
        write_fraction(outdir / f"SSE_chain_{safe}_fractions.csv", component_matrix, time_ns)
        chain_residues = [protein_residues[index] for index in indices]
        pd.DataFrame({
            "Chain_Local_Position": np.arange(1, len(indices) + 1),
            "Original_Residue": [f"{residue.name}{residue.resSeq}" for residue in chain_residues],
            "Helix_Persistence": (component_matrix == 2).mean(axis=0),
            "Sheet_Persistence": (component_matrix == 1).mean(axis=0),
            "Coil_Persistence": (component_matrix == 0).mean(axis=0),
        }).to_csv(outdir / f"Residue_SSE_Persistence_chain_{safe}.csv", index=False)
        print(f"Recovered SSE for {label}: {len(indices)} residues")
    print(f"Recovered DSSP/SSE for {matrix.shape[0]} frames and {matrix.shape[1]} protein residues")


if __name__ == "__main__":
    main()
