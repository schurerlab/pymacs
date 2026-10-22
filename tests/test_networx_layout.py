import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "3B_NETWORX.py"
SPEC = importlib.util.spec_from_file_location("networx", SCRIPT)
networx = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(networx)


def test_specific_edge_fraction_uses_distinct_trajectory_frames():
    framewise = pd.DataFrame({
        "Time_ns": [0.0, 0.0, 1.0, 2.0, 2.0],
        "Residue": ["ALA1", "ALA1", "ALA1", "ALA1", "GLY2"],
        "Type": ["Hydrophobic", "Hydrophobic", "Hydrophobic", "Ionic", "Ionic"],
        "LigAtomName": ["H1", "H1", "H1", "H1", "N2"],
    })

    edges = networx.summarize_interaction_edges(framewise, 0.0)
    hydrophobic = edges[(edges.Residue == "ALA1") & (edges.Type == "Hydrophobic")].iloc[0]
    ionic = edges[(edges.Residue == "ALA1") & (edges.Type == "Ionic")].iloc[0]

    assert hydrophobic.Frames == 2
    assert hydrophobic.Fraction == 2 / 3
    assert ionic.Fraction == 1 / 3


def test_guided_layout_never_stacks_a_dense_contact_hotspot():
    residues = [f"ALA{i}" for i in range(1, 25)]
    atom_positions = {"Br1": (0.0, 0.0), "C1": (1.0, 0.0)}
    positions = networx.auto_place_residues_guided(
        residues, {residue: ["Br1"] for residue in residues}, atom_positions, (0.5, 0.0)
    )
    coords = np.array([positions[residue] for residue in residues])
    pairwise = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(pairwise, np.inf)

    assert pairwise.min() >= 2 * networx.RESIDUE_RADIUS + 0.18 - 1e-8
