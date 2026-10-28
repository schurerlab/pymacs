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


def test_interactive_html_includes_the_persistence_switches(tmp_path):
    rows = []
    for frame in range(20):
        rows.append((frame, "TYR10", "Hydrophobic", "C1"))
        if frame < 19:
            rows.append((frame, "LYS11", "Ionic", "O1"))
        if frame < 18:
            rows.append((frame, "ASN12", "H-bond", "O1"))
        if frame < 17:
            rows.append((frame, "PHE13", "Hydrophobic", "C1"))
    framewise = pd.DataFrame(rows, columns=["Time_ns", "Residue", "Type", "LigAtomName"])

    output = networx.make_interactive_network_html(
        "LIG", framewise,
        {"C1": (0.0, 0.0), "O1": (1.0, 0.0)},
        {"C1": "C", "O1": "O"},
        [("C1", "O1", 1)], (0.5, 0.0), {}, tmp_path,
    )

    html = Path(output).read_text(encoding="utf-8")
    assert 'data-threshold="0.10"' in html
    assert 'data-threshold="0.50" aria-pressed="true"' in html
    assert 'data-threshold="0.80"' in html
    assert 'data-threshold="0.90"' in html
    assert "TYR10" in html
    assert "Hydrophobic" in html


def test_residue_badge_font_keeps_normal_labels_at_default_and_shrinks_long_ones():
    assert networx.fit_residue_badge_font_size("A", "LYS", "16", 23, 45) == 23
    assert networx.fit_residue_badge_font_size("CHAIN", "TRP", "123456", 23, 45) < 23
