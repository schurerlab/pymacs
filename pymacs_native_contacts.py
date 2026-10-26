"""Reference-defined native-contact retention (Q) helpers for PyMACS.

Q is the fraction of heavy-atom contacts present in a declared reference frame
that are retained in each sampled trajectory frame.  It is deliberately
reference-defined: it does not imply that the reference is experimental or a
biological native state.
"""

from __future__ import annotations

import numpy as np


def reference_contact_pairs(
    group_a_positions: np.ndarray,
    group_b_positions: np.ndarray,
    cutoff_angstrom: float,
) -> np.ndarray:
    """Return local heavy-atom index pairs within the reference cutoff."""
    a = np.asarray(group_a_positions, dtype=float)
    b = np.asarray(group_b_positions, dtype=float)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != 3 or b.shape[1] != 3:
        raise ValueError("Coordinates must have shape (n_atoms, 3).")
    if cutoff_angstrom <= 0:
        raise ValueError("cutoff_angstrom must be positive.")
    distances_sq = np.sum((a[:, None, :] - b[None, :, :]) ** 2, axis=-1)
    return np.argwhere(distances_sq <= cutoff_angstrom ** 2).astype(int, copy=False)


def contact_fraction(
    group_a_positions: np.ndarray,
    group_b_positions: np.ndarray,
    reference_pairs: np.ndarray,
    cutoff_angstrom: float,
) -> float:
    """Return the fraction of reference contact pairs retained in one frame."""
    pairs = np.asarray(reference_pairs, dtype=int)
    if pairs.size == 0:
        return float("nan")
    a = np.asarray(group_a_positions, dtype=float)
    b = np.asarray(group_b_positions, dtype=float)
    deltas = a[pairs[:, 0]] - b[pairs[:, 1]]
    retained = np.sum(np.einsum("ij,ij->i", deltas, deltas) <= cutoff_angstrom ** 2)
    return float(retained / len(pairs))


def reference_protein_ca_pairs(
    ca_positions: np.ndarray,
    residue_ids: np.ndarray,
    cutoff_angstrom: float = 8.0,
    minimum_sequence_separation: int = 3,
) -> np.ndarray:
    """Return nonlocal Cα contact pairs in a declared protein reference frame.

    ``minimum_sequence_separation=3`` excludes contacts between residues that
    are too close in primary sequence to describe tertiary structure.  Residue
    IDs, rather than array offsets, are used so numbering gaps and multiple
    chains retain their original identity.
    """
    positions = np.asarray(ca_positions, dtype=float)
    resids = np.asarray(residue_ids)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("Cα coordinates must have shape (n_atoms, 3).")
    if len(positions) != len(resids):
        raise ValueError("Cα coordinates and residue IDs must have the same length.")
    if cutoff_angstrom <= 0 or minimum_sequence_separation < 0:
        raise ValueError("Cutoff must be positive and sequence separation non-negative.")
    distances_sq = np.sum((positions[:, None, :] - positions[None, :, :]) ** 2, axis=-1)
    sequence_gap = np.abs(resids[:, None] - resids[None, :])
    return np.argwhere(
        np.triu(sequence_gap > minimum_sequence_separation, k=1)
        & (distances_sq <= cutoff_angstrom ** 2)
    ).astype(int, copy=False)


def intragroup_contact_fraction(
    positions: np.ndarray,
    reference_pairs: np.ndarray,
    cutoff_angstrom: float,
) -> float:
    """Return retained fraction for reference pairs within one atom group."""
    pairs = np.asarray(reference_pairs, dtype=int)
    if pairs.size == 0:
        return float("nan")
    coords = np.asarray(positions, dtype=float)
    deltas = coords[pairs[:, 0]] - coords[pairs[:, 1]]
    retained = np.sum(np.einsum("ij,ij->i", deltas, deltas) <= cutoff_angstrom ** 2)
    return float(retained / len(pairs))
