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
