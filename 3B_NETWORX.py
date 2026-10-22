#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NETWORX v3 — Journal-Grade Ligand–Protein Interaction Diagram
--------------------------------------------------------------
Implements:
 • Two-panel figure: (a) stacked bar summary, (b) oval network
 • Oval layout for residues
 • Curved edges from ligand atoms to residues
 • Edge labels at 75% distance (fractions)
 • Bond order rendering (1/2/3 parallel lines)
 • Full halogen / sulfur / phosphorus coloring
 • Hydrogen bond edges: thick dashed maroon
 • Hydrophobic edges: soft dashed light-green
 • Water bridges: dashed blue
 • Ionic edges omitted in network
 • Final figure saved at high DPI

Assumptions:
 • Run from MD system root (where Analysis_Results/ lives)
 • Required files:
     Analysis_Results/InteractionTypes_Summary.csv
     Analysis_Results/FilteredContacts_Framewise.csv
 • Ligand MOL2:
     <ligand>.cgenff.mol2
"""

import os
import math
import argparse
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib.patches import Patch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
import json
import shutil
import subprocess
import tempfile
from rdkit import Chem
from rdkit.Chem import rdDepictor
import re
import math

HALOGEN_ELEMENTS = {"F", "Cl", "Br", "I"}

def is_halogen(atom_name, atom_elements):
    """
    Return True if atom_name corresponds to a halogen (F, Cl, Br, I),
    allowing suffix digits (e.g. Br12, Cl3).
    """
    if atom_name not in atom_elements:
        return False

    elem = atom_elements[atom_name]  # e.g. "Br", "Cl"
    return elem in HALOGEN_ELEMENTS


def find_nearest_halogen(atom_positions, halogen_atoms, center):
    """
    If multiple halogens exist, pick the one *closest to center*.
    Otherwise return the single halogen available.
    """
    if not halogen_atoms:
        return None

    if len(halogen_atoms) == 1:
        return halogen_atoms[0]

    cx, cy = center
    best = None
    best_d = 9999

    for h in halogen_atoms:
        x, y = atom_positions[h]
        d = math.hypot(x - cx, y - cy)
        if d < best_d:
            best = h
            best_d = d

    return best




def _segments_properly_intersect(a, b, c, d):
    """True when two finite segments cross away from their endpoints."""
    def orient(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    o1, o2 = orient(a, b, c), orient(a, b, d)
    o3, o4 = orient(c, d, a), orient(c, d, b)
    return (o1 * o2 < -1e-8) and (o3 * o4 < -1e-8)


def _connector_crosses_skeleton(anchor, badge, anchor_atom, bonds, atom_positions):
    """Detect a connector that would visually cut through the ligand."""
    if not bonds:
        return False
    for a1, a2, _ in bonds:
        if anchor_atom in (a1, a2) or a1.startswith("H") or a2.startswith("H"):
            continue
        if a1 not in atom_positions or a2 not in atom_positions:
            continue
        if _segments_properly_intersect(anchor, badge, atom_positions[a1], atom_positions[a2]):
            return True
    return False


def auto_place_residues_guided(residues, residue_to_ligatom, atom_positions, center, bonds=None):
    """
    Place residues on an oval boundary, guided by the ligand atom they interact with.
    Ensures:
      - Residues sit closer to the ligand region they interact with
      - No residue overlaps another
      - No residue is allowed inside the molecular core
    """

    cx, cy = center

    # -----------------------------------------------------
    # 1. Determine ellipse radii based on ligand dimensions
    # -----------------------------------------------------
    xs = [p[0] for p in atom_positions.values()]
    ys = [p[1] for p in atom_positions.values()]

    ligand_w = (max(xs) - min(xs))
    ligand_h = (max(ys) - min(ys))

    # Base radii guided by ligand dimensions.  Do not use a fixed oval here:
    # a fixed oval is the reason crowded systems previously put residue labels
    # on top of one another.  The perimeter is expanded below until it can
    # physically accommodate every residue circle.
    # Keep the residue halo close enough that the ligand remains the visual
    # centrepiece.  The collision-expansion below, rather than a large fixed
    # multiplier, supplies any required clearance for dense systems.
    # A SID is a ligand drawing with residue annotations around it.  Keep
    # that ring close enough that the ligand occupies the centrepiece, while
    # collision relaxation still prevents label overlap.
    rx = ligand_w * 0.58
    ry = ligand_h * 0.66

    # Enforce sensible minimums so very compact ligand layouts don't
    # place residues on top of the ligand. These constants keep behaviour
    # local to the network builders and won't affect other graphs.
    rx = max(rx, ELLIPSE_RX * 0.9)
    ry = max(ry, ELLIPSE_RY * 0.9)

    # -----------------------------------------------------
    # 2. Compute target angle for each residue
    # -----------------------------------------------------
    angles = {}

    for res in residues:

        lig_atoms = residue_to_ligatom.get(res, [])
        if isinstance(lig_atoms, str):
            lig_atoms = [lig_atoms]

        # Default angle: evenly distributed if no known atom
        if not lig_atoms:
            idx = residues.index(res)
            angles[res] = 2*np.pi * idx / len(residues)
            continue

        # Compute average angle from center to its interacting ligand atoms
        angs = []
        for la in lig_atoms:
            if la not in atom_positions:
                continue
            xL, yL = atom_positions[la]
            ang = np.arctan2(yL - cy, xL - cx)
            angs.append(ang)

        if angs:
            angles[res] = np.mean(angs)
        else:
            idx = residues.index(res)
            angles[res] = 2*np.pi * idx / len(residues)

    # Sort once by the molecular region contacted.  Ordering around the
    # boundary is retained throughout relaxation, which prevents edge fans
    # from needlessly crossing the ligand drawing.
    sorted_res = sorted(residues, key=lambda r: (angles[r], residue_number(r), r))
    # Marker dimensions are specified in points while the placement operates
    # in ligand-coordinate units.  Empirical 2D clearance keeps the drawn
    # badges (and their labels) apart at publication DPI.
    min_dist = 2.35
    if len(sorted_res) == 1:
        return {sorted_res[0]: (cx + rx * math.cos(angles[sorted_res[0]]),
                                cy + ry * math.sin(angles[sorted_res[0]]))}

    # Preserve each residue's anchor angle.  The previous evenly-spaced-slot
    # approach prevented badge overlap but could move a contact halfway around
    # the ligand, producing long connectors that crossed the scaffold.  Here
    # only residues sharing a local hotspot are fanned apart.
    target = [angles[r] for r in sorted_res]
    for i in range(1, len(target)):
        while target[i] < target[i - 1]:
            target[i] += 2.0 * math.pi

    # Expand the halo only until adjacent badges can be separated.  A small
    # angular gap is then enforced in the original contact order.
    n = len(sorted_res)
    for _ in range(18):
        angular_gap = 2.0 * math.asin(min(0.99, min_dist / (2.0 * min(rx, ry))))
        placed = [target[0]]
        for desired in target[1:]:
            placed.append(max(desired, placed[-1] + angular_gap))
        if placed[-1] - placed[0] <= 2.0 * math.pi - angular_gap:
            break
        rx *= 1.12
        ry *= 1.12

    # Centre the local fan on the requested anchor directions.  This changes
    # neither order nor inter-badge separation, but minimizes total connector
    # length instead of distributing unused space around the full ellipse.
    shift = sum(t - p for t, p in zip(target, placed)) / n
    placed_angles = {res: placed[i] + shift for i, res in enumerate(sorted_res)}

    for _expand_pass in range(18):
        positions = {
            r: (cx + rx * math.cos(placed_angles[r]), cy + ry * math.sin(placed_angles[r]))
            for r in sorted_res
        }
        node_clear = RESIDUE_RADIUS + LIG_ATOM_RADIUS + 0.35
        ligand_clear = all(
            math.hypot(x - ax, y - ay) >= node_clear
            for x, y in positions.values() for ax, ay in atom_positions.values()
        )
        if ligand_clear:
            # Enforce line-of-sight to each actual contact atom.  For a
            # blocked direct connector, search nearby boundary locations and
            # choose the shortest visible one that also keeps badges apart.
            # This prevents the visually unacceptable case of a dashed edge
            # crossing a ligand bond merely to save a little distance.
            for res in sorted_res:
                lig_atoms = residue_to_ligatom.get(res, [])
                if isinstance(lig_atoms, str):
                    lig_atoms = [lig_atoms]
                anchor_atom = next((a for a in lig_atoms if a in atom_positions), None)
                if anchor_atom is None:
                    continue
                anchor = atom_positions[anchor_atom]
                if not _connector_crosses_skeleton(anchor, positions[res], anchor_atom, bonds, atom_positions):
                    continue
                base_angle = placed_angles[res]
                candidates = []
                for step in range(1, 25):
                    delta = step * math.pi / 36.0
                    candidates.extend((base_angle - delta, base_angle + delta))
                for candidate_angle in candidates:
                    candidate = (cx + rx * math.cos(candidate_angle),
                                 cy + ry * math.sin(candidate_angle))
                    if _connector_crosses_skeleton(anchor, candidate, anchor_atom, bonds, atom_positions):
                        continue
                    if any(other != res and math.dist(candidate, point) < min_dist
                           for other, point in positions.items()):
                        continue
                    positions[res] = candidate
                    break
            return positions
        rx *= 1.12
        ry *= 1.12

    return positions





def fix_edges_to_center_safe(residue_to_ligatom, atom_positions, atom_elements, center):
    """
    Only redirect ligand atoms that *remain unresolved* after
    hydrogen-parent mapping AND LP mapping.

    This prevents valid atoms like H5, H4, H6 from being replaced
    incorrectly by halogens. Also attempt similarity mapping for cases
    like 'N3' when possible.
    """

    def resolve_similar_atom(name, atom_positions, atom_elements, center=None):
        a = normalize_atom(name)
        if not a:
            return None
        if a in atom_positions:
            return a

        m = re.match(r"([A-Za-z]+)(\d*)", a)
        if not m:
            return None
        prefix, num = m.group(1), m.group(2)

        candidates = [nm for nm in atom_positions.keys()
                      if atom_elements.get(nm, "").upper() == prefix.upper() or nm.upper().startswith(prefix.upper())]
        if not candidates:
            candidates = [nm for nm in atom_positions.keys() if nm[0].upper() == prefix[0].upper()]

        if not candidates:
            return None
        if len(candidates) == 1:
            return candidates[0]

        if num:
            try:
                t = int(num)
                best, bestd = None, float("inf")
                for c in candidates:
                    dm = re.search(r"(\d+)$", c)
                    if dm:
                        v = int(dm.group(1))
                        d = abs(v - t)
                    else:
                        d = 10**6
                    if d < bestd:
                        best, bestd = c, d
                if best:
                    return best
            except Exception:
                pass

        if center is not None:
            cx, cy = center
            best, bestd = None, float("inf")
            for c in candidates:
                x, y = atom_positions[c]
                d = math.hypot(x - cx, y - cy)
                if d < bestd:
                    best, bestd = c, d
            return best

        return candidates[0]


    # Collect halogens
    halogens = [nm for nm, el in atom_elements.items() if el in {"F", "Cl", "Br", "I"}]

    if not halogens:
        print("⚠ No halogens available for fallback — attempting similarity mapping for unresolved atoms.")
        updated = residue_to_ligatom.copy()
        for res, ligatom in residue_to_ligatom.items():
            if not ligatom:
                continue
            ligatom_norm = normalize_atom(ligatom)
            if ligatom_norm in atom_positions:
                continue
            if ligatom_norm.startswith("H"):
                print(f"⚠ Unresolved hydrogen '{ligatom}' for {res} — keeping unresolved for now.")
                continue
            cand = resolve_similar_atom(ligatom_norm, atom_positions, atom_elements, center=center)
            if cand and cand in atom_positions:
                print(f"🔧 Mapping unresolved '{ligatom}' for {res} → similar atom '{cand}'")
                updated[res] = cand
            else:
                print(f"⚠ Could not map unresolved '{ligatom}' for {res}")
        return updated

    # fallback halogen path unchanged
    cx, cy = center
    halogen = min(
        halogens,
        key=lambda h: math.hypot(atom_positions[h][0] - cx, atom_positions[h][1] - cy)
    )

    updated = residue_to_ligatom.copy()
    for res, ligatom in residue_to_ligatom.items():
        if not ligatom:
            continue
        ligatom_norm = normalize_atom(ligatom)
        if ligatom_norm in atom_positions:
            continue
        if ligatom_norm.startswith("H"):
            print(f"⚠ Unresolved hydrogen '{ligatom}' for {res} — keeping center fallback.")
            continue
        print(f"🔧 Redirecting *UNRESOLVED* ligand atom '{ligatom}' for {res} → halogen {halogen}")
        updated[res] = halogen

    return updated


# --- small utility needed by several loaders (must be defined early) ---
def normalize_atom(a):
    """Clean ligand atom names (remove '.', spaces, etc.). Returns empty string for non-strings."""
    if not isinstance(a, str):
        return ""
    return "".join(c for c in a.strip() if c.isalnum())

def residue_number(resname: str) -> int:
    """Extract integer residue number from a label like 'LYS73' or 'ALA 123B'."""
    s = str(resname)
    digits = "".join(ch for ch in s if ch.isdigit())
    try:
        return int(digits) if digits else 0
    except Exception:
        return 0


def residue_badge_lines(residue_key: str):
    """Return (chain, residue code, number) for a three-line badge."""
    value = str(residue_key)
    chain, residue = (value.split("|", 1) + [""])[:2] if "|" in value else ("", value)
    match = re.search(r"([A-Za-z]{3})\s*(\d+[A-Za-z]?)", residue)
    if match:
        return chain or "–", match.group(1).upper(), match.group(2)
    return chain or "–", residue.upper(), ""

def _resolve_h_to_heavy(atom, parent_map, atom_positions, atom_elements, center, max_heavy_dist=2.5):
    """
    Map a possibly-hydrogen ligand atom name to the corresponding heavy atom.
    Strategy (in order):
      1) If atom is already a heavy atom in atom_positions -> return it
      2) If parent_map has direct mapping (H -> heavy) -> use it
      3) If atom starts with 'H', try prefix matching against parent_map keys
      4) If hydrogen has 2D coords, pick nearest heavy within max_heavy_dist
      5) Fallback: return heavy atom nearest to ligand center
    Returns heavy_atom_name or None.
    """
    a = normalize_atom(atom)
    if not a:
        return None

    # 1) already heavy
    if a in atom_positions and not a.startswith("H"):
        return a

    # 2) direct parent_map mapping
    if a in parent_map:
        parent = parent_map[a]
        if parent in atom_positions:
            return parent

    # 3) prefix match for hydrogens (H3 -> H31/H32)
    if a.startswith("H"):
        root = a.rstrip("0123456789")
        for h_mol2, heavy in parent_map.items():
            if h_mol2.startswith(root) and heavy in atom_positions:
                return heavy

    # 4) if H coord exists, pick nearest heavy within cutoff
    if a in atom_positions:
        hx, hy = atom_positions[a]
        best = None
        best_d = float("inf")
        for nm, (x, y) in atom_positions.items():
            # skip hydrogens
            if atom_elements.get(nm, "").upper() == "H" or nm.startswith("H"):
                continue
            d = math.hypot(x - hx, y - hy)
            if d < best_d:
                best = nm
                best_d = d
        if best and best_d <= max_heavy_dist:
            return best

    # 5) fallback: pick heavy atom nearest ligand center
    cx, cy = center
    best = None
    best_d = float("inf")
    for nm, (x, y) in atom_positions.items():
        if atom_elements.get(nm, "").upper() == "H" or nm.startswith("H"):
            continue
        d = math.hypot(x - cx, y - cy)
        if d < best_d:
            best = nm
            best_d = d
    return best
# --------------------------------------------------------------------
# Global paths
# --------------------------------------------------------------------
ANALYSIS_DIR = os.environ.get("PYMACS_ANALYSIS_DIR", "Analysis_Results")

# --------------------------------------------------------------------
# Color definitions
# --------------------------------------------------------------------

# Atom colors (Schrödinger-like)
ATOM_COLORS = {
    "C":  "#555555",   # dark gray
    "H":  "#D0D0D0",   # light gray
    "O":  "#E41A1C",   # red
    "N":  "#377EB8",   # blue
    "S":  "#FFD700",   # bright yellow
    "P":  "#FF69B4",   # hot pink

    "F":  "#7CFC00",   # light green
    "Cl": "#32CD32",   # lime green
    "Br": "#228B22",   # forest green
    "I":  "#800080",   # iodine purple
}

# Network edge colors
EDGE_HBOND = "#800000"    # maroon, dashed, thick
EDGE_HYDRO = "#2B8A3E"    # high-contrast green, dashed
EDGE_WATER = "#1E90FF"    # blue, dashed
    # FULL CONTACTS (for Panel C) — USE THE CORRECT FILE
# Bar plot colors
BAR_HBOND  = "#D62728"
BAR_HYDRO  = "#2CA02C"
BAR_IONIC  = "#FF7F0E"
BAR_WATER  = "#1F77B4"

# Radii for layout
# Geometry and visual radius of peripheral residue badges.  The clearance
# must accommodate their text, not just the marker face.
RESIDUE_RADIUS = 0.48
LIG_ATOM_RADIUS = 0.18   # radius of ligand atom circle (increased for visual clearance)

# Oval parameters (Option B)
# Collapsed network only
ELLIPSE_RX = 3.0
ELLIPSE_RY = 2.2

# Ligand 2D layout scaling factor (increase to spread ligand atom nodes)
# Publication SID mode gives the molecule visual priority.  Residues are
# placed around it afterwards, so enlarging this does not cause overlap.
LIGAND_LAYOUT_SCALE = 2.80

RESIDUE_CATEGORIES = {
    "hydrophobic": {"ALA","VAL","LEU","ILE","MET","PHE","TRP","TYR","PRO"},
    "positive":    {"LYS","ARG","HIS"},       # HIS treated as positive
    "negative":    {"ASP","GLU"},
    "polar":       {"SER","THR","ASN","GLN","CYS"},
}

RESIDUE_COLORS = {
    "hydrophobic": "#C9F7C5",  # soft green
    "positive":    "#C5F1FF",  # soft cyan
    "negative":    "#FFCBC5",  # soft salmon
    "polar":       "#E8E8E8",  # soft gray
    "peptide":     "#D8B4FE",  # soft violet
}

# --------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Generate journal-ready PROTAC interaction diagrams",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        "-l", "--ligand",
        type=str,
        help="Ligand residue name (e.g., A1D) or peptide group name."
    )

    p.add_argument(
        "-f", "--filter",
        dest="minfrac",
        type=float,
        default=None,
        help="Minimum contact fraction cutoff"
    )

    p.add_argument(
        "-o", "--output",
        type=str,
        help="Output base filename (without extension)"
    )

    p.add_argument(
        "--mode",
        choices=["auto", "ligand", "peptide"],
        default="auto",
        help="Auto-detect peptide mode from .chainmap.peptide.json, or force ligand/peptide rendering."
    )

    # NEW — allow node-only panel
    p.add_argument(
        "--no-edges",
        action="store_true",
        help="Render network diagram without edges (node-only mode)."
    )

    return p.parse_args()


def ask_if_missing(args):
    print("\n🧬 NETWORX — Two-Panel Journal Mode\n")

    chainmap = infer_mode_from_chainmap(".")

    if not args.ligand:
        if chainmap.get("is_peptide"):
            args.ligand = chainmap.get("group_name") or chainmap.get("name") or "PEPTIDE"
            print(f"🧬 Auto-detected peptide label from {chainmap['path']}: {args.ligand}")
        else:
            args.ligand = input("Ligand resname (e.g., A1D): ").strip()

    if args.minfrac is None:
        v = input("Minimum contact fraction [default 0.10]: ").strip()
        if v:
            v = v.replace(",", ".")
            args.minfrac = float(v)
        else:
            args.minfrac = 0.10

    # Keep output unset so we do not create a per-run subdirectory.
    if not args.output:
        args.output = None

    return args


# --------------------------------------------------------------------
# Ligand layout via MOL2 (2D drawing)
# --------------------------------------------------------------------
def mol2_atom_names_in_file_order(mol2_file):
    """Read MOL2 atom names without RDKit's PDB-name normalisation.

    CGenFF MOL2 files often call every hydrogen simply ``H``.  Their topology
    and contact CSV, however, use H1, H2, ... .  Giving duplicate hydrogens a
    stable sequential suffix restores the one-to-one mapping required for a
    faithful interaction diagram.
    """
    names = []
    in_atoms = False
    duplicate_counts = {}
    with open(mol2_file) as handle:
        for line in handle:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>"):
                if in_atoms:
                    break
                continue
            if not in_atoms or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            base = normalize_atom(parts[1])
            duplicate_counts[base] = duplicate_counts.get(base, 0) + 1
            # Preserve names such as C01 -> C1, but make repeated bare H
            # names addressable as H1, H2, ... just like the GROMACS topology.
            name = f"{base}{duplicate_counts[base]}" if base == "H" else base
            if name in names:
                name = f"{base}{duplicate_counts[base]}"
            names.append(name)
    return names


def openbabel_2d_coordinates(mol2_file, atom_names, atom_elements):
    """Return a conventional chemical 2D layout when Open Babel is present.

    The MOL/SDF atom block retains source atom order, so coordinates can be
    mapped back to the stable CGenFF names without relying on atom labels.
    ``None`` means that the portable RDKit layout should be retained.
    """
    obabel = shutil.which("obabel")
    if not obabel:
        return None
    try:
        with tempfile.TemporaryDirectory(prefix="pymacs-2d-") as tmpdir:
            mol_path = os.path.join(tmpdir, "ligand_2d.mol")
            completed = subprocess.run(
                [obabel, "-imol2", mol2_file, "-omol", "--gen2d", "-O", mol_path],
                capture_output=True, text=True, check=False,
            )
            if completed.returncode != 0 or not os.path.exists(mol_path):
                return None
            lines = open(mol_path, encoding="utf-8").read().splitlines()
            if len(lines) < 5:
                return None
            atom_count = int(lines[3][:3])
            if atom_count != len(atom_names) or len(lines) < 4 + atom_count:
                return None
            coords = {}
            for index, line in enumerate(lines[4:4 + atom_count]):
                fields = line.split()
                if len(fields) < 4:
                    return None
                name = atom_names[index]
                element = fields[3]
                if element != atom_elements.get(name):
                    return None
                coords[name] = (float(fields[0]), float(fields[1]))
            return coords
    except (OSError, ValueError):
        return None


def get_ligand_layout(ligand_resname: str):
    """
    Read <ligand_resname>.cgenff.mol2 and extract a 2D layout WITHOUT
    modifying or renaming atoms. This preserves CGenFF atom names (H1, H2, H3…),
    ensures correct hydrogen → heavy mapping, and avoids RDKit sanitization
    that breaks name consistency.

    Returns:
        atom_positions : {atom_name: (x, y)}
        atom_elements  : {atom_name: element_symbol}
        bonds          : [(atom1, atom2, order)]
        center         : (cx, cy)
        parent_map     : {Hname: heavy_atom_name}
    """

    mol2_file = f"{ligand_resname}.cgenff.mol2"
    if not os.path.exists(mol2_file):
        raise FileNotFoundError(f"MOL2 file not found: {mol2_file}")

    # ---------------------------------------------------------------
    # 1) Load MOL2 EXACTLY as-is.
    # ---------------------------------------------------------------
    # sanitize=False prevents RDKit from renaming atoms or changing valence
    mol = Chem.MolFromMol2File(mol2_file, sanitize=False, removeHs=False)
    mol_arom = restore_aromaticity(mol)

    if mol is None:
        raise ValueError(f"RDKit failed to load {mol2_file}")

    # ---------------------------------------------------------------
    # 2) Fix any missing atomic numbers (rare but possible)
    # ---------------------------------------------------------------
    pt = Chem.GetPeriodicTable()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            sym = atom.GetSymbol()
            atom.SetAtomicNum(pt.GetAtomicNumber(sym))

    # ---------------------------------------------------------------
    # 3) DO NOT SANITIZE — THIS BREAKS ATOM NAMES
    # ---------------------------------------------------------------
    # No SanitizeMol, No SetAromaticity, No Kekulize
    # These cause H13/H14 renaming, destroy mapping.

    # ---------------------------------------------------------------
    # 4) Generate 2D coordinates (safe & does NOT rename anything)
    # ---------------------------------------------------------------
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()

    # ---------------------------------------------------------------
    # 5) Extract atom positions exactly as named in MOL2
    # ---------------------------------------------------------------
    mol2_names = mol2_atom_names_in_file_order(mol2_file)
    if len(mol2_names) != mol.GetNumAtoms():
        raise ValueError(f"MOL2 atom table and RDKit molecule disagree for {mol2_file}")

    atom_positions = {}
    atom_elements = {}
    atom_names = []
    parent_map = {}

    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)

        # Use the raw MOL2 order/name, not RDKit's PDB monomer name.  The
        # latter silently turns duplicate H atoms into H13/H14/... .
        atom_name = mol2_names[i]

        atom_positions[atom_name] = (pos.x, pos.y)
        atom_elements[atom_name] = atom.GetSymbol()
        atom_names.append(atom_name)

    # Prefer a chemical depiction layout when available.  This gives rigid
    # rings, substituents, and heteroatom labels a conventional orientation;
    # RDKit's coordinate generation remains the dependency-free fallback.
    obabel_positions = openbabel_2d_coordinates(mol2_file, atom_names, atom_elements)
    if obabel_positions is not None:
        atom_positions = obabel_positions
        print("✓ Using Open Babel chemical 2D layout")

    # ---------------------------------------------------------------
    # 6) Bonds + parent map for hydrogens
    # ---------------------------------------------------------------
    # The MOL2 bond table is authoritative for names, numeric atom IDs, and
    # bond orders.  RDKit may reorder or partially reinterpret an unsanitized
    # CGenFF MOL2 (particularly heteroaromatics), which is unacceptable for a
    # publication diagram: it can attach real bond orders to the wrong atom.
    raw_bonds = parse_mol2_bonds(mol2_file)

    # Start every aromatic bond as a single graphical line.  We then add only
    # valence-safe companion strokes.  A global odd/even alternation is not
    # chemically valid for fused/heteroaromatic systems: it can make a
    # carbonyl carbon appear to have three double bonds.
    valence = {name: 0 for name in atom_elements}
    for a1, a2, raw_order in raw_bonds:
        contribution = 1 if raw_order == "AR" else int(raw_order)
        valence[a1] = valence.get(a1, 0) + contribution
        valence[a2] = valence.get(a2, 0) + contribution

    max_valence = {"C": 4, "N": 3, "O": 2, "S": 6, "P": 5,
                   "F": 1, "Cl": 1, "Br": 1, "I": 1}
    aromatic_double = set()
    for bond_index, (a1, a2, raw_order) in enumerate(raw_bonds):
        if raw_order != "AR":
            continue
        max1 = max_valence.get(atom_elements.get(a1, "C"), 4)
        max2 = max_valence.get(atom_elements.get(a2, "C"), 4)
        if valence.get(a1, 0) < max1 and valence.get(a2, 0) < max2:
            aromatic_double.add(bond_index)
            valence[a1] += 1
            valence[a2] += 1

    bonds = []
    for bond_index, (a1, a2, raw_order) in enumerate(raw_bonds):
        order = 2 if bond_index in aromatic_double else (1 if raw_order == "AR" else raw_order)
        bonds.append((a1, a2, order))

        # Build hydrogen → heavy-atom parent map from the same source table.
        sym1 = atom_elements.get(a1, "")
        sym2 = atom_elements.get(a2, "")
        if sym1 == "H" and sym2 != "H":
            parent_map[a1] = a2
        elif sym2 == "H" and sym1 != "H":
            parent_map[a2] = a1

    

    

    # ---------------------------------------------------------------
    # 7) Ligand center (for fallback)
    # ---------------------------------------------------------------
    xs = [p[0] for p in atom_positions.values()]
    ys = [p[1] for p in atom_positions.values()]
    center = (float(np.mean(xs)), float(np.mean(ys)))

    # ---------------------------------------------------------------
    # 8) Optionally scale ligand atom layout relative to center to spread nodes
    # ---------------------------------------------------------------
    if LIGAND_LAYOUT_SCALE != 1.0:
        cx, cy = center
        for nm, (x, y) in list(atom_positions.items()):
            dx = x - cx
            dy = y - cy
            atom_positions[nm] = (cx + dx * LIGAND_LAYOUT_SCALE, cy + dy * LIGAND_LAYOUT_SCALE)

    # ---------------------------------------------------------------
    # Debug
    # ---------------------------------------------------------------
    print("\n🔬 Ligand atom names from MOL2 (PRESERVED):")
    print(sorted(atom_positions.keys()))

    return atom_positions, atom_elements, bonds, center, parent_map, mol_arom, atom_names


def resolve_lp_mapping(ligand, residue_to_ligatom):
    """
    Replace LP atoms (LP1, LP2...) with nearest heavy atom based on ligand_init.pdb geometry.

    Requirements:
        ligand_init.pdb exists in working directory (e.g. A1D_ini.pdb)
    """
    import math

    pdb_file = f"{ligand}_ini.pdb"
    if not os.path.exists(pdb_file):
        print(f"[LP-MAP] No {pdb_file} found — skipping LP mapping.")
        return residue_to_ligatom

    # Load coordinates from PDB
    atom_coords = {}  # atom_name → (x,y,z)
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                nm = line[12:16].strip()   # atom name
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_coords[nm] = (x, y, z)

    # Identify heavy atoms (exclude hydrogens + LP atoms)
    heavy_atoms = [nm for nm in atom_coords if not nm.startswith("H") and not nm.startswith("LP")]

    def dist(a, b):
        return math.dist(atom_coords[a], atom_coords[b])

    updated_map = residue_to_ligatom.copy()

    for res, at in residue_to_ligatom.items():
        if not at or not at.startswith("LP"):
            continue

        if at not in atom_coords:
            continue  #

        # Find nearest heavy atom
        lp_coord = atom_coords[at]
        best = None
        best_d = 9999

        for hvy in heavy_atoms:
            d = dist(at, hvy)
            if d < best_d:
                best_d = d
                best = hvy

        if best:
            print(f"[LP-MAP] {res}: {at} → {best}  (distance {best_d:.2f} Å)")
            updated_map[res] = best

    return updated_map

def shorten_to_circle(x_center, y_center, x_out, y_out, r):
    """
    Shorten line from outer point (x_out, y_out) toward center (x_center, y_center)
    so that it ends at radius r from center.
    """
    dx = x_out - x_center
    dy = y_out - y_center
    d = math.hypot(dx, dy)
    if d == 0:
        return x_center, y_center
    scale = r / d
    return x_center + dx * scale, y_center + dy * scale


def _atom_node_color(elem: str):
    """Return the display color for a ligand atom element (node color).

    Carbon -> black, S -> yellow, P -> pink, halogens -> Miami green,
    O -> red, N -> blue, otherwise fall back to ATOM_COLORS mapping.
    """
    if not elem:
        return ATOM_COLORS.get("C", "#000000")
    if elem == "C":
        return "#000000"
    if elem == "S":
        return "#FFD700"
    if elem == "P":
        return "#FF69B4"
    if elem in ("F", "Cl", "Br", "I"):
        return "#007849"
    if elem == "O":
        return "#D72638"
    if elem == "N":
        return "#1E90FF"
    return ATOM_COLORS.get(elem, "#888")


def _draw_bicolor_arc(ax, x1, y1, x2, y2, rad, c1, c2, lw=2.0, linestyle='--', zorder=2, nseg=48):
    """Draw a curved edge from (x1,y1) to (x2,y2) with a color gradient from c1->c2.

    The curve is approximated as a quadratic Bezier using a single control
    point computed from `rad` (mimicking `arc3,rad=`). The curve is sampled
    and a LineCollection is used with per-segment colors interpolated between
    the two endpoint colors.
    """
    dx = x2 - x1
    dy = y2 - y1
    midx = 0.5 * (x1 + x2)
    midy = 0.5 * (y1 + y2)
    dist = math.hypot(dx, dy)
    if dist == 0:
        cp = (midx, midy)
    else:
        perp_x = -dy / dist
        perp_y = dx / dist
        cp = (midx + perp_x * rad * dist, midy + perp_y * rad * dist)

    ts = np.linspace(0.0, 1.0, nseg)
    pts = []
    for t in ts:
        x = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * cp[0] + t ** 2 * x2
        y = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * cp[1] + t ** 2 * y2
        pts.append((x, y))

    segments = [[pts[i], pts[i + 1]] for i in range(len(pts) - 1)]

    rgba1 = np.array(to_rgba(c1))
    rgba2 = np.array(to_rgba(c2))
    frac_ts = np.linspace(0.0, 1.0, len(segments))
    colors = [tuple(rgba1 * (1 - t) + rgba2 * t) for t in frac_ts]

    lc = LineCollection(segments, colors=colors, linewidths=lw, linestyles=linestyle, zorder=zorder, capstyle='round')
    ax.add_collection(lc)


def label_point_frac(x1, y1, x2, y2, f=0.75):
    """Return point at fraction f along line from (x1,y1) → (x2,y2)."""
    return (x1 + f * (x2 - x1), y1 + f * (y2 - y1))


def summarize_interaction_edges(df_full, min_fraction):
    """Return one chemically meaningful edge per residue/atom/type.

    The old renderer copied the residue-level contact fraction onto every
    framewise row.  That made a transient interaction look as persistent as
    the residue itself.  Here each percentage is the fraction of trajectory
    frames containing that specific residue--ligand atom--interaction type.
    """
    # The validated event schema is authoritative.  Legacy input remains
    # renderable for historical comparison, but never acquires a false
    # "validated" status merely by being read here.
    if "InteractionClass" in df_full.columns:
        status = df_full["ValidationStatus"].astype(str) if "ValidationStatus" in df_full else pd.Series("", index=df_full.index)
        df_full = df_full[status == "validated"].copy()
        df_full["Type"] = df_full["InteractionClass"]
        df_full["LigAtomName"] = df_full["LigandCanonicalAtom"]
    frame_column = "Frame" if "Frame" in df_full.columns else "Time_ns"
    required = {frame_column, "Residue", "Type", "LigAtomName"}
    missing = required.difference(df_full.columns)
    if missing:
        raise ValueError(f"InteractionTypes_Framewise.csv missing columns: {sorted(missing)}")
    # Occupancy must use every analyzed frame, including frames that produced
    # no interaction event.  Prefer 3A's explicit metadata; older event-only
    # CSVs are recovered from their regular sampled-frame grid rather than
    # silently using only frames that happen to contain an event.
    n_frames = None
    metadata_path = os.path.join(ANALYSIS_DIR, "InteractionFrameMetadata.json")
    try:
        with open(metadata_path, encoding="utf-8") as handle:
            n_frames = int(json.load(handle).get("total_analyzed_frames", 0)) or None
    except (OSError, ValueError, TypeError):
        pass
    if n_frames is None:
        frame_values = np.sort(df_full[frame_column].dropna().unique())
        if len(frame_values) >= 2:
            step = float(np.median(np.diff(frame_values)))
            inferred = int(round((float(frame_values[-1]) - float(frame_values[0])) / step)) + 1 if step > 0 else 0
            n_frames = max(len(frame_values), inferred)
        else:
            n_frames = len(frame_values)
    n_frames = max(1, int(n_frames))
    optional = ["ChainID"] if "ChainID" in df_full.columns else []
    work = df_full[list(required) + optional].copy()
    if "ChainID" in work.columns:
        chain = work["ChainID"].fillna("").astype(str).str.strip()
        # Chain is part of node identity, not merely its label: residues with
        # the same number on different PROTAC chains remain separate nodes.
        work["Residue"] = np.where(
            chain.ne(""), chain + "|" + work["Residue"].astype(str), work["Residue"]
        )
    work["LigAtomName"] = work["LigAtomName"].map(normalize_atom)
    work = work.dropna(subset=["Residue", "Type", "LigAtomName"])
    edges = (work.groupby(["Residue", "LigAtomName", "Type"], as_index=False)[frame_column]
                 .nunique()
                 .rename(columns={frame_column: "Frames"}))
    edges["Fraction"] = edges["Frames"] / float(n_frames)
    return edges[edges["Fraction"] >= min_fraction].sort_values(
        ["Fraction", "Residue", "LigAtomName"], ascending=[False, True, True]
    ).reset_index(drop=True)


def edge_style(interaction_type):
    """Visual convention used consistently by every ligand network panel."""
    value = str(interaction_type).lower()
    if "h-bond" in value or "hbond" in value:
        return EDGE_HBOND, 2.7, "--"
    if "hydrophobic" in value:
        return EDGE_HYDRO, 2.2, "--"
    if "water" in value:
        return EDGE_WATER, 2.3, "--"
    if "ionic" in value or "salt" in value:
        return "#8E5BBE", 2.5, "--"
    return "#777777", 1.8, "--"


def contact_residue_map(df_edges, fallback_map):
    """Use observed contact atoms to anchor residues around the ligand."""
    mapping = {}
    for residue, group in df_edges.groupby("Residue"):
        atoms = [a for a in group.sort_values("Fraction", ascending=False)["LigAtomName"].tolist() if a]
        mapping[residue] = list(dict.fromkeys(atoms))
    for residue, atom in fallback_map.items():
        mapping.setdefault(residue, [atom] if isinstance(atom, str) else atom)
    return mapping


def clear_edge_label_position(x1, y1, x2, y2, edge_index, occupied):
    """Choose a readable percentage-label position away from prior labels."""
    dx, dy = x2 - x1, y2 - y1
    length = max(math.hypot(dx, dy), 1e-6)
    px, py = -dy / length, dx / length
    # Alternating offsets keep bundles from the same atom legible.
    direction = 1 if edge_index % 2 == 0 else -1
    candidates = []
    for frac in (0.58, 0.68, 0.78, 0.88):
        bx, by = label_point_frac(x1, y1, x2, y2, frac)
        for offset in (0.0, 0.34, -0.34, 0.68, -0.68):
            candidates.append((bx + direction * offset * px, by + direction * offset * py))
    for x, y in candidates:
        if all(math.hypot(x - ox, y - oy) >= 0.62 for ox, oy in occupied):
            occupied.append((x, y))
            return x, y
    x, y = candidates[-1]
    occupied.append((x, y))
    return x, y


def arc3_label_point(x1, y1, x2, y2, rad, fraction=0.52):
    """Point on the quadratic approximation of Matplotlib's Arc3 connector."""
    mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
    dx, dy = x2 - x1, y2 - y1
    # Arc3's control point is a perpendicular displacement from the chord.
    cx, cy = mx + rad * dy, my - rad * dx
    t = fraction
    u = 1.0 - t
    return (u * u * x1 + 2.0 * u * t * cx + t * t * x2,
            u * u * y1 + 2.0 * u * t * cy + t * t * y2)


def classify_residue(res_name: str):
    # Extract the residue code from either TYR76 or A|TYR76.
    match = re.search(r"([A-Za-z]{3})(?=\s*\d)", str(res_name))
    code = match.group(1).upper() if match else "".join([c for c in str(res_name) if c.isalpha()]).upper()

    for cat, members in RESIDUE_CATEGORIES.items():
        if code in members:
            return cat
    return "polar"   # fallback


def auto_place_residues_no_overlap(residue_list, atom_positions,
                                   start_rx=4.2, start_ry=2.3,
                                   scale_factor=1.12):
    """
    Expands ellipse radii until no residue node circle overlaps ligand atoms.
    Returns:
        res_positions, final_rx, final_ry
    """
    rx, ry = start_rx, start_ry

    while True:
        # First: try layout
        n = len(residue_list)
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        res_positions = {}
        for res, t in zip(residue_list, angles):
            x = rx * math.cos(t)
            y = ry * math.sin(t)
            res_positions[res] = (x, y)

        # Check for overlaps
        overlap = False
        for (ax, ay) in atom_positions.values():
            for (rx0, ry0) in res_positions.values():
                dist = math.hypot(ax - rx0, ay - ry0)
                if dist < (RESIDUE_RADIUS + LIG_ATOM_RADIUS + 0.2):
                    overlap = True
                    break
            if overlap:
                break

        if not overlap:
            return res_positions, rx, ry

        # Expand ellipse
        rx *= scale_factor
        ry *= scale_factor



def extract_mol2_bond_orders(mol2_file):
    """
    Parse the MOL2 file directly and extract:
        atom_name → atom_name → bond_order (1,2,3)
    Aromatic bonds ('ar') are treated as alternating 1/2 bond pattern.
    """

    bonds = []
    in_bond_section = False

    with open(mol2_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("@<TRIPOS>BOND"):
                in_bond_section = True
                continue
            if line.startswith("@<TRIPOS>") and in_bond_section:
                break

            if in_bond_section and line:
                parts = line.split()
                if len(parts) < 4:
                    continue

                _, a1, a2, order_raw = parts[:4]

                # convert bond order
                if order_raw == "1":
                    order = 1
                elif order_raw == "2":
                    order = 2
                elif order_raw == "3":
                    order = 3
                elif order_raw in ("ar", "4"):
                    # assign aromatic alternating pattern later
                    order = "AROMATIC"
                else:
                    order = 1

                bonds.append((a1, a2, order))

    return bonds



def make_residue_legend(ax, interaction_types=None, residue_categories=None):
    """
    Unified legend for Panel C — identical style to Panel B.
    Includes:
      • Edge types (H-bond, Hydrophobic, Water bridge, Ionic)
      • Residue category colors (Hydrophobic, Positive, Negative, Polar)
    """

    style_by_type = {
        "H-bond": (EDGE_HBOND, 2.8), "Hydrophobic": (EDGE_HYDRO, 1.8),
        "Water bridge": (EDGE_WATER, 2.1), "Ionic": ("#C39BD3", 2.3),
    }
    interaction_types = interaction_types or list(style_by_type)
    edge_handles = [
        Line2D([0], [0], color=style_by_type[label][0], linestyle="--",
               linewidth=style_by_type[label][1], label=label)
        for label in style_by_type if label in interaction_types
    ]
    residue_categories = residue_categories or ["hydrophobic", "positive", "negative", "polar"]
    residue_handles = [
        Patch(facecolor=RESIDUE_COLORS[category], edgecolor="black", label=category.title())
        for category in ["hydrophobic", "positive", "negative", "polar"] if category in residue_categories
    ]

    handles = edge_handles + residue_handles

    ax.legend(
        handles=handles,
        fontsize=9,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=min(4, max(1, len(handles)))
    )



def find_active_residues(df_full, minfrac):
    """
    Return a sorted list of residues that have ≥1 interaction
    with Fraction >= minfrac.
    """
    active = (
        df_full[df_full["Fraction"] >= minfrac]["Residue"]
        .unique()
        .tolist()
    )
    active = sorted(active, key=residue_number)
    return active




def safe_kekulize_and_get_rings(mol, atom_names):
    """
    Try to kekulize the molecule and return:
        - valid_rings: list of rings (list of atom names)
        - bad_atoms: atoms RDKit said were not kekulizable

    Only rings fully composed of OK atoms will be returned.
    """
    m2 = Chem.Mol(mol)  # copy
    bad_atoms = set()

    try:
        Chem.SanitizeMol(
            m2,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                        Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        )
        print("✓ Kekulization succeeded — ring detection active.")
    except Exception as e:
        msg = str(e)
        print("⚠ Kekulization failed:", msg)

        # Extract atoms RDKit complains about
        # Format usually: "Can't kekulize mol. Unkekulized atoms: 6 7 8 9 10"
        m = re.search(r"atoms:\s*([\d\s]+)", msg)
        if m:
            idxs = m.group(1).split()
            bad_atoms = {atom_names[int(i)].strip() for i in idxs}
            print("⚠ Problematic atoms:", bad_atoms)

        # Continue but with no valid rings
        return [], bad_atoms

    # If kekulized: get rings
    ri = m2.GetRingInfo().AtomRings()

    valid_rings = []
    for ring in ri:
        # Convert indices → atom names
        ring_atoms = [atom_names[i] for i in ring]

        # Exclude rings containing bad atoms
        if any(a in bad_atoms for a in ring_atoms):
            continue

        valid_rings.append(ring_atoms)

    return valid_rings, bad_atoms



def fix_kekulization_errors(mol, bonds, atom_names):
    """
    Post-processing correction that:

    • Attempts a safe kekulization on a COPY of the mol
    • Restores proper bond orders (1/2/3) for aromatic / conjugated rings
    • Leaves atom names, indices, coordinates untouched
    • Updates the existing 'bonds' list IN PLACE

    Returns: corrected_bonds (new list of (atom1, atom2, order))
    """

    # Always work on a copy — do NOT touch original
    m2 = Chem.Mol(mol)

    # --- Attempt sanitization to assign proper bond orders ---
    try:
        Chem.SanitizeMol(
            m2,
            sanitizeOps = Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                          Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        )
        print("✓ Kekulization successful.")
    except Exception as e:
        print("⚠ Kekulization failed softly:", e)
        # Fallback: only aromaticity perception
        try:
            Chem.SanitizeMol(m2, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            print("✓ Aromaticity restored without kekulization.")
        except Exception as e2:
            print("⚠ Aromaticity restoration failed:", e2)
            return bonds   # keep original

    # --- Extract corrected orders ---
    corrected = []
    for b in m2.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()

        a1 = atom_names[i]
        a2 = atom_names[j]

        bt = b.GetBondType()
        if bt == Chem.BondType.SINGLE:
            order = 1
        elif bt == Chem.BondType.DOUBLE:
            order = 2
        elif bt == Chem.BondType.TRIPLE:
            order = 3
        else:
            # AROMATIC now has explicit Kekulé form → either 1 or 2
            order = 1

        corrected.append((a1, a2, order))

    print("✓ Bond orders updated using kekulized molecule.")
    return corrected





def draw_curved_edge(ax, x1, y1, x2, y2,
                     color="#000", lw=1.8, style="solid",
                     curvature=0.25, zorder=1):
    """
    Draws a quadratic Bezier curved edge between two points.
    Identical behavior to Panel B.
    """

    # Midpoint
    mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0

    # Perpendicular offset for curvature
    dx, dy = x2 - x1, y2 - y1
    nx, ny = -dy, dx
    norm = math.sqrt(nx*nx + ny*ny) or 1.0
    nx /= norm; ny /= norm

    cx = mx + curvature * nx
    cy = my + curvature * ny

    path = matplotlib.path.Path([
        (x1, y1),
        (cx, cy),
        (x2, y2)
    ], [matplotlib.path.Path.MOVETO,
        matplotlib.path.Path.CURVE3,
        matplotlib.path.Path.CURVE3])

    patch = matplotlib.patches.PathPatch(
        path,
        edgecolor=color,
        linewidth=lw,
        linestyle=style,
        fill=False,
        zorder=zorder
    )

    ax.add_patch(patch)



from rdkit.Chem.rdchem import BondType





def restore_aromaticity(mol):
    """
    Safely restore aromaticity perception *without* touching atom names,
    coordinates, or bond orders. Sanitization steps that would modify
    the molecule (e.g., kekulization, valence checking) are skipped.
    """

    params = Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | \
             Chem.SanitizeFlags.SANITIZE_KEKULIZE

    m2 = Chem.Mol(mol)  # shallow copy

    try:
        Chem.SanitizeMol(m2, sanitizeOps=params)
    except:
        # KEKULIZE fails for some CGenFF heteroaromatics — run only aromaticity
        try:
            Chem.SanitizeMol(
                m2,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            )
        except Exception as e:
            print("⚠ Aromaticity restoration failed:", e)

    return m2



def parse_mol2_bonds(mol2_path):
    """
    Reads @<TRIPOS>BOND entries from a MOL2 file and returns:
        [(atom_name1, atom_name2, bond_order), ...]
    where bond_order is integer 1/2/3 or "AR" for aromatic.
    """

    bonds = []
    atom_ids = []
    raw_bonds = []
    reading_atoms = False
    reading_bonds = False

    with open(mol2_path, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("@<TRIPOS>ATOM"):
                reading_atoms = True
                reading_bonds = False
                continue

            if line.startswith("@<TRIPOS>BOND"):
                reading_atoms = False
                reading_bonds = True
                continue

            if line.startswith("@<TRIPOS>"):
                reading_atoms = False
                reading_bonds = False
                continue

            # Parse atom block
            if reading_atoms and line:
                parts = line.split()
                # MOL2 atoms: index name x y z type subst subst_id charge
                idx = int(parts[0])
                atom_ids.append(idx)

            # Parse bond block
            if reading_bonds and line:
                parts = line.split()
                if len(parts) < 4:
                    continue

                _, a1_idx, a2_idx, order_raw = parts[:4]
                a1_idx = int(a1_idx)
                a2_idx = int(a2_idx)

                raw_bonds.append((a1_idx, a2_idx, order_raw))

    # Atom labels in CGenFF MOL2 files are often duplicated (for example,
    # several carbons are simply named ``C``).  Use exactly the same stable,
    # de-duplicated names used for coordinates, keyed by the source numeric
    # atom ID.  Using raw labels here was the cause of bond endpoints drifting
    # onto the wrong carbon in the rendered structure.
    names = mol2_atom_names_in_file_order(mol2_path)
    if len(atom_ids) != len(names):
        raise ValueError(f"MOL2 atom table could not be reconciled for {mol2_path}")
    atoms = dict(zip(atom_ids, names))

    for a1_idx, a2_idx, order_raw in raw_bonds:
        atom1 = atoms[a1_idx]
        atom2 = atoms[a2_idx]
        if order_raw == "1":
            order = 1
        elif order_raw == "2":
            order = 2
        elif order_raw == "3":
            order = 3
        elif order_raw in ("ar", "4"):
            order = "AR"
        else:
            order = 1
        bonds.append((atom1, atom2, order))

    return bonds




def normalize_bond_endpoint(a, atom_names):
    """
    Convert RDKit index -> atom name, or return atom name directly.
    """
    # Case 1: integer index
    if isinstance(a, int):
        return atom_names[a]

    # Case 2: string atom name (MOL2)
    if isinstance(a, str):
        return a

    print(f"[WARN] Unexpected bond endpoint type: {a}")
    return None




def find_aromatic_rings(mol):
    """
    Return rings where all atoms are aromatic.
    """
    aromatic_rings = []
    ri = mol.GetRingInfo().AtomRings()

    for ring in ri:
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings.append(ring)

    return aromatic_rings



def ring_center(positions, atom_indices):
    xs = [positions[i][0] for i in atom_indices]
    ys = [positions[i][1] for i in atom_indices]
    return (np.mean(xs), np.mean(ys))


def draw_bond(ax, x1, y1, x2, y2, order, elem1="C", elem2="C", a1_name=None, a2_name=None):
    """
    Draw a single, double, or triple bond with correct parallel offsets.
    Quiet, robust: handle None names gracefully and only draw valid geometry.
    """
    dx = x2 - x1
    dy = y2 - y1
    L = math.hypot(dx, dy)
    if L == 0:
        return

    # Unit / perpendicular
    ux = dx / L
    uy = dy / L
    px = -uy
    py = ux

    # Default styling
    base_offset = 0.18
    # Draw ligand bonds below ligand atoms: use low bond_zorder
    bond_zorder = 1.0
    lw = 3.6

    # determine endpoint colors from element types
    c1 = _atom_node_color(elem1)
    c2 = _atom_node_color(elem2)

    def plot_split_line(xa, ya, xb, yb, col_a, col_b, linewidth=lw, z=bond_zorder):
        # split the line exactly at midpoint and draw two colored halves
        mx = 0.5 * (xa + xb)
        my = 0.5 * (ya + yb)
        ax.plot([xa, mx], [ya, my], color=col_a, linewidth=linewidth, zorder=z)
        ax.plot([mx, xb], [my, yb], color=col_b, linewidth=linewidth, zorder=z)

    # For triple bonds: increase spacing and make outer lines more prominent
    if order == 3:
        offset = 0.34
        outer_lw = 3.4
        mid_lw = 3.6
        # middle line
        plot_split_line(x1, y1, x2, y2, c1, c2, linewidth=mid_lw, z=bond_zorder)
        # outer parallel lines
        plot_split_line(x1 + px*offset, y1 + py*offset, x2 + px*offset, y2 + py*offset, c1, c2, linewidth=outer_lw, z=bond_zorder)
        plot_split_line(x1 - px*offset, y1 - py*offset, x2 - px*offset, y2 - py*offset, c1, c2, linewidth=outer_lw, z=bond_zorder)
        return

    # Double bond
    if order == 2:
        offset = base_offset
        plot_split_line(x1 + px*offset, y1 + py*offset, x2 + px*offset, y2 + py*offset, c1, c2, linewidth=lw, z=bond_zorder)
        plot_split_line(x1 - px*offset, y1 - py*offset, x2 - px*offset, y2 - py*offset, c1, c2, linewidth=lw, z=bond_zorder)
        return

    # Single bond
    if order == 1:
        plot_split_line(x1, y1, x2, y2, c1, c2, linewidth=lw, z=bond_zorder)
        return

    # fallback
    plot_split_line(x1, y1, x2, y2, c1, c2, linewidth=lw, z=bond_zorder)


def draw_sid_bond(ax, x1, y1, x2, y2, order):
    """Draw a conventional 2D chemical bond for the SID centrepiece.

    Atom spheres work for a network graph but obscure the chemical silhouette.
    This deliberately follows the standard structure-drawing convention:
    carbon is implicit and bond order carries the visual information.
    """
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy)
    if length < 1e-8:
        return
    px, py = -dy / length, dx / length
    # Deliberately heavier than a default molecular drawing: the SID is
    # commonly viewed beside a contact plot or after PDF downsampling.
    color, linewidth = "#202020", 2.15
    # A conventional bond is a main line plus a short, inset companion line
    # for double/triple order.  Two full-width heavy rails created the black
    # knot that made aromatic systems unreadable in the prior render.
    ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth,
            solid_capstyle="round", zorder=3)
    inset = min(0.24, length * 0.20)
    ux, uy = dx / length, dy / length
    if order == 2:
        ax.plot([x1 + ux * inset + px * 0.40, x2 - ux * inset + px * 0.40],
                [y1 + uy * inset + py * 0.40, y2 - uy * inset + py * 0.40],
                color=color, linewidth=1.70, solid_capstyle="round", zorder=3)
    elif order == 3:
        for offset in (-0.40, 0.40):
            ax.plot([x1 + ux * inset + px * offset, x2 - ux * inset + px * offset],
                    [y1 + uy * inset + py * offset, y2 - uy * inset + py * offset],
                    color=color, linewidth=1.55, solid_capstyle="round", zorder=3)


def draw_sid_atom_label(ax, element, x, y):
    """Draw only explicit non-carbon atoms, as in a 2D ligand diagram."""
    colors = {"O": "#E53935", "N": "#1565C0", "S": "#C69200", "P": "#D81B60"}
    color = colors.get(element, "#087F5B" if element in HALOGEN_ELEMENTS else "#303030")
    ax.text(x, y, element, color=color, fontsize=16, weight="medium",
            ha="center", va="center", zorder=7,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.97))



# --------------------------------------------------------------------
# Panel A — stacked bar plot
# --------------------------------------------------------------------
def draw_contact_bar_panel(ax, df, total_frames, minfrac):
    """
    Draw stacked bar plot of interaction fractions (Panel a).
    Uses InteractionTypes_Summary.csv.
    """
    df_filt = df[df["Fraction"] >= minfrac].copy()
    if df_filt.empty:
        raise ValueError(f"No residues exceed minfrac={minfrac} for bar plot")

    df_filt = df_filt.sort_values("Residue", key=lambda s: s.map(residue_number))

    # Identify interaction columns present
    type_cols = []
    if "H-bond" in df_filt.columns:
        type_cols.append("H-bond")
    if "Hydrophobic" in df_filt.columns:
        type_cols.append("Hydrophobic")
    if "Ionic" in df_filt.columns:
        type_cols.append("Ionic")
    if "Water bridge" in df_filt.columns:
        type_cols.append("Water bridge")

    if not type_cols:
        raise ValueError("No known interaction type columns found for bar plot")

    frac_data = {}
    for col in type_cols:
        frac_data[col] = df_filt[col].fillna(0) / float(total_frames)

    residues = df_filt["Residue"].tolist()
    x = np.arange(len(residues))
    bottom = np.zeros(len(residues))

    color_map = {
        "H-bond":      (BAR_HBOND, "H-bonds"),
        "Hydrophobic": (BAR_HYDRO, "Hydrophobic"),
        "Ionic":       (BAR_IONIC, "Ionic"),
        "Water bridge": (BAR_WATER, "Water bridges"),
    }

    for col in type_cols:
        frac_vals = frac_data[col].values
        c, label = color_map.get(col, ("#999999", col))
        ax.bar(x, frac_vals, bottom=bottom,
               color=c, edgecolor="black", linewidth=0.6, label=label)
        bottom += frac_vals

    ax.set_xticks(x)
    ax.set_xticklabels(residues, rotation=60, ha="right", fontsize=7.5)
    ax.set_ylabel("Interaction fraction", fontsize=9)
    ax.set_title("Protein–Ligand Contacts", fontsize=12, weight="bold", pad=6)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend on the right, small
    handles = [
        Patch(facecolor=BAR_HBOND, edgecolor="black", label="H-bonds"),
        Patch(facecolor=BAR_HYDRO, edgecolor="black", label="Hydrophobic"),
        Patch(facecolor=BAR_IONIC, edgecolor="black", label="Ionic"),
        Patch(facecolor=BAR_WATER, edgecolor="black", label="Water bridges"),
    ]

    # Compact legend placed slightly left of center above the title so it doesn't overlap.
    # Move x a small amount (0.44) to shift the legend a few pixels left in axes coords.
    ax.legend(
        handles=handles,
        fontsize=7.5,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.84, 1.12),   # small left shift from center
        ncol=2,
        handletextpad=0.2,
        columnspacing=0.6,
        borderaxespad=0.2
    )

# --------------------------------------------------------------------
# Panel B — network diagram
# --------------------------------------------------------------------
def draw_network_panel(ax, ligand, df_filt, residue_to_ligatom):
    """
    Schrödinger-style ligand–protein interaction network.
    Multi-edge + hydrogen-safe + LP-mapping + aromatic rings preserved.

    Features:
      ✔ Multi-atom edges per residue (NO COLLAPSE)
      ✔ Hydrogen → heavy-atom remapping (prefix logic)
      ✔ LP1/LP2 → heavy-atom fallback + halogen fallback
      ✔ No center-collapsed edges unless truly unresolved
      ✔ Multiple arcs per residue with alternating curvature
      ✔ Correct % labels for EACH interaction
    """

    # ---------------------------------------------------------
    # 1. Load ligand layout (atom positions, bonds, parents)
    # ---------------------------------------------------------
    atom_positions, atom_elements, bonds, center, parent_map, mol, atom_names = \
        get_ligand_layout(ligand)
    rings, bad_atoms = safe_kekulize_and_get_rings(mol, atom_names)

    print("Valid rings:", rings)
    print("Skipped (bad) atoms:", bad_atoms)

    # ---------------------------------------------------------
    # Optional safety: if ligand atom layout is crowded (many atoms
    # overlapping visually), gently scale ligand atom positions away
    # from the ligand center so atoms/bonds and residue edges don't
    # occlude each other. This preserves original MOL2 coords but
    # improves figure readability when RDKit produced compact layouts.
    # ---------------------------------------------------------
    try:
        # consider only non-hydrogen atoms for crowding metric
        non_h_coords = [pos for nm, pos in atom_positions.items() if not nm.startswith("H")]
        if len(non_h_coords) > 1:
            min_d = float("inf")
            for i in range(len(non_h_coords)):
                xi, yi = non_h_coords[i]
                for j in range(i+1, len(non_h_coords)):
                    xj, yj = non_h_coords[j]
                    d = math.hypot(xi - xj, yi - yj)
                    if d < min_d:
                        min_d = d
            # if minimum inter-atom distance is small, scale outward
            if min_d < 0.9:  # threshold (Angstrom-like units from RDKit 2D coords)
                cx, cy = center
                # compute a conservative scale: either module default or enough
                # to push min_d to ~1.1
                safe_scale = max(LIGAND_LAYOUT_SCALE, 1.1 / max(min_d, 1e-6))
                for nm, (x, y) in list(atom_positions.items()):
                    dx = x - cx
                    dy = y - cy
                    atom_positions[nm] = (cx + dx * safe_scale, cy + dy * safe_scale)
    except Exception:
        # non-fatal: if something goes wrong here, keep original positions
        pass


    # ---------------------------------------------------------
    # 2. Fix LP atoms → nearest heavy atoms
    # ---------------------------------------------------------
    residue_to_ligatom = resolve_lp_mapping(ligand, residue_to_ligatom)

    # ---------------------------------------------------------
    # 3. Fallback unresolved atoms → halogen closest to center
    # ---------------------------------------------------------
    residue_to_ligatom = fix_edges_to_center_safe(
        residue_to_ligatom, atom_positions, atom_elements, center
    )

    # Prepare canvas
    ax.set_aspect("equal")
    ax.axis("off")

    # Draw ligand bonds earlier/later no longer required — zorder handles stacking.
    # ---------------------------------------------------------
    # 4. Draw ligand bonds (1/2/3-line rendering)
    # ---------------------------------------------------------
    for a1_idx, a2_idx, order in bonds:

        a1 = normalize_bond_endpoint(a1_idx, atom_names)
        a2 = normalize_bond_endpoint(a2_idx, atom_names)

        if a1 is None or a2 is None:
            continue

        # Skip hydrogens
        if a1.startswith("H") or a2.startswith("H"):
            continue

        x1, y1 = atom_positions[a1]
        x2, y2 = atom_positions[a2]

        elem1 = atom_elements.get(a1, "C")
        elem2 = atom_elements.get(a2, "C")

        draw_bond(ax, x1, y1, x2, y2, order, elem1, elem2)



    # ---------------------------------------------------------
    # 5. Aromatic rings
    # ---------------------------------------------------------
    # try:
    #     aromatic_rings = find_aromatic_rings(mol)
    #     for ring in aromatic_rings:
    #         ring_names = [atom_names[i] for i in ring]
    #         coords = [atom_positions[n] for n in ring_names]

    #         cx = np.mean([c[0] for c in coords])
    #         cy = np.mean([c[1] for c in coords])

    #         distances = [
    #             math.hypot(coords[i][0] - coords[(i+1) % len(coords)][0],
    #                        coords[i][1] - coords[(i+1) % len(coords)][1])
    #             for i in range(len(coords))
    #         ]
    #         r = np.mean(distances) * 0.38

    #         ax.add_patch(
    #             plt.Circle((cx, cy), r, fill=False, linewidth=3,
    #                        color="#555", zorder=1.5)
    #         )
    # except Exception as e:
    #     print("Aromatic ring drawing failed:", e)

    # ---------------------------------------------------------
    # 6. Draw ligand atom spheres (increase zorder so atoms are on top)
    # ---------------------------------------------------------
    for nm, (x, y) in atom_positions.items():

        # 🔥 HIDE ALL HYDROGEN ATOM NODES
        if nm.startswith("H"):
            continue

        elem = atom_elements.get(nm, "C")
        # Node color: carbon -> black, otherwise element-specific highlights
        if elem == "C":
            node_col = "#000000"
        elif elem == "S":
            node_col = "#FFD700"  # sulfur -> yellow
        elif elem == "P":
            node_col = "#FF69B4"  # phosphorus -> pink
        elif elem in ("F", "Cl", "Br", "I"):
            node_col = "#007849"  # halogens -> Miami green
        elif elem == "O":
            node_col = "#D72638"  # oxygen -> fire/blood red
        elif elem == "N":
            node_col = "#1E90FF"  # nitrogen -> cool blue
        else:
            node_col = ATOM_COLORS.get(elem, "#888")

        # Label: show element symbol for non-carbon heavy atoms.
        if elem == "C":
            label = None
        else:
            label = elem

        # Label color scheme: phosphorus & oxygen -> Miami green, others -> Miami orange
        # All labels should be white
        label_color = "white"
    # atom spheres must sit above bonds/edges
        ax.scatter(
            x, y, s=900, color=node_col,
            edgecolors="black", linewidth=2.2,
            zorder=4
        )
        # add element label for non-carbon heavy atoms
        if label:
            ax.text(x, y, label, fontsize=12, weight="bold",
                    ha="center", va="center", color=label_color, zorder=5)
    # ---------------------------------------------------------
    # 7. Place residues around ligand with automatic no-overlap
    # ---------------------------------------------------------
    residues = df_filt["Residue"].tolist()
    res_positions = auto_place_residues_guided(
        residues, residue_to_ligatom, atom_positions, center
    )

    # ---------------------------------------------------------
    # 8. Draw EDGES — MULTI-LIGATOM SUPPORT
    # ---------------------------------------------------------
    drawn_residues = set()   # collect residues that actually get edges

    for e_idx, (_, row) in enumerate(df_filt.iterrows()):

        res = row["Residue"]
        if res not in res_positions:
            continue

        xR, yR = res_positions[res]

        # Multiple ligand atoms are allowed
        lig_atoms = residue_to_ligatom.get(res, [])
        if isinstance(lig_atoms, str):
            lig_atoms = [lig_atoms]

        # ---------------------------
        # Determine interaction type
        # ---------------------------
        has_hbond = bool(row.get("HasHbond", False))
        hydro = row.get("Hydrophobic", 0) or 0
        water = row.get("Water bridge", 0) or 0
        ionic = row.get("Ionic", 0) or 0

        if has_hbond:
            col, lw = EDGE_HBOND, 2.8
        elif hydro > 0:
            col, lw = EDGE_HYDRO, 1.6
        elif water > 0:
            col, lw = EDGE_WATER, 2.0
        elif ionic > 0:
            col, lw = "#C39BD3", 2.2
        else:
            continue

        # Draw each atom → residue interaction separately
        for lig_atom in lig_atoms:

            # ---------------------------
            # Fix ligand atom coordinate — map hydrogens to heavy atoms
            # ---------------------------
            heavy = _resolve_h_to_heavy(lig_atom, parent_map, atom_positions, atom_elements, center)
            if heavy and heavy in atom_positions:
                xL, yL = atom_positions[heavy]
            else:
                # Unresolved ligand atom — skip this edge to avoid stray center edges
                print(f"⚠ Skipping unresolved ligand atom '{lig_atom}' for residue {res}")
                continue

            # ---------------------------
            # Compute arc geometry
            # ---------------------------
            sx, sy = shorten_to_circle(xL, yL, xR, yR, LIG_ATOM_RADIUS)
            ex, ey = shorten_to_circle(xR, yR, xL, yL, RESIDUE_RADIUS)

            rad = 0.18 if e_idx % 2 == 0 else -0.18

            # draw bicolor curved edge: ligand atom color -> residue node color
            lig_elem = atom_elements.get(heavy, 'C')
            lig_node_col = _atom_node_color(lig_elem)
            res_node_col = RESIDUE_COLORS.get(classify_residue(res), '#999999')
            _draw_bicolor_arc(ax, sx, sy, ex, ey, rad, lig_node_col, res_node_col, lw=lw, linestyle='--', zorder=2)

            drawn_residues.add(res)   # mark residue as having an edge

            # Label (%) — compute pct from the row BEFORE using it
            pct = row.get("Fraction", 0.0) * 100.0
            lx, ly = label_point_frac(sx, sy, ex, ey, 0.75)
            ax.text(
                lx, ly, f"{pct:.0f}%",
                fontsize=9, weight="bold",
                ha="center", va="center", color=col,
                zorder=3
            )

    # ---------------------------------------------------------
    # 8b. Re-draw ligand atom markers with higher zorder so ligand nodes are visually in front of edges.
    for nm, (x, y) in atom_positions.items():
        if nm.startswith("H"):
            continue
        elem = atom_elements.get(nm, "C")
        UM_ORANGE = "#FF6A00"
        UM_GREEN = "#007849"

        # choose node color as above (element-specific overrides)
        if elem == "C":
            node_col = "#000000"
        elif elem == "S":
            node_col = "#FFD700"
        elif elem == "P":
            node_col = "#FF69B4"
        elif elem in ("F", "Cl", "Br", "I"):
            node_col = "#007849"
        elif elem == "O":
            node_col = "#D72638"
        elif elem == "N":
            node_col = "#1E90FF"
        else:
            node_col = ATOM_COLORS.get(elem, "#888")

        if elem == "C":
            label = None
        else:
            label = elem

        # all labels white
        label_color = "white"

        ax.scatter(
            x, y, s=360, color=node_col,
            edgecolors="black", linewidth=1.8,
            zorder=5
        )
        if label:
            ax.text(x, y, label, fontsize=9, weight="bold",
                    ha="center", va="center", color=label_color, zorder=6)

    # ---------------------------------------------------------
    # 9. Draw residue circles + labels ONLY for residues with edges
    # ---------------------------------------------------------
    for res in sorted(drawn_residues, key=residue_number):
        x, y = res_positions[res]
        cat = classify_residue(res)
        col = RESIDUE_COLORS[cat]

        # residue nodes should be above everything else
        ax.scatter(
            x, y, s=3000, color=col,
            edgecolors="black", linewidth=2.0,
            zorder=6
        )
        ax.text(x, y, res,
                fontsize=11, weight="bold",
                ha="center", va="center",
                zorder=7)



def infer_mode_from_chainmap(run_dir="."):
    candidates = [
        os.path.join(run_dir, ".chainmap.peptide.json"),
        os.path.join(run_dir, "Analysis_Results", ".chainmap.peptide.json"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            with open(path, "r") as fh:
                data = json.load(fh)
            chain_id = str(data.get("chain_id", "")).strip()
            if not chain_id:
                raise ValueError(f"Peptide chainmap exists but missing chain_id: {path}")
            return {
                "is_peptide": True,
                "path": path,
                "chain_id": chain_id,
                "name": str(data.get("name", "Peptide")).strip(),
                "group_name": str(data.get("group_name", "PEPTIDE")).strip(),
                "atom_range": str(data.get("atom_range", "")).strip(),
            }
    return {"is_peptide": False}


def candidate_topology_files_networx():
    return [
        "Final_Trajectory.pdb",
        "binding_pocket_only.pdb",
        "protein.pdb",
    ] + sorted([f for f in os.listdir('.') if f.lower().endswith('.pdb')])


def parse_pdb_atom_table_networx(pdb_path):
    rows = []
    idx0 = 0
    with open(pdb_path, "r") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip()
            resid_raw = line[22:26].strip()
            try:
                resid = int(re.sub(r"[^\d\-]", "", resid_raw))
            except Exception:
                resid = None
            rows.append({
                "index0": idx0,
                "atom_name": atom_name,
                "resname": resname,
                "chain": chain,
                "resid": resid,
                "reslabel": f"{resname}{resid}" if resid is not None else f"{resname}{idx0}",
            })
            idx0 += 1
    return pd.DataFrame(rows)


def get_topology_atom_table_networx():
    seen = set()
    for path in candidate_topology_files_networx():
        if path in seen:
            continue
        seen.add(path)
        if os.path.exists(path):
            df = parse_pdb_atom_table_networx(path)
            if not df.empty:
                print(f"🔎 Using topology for peptide atom mapping: {path}")
                return df, path
    raise FileNotFoundError("Could not find a suitable PDB topology file to map peptide atom indices.")


def parse_atom_range_networx(atom_range):
    if not atom_range or "-" not in atom_range:
        raise ValueError(f"Invalid atom_range in chainmap: {atom_range!r}")
    s, e = atom_range.split("-", 1)
    return int(s), int(e)


def build_peptide_index_to_residue_map_networx(chainmap):
    topo_df, topo_path = get_topology_atom_table_networx()
    start1, end1 = parse_atom_range_networx(chainmap["atom_range"])
    start0 = start1 - 1
    end0 = end1 - 1

    sub = topo_df[(topo_df["index0"] >= start0) & (topo_df["index0"] <= end0)].copy()
    if sub.empty:
        raise ValueError(f"No topology atoms overlap peptide atom_range {chainmap['atom_range']} using {topo_path}")

    if chainmap.get("chain_id"):
        sub_chain = sub[sub["chain"] == chainmap["chain_id"]]
        if not sub_chain.empty:
            sub = sub_chain

    return {int(row.index0): row.reslabel for row in sub.itertuples(index=False)}


def infer_peptide_pair_summary_networx(df_full, chainmap):
    if "LigAtomIndex" not in df_full.columns:
        raise ValueError("InteractionTypes_Framewise.csv missing LigAtomIndex, needed for peptide mode.")

    idx_to_pep_res = build_peptide_index_to_residue_map_networx(chainmap)
    dfx = df_full.copy()
    dfx["LigAtomIndex"] = pd.to_numeric(dfx["LigAtomIndex"], errors="coerce")
    dfx = dfx.dropna(subset=["LigAtomIndex"]).copy()
    dfx["LigAtomIndex"] = dfx["LigAtomIndex"].astype(int)
    dfx["PeptideResidue"] = dfx["LigAtomIndex"].map(idx_to_pep_res)
    dfx = dfx.dropna(subset=["PeptideResidue"]).copy()

    if dfx.empty:
        raise ValueError("No framewise rows could be mapped from LigAtomIndex to peptide residues.")

    n_frames = int(dfx["Time_ns"].nunique())
    pair_summary = (dfx.groupby(["PeptideResidue", "Residue", "Type"])["Time_ns"]
                    .nunique()
                    .reset_index(name="Frames"))
    pair_summary["Fraction"] = pair_summary["Frames"] / float(n_frames)

    pair_frac = (dfx.groupby(["PeptideResidue", "Residue"])["Time_ns"]
                 .nunique()
                 .reset_index(name="Frames"))
    pair_frac["Fraction"] = pair_frac["Frames"] / float(n_frames)
    return dfx, pair_summary, pair_frac, n_frames


def draw_peptide_network_networx(ax, pair_summary, pair_frac, minfrac, no_edges=False, title=None):
    pair_frac = pair_frac[pair_frac["Fraction"] >= minfrac].copy()
    if pair_frac.empty:
        raise ValueError(f"No peptide-residue pairs exceed minfrac={minfrac}")

    peptide_res = sorted(pair_frac["PeptideResidue"].unique(), key=residue_number)
    receptor_res = sorted(pair_frac["Residue"].unique(), key=residue_number)

    pep_weight = pair_frac.groupby("PeptideResidue")["Fraction"].sum().to_dict()
    rec_weight = pair_frac.groupby("Residue")["Fraction"].sum().to_dict()
    peptide_res = sorted(peptide_res, key=lambda r: (-pep_weight.get(r, 0), residue_number(r)))
    receptor_res = sorted(receptor_res, key=lambda r: (-rec_weight.get(r, 0), residue_number(r)))

    ax.set_aspect("equal")
    ax.axis("off")

    left_x, right_x = -5.5, 5.5
    y_left = np.linspace(8.5, -8.5, len(peptide_res)) if peptide_res else []
    y_right = np.linspace(8.5, -8.5, len(receptor_res)) if receptor_res else []
    pep_pos = {res: (left_x, y) for res, y in zip(peptide_res, y_left)}
    rec_pos = {res: (right_x, y) for res, y in zip(receptor_res, y_right)}

    if not no_edges:
        for i, row in pair_summary.iterrows():
            pres, rres = row["PeptideResidue"], row["Residue"]
            frac = float(row["Fraction"])
            if frac < minfrac or pres not in pep_pos or rres not in rec_pos:
                continue
            x1, y1 = pep_pos[pres]
            x2, y2 = rec_pos[rres]
            itype = str(row["Type"])
            if "H-bond" in itype:
                col, lw, style = EDGE_HBOND, 2.4, "--"
            elif "Hydrophobic" in itype:
                col, lw, style = EDGE_HYDRO, 1.8, "--"
            elif "Water" in itype:
                col, lw, style = EDGE_WATER, 1.8, "--"
            elif "Ionic" in itype:
                col, lw, style = "#C39BD3", 2.1, "--"
            else:
                col, lw, style = "#999999", 1.4, "--"
            rad = 0.12 if (i % 2 == 0) else -0.12
            edge = FancyArrowPatch((x1 + 0.55, y1), (x2 - 0.55, y2),
                                   arrowstyle="-", connectionstyle=f"arc3,rad={rad}",
                                   linewidth=lw, linestyle=style, color=col, alpha=0.80, zorder=1)
            ax.add_patch(edge)

        for _, row in pair_frac.iterrows():
            pres, rres = row["PeptideResidue"], row["Residue"]
            frac = float(row["Fraction"])
            if pres not in pep_pos or rres not in rec_pos:
                continue
            x1, y1 = pep_pos[pres]
            x2, y2 = rec_pos[rres]
            lx, ly = label_point_frac(x1 + 0.55, y1, x2 - 0.55, y2, 0.55)
            ax.text(lx, ly, f"{frac*100:.0f}%", fontsize=7.5, color="black",
                    ha="center", va="center", zorder=2,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))

    for res, (x, y) in pep_pos.items():
        ax.scatter(x, y, s=2500, color=RESIDUE_COLORS["peptide"], edgecolors="black", linewidth=2.0, zorder=3)
        ax.text(x, y, res, fontsize=9.5, weight="bold", ha="center", va="center", zorder=4)

    for res, (x, y) in rec_pos.items():
        cat = classify_residue(res)
        col = RESIDUE_COLORS.get(cat, "#E8E8E8")
        ax.scatter(x, y, s=2500, color=col, edgecolors="black", linewidth=2.0, zorder=3)
        ax.text(x, y, res, fontsize=9.5, weight="bold", ha="center", va="center", zorder=4)

    ax.text(left_x, 10.0, "Peptide residues", fontsize=12, weight="bold", ha="center")
    ax.text(right_x, 10.0, "Protein residues", fontsize=12, weight="bold", ha="center")
    ax.set_title(title or "Peptide–Protein Residue Interaction Network", fontsize=13, weight="bold")

    handles = [
        Line2D([0],[0], color=EDGE_HBOND, linestyle="--", linewidth=2.4, label="H-bond"),
        Line2D([0],[0], color=EDGE_HYDRO, linestyle="--", linewidth=1.8, label="Hydrophobic"),
        Line2D([0],[0], color=EDGE_WATER, linestyle="--", linewidth=1.8, label="Water bridge"),
        Line2D([0],[0], color="#C39BD3", linestyle="--", linewidth=2.1, label="Ionic"),
        Patch(facecolor=RESIDUE_COLORS["peptide"], edgecolor="black", label="Peptide residue"),
        Patch(facecolor=RESIDUE_COLORS["hydrophobic"], edgecolor="black", label="Hydrophobic residue"),
        Patch(facecolor=RESIDUE_COLORS["positive"], edgecolor="black", label="Positive residue"),
        Patch(facecolor=RESIDUE_COLORS["negative"], edgecolor="black", label="Negative residue"),
        Patch(facecolor=RESIDUE_COLORS["polar"], edgecolor="black", label="Polar residue"),
    ]
    ax.legend(handles=handles, fontsize=8, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.02), ncol=5)


def save_peptide_network_networx(name, chainmap, df_summary, minfrac, no_edges=False, outdir=None):
    df_full = pd.read_csv(os.path.join(ANALYSIS_DIR, "InteractionTypes_Framewise.csv"))
    dfx, pair_summary, pair_frac, n_frames = infer_peptide_pair_summary_networx(df_full, chainmap)
    if outdir is None:
        outdir = os.path.join(ANALYSIS_DIR, "NETWORX")
    path = os.path.join(outdir, f"{name}_peptide_network.png")
    fig, ax = plt.subplots(figsize=(14, 11))
    draw_peptide_network_networx(ax, pair_summary, pair_frac, minfrac, no_edges=no_edges,
                         title=f"{chainmap.get('name', 'Peptide')}–Protein Residue Interaction Network")
    fig.savefig(path, dpi=700, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Peptide network saved → {path}")
    return path

def build_two_panel_figure(ligand: str, minfrac: float, output: str, no_edges: bool=False, mode: str="auto"):
    """
    Build figures. Auto-detect peptide mode from .chainmap.peptide.json unless forced.
    Keeps the original ligand NETWORX functionality intact, and swaps to a peptide
    residue-network only when the peptide sidecar exists.
    """
    outdir = os.path.join(ANALYSIS_DIR, "NETWORX")
    os.makedirs(outdir, exist_ok=True)

    df_summary, total_frames = load_interaction_summary()
    chainmap = infer_mode_from_chainmap(".")
    is_peptide = chainmap.get("is_peptide", False)
    if mode == "ligand":
        is_peptide = False
    elif mode == "peptide":
        is_peptide = True
        if not chainmap.get("is_peptide", False):
            raise RuntimeError("Forced peptide mode, but .chainmap.peptide.json was not found.")

    bar_out = make_bar_figure(ligand, df_summary, total_frames, minfrac, outdir, dpi=900)

    if is_peptide:
        print(f"\n🧬 Auto-detected peptide mode from {chainmap['path']}")
        print(f"   • chain_id   = {chainmap['chain_id']}")
        print(f"   • name       = {chainmap['name']}")
        print(f"   • group_name = {chainmap['group_name']}")
        print(f"   • atom_range = {chainmap['atom_range']}")

        net_out = save_peptide_network_networx(ligand, chainmap, df_summary, minfrac, no_edges=no_edges, outdir=outdir)
        exp_out = net_out
        clean_out = net_out
    else:
        residue_to_ligatom = load_residue_ligand_atom_map()
        atom_positions, atom_elements, bonds, center, parent_map, mol, atom_names =             get_ligand_layout(ligand)
        residue_to_ligatom_fixed = resolve_lp_mapping(ligand, residue_to_ligatom)
        residue_to_ligatom_fixed = fix_edges_to_center_safe(
            residue_to_ligatom_fixed, atom_positions, atom_elements, center
        )
        n_fixed = sum(1 for r in residue_to_ligatom if residue_to_ligatom.get(r) != residue_to_ligatom_fixed.get(r))
        print(f"\n🔁 Applied ligand mapping fixes — {n_fixed} residues remapped (persisted for all panels)\n")
        residue_to_ligatom = residue_to_ligatom_fixed

        df_filt = df_summary[df_summary["Fraction"] >= minfrac].copy()
        if df_filt.empty:
            raise ValueError(f"No residues exceed minfrac={minfrac}")

        framewise_path = os.path.join(ANALYSIS_DIR, "ValidatedInteractionEvents.csv")
        if os.path.exists(framewise_path):
            df_full = pd.read_csv(framewise_path)
            status = df_full["ValidationStatus"].astype(str) if "ValidationStatus" in df_full else pd.Series("", index=df_full.index)
            df_full = df_full[status == "validated"].copy()
            df_full["Type"] = df_full["InteractionClass"]
            df_full["LigAtomName"] = df_full["LigandCanonicalAtom"]
        else:
            framewise_path = os.path.join(ANALYSIS_DIR, "InteractionTypes_Framewise.csv")
            print("⚠️ Rendering legacy InteractionTypes framewise data; labels are unvalidated.")
            df_full = pd.read_csv(framewise_path)
        required_cols = ["Time_ns", "Residue", "Type", "LigAtomName"]
        for col in required_cols:
            if col not in df_full.columns:
                raise ValueError(f"InteractionTypes_Framewise.csv missing required column: {col}")
        df_full["LigAtomName"] = df_full["LigAtomName"].apply(normalize_atom)

        n_frames = df_full["Time_ns"].nunique()
        res_frame_counts = df_full.groupby("Residue")["Time_ns"].nunique()
        frac_map = (res_frame_counts / float(n_frames)).to_dict()
        df_full["Fraction"] = df_full["Residue"].map(frac_map).fillna(0.0)

        global MIN_FRACTION
        MIN_FRACTION = minfrac

        net_out = make_collapsed_network_figure(ligand, df_filt, residue_to_ligatom, outdir, dpi=900)
        rings, bad_atoms = safe_kekulize_and_get_rings(mol, atom_names)
        print("Valid rings:", rings)
        print("Skipped (bad) atoms:", bad_atoms)

        exp_out = make_expanded_network_figure(
            ligand,
            df_full=df_full,
            residue_to_ligatom=residue_to_ligatom,
            atom_positions=atom_positions,
            atom_elements=atom_elements,
            bonds=bonds,
            mol=mol,
            atom_names=atom_names,
            center=center,
            parent_map=parent_map,
            outdir=outdir,
            show_edges=not no_edges,
            dpi=900
        )

        clean_out = os.path.join(outdir, f"{ligand}_expanded_clean.png")
        make_expanded_network_figure(
            ligand,
            df_full=df_full,
            residue_to_ligatom=residue_to_ligatom,
            atom_positions=atom_positions,
            atom_elements=atom_elements,
            bonds=bonds,
            mol=mol,
            atom_names=atom_names,
            center=center,
            parent_map=parent_map,
            outdir=outdir,
            show_edges="single" if not no_edges else False,
            dpi=900,
            out_name=clean_out,
        )
        print(f"  ✓ Clean expanded network saved → {clean_out}")

    combined_out = os.path.join(outdir, f"{ligand}_combined.pdf")
    figD = plt.figure(figsize=(12, 16))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.8, 1.2], hspace=0.25)
    axD1 = figD.add_subplot(gs[0])
    axD1.imshow(plt.imread(bar_out))
    axD1.axis("off")
    axD2 = figD.add_subplot(gs[1])
    axD2.imshow(plt.imread(exp_out))
    axD2.axis("off")
    figD.savefig(combined_out, dpi=600)
    plt.close(figD)

    print("\n===================================================")
    print(f"  ✓ Bar plot saved → {bar_out}")
    print(f"  ✓ Collapsed network saved → {net_out}")
    print(f"  ✓ Expanded network saved → {exp_out}")
    print(f"  ✓ Combined PDF saved → {combined_out}")
    print("===================================================")

# --- Minimal loader + figure-builder fallbacks (keeps script self-contained) ---
def load_interaction_summary():
    """
    Simple loader for Analysis_Results/InteractionTypes_Summary.csv.
    Returns (df, total_frames)
    """
    validated_path = os.path.join(ANALYSIS_DIR, "ValidatedInteractionEvents.csv")
    if os.path.exists(validated_path):
        events = pd.read_csv(validated_path)
        status = events["ValidationStatus"].astype(str) if "ValidationStatus" in events else pd.Series("", index=events.index)
        events = events[status == "validated"]
        events = events[events["InteractionClass"].isin(["H-bond", "Hydrophobic", "Ionic"])]
        total_frames = max(1, int(pd.read_csv(validated_path)["Frame"].nunique()))
        if events.empty:
            return pd.DataFrame(columns=["Residue", "H-bond", "Hydrophobic", "Ionic", "Total", "Fraction", "HasHbond"]), total_frames
        unique = events.drop_duplicates(["ResidueID", "Residue", "InteractionClass", "Frame"])
        df = unique.groupby(["Residue", "InteractionClass"])["Frame"].nunique().unstack(fill_value=0).reset_index()
        for column in ["H-bond", "Hydrophobic", "Ionic"]:
            if column not in df:
                df[column] = 0
        df["Total"] = unique.drop_duplicates(["ResidueID", "Residue", "Frame"]).groupby("Residue")["Frame"].nunique().reindex(df["Residue"]).to_numpy()
        df["Fraction"] = df["Total"] / float(total_frames)
        df["HasHbond"] = df["H-bond"] > 0
        return df, total_frames

    fpath = os.path.join(ANALYSIS_DIR, "InteractionTypes_Summary.csv")
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"{fpath} not found")

    df = pd.read_csv(fpath)

    if "Residue" not in df.columns:
        raise ValueError("InteractionTypes_Summary.csv missing 'Residue' column")

    # Coerce numeric columns where possible
    for c in df.columns:
        if c == "Residue":
            continue
        try:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        except Exception:
            pass

    if "Total" not in df.columns:
        raise ValueError("InteractionTypes_Summary.csv must contain 'Total' column")

    print("⚠️ Using legacy InteractionTypes CSV: chemistry labels are unvalidated historical diagnostics.")
    total_frames = max(1, int(df["Total"].max()))
    df["Fraction"] = df["Total"] / float(total_frames)

    # boolean convenience flag for H-bonds
    if "H-bond" in df.columns:
        df["HasHbond"] = df["H-bond"].fillna(0) > 0
    else:
        df["HasHbond"] = False

    return df, total_frames


def load_residue_ligand_atom_map():
    """
    Build residue -> ligand-atom map by preferring FilteredContacts over full framewise file.
    Returns mapping dict.
    """
    f_validated = os.path.join(ANALYSIS_DIR, "ValidatedInteractionEvents.csv")
    if os.path.exists(f_validated):
        events = pd.read_csv(f_validated)
        status = events["ValidationStatus"].astype(str) if "ValidationStatus" in events else pd.Series("", index=events.index)
        events = events[status == "validated"]
        events = events[events["InteractionClass"].isin(["H-bond", "Hydrophobic", "Ionic"])]
        if not events.empty:
            return (events.groupby("Residue")["LigandCanonicalAtom"].agg(lambda x: x.value_counts().idxmax())
                    .apply(normalize_atom).to_dict())
    f_filt = os.path.join(ANALYSIS_DIR, "FilteredContacts_Framewise.csv")
    f_full = os.path.join(ANALYSIS_DIR, "InteractionTypes_Framewise.csv")

    if not os.path.exists(f_filt):
        raise FileNotFoundError(f"{f_filt} not found")
    if not os.path.exists(f_full):
        raise FileNotFoundError(f"{f_full} not found")

    df_filt = pd.read_csv(f_filt)
    df_full = pd.read_csv(f_full)

    # safe normalization
    if "LigAtomName" in df_filt.columns:
        df_filt["LigAtomName"] = df_filt["LigAtomName"].astype(str)
    else:
        df_filt["LigAtomName"] = ""

    if "LigAtomName" in df_full.columns:
        df_full["LigAtomName"] = df_full["LigAtomName"].astype(str)
    else:
        df_full["LigAtomName"] = ""

    # most frequent atom per residue
    if "Residue" in df_filt.columns and not df_filt.empty:
        filt_map = (
            df_filt.groupby("Residue")["LigAtomName"]
            .agg(lambda x: x.value_counts().idxmax())
            .apply(normalize_atom)
            .to_dict()
        )
    else:
        filt_map = {}

    if "Residue" in df_full.columns and not df_full.empty:
        full_map = (
            df_full.groupby("Residue")["LigAtomName"]
            .agg(lambda x: x.value_counts().idxmax())
            .apply(normalize_atom)
            .to_dict()
        )
    else:
        full_map = {}

    merged = full_map.copy()
    merged.update(filt_map)

    # minimal logging
    print("\n🔗 Residue → LigAtom (merged):")
    for k, v in sorted(merged.items()):
        print(f"  {k} → {v}")

    return merged


def make_bar_figure(ligand, df_summary, total_frames, minfrac, outdir, dpi=600):
    """Create and save Panel A (stacked bar). Returns path."""
    fig, ax = plt.subplots(figsize=(10, 4))
    draw_contact_bar_panel(ax, df_summary, total_frames, minfrac)
    path = os.path.join(outdir, f"{ligand}_bar.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Bar panel saved → {path}")
    return path


def draw_expanded_network_panel(ax, ligand, df_full, residue_to_ligatom,
                                atom_positions, atom_elements, bonds, mol,
                                atom_names, center, parent_map, show_edges=True,
                                show_edge_labels=True, one_edge_per_residue=False):
    """Draw a structure-preserving 2D ligand interaction diagram.

    This is deliberately a chemical diagram first and a network second: the
    ligand bond graph remains intact in the centre, while residues are placed
    on the side of the molecule containing their observed contact atom.
    """
    ax.set_aspect("equal")
    ax.axis("off")

    # Chemical skeleton first, so interaction lines never hide bond order.
    # This is intentionally a structure drawing (implicit carbon), not a
    # ball-and-stick/network depiction.
    for a1_idx, a2_idx, order in bonds:
        a1 = normalize_bond_endpoint(a1_idx, atom_names)
        a2 = normalize_bond_endpoint(a2_idx, atom_names)
        if not a1 or not a2 or a1.startswith("H") or a2.startswith("H"):
            continue
        draw_sid_bond(ax, *atom_positions[a1], *atom_positions[a2], order)

    df_edges = summarize_interaction_edges(df_full, MIN_FRACTION)
    if df_edges.empty:
        show_edges = False
        residues = []
        res_positions = {}
    else:
        # Resolve virtual sites and hydrogens before placement.  This keeps
        # contact edges at the chemically associated heavy atom, never at the
        # ligand centre.
        raw_map = {str(i): name for i, name in enumerate(df_edges["LigAtomName"])}
        resolved_map = fix_edges_to_center_safe(raw_map, atom_positions, atom_elements, center)
        df_edges = df_edges.copy()
        df_edges["DisplayAtom"] = [resolved_map[str(i)] for i in range(len(df_edges))]
        df_edges["DisplayAtom"] = [
            _resolve_h_to_heavy(atom, parent_map, atom_positions, atom_elements, center)
            for atom in df_edges["DisplayAtom"]
        ]
        df_edges = df_edges[df_edges["DisplayAtom"].isin(atom_positions)].copy()
        # One residue--atom anchor gets one unambiguous line style.  A real
        # trajectory can contain several event classes over time, but mixing
        # their colours on the same line implies a non-existent hybrid type.
        priority = {"H-bond": 3, "Ionic": 2, "Hydrophobic": 1, "Water bridge": 0}
        df_edges["_priority"] = df_edges["Type"].map(priority).fillna(-1)
        df_edges = (df_edges.sort_values(["Fraction", "_priority"], ascending=[False, False])
                    .drop_duplicates(["Residue", "DisplayAtom"], keep="first")
                    .drop(columns="_priority"))
        if one_edge_per_residue:
            df_edges = (df_edges.sort_values("Fraction", ascending=False)
                        .drop_duplicates(["Residue"], keep="first"))
        residue_atoms = {
            residue: list(dict.fromkeys(group.sort_values("Fraction", ascending=False)["DisplayAtom"].tolist()))
            for residue, group in df_edges.groupby("Residue")
        }
        for residue, atom in residue_to_ligatom.items():
            residue_atoms.setdefault(residue, [atom] if isinstance(atom, str) else atom)
        residues = sorted(df_edges["Residue"].unique(), key=residue_number)
        res_positions = auto_place_residues_guided(
            residues, residue_atoms, atom_positions, center, bonds=bonds
        )

    # Interaction edges and percentages.  A distinct curvature per bundle
    # makes parallel contacts traceable without sacrificing the molecular core.
    if show_edges and not df_edges.empty:
        bundle_order = {}
        occupied_labels = []
        for edge_index, (_, row) in enumerate(df_edges.iterrows()):
            residue, atom = row["Residue"], row["DisplayAtom"]
            if residue not in res_positions or atom not in atom_positions:
                continue
            x_lig, y_lig = atom_positions[atom]
            x_res, y_res = res_positions[residue]
            sx, sy = shorten_to_circle(x_lig, y_lig, x_res, y_res, LIG_ATOM_RADIUS)
            ex, ey = shorten_to_circle(x_res, y_res, x_lig, y_lig, RESIDUE_RADIUS)
            # Contacts sharing one ligand atom should leave that hotspot as a
            # small fan.  Grouping by residue previously reset the curvature
            # for every edge, stacking coincident connectors on top of one
            # another.
            key = atom
            ordinal = bundle_order.get(key, 0)
            bundle_order[key] = ordinal + 1
            # Give each connector a visible but compact bow.  The placement
            # has already established a clear line of sight; this curvature
            # makes closely related contacts traceable back to their exact
            # atom without reintroducing long cross-scaffold routes.
            rad = (0.12 + 0.045 * (ordinal // 2)) * (1 if ordinal % 2 == 0 else -1)
            color, linewidth, style = edge_style(row["Type"])
            ax.add_patch(FancyArrowPatch(
                (sx, sy), (ex, ey), arrowstyle="-", connectionstyle=f"arc3,rad={rad}",
                linewidth=linewidth, linestyle=style, color=color, alpha=0.86, zorder=2,
            ))
            if show_edge_labels:
                # Put the occupancy on the connector itself.  An opaque white
                # box masks a short stretch of dashes, creating an intentional
                # inline gap rather than a detached, potentially ambiguous
                # number near a bundle of contacts.
                # Try separated positions along the same curved connector.
                # Hotspots commonly have several lines leaving one atom; a
                # single fixed midpoint makes their occupancy labels overlap.
                preferred_fractions = (
                    [0.52, 0.72, 0.34, 0.86, 0.44, 0.62]
                    if ordinal % 2 == 0 else
                    [0.70, 0.36, 0.56, 0.84, 0.46, 0.64]
                )
                chosen = None
                for label_fraction in preferred_fractions:
                    candidate = arc3_label_point(sx, sy, ex, ey, rad, label_fraction)
                    if all(math.dist(candidate, placed) >= 0.95 for placed in occupied_labels):
                        chosen = candidate
                        break
                if chosen is None:
                    # Dense diagrams remain deterministic even when every
                    # ideal point is occupied: use the farthest available
                    # inline point rather than detaching the label.
                    chosen = max(
                        (arc3_label_point(sx, sy, ex, ey, rad, t) for t in preferred_fractions),
                        key=lambda point: min((math.dist(point, placed) for placed in occupied_labels), default=float("inf")),
                    )
                lx, ly = chosen
                occupied_labels.append(chosen)
                ax.text(lx, ly, f"{float(row['Fraction']) * 100:.0f}%", fontsize=8.6,
                        weight="bold", color=color, ha="center", va="center", zorder=4,
                        bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="none", alpha=1.0))

    # Draw only explicit heteroatoms over the conventional chemical skeleton.
    # Carbon vertices remain implicit, which is the key visual difference
    # between a molecular interaction diagram and the previous node graph.
    for name, (x, y) in atom_positions.items():
        if name.startswith("H"):
            continue
        element = atom_elements.get(name, "C")
        if element != "C":
            draw_sid_atom_label(ax, element, x, y)

    for residue in residues:
        x, y = res_positions[residue]
        chain, code, number = residue_badge_lines(residue)
        # Peripheral residues must remain annotations to the ligand, not the
        # visual centre of the diagram.  Keep their actual marker size in
        # sync with the collision radius used by the placement routine.
        ax.scatter(x, y, s=1600, color=RESIDUE_COLORS[classify_residue(residue)],
                   edgecolors="black", linewidth=1.45, zorder=6)
        ax.text(x, y, f"{chain}\n{code}\n{number}", fontsize=9.2, weight="bold",
                linespacing=0.90, ha="center", va="center", zorder=7)

    if show_edges:
        make_residue_legend(
            ax,
            interaction_types=df_edges["Type"].unique().tolist(),
            residue_categories=[classify_residue(residue) for residue in residues],
        )
    ax.relim()
    ax.autoscale_view()
    ax.margins(0.16)


def make_collapsed_network_figure(ligand, df_filt, residue_to_ligatom, outdir, dpi=600):
    """Create and save Panel B (collapsed network). Returns path."""
    fig, ax = plt.subplots(figsize=(12, 10))
    draw_network_panel(ax, ligand, df_filt, residue_to_ligatom)
    path = os.path.join(outdir, f"{ligand}_network.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Collapsed network saved → {path}")
    return path


def make_expanded_network_figure(ligand, df_full, residue_to_ligatom,
    atom_positions, atom_elements, bonds,
    mol, atom_names, center, parent_map,
    outdir, show_edges=True, dpi=600, out_name=None):
    """Create and save Panel C (expanded network). Returns path.

    show_edges may be:
      - True: full expanded rendering
      - False: nodes-only (no edges)
      - "single": simplified fallback that draws one representative edge per residue
    """
    fig, ax = plt.subplots(figsize=(12, 10))

    # Use the structure-preserving renderer for every normal display mode.
    # ``single`` means one strongest edge per residue, not the historical
    # ball-and-stick fallback.
    panel_fn = globals().get("draw_expanded_network_panel", None)
    if callable(panel_fn):
        panel_fn(
            ax,
            ligand=ligand,
            df_full=df_full,
            residue_to_ligatom=residue_to_ligatom,
            atom_positions=atom_positions,
            atom_elements=atom_elements,
            bonds=bonds,
            mol=mol,
            atom_names=atom_names,
            center=center,
            parent_map=parent_map,
            show_edges=bool(show_edges),
            show_edge_labels=bool(show_edges),
            one_edge_per_residue=(show_edges == "single"),
        )
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network saved → {path}")
        return path

    # Fallback branch (either no panel_fn, or show_edges == "single")
    ax.set_aspect("equal")
    ax.axis("off")

    # 1) ligand atoms
    for nm, (x, y) in atom_positions.items():
        if nm.startswith("H"):
            continue
        elem = atom_elements.get(nm, "C")
        # node color mapping (carbon -> black)
        if elem == "C":
            node_col = "#000000"
        elif elem == "S":
            node_col = "#FFD700"
        elif elem == "P":
            node_col = "#FF69B4"
        elif elem in ("F", "Cl", "Br", "I"):
            node_col = "#007849"
        elif elem == "O":
            node_col = "#D72638"
        elif elem == "N":
            node_col = "#1E90FF"
        else:
            node_col = ATOM_COLORS.get(elem, "#888")

        # label presence & color (all labels white)
        if elem == "C":
            label = None
        else:
            label = elem
        label_color = "white"

        ax.scatter(x, y, s=900, color=node_col, edgecolors="black", linewidth=2.2, zorder=3)
        if label:
            ax.text(x, y, label, fontsize=11, weight="bold", ha="center", va="center", color=label_color, zorder=4)

    # 2) ligand bonds
    for a1_idx, a2_idx, order in bonds:
        a1 = normalize_bond_endpoint(a1_idx, atom_names)
        a2 = normalize_bond_endpoint(a2_idx, atom_names)
        if a1 is None or a2 is None:
            continue
        if a1.startswith("H") or a2.startswith("H"):
            continue
        x1, y1 = atom_positions[a1]
        x2, y2 = atom_positions[a2]
        draw_bond(ax, x1, y1, x2, y2, order, atom_elements.get(a1, "C"), atom_elements.get(a2, "C"))

    # 3) aggregated edge table
    df_edges = df_full.copy()
    df_edges["LigAtomName"] = df_edges["LigAtomName"].apply(normalize_atom)
    df_edges = df_edges[(df_edges["Fraction"] >= MIN_FRACTION)]

    if df_edges.empty:
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network (empty) saved → {path}")
        return path

    # collapse to unique Residue,LigAtom,Type with max fraction
    df_edges = df_edges.groupby(["Residue", "LigAtomName", "Type"], as_index=False)["Fraction"].max()

    residues_with_edges = sorted(df_edges["Residue"].unique(), key=residue_number)
    if not residues_with_edges:
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network (no residues) saved → {path}")
        return path

    residue_pos, rx, ry = auto_place_residues_no_overlap(residues_with_edges, atom_positions)

    # draw residue nodes
    for res in residues_with_edges:
        x, y = residue_pos[res]
        cat = classify_residue(res)
        col = RESIDUE_COLORS[cat]
        ax.scatter(x, y, s=3000, color=col, edgecolors="black", linewidth=2.0, zorder=6)
        ax.text(x, y, res, fontsize=11, weight="bold", ha="center", va="center", zorder=7)

    # nodes-only requested -> finish here
    if show_edges is False:
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network (nodes-only) saved → {path}")
        return path

    # single-edge-per-residue mode
    if show_edges == "single":
        idx = df_edges.groupby("Residue")["Fraction"].idxmax()
        best = df_edges.loc[idx].reset_index(drop=True)
        for e_idx, row in best.iterrows():
            res = row["Residue"]
            if res not in residue_pos:
                continue
            xR, yR = residue_pos[res]
            raw_atom = residue_to_ligatom.get(res, row["LigAtomName"])
            atom = normalize_atom(raw_atom)

            if atom in parent_map and parent_map.get(atom) in atom_positions:
                atom = parent_map[atom]
            if atom not in atom_positions and atom.startswith("H"):
                root = atom.rstrip("0123456789")
                for h_mol2, heavy in parent_map.items():
                    if h_mol2.startswith(root) and heavy in atom_positions:
                        atom = heavy
                        break

            if atom and atom in atom_positions:
                xL, yL = atom_positions[atom]
            else:
                print(f"⚠ Skipping unresolved ligand atom '{raw_atom}' for residue {res}")
                continue

            sx, sy = shorten_to_circle(xL, yL, xR, yR, LIG_ATOM_RADIUS)
            ex, ey = shorten_to_circle(xR, yR, xL, yL, RESIDUE_RADIUS)

            itype = str(row["Type"])
            if "H-bond" in itype:
                col, lw = EDGE_HBOND, 2.5
            elif "Hydrophobic" in itype:
                col, lw = EDGE_HYDRO, 2.0
            elif "Water" in itype:
                col, lw = EDGE_WATER, 2.0
            elif "Ionic" in itype:
                col, lw = "#C39BD3", 2.3
            else:
                col, lw = "#999999", 1.8

            rad = 0.18 if (e_idx % 2 == 0) else -0.18
            # draw curved edge colored by interaction type (match legend)
            edge_patch = FancyArrowPatch((sx, sy), (ex, ey), connectionstyle=f"arc3,rad={rad}", arrowstyle="-",
                                        linestyle="--", linewidth=lw, color=col, zorder=2)
            ax.add_patch(edge_patch)

            pct = float(row["Fraction"]) * 100.0
            lx, ly = label_point_frac(sx, sy, ex, ey, 0.75)
            ax.text(lx, ly, f"{pct:.0f}%", fontsize=8, weight="bold", ha="center", va="center", color=col, zorder=5)

        make_residue_legend(ax)
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network (single-edge) saved → {path}")
        return path

    # full fallback: draw all edges
    for e_idx, row in df_edges.iterrows():
        res = row["Residue"]
        if res not in residue_pos:
            continue
        xR, yR = residue_pos[res]

        raw_atom = residue_to_ligatom.get(res, row["LigAtomName"])
        atom = normalize_atom(raw_atom)

        if atom in parent_map and parent_map.get(atom) in atom_positions:
            atom = parent_map[atom]
        if atom not in atom_positions and atom.startswith("H"):
            root = atom.rstrip("0123456789")
            for h_mol2, heavy in parent_map.items():
                if h_mol2.startswith(root) and heavy in atom_positions:
                    atom = heavy
                    break

        if atom and atom in atom_positions:
            xL, yL = atom_positions[atom]
        else:
            print(f"⚠ Skipping unresolved ligand atom '{raw_atom}' for residue {res}")
            continue

        sx, sy = shorten_to_circle(xL, yL, xR, yR, LIG_ATOM_RADIUS)
        ex, ey = shorten_to_circle(xR, yR, xL, yL, RESIDUE_RADIUS)

        itype = str(row["Type"])
        if "H-bond" in itype:
            col, lw = EDGE_HBOND, 2.5
        elif "Hydrophobic" in itype:
            col, lw = EDGE_HYDRO, 2.0
        elif "Water" in itype:
            col, lw = EDGE_WATER, 2.0
        elif "Ionic" in itype:
            col, lw = "#C39BD3", 2.3
        else:
            col, lw = "#999999", 1.8

        rad = 0.18 if (e_idx % 2 == 0) else -0.18
        edge = FancyArrowPatch((sx, sy), (ex, ey), connectionstyle=f"arc3,rad={rad}", arrowstyle="-",
                        linestyle="--", linewidth=lw, color=col, zorder=2)
        ax.add_patch(edge)

        pct = float(row["Fraction"]) * 100.0
        lx, ly = label_point_frac(sx, sy, ex, ey, 0.75)
        ax.text(lx, ly, f"{pct:.0f}%", fontsize=8, weight="bold", ha="center", va="center", color=col, zorder=5)

    make_residue_legend(ax)

    if out_name:
        path = out_name
    else:
        path = os.path.join(outdir, f"{ligand}_expanded.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Expanded network saved → {path}")
    return path
# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    chainmap = infer_mode_from_chainmap(".")
    if (not args.ligand) or args.minfrac is None or not args.output:
        args = ask_if_missing(args)

    auto_name = args.ligand
    if not auto_name:
        if chainmap.get("is_peptide"):
            auto_name = chainmap.get("group_name", "PEPTIDE")
        else:
            auto_name = os.path.basename(os.getcwd())

    build_two_panel_figure(auto_name, args.minfrac, args.output, no_edges=args.no_edges, mode=args.mode)
