#!/usr/bin/env python3
"""Place the fully hydrogenated compound-32 MOL2 in the docked XTL pose."""
from pathlib import Path
import sys
import numpy as np

PDB = Path(sys.argv[1])
MOL2 = Path(sys.argv[2])
OUT = Path(sys.argv[3])

# Map MOL2 heavy-atom names to the corresponding docked XTL names.
NAME_MAP = {
    "C01": "C1", "C02": "C2", "C04": "C3", "C05": "C4",
    "C06": "C5", "C07": "C6", "C08": "C7", "C09": "C8",
    "N03": "N9", "N10": "N10", "N11": "N11", "BR12": "Br12",
}
HYDROGEN_NAMES = {13: "H1", 14: "H2", 15: "H3", 16: "H4", 17: "H5", 18: "H6"}


def parse_mol2(path):
    atoms, bonds = [], []
    section = None
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line[9:]
            continue
        if section == "ATOM" and line.strip():
            fields = line.split()
            atoms.append((int(fields[0]), fields[1], np.array(list(map(float, fields[2:5]))), fields[5]))
        elif section == "BOND" and line.strip():
            fields = line.split()
            bonds.append((int(fields[1]), int(fields[2])))
    return atoms, bonds


def docked_xtl(path):
    coords = {}
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "XTL":
            coords[line[12:16].strip()] = np.array([
                float(line[30:38]), float(line[38:46]), float(line[46:54])
            ])
    return coords


atoms, bonds = parse_mol2(MOL2)
target = docked_xtl(PDB)
source_by_name = {name: xyz for _, name, xyz, _ in atoms}
src = np.array([source_by_name[source] for source in NAME_MAP])
dst = np.array([target[target_name] for target_name in NAME_MAP.values()])

# Row-vector Kabsch alignment: source @ rot + shift = target.
src_centroid, dst_centroid = src.mean(axis=0), dst.mean(axis=0)
u, _, vt = np.linalg.svd((src - src_centroid).T @ (dst - dst_centroid))
rot = u @ vt
if np.linalg.det(rot) < 0:
    u[:, -1] *= -1
    rot = u @ vt
shift = dst_centroid - src_centroid @ rot
aligned = {idx: xyz @ rot + shift for idx, _, xyz, _ in atoms}
rmsd = np.sqrt(np.mean(np.sum((src @ rot + shift - dst) ** 2, axis=1)))

# Preserve the exact docked heavy-atom coordinates; use transformed MOL2 coordinates for H.
serial0 = min(
    int(line[6:11]) for line in PDB.read_text().splitlines()
    if line.startswith("HETATM") and line[17:20].strip() == "XTL"
)
ligand_lines = []
for idx, source_name, xyz, atom_type in atoms:
    name = NAME_MAP.get(source_name, HYDROGEN_NAMES.get(idx, source_name))
    xyz_out = target[name] if source_name in NAME_MAP else aligned[idx]
    element = "H" if atom_type.upper().startswith("H") else ("BR" if source_name == "BR12" else source_name[0])
    ligand_lines.append(
        f"HETATM{serial0 + idx - 1:5d} {name:>4s} XTL X   1    "
        f"{xyz_out[0]:8.3f}{xyz_out[1]:8.3f}{xyz_out[2]:8.3f}  1.00  0.00          {element:>2s}\n"
    )

# Replace the 12 heavy-atom XTL records with all 18 ligand atoms.
out, inserted = [], False
for line in PDB.read_text().splitlines(keepends=True):
    is_xtl = line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "XTL"
    if is_xtl:
        if not inserted:
            out.extend(ligand_lines)
            inserted = True
        continue
    out.append(line)
OUT.write_text("".join(out))

# Validate ligand bond geometry after placing its H atoms.
atom_xyz = {idx: (target[NAME_MAP[name]] if name in NAME_MAP else aligned[idx]) for idx, name, _, _ in atoms}
source_xyz = {idx: xyz for idx, _, xyz, _ in atoms}
errors = [abs(np.linalg.norm(atom_xyz[a] - atom_xyz[b]) - np.linalg.norm(source_xyz[a] - source_xyz[b])) for a, b in bonds]
print(f"heavy_atom_fit_rmsd_A={rmsd:.4f}")
print(f"max_bond_length_change_A={max(errors):.4f}")
print(f"wrote={OUT}")
