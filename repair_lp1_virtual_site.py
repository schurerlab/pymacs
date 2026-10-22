#!/usr/bin/env python3
"""Initialize a CGenFF type-2 LP virtual site in a GROMACS .gro file.

For virtual_sites2 function 2, the final parameter is a distance (nm):
the site is placed that distance from the first parent, in the direction of
the second parent.  A negative value therefore puts LP1 beyond the first
parent, opposite the second.  GROMACS will reconstruct this site during MD;
this makes its *input* coordinate internally consistent for preprocessing.
"""
import math
import sys
from pathlib import Path

source, target = map(Path, sys.argv[1:3])
lines = source.read_text().splitlines()
natoms = int(lines[1].strip())
atoms = lines[2 : 2 + natoms]

def fields(line):
    # GRO columns: residue number (0:5), residue name (5:10), atom name (10:15).
    return line[5:10].strip(), line[10:15].strip(), [float(line[i:i + 8]) for i in (20, 28, 36)]

selected = [(i, *fields(line)) for i, line in enumerate(atoms)]
lookup = {name: (i, xyz) for i, res, name, xyz in selected if res == "XTL"}
br_i, br = lookup["Br12"]
c7_i, c7 = lookup["C7"]
lp_i, _ = lookup["LP1"]

vector = [c7[k] - br[k] for k in range(3)]
length = math.sqrt(sum(component * component for component in vector))
distance_nm = -0.189
if length == 0:
    raise ValueError("Br12 and C7 have identical coordinates")
lp = [br[k] + distance_nm * vector[k] / length for k in range(3)]
old = atoms[lp_i]
atoms[lp_i] = old[:20] + "".join(f"{value:8.3f}" for value in lp) + old[44:]
target.write_text("\n".join([lines[0], lines[1], *atoms, *lines[2 + natoms :]]) + "\n")
print(f"LP1: ({lp[0]:.3f}, {lp[1]:.3f}, {lp[2]:.3f}) nm; Br12-C7 = {length:.3f} nm")
