#!/usr/bin/env python3
"""Recover a ligand MOL2 from a published final-trajectory PDB.

This utility is deliberately for *display/analysis recovery*, not for changing
an MD system.  It extracts one HETATM residue from ``Final_Trajectory.pdb``,
uses Open Babel to write ``<RES>.cgenff.mol2`` with the source atom order, and
can optionally send that MOL2 to the authenticated Kyle CGenFF Funnel to
recover the accompanying ``.str`` file.

CGenFF parameterizes an already-prepared MOL2; it does not create the MOL2.
Keeping those two steps explicit avoids implying that inferred connectivity is
validated topology data.  Review the generated MOL2 before using it for a new
simulation.  It is sufficient for the NETWORX interaction explorer, which
uses it only for atom names, elements, bonds, and 2-D drawing.
"""

from __future__ import annotations

import argparse
import collections
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable


WATER_AND_COMMON_IONS = {
    "HOH", "WAT", "SOL", "TIP3", "TIP3P", "NA", "CL", "K", "CA", "MG",
    "ZN", "MN", "FE", "CU", "CS", "BR", "IOD",
}


def hetatm_residues(pdb_path: Path) -> dict[tuple[str, str, str], list[str]]:
    """Return HETATM lines keyed by ``(resname, chain, resid)``.

    Alternate locations other than blank/A are skipped, which prevents the
    same atom being sent twice to Open Babel.
    """
    groups: dict[tuple[str, str, str], list[str]] = collections.defaultdict(list)
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        altloc = line[16:17]
        if altloc not in {" ", "A"}:
            continue
        resname = line[17:20].strip().upper()
        chain = line[21:22].strip() or "_"
        resid = (line[22:26].strip() + line[26:27].strip()).strip() or "0"
        if resname:
            groups[(resname, chain, resid)].append(line)
    return dict(groups)


def select_ligand(groups: dict[tuple[str, str, str], list[str]], code: str | None) -> tuple[tuple[str, str, str], list[str]]:
    """Select the requested residue or the largest non-solvent HETATM group."""
    if code:
        requested = code.upper().strip()
        candidates = [(key, lines) for key, lines in groups.items() if key[0] == requested]
        if not candidates:
            raise ValueError(f"No HETATM residue named {requested!r} was found.")
        return max(candidates, key=lambda item: len(item[1]))
    candidates = [(key, lines) for key, lines in groups.items() if key[0] not in WATER_AND_COMMON_IONS]
    if not candidates:
        raise ValueError("No non-solvent HETATM residue was found; use --ligand-code if appropriate.")
    return max(candidates, key=lambda item: len(item[1]))


def write_extracted_pdb(lines: Iterable[str], destination: Path) -> None:
    destination.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def make_mol2(extracted_pdb: Path, mol2_path: Path, ligand_code: str) -> None:
    obabel = shutil.which("obabel")
    if not obabel:
        raise RuntimeError("Open Babel (obabel) is required to create the recovery MOL2.")
    completed = subprocess.run(
        [obabel, "-ipdb", str(extracted_pdb), "-omol2", "-O", str(mol2_path), "--title", ligand_code],
        text=True, capture_output=True, check=False, timeout=120,
    )
    if completed.returncode or not mol2_path.is_file() or not mol2_path.stat().st_size:
        detail = (completed.stderr or completed.stdout or "Open Babel failed").strip()
        raise RuntimeError(f"Open Babel could not create {mol2_path.name}: {detail[-1000:]}")
    if "@<TRIPOS>MOLECULE" not in mol2_path.read_text(encoding="utf-8", errors="replace").upper():
        raise RuntimeError(f"Open Babel created an invalid MOL2: {mol2_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Published result root containing Final_Trajectory.pdb.")
    parser.add_argument("--ligand-code", help="Three/four-character HETATM residue name. Defaults to largest non-solvent HETATM residue.")
    parser.add_argument("--pdb", type=Path, help="Override Final_Trajectory.pdb path.")
    parser.add_argument("--parameterize", action="store_true", help="Also request a matching CGenFF .str from the Kyle Funnel.")
    parser.add_argument("--funnel-url", default="https://kyle.rove-vernier.ts.net/cgenff", help="Kyle CGenFF Funnel base URL.")
    parser.add_argument("--token-file", type=Path, default=Path("~/.config/pymacs/cgenff-funnel.token").expanduser())
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.expanduser().resolve()
    pdb = (args.pdb.expanduser().resolve() if args.pdb else next(iter(root.rglob("Final_Trajectory.pdb")), None))
    if pdb is None or not pdb.is_file():
        raise RuntimeError(f"No Final_Trajectory.pdb found below {root}")
    groups = hetatm_residues(pdb)
    key, lines = select_ligand(groups, args.ligand_code)
    code = key[0]
    if len(code) > 4:
        raise RuntimeError(f"Ligand residue {code!r} is too long for CGenFF; pass a 1-4 character --ligand-code.")

    output_dir = pdb.parent
    extracted = output_dir / f"{code}_recovered_from_final_trajectory.pdb"
    mol2 = output_dir / f"{code}.cgenff.mol2"
    write_extracted_pdb(lines, extracted)
    make_mol2(extracted, mol2, code)
    print(f"Recovered {code} from chain {key[1]}, residue {key[2]} ({len(lines)} HETATM records)")
    print(f"MOL2: {mol2}")

    if args.parameterize:
        from pymacs_cgenff_funnel import parameterize_mol2
        stream, service_log = parameterize_mol2(
            endpoint=args.funnel_url, token_file=args.token_file, ligand_code=code, mol2_path=mol2,
        )
        stream_path = output_dir / f"{code}.str"
        stream_path.write_text(stream, encoding="utf-8")
        (output_dir / f"{code}.cgenff-recovery.log").write_text(service_log, encoding="utf-8")
        print(f"CGenFF stream: {stream_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
