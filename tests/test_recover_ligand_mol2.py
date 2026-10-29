import tempfile
import unittest
from pathlib import Path

from pymacs_recover_ligand_mol2 import hetatm_residues, select_ligand, write_extracted_pdb


def hetatm(serial, atom, residue, chain, resid, element):
    return f"HETATM{serial:5d} {atom:<4} {residue:>3} {chain}{resid:4d}    {serial:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2}"


class RecoverLigandMol2Tests(unittest.TestCase):
    def test_selects_largest_non_solvent_hetatm_residue(self):
        with tempfile.TemporaryDirectory() as temporary:
            pdb = Path(temporary) / "Final_Trajectory.pdb"
            pdb.write_text("\n".join([
                hetatm(1, "O", "HOH", "A", 1, "O"),
                hetatm(2, "C1", "LIG", "B", 3, "C"),
                hetatm(3, "O1", "LIG", "B", 3, "O"),
                hetatm(4, "N1", "LIG", "B", 3, "N"),
                hetatm(5, "C", "COF", "A", 2, "C"),
            ]) + "\n", encoding="utf-8")
            key, lines = select_ligand(hetatm_residues(pdb), None)
            self.assertEqual(key, ("LIG", "B", "3"))
            self.assertEqual(len(lines), 3)

    def test_extract_writes_only_selected_hetatoms(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "ligand.pdb"
            write_extracted_pdb([hetatm(1, "C1", "LIG", "A", 1, "C")], output)
            self.assertEqual(output.read_text(encoding="utf-8").splitlines()[-1], "END")

