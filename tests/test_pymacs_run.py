import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pymacs_run


REPO = Path(__file__).resolve().parents[1]


def write_config(folder: Path):
    config = {
        "input": "complex.pdb",
        "setup": {
            "mode": "protein",
            "chain_map": "A:Receptor,B:Peptide",
            "keep_chains": "A,B",
            "remove_input_waters": True,
            "remove_input_ions": True,
            "box_type": "dodecahedron",
            "box_distance_nm": 1.0,
        },
        "simulation": {"mode": "peptide", "length_ns": 0.25, "threads": 16, "compute": "GPU"},
    }
    (folder / "pymacs_run.json").write_text(json.dumps(config), encoding="utf-8")
    return config


class PyMACSRunTests(unittest.TestCase):
    def test_inspect_pdb_reports_chains_and_hetero(self):
        with tempfile.TemporaryDirectory() as temporary:
            structure = Path(temporary) / "complex.pdb"
            structure.write_text(
                "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n"
                "ATOM      2  N   GLY B   1      0.000   0.000   0.000  1.00  0.00           N\n"
                "HETATM    3  C1  LIG C   1      0.000   0.000   0.000  1.00  0.00           C\n"
                "HETATM    4  O   HOH D   1      0.000   0.000   0.000  1.00  0.00           O\n",
                encoding="utf-8",
            )
            self.assertEqual(pymacs_run.inspect_pdb(structure), (["A", "B"], ["LIG"]))

    def test_commands_and_lsf_render_for_peptide_folder(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / "complex.pdb").write_text("ATOM\n", encoding="utf-8")
            (folder / "1_AutomateGromacs_MPI.py").write_text("", encoding="utf-8")
            for name in ("topol.top", "solv_ions.gro", "em.tpr"):
                (folder / name).write_text("ok\n", encoding="utf-8")
            config = write_config(folder)
            pymacs_run.validate_config(config, folder)

            setup = pymacs_run.setup_command(folder, config)
            self.assertIn("protein", setup)
            self.assertIn("--chain-map", setup)
            self.assertIn("--keep-chains", setup)
            command = pymacs_run.simulation_command(config)
            self.assertEqual(command[command.index("--mode") + 1], "peptide")
            self.assertIn("--external-mpi", command)

            rendered = pymacs_run.render_lsf(folder, config, "triton", resume=False)
            contents = rendered.read_text(encoding="utf-8")
            self.assertIn("#BSUB -q normal", contents)
            self.assertIn("gromacs/2025.1-gcc-13.4.0-6vq7xfo", contents)
            self.assertIn("--mode peptide", contents)
            self.assertIn("--ntomp 16", contents)

    def test_cli_dry_run_submission(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / "complex.pdb").write_text("ATOM\n", encoding="utf-8")
            for name in ("topol.top", "solv_ions.gro", "em.tpr"):
                (folder / name).write_text("ok\n", encoding="utf-8")
            write_config(folder)
            result = subprocess.run(
                [sys.executable, str(REPO / "pymacs_run.py"), "submit", "--dry-run"],
                cwd=folder,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('#BSUB -gpu "num=1"', result.stdout)


if __name__ == "__main__":
    unittest.main()
