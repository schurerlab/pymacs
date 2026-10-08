import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

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

    def test_inspect_cif_reports_chains_and_hetero(self):
        with tempfile.TemporaryDirectory() as temporary:
            structure = Path(temporary) / "complex.cif"
            structure.write_text(
                "data_test\nloop_\n_atom_site.group_PDB\n_atom_site.label_comp_id\n"
                "_atom_site.label_asym_id\n_atom_site.auth_asym_id\nATOM ALA A A\n"
                "ATOM GLY B B\nHETATM LIG C C\nHETATM HOH D D\n#\n",
                encoding="utf-8",
            )
            self.assertEqual(pymacs_run.inspect_cif(structure), (["A", "B"], ["LIG"]))

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
            self.assertIn("#BSUB -P brd", contents)
            self.assertIn("#BSUB -L /bin/bash", contents)
            self.assertIn("source /etc/bashrc", contents)
            self.assertIn(str(folder.resolve()), contents)
            self.assertNotIn("source $HOME/miniforge3", contents)
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

    def test_cli_validate_prints_the_next_setup_command(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / "complex.pdb").write_text("ATOM\n", encoding="utf-8")
            write_config(folder)
            result = subprocess.run(
                [sys.executable, str(REPO / "pymacs_run.py"), "validate"],
                cwd=folder,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Next: run PyMACS Step 1 setup:", result.stdout)
            self.assertIn("python pymacs_run.py setup", result.stdout)

    def test_ligand_simulation_command_carries_ligand_flags(self):
        config = {
            "setup": {"ligand": "LIG", "cofactors": "HEM"},
            "simulation": {"mode": "ligand", "length_ns": 1, "threads": 16, "compute": "GPU"},
        }
        command = pymacs_run.simulation_command(config)
        self.assertEqual(command[command.index("--ligand") + 1], "LIG")
        self.assertEqual(command[command.index("--cofactors") + 1], "HEM")

    def test_funnel_ligand_setup_command_keeps_token_out_of_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / "complex.pdb").write_text("ATOM\n", encoding="utf-8")
            (folder / "1_AutomateGromacs_MPI.py").write_text("", encoding="utf-8")
            config = {
                "input": "complex.pdb",
                "setup": {
                    "mode": "ligand", "ligand": "DR7", "box_type": "dodecahedron", "box_distance_nm": 1.0,
                    "cgenff_backend": "funnel", "cgenff_funnel_url": "https://kyle.example.ts.net",
                    "cgenff_token_file": "~/.config/pymacs/cgenff-funnel.token",
                },
                "simulation": {"mode": "ligand", "length_ns": 1, "threads": 16, "compute": "GPU"},
            }
            pymacs_run.validate_config(config, folder)
            command = pymacs_run.setup_command(folder, config)
            self.assertIn("--cgenff-funnel-url", command)
            self.assertIn("--cgenff-token-file", command)
            self.assertNotIn("token-value", " ".join(command))

    def test_configure_option_five_writes_funnel_backend(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / "complex.pdb").write_text(
                "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n",
                encoding="utf-8",
            )
            answers = iter(["1", "Receptor", "DR7", "", "1.0", "10", "https://kyle.example.ts.net", "~/.config/pymacs/cgenff-funnel.token"])
            with patch("pymacs_run.ask", side_effect=lambda *_args, **_kwargs: next(answers)), \
                 patch("pymacs_run.ask_choice", side_effect=["ligand_funnel", "dodecahedron"]), \
                 patch("pymacs_run.ask_yes_no", side_effect=[False, True, True, True]):
                pymacs_run.configure(folder)
            config = json.loads((folder / "pymacs_run.json").read_text(encoding="utf-8"))
            self.assertEqual(config["setup"]["cgenff_backend"], "funnel")
            self.assertEqual(config["setup"]["mode"], "ligand")
            self.assertEqual(config["simulation"]["mode"], "ligand")


if __name__ == "__main__":
    unittest.main()
