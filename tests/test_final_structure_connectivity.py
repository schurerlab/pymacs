import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "5_PyMACS_Result_Dashboard.py"

def load_connectivity_helper():
    """Load only the dependency-free helper from the dashboard source."""
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    function = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "mol2_with_final_pdb_coordinates"
    )
    module = ast.Module(
        body=[ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0), function],
        type_ignores=[],
    )
    namespace = {"Path": Path}
    exec(compile(ast.fix_missing_locations(module), str(SCRIPT), "exec"), namespace)
    return namespace["mol2_with_final_pdb_coordinates"]


def pdb_line(atom: str, resname: str, x: float, y: float, z: float) -> str:
    return f"HETATM    1 {atom:>4} {resname:>3} B   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C\n"


def test_rebuilds_cgenff_mol2_with_final_pdb_coordinates(tmp_path):
    rebuild = load_connectivity_helper()
    mol2 = tmp_path / "LIG.cgenff.mol2"
    mol2.write_text(
        "@<TRIPOS>MOLECULE\nLIG\n 2 1 0 0 0\nSMALL\nGASTEIGER\n\n"
        "@<TRIPOS>ATOM\n"
        "1 C01 1.0 2.0 3.0 C.3 1 LIG 0.0\n"
        "2 N02 2.0 3.0 4.0 N.3 1 LIG 0.0\n"
        "\n"
        "@<TRIPOS>BOND\n1 1 2 1\n",
        encoding="utf-8",
    )
    pdb = tmp_path / "Final_Trajectory.pdb"
    pdb.write_text(
        pdb_line("C01", "LIG", 11.0, 12.0, 13.0)
        + pdb_line("N02", "LIG", 21.0, 22.0, 23.0)
        + pdb_line("LP1", "LIG", 0.0, 0.0, 0.0),
        encoding="utf-8",
    )

    rebuilt, resname = rebuild(pdb, mol2, "")

    assert resname == "LIG"
    assert rebuilt is not None
    assert "11.0000    12.0000    13.0000" in rebuilt
    assert "21.0000    22.0000    23.0000" in rebuilt
    assert "@<TRIPOS>BOND\n1 1 2 1" in rebuilt


def test_rebuild_requires_a_complete_atom_name_match(tmp_path):
    rebuild = load_connectivity_helper()
    mol2 = tmp_path / "LIG.cgenff.mol2"
    mol2.write_text(
        "@<TRIPOS>ATOM\n1 C01 0 0 0 C.3 1 LIG 0\n2 N02 0 0 0 N.3 1 LIG 0\n@<TRIPOS>BOND\n1 1 2 1\n",
        encoding="utf-8",
    )
    pdb = tmp_path / "Final_Trajectory.pdb"
    pdb.write_text(pdb_line("C01", "LIG", 1.0, 2.0, 3.0), encoding="utf-8")

    rebuilt, resname = rebuild(pdb, mol2, "LIG")

    assert rebuilt is None
    assert resname == ""
