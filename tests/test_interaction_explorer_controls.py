from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_ligand_explorer_offers_full_persistence_range_without_suppression():
    source = (ROOT / "3B_NETWORX.py").read_text(encoding="utf-8")
    assert "range(10, 100, 10)" in source
    assert "default_threshold = 0.50" in source
    assert "must never disappear" in source


def test_interface_explorer_keeps_chain_chemistry_modes_mutually_exclusive():
    source = (ROOT / "3C_Interface_RIN.py").read_text(encoding="utf-8")
    assert "range(10 if has_peptide else 50, 100, 10)" in source
    assert "colourMode='chain'" in source
    assert "Color ${esc(p.left)} by chemistry" in source
    assert "Color ${esc(p.right)} by chemistry" in source
    assert "Exactly one chain may use residue-chemistry colours at a time" in source
