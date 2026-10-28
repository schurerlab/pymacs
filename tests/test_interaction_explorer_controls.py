from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_ligand_explorer_offers_full_persistence_range_without_suppression():
    source = (ROOT / "3B_NETWORX.py").read_text(encoding="utf-8")
    assert "range(10, 100, 10)" in source
    assert "default_threshold = 0.50" in source
    assert "must never disappear" in source
    assert "760.0 / span_x" in source
    assert "Hydrogens/virtual sites are not drawn" in source
    assert 'id="drag-toggle"' not in source
    assert 'id="label-size"' in source
    assert 'type="range"' in source
    assert "labelSize.addEventListener('input'" in source
    assert 'min="20" max="64" value="26"' in source
    assert "node.setPointerCapture" in source
    assert "event.detail>=2" in source
    assert "labelSize.oninput" in source
    assert "(?:const EMPTY_MESSAGE=.*?;)?const EDGES" in source
    assert 'id="export-svg"' in source
    assert "Changing the retained threshold resets your node arrangement" in source
    assert "node.addEventListener('dblclick'" in source
    assert "ring_inset_center" in source
    assert "move-hint" in source


def test_dashboard_embed_allows_export_and_has_room_for_the_full_explorer():
    source = (ROOT / "5_PyMACS_Result_Dashboard.py").read_text(encoding="utf-8")
    assert 'sandbox="allow-scripts allow-downloads"' in source
    assert "height:1320px" in source


def test_interface_explorer_keeps_chain_chemistry_modes_mutually_exclusive():
    source = (ROOT / "3C_Interface_RIN.py").read_text(encoding="utf-8")
    assert "range(10 if has_peptide else 50, 100, 10)" in source
    assert "colourMode='chain'" in source
    assert "Color ${esc(p.left)} by chemistry" in source
    assert "Color ${esc(p.right)} by chemistry" in source
    assert "Exactly one chain may use residue-chemistry colours at a time" in source
    assert "radiusFor=side=>isPeptide(side)?40:proteinRadius" in source
