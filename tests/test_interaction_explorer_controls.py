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
    assert "edgeControls={}" in source
    assert "kind:'edge'" in source
    assert "updateVisibleEdge(edge)" in source
    assert "event.detail>=2" in source
    assert "labelSize.oninput" in source
    assert "(?:const EMPTY_MESSAGE=.*?;)?const EDGES" in source
    assert "sourceNodes=[...svg.querySelectorAll('.residue-node')]" in source
    assert "exportedCircle.style.setProperty('fill',style.fill,'important')" in source
    assert "sourceLabels=[...svg.querySelectorAll('.edge-label')]" in source
    assert "atomColours={N:'#58b9ff',O:'#ff6478',S:'#e5c400',P:'#ff927b'}" in source
    assert "labelTrims={F:11,Cl:18,Br:18,I:8,P:12}" in source
    assert "firstTrim=Math.min(labelTrims[line.dataset.firstElement]||0" in source
    assert "copy.querySelectorAll('.atom-N,.atom-O,.atom-S,.atom-C').forEach(atom=>atom.remove())" in source
    assert 'data-first-element' in source
    assert "Changing the retained threshold resets your node arrangement. Continue?" not in source
    assert 'id="export-svg"' in source
    assert "Choosing a different retained threshold resets these temporary edits" in source
    assert "svg.addEventListener('dblclick'" in source
    assert "ring_inset_center" in source
    assert 'id="move-info"' in source
    assert 'id="move-popover"' in source
    assert 'id="export-legend-horizontal"' in source
    assert 'id="export-legend-vertical"' in source
    assert "function downloadLegend(orientation)" in source
    assert "Double-click anywhere on the diagram to clear" in source
    assert "svg.addEventListener('dblclick',event=>{event.preventDefault();event.stopPropagation();select(null)})" in source
    assert "label_clearance = 0.0 if element == \"C\"" in source
    assert '"atom_clearance": round(label_clearance, 2)' in source
    assert "sourceX=atomX+toTargetX/toTargetLength*clearance" in source


def test_dashboard_embed_allows_export_and_has_room_for_the_full_explorer():
    source = (ROOT / "5_PyMACS_Result_Dashboard.py").read_text(encoding="utf-8")
    assert 'sandbox="allow-scripts allow-downloads"' in source
    assert "height:1320px" in source
    assert 'if chain else \'\'' in source
    assert "if(mode==='complex'){{pdbModel.setStyle" in source
    assert "viewer.zoomTo()}}if(mode==='chain'" in source
    assert "data-final-ligand-b64" in source
    assert "mol2_with_final_pdb_coordinates" in source
    assert "metadata_ligand" in source
    assert "or final_ligand_payload" in source
    assert "singleBonds:false" in source
    assert "aromaticStyle:'solid'" in source
    assert "def ligand_scorecard_html" in source
    assert "def contact_lifetime_plot" in source
    assert "def contact_turnover_plot" in source
    assert "def contact_fingerprint_html" in source
    assert "def chemistry_composition_plot" in source
    assert "def residue_chemistry_heatmap" in source
    assert "def representative_snapshots_html" in source
    assert '"Binding scorecard"' in source
    assert '"Ligand contact-state timeline"' in source
    assert '"Residue contact fingerprint"' in source


def test_legacy_ligand_evidence_is_not_hidden_by_missing_metadata():
    source = (ROOT / "5_PyMACS_Result_Dashboard.py").read_text(encoding="utf-8")
    assert "has_ligand_evidence" in source
    assert "Older published runs occasionally lack ligand metadata" in source
    assert 'any(analysis.glob("*_contact_frequency_filtered.csv"))' in source
    assert 'str.extract(r"(\\d+)", expand=False)' in source
    assert '"Most continuously proximal"' in source
    assert '"Strongest typed interaction"' in source
    assert "Why the two interaction occupancies can differ" in source
    assert "Reuse the resolved edge records published by 3B's interaction explorer" in source
    assert 'const EDGES=(\\[.*?\\]);\\s*const ' in source
    assert "def interaction_explorer_summary" in source
    assert 'Strongest observed edge:\\s*([0-9.]+)%' in source
    assert "no exact residue–atom/type pair reaches the 10% explorer threshold" in source


def test_ligand_evidence_promotes_a_stale_protein_context():
    source = (ROOT / "5_PyMACS_Result_Dashboard.py").read_text(encoding="utf-8")
    assert "def promote_ligand_context_when_evidenced" in source
    assert 'context.partner_kind not in {"ligand", "peptide"}' in source
    assert '"ligand", ligand_label_from_evidence' in source
    assert "context, has_ligand_evidence = promote_ligand_context_when_evidenced(" in source


def test_legacy_ligand_receptor_name_uses_historical_source_before_generic_protein():
    source = (ROOT / "5_PyMACS_Result_Dashboard.py").read_text(encoding="utf-8")
    assert "def legacy_receptor_label" in source
    assert 'run.get("historical_source")' in source
    assert 're.fullmatch(r"([A-Za-z][A-Za-z0-9-]{1,30})_[0-9][A-Za-z0-9]{3}"' in source
    assert 'return f"Receptor (chain {chain})"' in source
    assert "or legacy_receptor_label(root)" in source


def test_interface_explorer_keeps_chain_chemistry_modes_mutually_exclusive():
    source = (ROOT / "3C_Interface_RIN.py").read_text(encoding="utf-8")
    assert "range(10 if has_peptide else 50, 100, 10)" in source
    assert "colourMode='chain'" in source
    assert "Color ${esc(p.left)} by chemistry" in source
    assert "Color ${esc(p.right)} by chemistry" in source
    assert "Exactly one chain may use residue-chemistry colours at a time" in source
    assert "radiusFor=side=>isPeptide(side)?40:proteinRadius" in source
