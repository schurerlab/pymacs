#!/usr/bin/env python3
"""Build a self-contained interactive PyMACS result dossier for one MD run.

This intentionally complements (rather than replaces) 00_Pymacs_X_analysis.py.
The latter is a comparative apo/ligand dashboard; this program is the automatic
per-system report that follows a normal 3A/PROTAC analysis job.
"""
from __future__ import annotations

import argparse
import base64
import datetime as dt
import json
import os
import re
import subprocess
import textwrap
from dataclasses import dataclass
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot


BG, PANEL, GRID, TEXT, MUTED = "#05060a", "#0b1020", "#1c2a4a", "#cfe8ff", "#7aa6d9"
CYAN, GREEN, MAGENTA, ORANGE = "#35d0ff", "#18ff8b", "#ff2bd6", "#ff8a3d"
# A browser report should open immediately even when its source CSV contains a
# full 500 ns trajectory at sub-ps sampling.  The CSV download remains exact;
# this only limits the points drawn in the browser.
MAX_RENDER_POINTS = 2000
MAX_HEATMAP_TIMEPOINTS = 250
MAX_EVENT_POINTS = 12000
MAX_TABLE_ROWS = 50000

# The analysis pipeline intentionally writes both general-purpose calculations
# and specialised diagnostics.  Some of their traces describe essentially the
# same question at different selections (for example Protein_RMSD,
# FullComplex_RMSD and RMSD_<chain>).  They all remain in the run directory,
# but the browser gallery needs one canonical reading path instead of a wall
# of near-duplicates.  This list is also written as a small, machine-readable
# contract for the portal publisher.
RESULTS_MANIFEST_NAME = "portal_results_manifest.json"
CANONICAL_RESULT_ORDER = (
    "DSSP_encoded_numeric.csv",
    "Protein_RMSD.csv",
    "RMSF_*_per_residue.csv",
    "Radius_of_Gyration_Protein_Ligand_Complex.csv",
    "SSE_Fractions.csv",
    "*_Ligand_RMSD.csv",
    "*_Ligand_RMSF.csv",
    "*_contact_frequency_filtered.csv",
    "*_contact_map_residue_time.csv",
    "InteractionTypes_Framewise.csv",
    "native_contact_q.csv",
    "protein_native_contact_q.csv",
)


@dataclass(frozen=True)
class PartnerContext:
    """One report vocabulary for ligand, peptide, chain-interface, and apo runs."""
    partner_kind: str
    partner_label: str
    partner_chain_ids: tuple[str, ...]
    interface_pairs: tuple[tuple[str, str], ...]
    system_class: str

    @property
    def is_interface(self) -> bool:
        return self.partner_kind in {"peptide", "protein_chain"}

    @property
    def contact_tab_label(self) -> str:
        return {"peptide": "Peptide–protein interface", "protein_chain": "Chain interfaces", "none": "Contacts"}.get(self.partner_kind, "Contacts & binding")

    @property
    def relationship_label(self) -> str:
        return {"ligand": "ligand–protein", "peptide": "peptide–protein", "protein_chain": "chain-pair", "none": "intramolecular"}[self.partner_kind]


def atomindex_components(root: Path) -> list[dict[str, str]]:
    """Read component names and simulation chains from atomIndex's provenance comment."""
    entries: list[dict[str, str]] = []
    try:
        lines = (root / "atomIndex.txt").read_text(encoding="utf-8").splitlines()
    except OSError:
        return entries
    for raw in lines:
        body, _, comment = raw.partition("#")
        match = re.match(r"^\s*(\S+)\s+\d+\s*-\s*\d+\s*$", body)
        if not match:
            continue
        metadata = dict(re.findall(r"(\w+)=([^\s#]+)", comment))
        entries.append({"label": match.group(1), **metadata})
    return entries


def component_display_labels(root: Path, context: PartnerContext) -> list[str]:
    labels = []
    for entry in atomindex_components(root):
        chain = entry.get("current_chain") or entry.get("source_chain") or ""
        is_partner = context.partner_kind == "peptide" and chain in context.partner_chain_ids
        role = "peptide" if is_partner else "receptor" if context.partner_kind == "peptide" else entry.get("chain_type", "protein")
        # Older runs can have a valid component name but no explicit chain
        # provenance. Never render a fictitious-looking “chain ?” label.
        labels.append(f"{entry['label']} ({role}{f', chain {chain}' if chain else ''})")
    return labels


def partner_context(root: Path, analysis: Path, setup: dict, metadata: dict,
                    simulation: dict | None = None) -> PartnerContext:
    """Infer the report partner from authoritative run/interface metadata, not filenames alone."""
    chain_ids: list[str] = []
    chain_types: dict[str, str] = {}
    peptide_chain_ids: list[str] = []
    pairs: list[tuple[str, str]] = []
    manifest_path = analysis / "Interface_RIN" / "interface_analysis_manifest.json"
    try:
        interface_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for entry in interface_manifest.get("chain_metadata", []):
            label = str(entry.get("display_label") or entry.get("current_chain") or "").strip()
            chain_id = str(entry.get("current_chain") or entry.get("chain_id") or "").strip()
            if label:
                chain_types[label] = str(entry.get("chain_type") or "").casefold()
            if chain_id and chain_id not in chain_ids:
                chain_ids.append(chain_id)
            if chain_id and "peptide" in str(entry.get("chain_type") or "").casefold():
                peptide_chain_ids.append(chain_id)
    except (OSError, ValueError, TypeError):
        pass
    summary_path = analysis / "Interface_RIN" / "interface_chain_pair_summary.csv"
    try:
        summary = pd.read_csv(summary_path, usecols=["ChainA_Label", "ChainB_Label"])
        pairs = list(dict.fromkeys((str(row.ChainA_Label), str(row.ChainB_Label)) for row in summary.itertuples(index=False)))
    except (OSError, ValueError, pd.errors.ParserError):
        pass
    configured = []
    for entry in str(setup.get("chain_map") or "").split(","):
        chain, separator, label = entry.strip().partition(":")
        if chain and separator:
            configured.append((chain.strip(), label.strip()))
    if not chain_ids:
        chain_ids = [chain for chain, _ in configured]
    atomindex_entries = atomindex_components(root)
    if not chain_ids:
        chain_ids = [entry["current_chain"] for entry in atomindex_entries if entry.get("current_chain")]
    peptide_identity: dict = {}
    for identity_path in (analysis / "portal_peptide_identity.json", root / ".chainmap.peptide.json"):
        try:
            candidate = json.loads(identity_path.read_text(encoding="utf-8"))
            if isinstance(candidate, dict):
                peptide_identity = candidate
                break
        except (OSError, ValueError, TypeError):
            continue
    identity_chain = str(peptide_identity.get("chain_id") or peptide_identity.get("chain") or "").strip()
    identity_name = str(peptide_identity.get("name") or peptide_identity.get("group_name") or "").strip()
    if identity_chain and identity_chain not in peptide_chain_ids:
        peptide_chain_ids.append(identity_chain)
    has_peptide_metrics = any((analysis / filename).is_file() for filename in (
        "Peptide_RMSD.csv", "Peptide_RMSF.csv", "peptide_native_contact_q.csv",
    ))
    peptide_pairs = [pair for pair in pairs if any("peptide" in chain_types.get(label, "") for label in pair)]
    ligand = str(metadata.get("ligand_name") or metadata.get("ligand_residue") or setup.get("ligand") or "").strip()
    mode = str((simulation or {}).get("mode") or setup.get("mode") or "").casefold()
    if ligand and ligand.casefold() not in {"none", "no ligand label", "proteinonly"}:
        return PartnerContext("ligand", ligand, tuple(chain_ids), tuple(pairs), "ligand_complex")
    atomindex_peptide = next((entry for entry in atomindex_entries
                              if entry.get("current_chain") in peptide_chain_ids or entry.get("chain_type", "").casefold() == "peptide"), None)
    if peptide_pairs or mode == "peptide" or identity_name or has_peptide_metrics:
        peptide_labels = [label for pair in peptide_pairs for label in pair if "peptide" in chain_types.get(label, "")]
        label = identity_name or (atomindex_peptide or {}).get("label") or next(iter(peptide_labels), "Peptide")
        return PartnerContext("peptide", label, tuple(peptide_chain_ids or chain_ids), tuple(pairs), "peptide_protein_complex")
    if pairs or len(chain_ids) > 1:
        label = " / ".join(chain_ids) if chain_ids else "Chain interface"
        return PartnerContext("protein_chain", label, tuple(chain_ids), tuple(pairs), "protein_complex")
    label = "/".join(chain_ids) or "Protein"
    return PartnerContext("none", label, tuple(chain_ids), tuple(), "apo_protein")


def existing_result_matches(analysis: Path, pattern: str) -> list[Path]:
    """Return deterministic matches for a canonical result pattern."""
    return sorted(path for path in analysis.glob(pattern) if path.is_file())


def build_results_manifest(root: Path, analysis: Path, output: Path,
                           system_name: str, ligand_name: str) -> dict:
    """Describe the standard, presentation-safe outputs of a completed run.

    This deliberately separates the canonical gallery from supplementary
    files.  No science is discarded: supplementary diagnostics remain in the
    analysis directory and data bundle, while the gallery avoids showing three
    subtly different RMSDs or atom- and residue-level RMSFs as if they were
    interchangeable measurements.
    """
    canonical: list[str] = []
    for pattern in CANONICAL_RESULT_ORDER:
        for path in existing_result_matches(analysis, pattern):
            relative = path.relative_to(root).as_posix()
            if relative not in canonical:
                canonical.append(relative)
    report_relative = output.relative_to(root).as_posix()
    ligand_svg = analysis / "Ligand_2D.svg"
    if ligand_svg.is_file():
        canonical.append(ligand_svg.relative_to(root).as_posix())
    interface_explorer = analysis / "Interface_RIN" / "Interface_RIN_Explorer.html"
    if interface_explorer.is_file():
        canonical.append(interface_explorer.relative_to(root).as_posix())
    canonical.append(report_relative)

    all_outputs = sorted(
        path.relative_to(root).as_posix()
        for path in analysis.rglob("*")
        if path.is_file() and path.name not in {RESULTS_MANIFEST_NAME}
    )
    supplementary = [relative for relative in all_outputs if relative not in canonical]
    required = {
        "global_stability": "Analysis_Results/Protein_RMSD.csv",
        "secondary_structure": "Analysis_Results/DSSP_encoded_numeric.csv",
        "interactive_report": report_relative,
    }
    optional = {
        "ligand_pose": next((path for path in canonical if path.endswith("_Ligand_RMSD.csv")), None),
        "ligand_contacts": next((path for path in canonical if path.endswith("InteractionTypes_Framewise.csv")), None),
        "native_contact_q": next((path for path in canonical if path.endswith("native_contact_q.csv")), None),
        "interface_rin": next((path for path in canonical if path.endswith("Interface_RIN_Explorer.html")), None),
    }
    return {
        "schema_version": 1,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "system_name": system_name,
        "ligand_name": ligand_name,
        "analysis_complete": all(Path(root / relative).is_file() for relative in required.values()),
        "report": report_relative,
        "canonical_assets": canonical,
        "supplementary_assets": supplementary,
        "required_outputs": required,
        "optional_outputs": optional,
        "notes": {
            "canonical_policy": "The interactive gallery shows one primary trace per scientific question; supplementary diagnostics remain downloadable.",
            "native_contact_q": "Optional until native-contact analysis has completed.",
        },
    }


def canonical_gallery_csvs(analysis: Path) -> set[Path]:
    """CSV inputs shown in the main gallery; all other CSVs remain exports."""
    selected: set[Path] = set()
    for pattern in CANONICAL_RESULT_ORDER:
        selected.update(existing_result_matches(analysis, pattern))
    # Matrix outputs are rendered as bespoke heatmaps rather than XY cards.
    return {
        path for path in selected
        if path.name not in {"DSSP_encoded_numeric.csv", "InteractionTypes_Framewise.csv"}
        and not path.name.endswith("_contact_map_residue_time.csv")
    }


GRAPH_GUIDES = {
    "protein_rmsd": (
        "Global protein stability",
        "Protein RMSD measures backbone movement relative to the starting reference after alignment.",
        "A short early adjustment followed by a plateau is common. Persistent upward drift or large late spikes merits inspection; RMSD alone is not a binding-affinity measurement.",
    ),
    "ligand_rmsd": (
        "Ligand pose stability",
        "Ligand RMSD tracks how much the ligand changes position or conformation from its starting bound pose.",
        "A stable, low-to-moderate plateau together with persistent contacts supports retention of the modeled pose. A steadily rising trace can signal reorientation or possible departure from the pocket.",
    ),
    "ligand_rmsf": (
        "Ligand flexibility",
        "Ligand RMSF reports how much each ligand atom or position fluctuates during the trajectory.",
        "Mobile terminal groups are normal. Broadly high values across the ligand are more consistent with a flexible or changing pose; interpret alongside ligand RMSD and contacts.",
    ),
    "protein_rmsf": (
        "Local protein flexibility",
        "Protein RMSF shows which residues move most across the simulation.",
        "Flexible loops and termini often peak. Unexpected peaks in the structured core or binding pocket deserve a closer look, especially if they coincide with contact loss.",
    ),
    "gyration": (
        "Protein compactness",
        "Radius of gyration (Rg) summarizes how spread out the protein/complex is around its center of mass.",
        "A flat trace supports stable compactness. A sustained rise can reflect opening or expansion; a sustained fall can reflect compaction. Neither is by itself proof of a favorable or unfavorable ligand.",
    ),
    "contact_frequency": (
        "Contact persistence",
        "This reports how frequently ligand-contacting residues are observed across sampled frames.",
        "Persistent high-contact residues identify a stable interaction fingerprint. Contacts that appear and disappear can reflect normal local dynamics or a changing pose.",
    ),
    "interaction": (
        "Interaction chemistry",
        "Hydrogen-bond, ionic, and hydrophobic contacts summarize the chemical character of ligand–protein engagement.",
        "A coherent, persistent mixture is often more informative than a raw contact count. Categories describe the modeled trajectory; they are not an affinity score.",
    ),
    "sse": (
        "Secondary-structure evolution",
        "Each row follows a residue’s secondary structure through time: coil, sheet, or helix.",
        "Long continuous blocks indicate retained secondary structure. Repeated color changes identify local transitions. This is a structural-stability readout, not a binding measurement.",
    ),
    "sse_fraction": (
        "Secondary-structure fractions",
        "This is not a yes/no trace. At each time point it reports the fraction of residues assigned as helix, sheet, or coil, so a value such as 0.25 means 25% of the selected protein residues are helical.",
        "Stable fractions support retention of the overall fold. A sustained change in one fraction together with changes in the residue-by-time SSE map is stronger evidence for a structural transition than small frame-to-frame noise.",
    ),
    "q": (
        "Native-contact retention (Q)",
        "Q is the fraction of reference contacts retained from the first analyzed frame. Depending on the selection, it reports a ligand/partner interface or protein-fold contact retention.",
        "Higher, sustained Q means the stated reference interface or fold is retained. Q evaluates persistence relative to that model; it does not establish experimental binding affinity.",
    ),
    "generic": (
        "Analysis readout",
        "This plot is computed directly from the local trajectory analysis output shown in its axes.",
        "Read the overall trend first, then inspect large transitions or late-time changes. Cross-check with the companion stability and interaction panels before drawing a conclusion.",
    ),
}


def graph_guide_key(name: str) -> str:
    key = name.lower()
    if "sse" in key and "fraction" in key:
        return "sse_fraction"
    if "dssp" in key or "secondary" in key or "sse" in key or "helix" in key or "sheet" in key:
        return "sse"
    if "interaction" in key or "contact" in key or "persistence" in key:
        return "interaction" if "interaction" in key else "contact_frequency"
    if "ligand" in key and "rmsd" in key:
        return "ligand_rmsd"
    if "ligand" in key and "rmsf" in key:
        return "ligand_rmsf"
    if "rmsd" in key:
        return "protein_rmsd"
    if "rmsf" in key or "flexibility" in key:
        return "protein_rmsf"
    if "gyration" in key or key.startswith("rg"):
        return "gyration"
    if "native" in key or key == "qx":
        return "q"
    return "generic"


def guide_html(key: str) -> str:
    title, meaning, expected = GRAPH_GUIDES.get(key, GRAPH_GUIDES["generic"])
    return (
        '<details class="guide"><summary>ⓘ What this means</summary>'
        f'<p><b>{escape(title)}.</b> {escape(meaning)}</p>'
        f'<p><b>What to look for:</b> {escape(expected)}</p></details>'
    )


def ligand_inputs(root: Path, ligand_name: str) -> tuple[Path | None, Path | None]:
    """Prefer the parameterized ligand file and its perceived SMILES sidecar."""
    code = "".join(ch for ch in ligand_name if ch.isalnum())[:12]
    candidates = []
    for suffix in ("*.cgenff.mol2", "*.mol2", "*.sdf", "*.mol"):
        candidates.extend(sorted(root.glob(suffix)))
    if code:
        preferred = [path for path in candidates if path.name.upper().startswith(code.upper())]
        candidates = preferred + [path for path in candidates if path not in preferred]
    structure = candidates[0] if candidates else None
    smiles = None
    if structure:
        candidate = structure.with_name(structure.name.split(".", 1)[0] + ".perceived.smi")
        if candidate.is_file():
            smiles = candidate
    if smiles is None:
        sidecars = sorted(root.glob("*.perceived.smi"))
        smiles = sidecars[0] if sidecars else None
    return structure, smiles


def mol2_with_final_pdb_coordinates(
    final_pdb: Path, cgenff_mol2: Path | None, preferred_resname: str = "",
) -> tuple[str | None, str]:
    """Build a display MOL2 with final MD coordinates and CGenFF connectivity.

    GROMACS final PDB snapshots generally omit ``CONECT`` records.  Browsers
    therefore guess small-molecule bonds from distances, which loses aromatic
    bond orders and can make a ligand look chemically wrong.  CGenFF MOL2
    files contain the authoritative atom names and bonds, but their coordinates
    are pre-MD.  This joins those two sources without altering either input.

    Returns ``(mol2_text, residue_name)``.  The CHARMM LP1 virtual site is
    intentionally ignored: it is not a chemical atom and is not part of the
    MOL2 topology.
    """
    if cgenff_mol2 is None or not cgenff_mol2.is_file() or not final_pdb.is_file():
        return None, ""
    try:
        mol2_lines = cgenff_mol2.read_text(encoding="utf-8", errors="ignore").splitlines()
        pdb_lines = final_pdb.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return None, ""

    atom_start = next((index + 1 for index, line in enumerate(mol2_lines) if line.strip() == "@<TRIPOS>ATOM"), None)
    bond_start = next((index for index, line in enumerate(mol2_lines) if line.strip() == "@<TRIPOS>BOND"), None)
    if atom_start is None or bond_start is None or atom_start >= bond_start:
        return None, ""
    atom_records: list[tuple[int, list[str]]] = []
    for index in range(atom_start, bond_start):
        fields = mol2_lines[index].split()
        if not fields:
            continue
        if len(fields) < 6:
            return None, ""
        atom_records.append((index, fields))
    if not atom_records:
        return None, ""
    mol2_names = {fields[1].upper() for _, fields in atom_records}

    # Collect every HETATM residue, then select the requested residue if it is
    # complete.  If metadata has no useful ligand code, select the residue that
    # covers the greatest number of CGenFF atom names (normally its 3-letter
    # PDB resname, such as LIG).  This is deliberately a strict all-atom match.
    pdb_by_resname: dict[str, dict[str, tuple[float, float, float]]] = {}
    for line in pdb_lines:
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        resname = line[17:20].strip().upper()
        atom_name = line[12:16].strip().upper()
        if not resname or not atom_name:
            continue
        try:
            coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        pdb_by_resname.setdefault(resname, {})[atom_name] = coords
    preferred = "".join(ch for ch in preferred_resname if ch.isalnum()).upper()[:4]
    candidates = sorted(
        pdb_by_resname.items(),
        key=lambda item: (
            item[0] == preferred,
            len(mol2_names & set(item[1])),
            -len(item[1]),
        ),
        reverse=True,
    )
    for resname, coordinates in candidates:
        if not mol2_names.issubset(coordinates):
            continue
        rebuilt = list(mol2_lines)
        for line_index, fields in atom_records:
            x, y, z = coordinates[fields[1].upper()]
            rebuilt[line_index] = "{:<7} {:<8} {:>10.4f} {:>10.4f} {:>10.4f} {}".format(
                fields[0], fields[1], x, y, z, " ".join(fields[5:]),
            )
        return "\n".join(rebuilt) + "\n", resname
    return None, ""


def ligand_svg(structure: Path | None, smiles: Path | None) -> str | None:
    """Create a themed 2D SVG through the already-installed cgenff RDKit env."""
    if structure is None and smiles is None:
        return None
    python = Path.home() / "miniforge3" / "envs" / "cgenff" / "bin" / "python"
    if not python.is_file():
        return None
    source = textwrap.dedent("""
        import sys
        from rdkit import Chem
        from rdkit.Chem import AllChem, Draw
        structure, smiles = sys.argv[1:3]
        mol = None
        if smiles:
            try:
                line = open(smiles, encoding='utf-8', errors='ignore').readline().strip().split()
                mol = Chem.MolFromSmiles(line[0]) if line else None
            except Exception:
                pass
        if mol is None and structure:
            ext = structure.rsplit('.', 1)[-1].lower()
            mol = Chem.MolFromMol2File(structure, sanitize=True, removeHs=False) if ext == 'mol2' else Chem.MolFromMolFile(structure, sanitize=True, removeHs=False)
        if mol is None: raise SystemExit(2)
        AllChem.Compute2DCoords(mol)
        drawer = Draw.MolDraw2DSVG(680, 390)
        opts = drawer.drawOptions(); opts.clearBackground = False; opts.bondLineWidth = 2.6; opts.padding = .08
        opts.updateAtomPalette({1:(.62,.72,.84),6:(1.,.90,.20),7:(.21,.82,1.),8:(1.,.17,.17),9:(.49,1.,.61),15:(1.,.70,.28),16:(1.,.85,.30),17:(0.,1.,.40),35:(1.,.54,.24),53:(.83,.42,1.)})
        drawer.DrawMolecule(mol); drawer.FinishDrawing(); print(drawer.GetDrawingText())
    """)
    try:
        completed = subprocess.run(
            [str(python), "-c", source, str(structure or ""), str(smiles or "")],
            text=True, capture_output=True, timeout=45, check=True,
        )
        svg = completed.stdout.strip()
        return svg if "<svg" in svg else None
    except (OSError, subprocess.SubprocessError):
        return None


def ligand_html(root: Path, ligand_name: str) -> str:
    structure, smiles = ligand_inputs(root, ligand_name)
    svg = ligand_svg(structure, smiles)
    molecule_payload = ""
    molecule_format = ""
    if structure:
        try:
            molecule_payload = base64.b64encode(structure.read_bytes()).decode("ascii")
            molecule_format = {".mol2": "mol2", ".sdf": "sdf", ".mol": "sdf"}.get(structure.suffix.lower(), "mol2")
        except OSError:
            pass
    if not svg and not molecule_payload:
        return '<section class="ligand-card"><div class="eyebrow">Ligand identity</div><h2>Ligand structure</h2><p class="muted">A ligand chemistry file was not published with this legacy run. The 3D complex viewer and interaction panels remain available below.</p></section>'
    source_label = escape(structure.name if structure else "submitted ligand")
    two_d = svg or '<div class="missing-visual">2D depiction unavailable for this ligand file.</div>'
    viewer = (f'<div class="ligand-3d js-ligand-3d" data-format="{escape(molecule_format)}" data-structure-b64="{molecule_payload}"></div>' if molecule_payload else '<div class="missing-visual">3D ligand payload unavailable.</div>')
    smiles_text = ""
    if smiles and smiles.is_file():
        try:
            smiles_text = smiles.read_text(encoding="utf-8", errors="ignore").splitlines()[0].split()[0]
        except (OSError, IndexError):
            pass
    identifier = f'<details class="guide"><summary>Chemical identifier</summary><code>{escape(smiles_text)}</code></details>' if smiles_text else ""
    return f'''<section class="ligand-card"><div class="eyebrow">Ligand identity</div><h2>{escape(ligand_name)}</h2><p class="muted">2D chemical identity and the parameterized 3D geometry used for this run.</p><div class="ligand-grid"><div class="ligand-2d">{two_d}</div><div><div class="viewer-label">3D ligand · drag to rotate · scroll to zoom</div>{viewer}{identifier}</div></div><p class="asset-note">Source: {source_label}</p></section>'''


def interactive_networx_html(analysis: Path) -> str:
    """Embed the explorer source directly; never navigate to an attachment."""
    explorers = sorted((analysis / "NETWORX").glob("*_interactive.html"))
    if not explorers:
        return ""
    explorer = explorers[0]
    try:
        document = explorer.read_text(encoding="utf-8")
    except OSError:
        return ""
    # srcdoc makes the explorer part of this report.  Its CSS/IDs/scripts are
    # deliberately isolated, so the Overview and Contacts copies can coexist.
    return f'''<section class="gallery-section"><div class="section-head"><div><div class="eyebrow">Persistent interaction explorer</div><h2>Compare retained interactions</h2><p>Choose the displayed occupancy threshold to inspect the chemically classified residue–ligand-atom interactions. The larger embedded canvas supports node rearrangement and SVG export without opening a separate page.</p></div></div><iframe title="Interactive ligand interaction explorer" srcdoc="{escape(document, quote=True)}" sandbox="allow-scripts allow-downloads" style="display:block;width:100%;height:1320px;border:1px solid #285a79;border-radius:16px;background:#05070d" loading="lazy"></iframe></section>'''


def interface_rin_context(analysis: Path) -> dict:
    """Read the generated RIN explorer and retain one independently viewable panel per chain pair."""
    interface_dir = analysis / "Interface_RIN"
    explorer = interface_dir / "Interface_RIN_Explorer.html"
    if not explorer.is_file():
        return {"overview": "", "protein_detail": "", "peptide_detail": ""}
    try:
        document = explorer.read_text(encoding="utf-8")
        summary = pd.read_csv(interface_dir / "interface_chain_pair_summary.csv")
    except (OSError, ValueError, pd.errors.ParserError):
        return {"overview": "", "protein_detail": "", "peptide_detail": ""}
    required = {"ChainA_Label", "ChainB_Label"}
    if summary.empty or not required.issubset(summary.columns):
        return {"overview": "", "protein_detail": "", "peptide_detail": ""}

    chain_types: dict[str, str] = {}
    try:
        manifest = json.loads((interface_dir / "interface_analysis_manifest.json").read_text(encoding="utf-8"))
        for entry in manifest.get("chain_metadata", []):
            chain_types[str(entry.get("display_label", ""))] = str(entry.get("chain_type", ""))
    except (OSError, ValueError, TypeError):
        pass
    pairs = []
    for row in summary.drop_duplicates(["ChainA_Label", "ChainB_Label"]).itertuples(index=False):
        left, right = str(row.ChainA_Label), str(row.ChainB_Label)
        pair_id = f"{left}__{right}"
        types = {chain_types.get(left, "").casefold(), chain_types.get(right, "").casefold()}
        peptide = any("peptide" in item for item in types)
        pairs.append((pair_id, left, right, peptide))
    if not pairs:
        return {"overview": "", "protein_detail": "", "peptide_detail": ""}

    def embedded(pair_id: str, title: str, *, compact: bool) -> str:
        # The RIN exporter deliberately keeps every chain pair in one source
        # file.  Set its initial pair before embedding so a dossier panel is a
        # readable A:B view rather than a crowded all-chain graph.
        focused = document.replace(
            "let pair=DATA.pairs[0]?.id,threshold=.90,selected=null;",
            f"let pair={json.dumps(pair_id)},threshold=.90,selected=null;",
            1,
        )
        if compact:
            focused = focused.replace(".controls{display:flex", "#pairs{display:none}.controls{display:flex", 1)
        return (
            f'<iframe title="{escape(title)}" srcdoc="{escape(focused, quote=True)}" '
            'sandbox="allow-scripts" style="display:block;width:100%;height:880px;border:1px solid #285a79;'
            'border-radius:16px;background:#05070d" loading="lazy"></iframe>'
        )

    first_id, first_left, first_right, first_peptide = pairs[0]
    first_kind = "Peptide–protein" if first_peptide else "Protein–protein"
    overview = f'''<section class="gallery-section"><div class="section-head"><div><div class="eyebrow">Interface overview</div><h2>{escape(first_kind)} · {escape(first_left)} ↔ {escape(first_right)}</h2><p>This is the first analyzed interface. Each chain pair has its own panel in the Interface networks tab, avoiding a single unreadable all-protein graph.</p></div></div>{embedded(first_id, f'{first_left} to {first_right} interface explorer', compact=True)}</section>'''
    protein_panels, peptide_panels = [], []
    for index, (pair_id, left, right, peptide) in enumerate(pairs):
        kind = "Peptide–protein interface" if peptide else "Protein–protein interface"
        note = (
            "The complete peptide sequence is shown below the RIN; highlighted residues contact the protein at the active threshold."
            if peptide else
            "This pair is isolated from other chains so its persistent interface chemistry stays readable."
        )
        panel = (
            f'<details class="interface-pair"{" open" if index == 0 else ""}><summary>{escape(left)} ↔ {escape(right)} <span>{escape(kind)}</span></summary>'
            f"<p class=\"muted\">{escape(note)}</p>{embedded(pair_id, f'{left} to {right} interface explorer', compact=True)}</details>"
        )
        (peptide_panels if peptide else protein_panels).append(panel)

    def detail(title: str, eyebrow: str, text: str, panels: list[str]) -> str:
        if not panels:
            return ""
        return f'''<section class="gallery-section interface-gallery"><div class="section-head"><div><div class="eyebrow">{escape(eyebrow)}</div><h2>{escape(title)}</h2><p>{escape(text)}</p></div></div>{''.join(panels)}</section>'''

    return {
        "overview": overview,
        "protein_detail": detail("One network per protein pair", "Protein–protein interfaces", "Open a chain-pair container to inspect its contact persistence and validated interaction types without mixing other proteins into the graph.", protein_panels),
        "peptide_detail": detail("Peptide sequence and protein interface", "Peptide–protein interfaces", "Each panel keeps the full peptide sequence together, highlights contacting residues, and shows only its selected protein interface.", peptide_panels),
    }


def interface_evidence_html(analysis: Path, context: PartnerContext, rin: dict) -> str:
    """Make interface-network availability obvious; never silently hide a peptide relationship."""
    if not context.is_interface:
        return ""
    interface_dir = analysis / "Interface_RIN"
    explorer = interface_dir / "Interface_RIN_Explorer.html"
    has_maps = any((interface_dir / filename).is_file() for filename in (
        "Interface_RIN_Network.png", "Interface_ChainLevel_Network.png", "Interface_ResiduePair_ContactHeatmap.png",
    ))
    if explorer.is_file() or has_maps:
        previews = []
        for filename, label in (
            ("Interface_RIN_Network.png", "Residue interaction-node map"),
            ("Interface_ChainLevel_Network.png", "Chain-level interface map"),
            ("Interface_ResiduePair_ContactHeatmap.png", "Residue-pair occupancy heatmap"),
        ):
            path = interface_dir / filename
            try:
                encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            except OSError:
                continue
            previews.append(
                f'<figure class="interface-preview"><img alt="{escape(label)}" src="data:image/png;base64,{encoded}"><figcaption>{escape(label)}</figcaption></figure>'
            )
        static_maps = f'<div class="interface-preview-grid">{"".join(previews)}</div>' if previews else ""
        return (
            '<section class="gallery-section"><div class="section-head"><div><div class="eyebrow">Interface evidence</div>'
            f'<h2>{escape(context.relationship_label.title())} residue maps</h2>'
            f'<p>{"Start with the interactive residue-node network, then" if explorer.is_file() else "Use the published node maps, then"} use the occupancy heatmap to confirm which residue pairs persist. '
            'Each selected chain pair remains separate so unrelated interfaces are never merged.</p></div></div>'
            f'{static_maps}</section>'
        )
    expected = "Interface_RIN/Interface_RIN_Explorer.html and the interface residue-pair tables"
    return (
        '<section class="gallery-section"><div class="section-head"><div><div class="eyebrow">Interface status</div>'
        f'<h2>{escape(context.relationship_label.title())} analysis was not published</h2>'
        f'<p>The run is recognized as a {escape(context.system_class.replace("_", " "))}, but its residue-node maps are absent. '
        f'Expected analysis assets: <code>{escape(expected)}</code>. Re-run the interface RIN stage for the selected peptide/receptor chain pair; '
        'do not interpret this as evidence of no interaction.</p></div></div></section>'
    )


def final_structure_html(root: Path, context: PartnerContext) -> str:
    """Embed the final complex or apo PDB in the report's 3Dmol viewer."""
    candidates = (root / "Final_Trajectory.pdb", root / "binding_pocket_only.pdb")
    structure = next((path for path in candidates if path.is_file()), None)
    if structure is None:
        return (
            '<section class="structure-card"><div class="eyebrow">Final structure</div>'
            '<h2>3D complex / apo viewer</h2><p class="muted">No final PDB was published for this run.</p></section>'
        )
    try:
        payload = base64.b64encode(structure.read_bytes()).decode("ascii")
    except OSError:
        return ""
    ligand = "".join(ch for ch in context.partner_label if ch.isalnum())[:4] if context.partner_kind == "ligand" else ""
    final_ligand_payload = ""
    pdb_ligand_resname = ligand
    if context.partner_kind == "ligand":
        # Keep the final MD pose from the PDB while recovering the explicit
        # CGenFF bond graph.  If the metadata only says "Ligand", the strict
        # atom-name match below finds the actual 3-letter PDB residue (LIG is
        # the common fallback) without guessing a chemical identity.
        cgenff_mol2, _ = ligand_inputs(root, ligand)
        final_ligand, matched_resname = mol2_with_final_pdb_coordinates(structure, cgenff_mol2, ligand)
        if final_ligand:
            final_ligand_payload = base64.b64encode(final_ligand.encode("utf-8")).decode("ascii")
            pdb_ligand_resname = matched_resname
    chain_buttons = "".join(f'<button type="button" data-structure-view="chain" data-chain-id="{escape(chain_id)}">Chain {escape(chain_id)}</button>' for chain_id in context.partner_chain_ids)
    if context.partner_kind == "ligand":
        focused_button = '<button type="button" data-structure-view="pocket">Ligand + pocket</button>'
    elif context.partner_kind == "peptide":
        focused_button = '<button type="button" data-structure-view="partner">Peptide</button><button type="button" data-structure-view="interface">Peptide interface</button>'
    elif context.partner_kind == "protein_chain":
        focused_button = '<button type="button" data-structure-view="interface">Selected interface</button>'
    else:
        focused_button = ''
    protein_label = "Protein/receptor" if context.partner_kind == "peptide" else "Protein"
    bond_note = " Ligand bonds use the matching CGenFF topology." if final_ligand_payload else ""
    return f'''<section class="structure-card"><div class="section-head"><div><div class="eyebrow">Final structure</div><h2>3D complex / apo viewer</h2><p>Drag to rotate, scroll to zoom, and choose a focused view. The model is the final published PDB, not a representative image.{bond_note}</p></div></div><div class="structure-toolbar" role="group" aria-label="Structure display"><button type="button" data-structure-view="complex">Full complex</button><button type="button" data-structure-view="protein">{protein_label}</button>{chain_buttons}{focused_button}</div><div class="complex-3d js-complex-3d" data-ligand="{escape(ligand)}" data-pdb-ligand="{escape(pdb_ligand_resname)}" data-final-ligand-b64="{final_ligand_payload}" data-partner-kind="{escape(context.partner_kind)}" data-partner-chains="{escape(','.join(context.partner_chain_ids))}" data-structure-b64="{payload}"></div><p class="asset-note">Source: {escape(structure.name)} · final published coordinate snapshot.</p></section>'''


def read_metadata(root: Path) -> dict:
    try:
        return json.loads((root / "pymacs_run.json").read_text()).get("metadata", {})
    except (OSError, ValueError, TypeError):
        return {}


def numeric_plot(csv_file: Path, color: str) -> str | None:
    """Render an ordinary two-column analysis CSV as an interactive Plotly card."""
    # Analysis directories also contain intermediate/framewise/contact tables.
    # Those are valuable downloads, but are not automatically meaningful as a
    # simple XY trace.  Publishing them as empty cards is worse than omitting
    # them until their dedicated visualisation is available.
    stem = csv_file.stem.lower()
    unsupported = ("framewise", "filteredcontacts", "interactiontypes", "contactpersistence", "contact_map", "dssp_raw", "dssp_encoded", "dssp_chain", "fullcomplex")
    if any(token in stem for token in unsupported):
        return None
    try:
        frame = pd.read_csv(csv_file)
    except Exception:
        return None
    # ``DSSP_chain_B.csv`` can be an empty placeholder for a single-chain
    # system.  Its derived fraction table is 100% coil by construction, not a
    # statement about a ligand or a real second protein chain.  Do not publish
    # that misleading flat panel.
    if "sse_chain_" in stem and {"Helix", "Sheet", "Coil"}.issubset(frame.columns):
        helix = pd.to_numeric(frame["Helix"], errors="coerce").fillna(0).to_numpy()
        sheet = pd.to_numeric(frame["Sheet"], errors="coerce").fillna(0).to_numpy()
        coil = pd.to_numeric(frame["Coil"], errors="coerce").fillna(0).to_numpy()
        if len(helix) and np.allclose(helix, 0) and np.allclose(sheet, 0) and np.allclose(coil, 1):
            return None
    # Secondary-structure fraction files contain three categorical fractions;
    # plot all three with solid legend colors and a fixed 0–1 axis.  These are
    # population fractions (therefore legitimately intermediate values), not
    # per-residue categorical assignments.
    if "sse" in stem and "fraction" in stem and {"Time_ns", "Helix", "Sheet", "Coil"}.issubset(frame.columns):
        x = pd.to_numeric(frame["Time_ns"], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(x)
        if int(valid.sum()) < 3:
            return None
        if int(valid.sum()) > MAX_RENDER_POINTS:
            keep = np.linspace(0, int(valid.sum()) - 1, MAX_RENDER_POINTS, dtype=int)
            x = x[valid][keep]
            index = np.flatnonzero(valid)[keep]
        else:
            x = x[valid]
            index = np.flatnonzero(valid)
        figure = go.Figure()
        for column, trace_color in (("Helix", "#d32f2f"), ("Sheet", "#1e88e5"), ("Coil", "#ffd700")):
            y = pd.to_numeric(frame[column], errors="coerce").fillna(0).to_numpy(dtype=float)[index]
            figure.add_trace(go.Scatter(
                x=x, y=y, mode="lines", name=column,
                line={"color": trace_color, "width": 2},
                hovertemplate=f"Time=%{{x:.3f}} ns<br>{column}=%{{y:.3f}}<extra></extra>",
            ))
        figure.update_layout(
            title=csv_file.stem.replace("_", " "), paper_bgcolor=PANEL, plot_bgcolor=PANEL,
            font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 55, "r": 20, "t": 55, "b": 48},
            xaxis={"title": "Time (ns)", "gridcolor": GRID, "zerolinecolor": GRID},
            yaxis={"title": "Fraction of residues", "range": [0, 1], "gridcolor": GRID, "zerolinecolor": GRID},
            hovermode="x unified", height=360, legend={"orientation": "h", "y": 1.13},
        )
        return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})
    numeric = [column for column in frame.columns if pd.api.types.is_numeric_dtype(frame[column])]
    if not numeric:
        return None
    x_name = numeric[0]
    y_name = numeric[1] if len(numeric) > 1 else numeric[0]
    if len(numeric) == 1:
        x = list(range(len(frame)))
        x_name = "Frame"
    else:
        x = pd.to_numeric(frame[x_name], errors="coerce")
    y = pd.to_numeric(frame[y_name], errors="coerce")
    mask = np.isfinite(np.asarray(x, dtype=float)) & np.isfinite(np.asarray(y, dtype=float))
    if int(mask.sum()) < 3:
        return None
    x = np.asarray(x, dtype=float)[mask]
    y = np.asarray(y, dtype=float)[mask]
    if len(x) > MAX_RENDER_POINTS:
        # Evenly retain the complete time range, including the first and last
        # frame.  This avoids shipping hundreds of thousands of SVG points per
        # card while preserving the visible trajectory trend and hover data.
        keep = np.linspace(0, len(x) - 1, MAX_RENDER_POINTS, dtype=int)
        x, y = x[keep], y[keep]
    # Do not use Scattergl here.  A rich single-system report commonly has
    # 10–20 plots, and Chrome reclaims old WebGL canvases once its context
    # limit is reached; that makes earlier graphs appear and then go blank.
    # SVG Scatter remains interactive and is stable for the downsampled
    # analysis series exported by PyMACS.
    figure = go.Figure(go.Scatter(x=x, y=y, mode="lines", line={"color": color, "width": 2}, hovertemplate=f"{escape(x_name)}=%{{x}}<br>{escape(y_name)}=%{{y:.4g}}<extra></extra>"))
    figure.update_layout(
        title=csv_file.stem.replace("_", " "), paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 55, "r": 20, "t": 55, "b": 48},
        xaxis={"title": x_name, "gridcolor": GRID, "zerolinecolor": GRID},
        yaxis={"title": y_name, "gridcolor": GRID, "zerolinecolor": GRID},
        hovermode="x unified", height=360,
    )
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def _read_table(csv_file: Path) -> pd.DataFrame | None:
    try:
        # Frame-wise contact tables can contain millions of rows.  A bounded
        # scientific sample keeps the dashboard responsive; the exact CSV is
        # still linked from every card for full-resolution analysis.
        table = pd.read_csv(csv_file, nrows=MAX_TABLE_ROWS)
        return table if not table.empty else None
    except Exception:
        return None


def _sample_table(table: pd.DataFrame, maximum: int = MAX_RENDER_POINTS) -> pd.DataFrame:
    if len(table) <= maximum:
        return table
    return table.iloc[np.linspace(0, len(table) - 1, maximum, dtype=int)]


def multiseries_plot(csv_file: Path, title: str, x_column: str, y_columns: list[str],
                     y_title: str, colors: list[str] | None = None) -> str | None:
    """Render selected scientific columns as hoverable lines, never raw table rows."""
    table = _read_table(csv_file)
    if table is None or x_column not in table or not any(column in table for column in y_columns):
        return None
    table = table.copy()
    table[x_column] = pd.to_numeric(table[x_column], errors="coerce")
    table = table.dropna(subset=[x_column]).sort_values(x_column)
    table = _sample_table(table)
    if len(table) < 2:
        return None
    palette = colors or [CYAN, GREEN, MAGENTA, ORANGE, "#8b5cf6", "#facc15"]
    figure = go.Figure()
    for index, column in enumerate(y_columns):
        if column not in table:
            continue
        values = pd.to_numeric(table[column], errors="coerce")
        if not values.notna().any():
            continue
        figure.add_trace(go.Scatter(
            x=table[x_column], y=values, mode="lines", name=column.replace("_", " "),
            line={"color": palette[index % len(palette)], "width": 2},
            hovertemplate=f"{escape(x_column)}=%{{x:.4g}}<br>{escape(column)}=%{{y:.4g}}<extra></extra>",
        ))
    if not figure.data:
        return None
    figure.update_layout(
        title=title, paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"},
        margin={"l": 60, "r": 20, "t": 55, "b": 52}, hovermode="x unified", height=380,
        xaxis={"title": x_column.replace("_", " "), "gridcolor": GRID},
        yaxis={"title": y_title, "gridcolor": GRID},
        legend={"orientation": "h", "y": 1.13},
    )
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def interface_rin_plots(analysis: Path) -> list[tuple[str, str, Path, str]]:
    """Render each chain pair independently so multi-interface systems remain readable."""
    interface_dir = analysis / "Interface_RIN"
    time_path = interface_dir / "interface_contact_timeseries.csv"
    pair_path = interface_dir / "interface_residue_pair_contacts.csv"
    node_path = interface_dir / "interface_node_summary.csv"
    results: list[tuple[str, str, Path, str]] = []
    try:
        time_data = pd.read_csv(time_path)
    except (OSError, ValueError, pd.errors.ParserError):
        time_data = pd.DataFrame()
    if {"ChainA_Label", "ChainB_Label", "Time_ns"}.issubset(time_data.columns):
        for (left, right), table in time_data.groupby(["ChainA_Label", "ChainB_Label"], sort=False):
            table = _sample_table(table.sort_values("Time_ns"))
            for column, title, unit in (
                ("ContactCount", "interface contacts over time", "Residue-pair contacts"),
                ("MinDistance_A", "closest interface distance over time", "Closest distance (Å)"),
            ):
                if column not in table:
                    continue
                figure = go.Figure(go.Scatter(x=table["Time_ns"], y=pd.to_numeric(table[column], errors="coerce"), mode="lines", line={"color": CYAN if column == "ContactCount" else ORANGE, "width": 2}, hovertemplate=f"Time=%{{x:.3f}} ns<br>{escape(column)}=%{{y:.4g}}<extra></extra>"))
                figure.update_layout(title=f"{left} ↔ {right}: {title}", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 52}, height=360, xaxis={"title": "Time (ns)", "gridcolor": GRID}, yaxis={"title": unit, "gridcolor": GRID})
                results.append(("contacts" if column == "ContactCount" else "distances", f"{left} ↔ {right}: {title}", time_path, plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})))
    try:
        pairs = pd.read_csv(pair_path)
    except (OSError, ValueError, pd.errors.ParserError):
        pairs = pd.DataFrame()
    if {"ChainA_Label", "ChainB_Label", "Pair_ID", "ContactFraction"}.issubset(pairs.columns):
        for (left, right), table in pairs.groupby(["ChainA_Label", "ChainB_Label"], sort=False):
            table = table.nlargest(30, "ContactFraction").sort_values("ContactFraction")
            figure = go.Figure(go.Bar(y=table["Pair_ID"], x=table["ContactFraction"], orientation="h", marker_color=GREEN, hovertemplate="Residue pair=%{y}<br>Occupancy=%{x:.1%}<extra></extra>"))
            figure.update_layout(title=f"{left} ↔ {right}: top residue-pair persistence", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 160, "r": 20, "t": 55, "b": 52}, height=max(360, 30 * len(table)), xaxis={"title": "Contact occupancy", "range": [0, 1], "gridcolor": GRID}, yaxis={"title": "Residue pair", "gridcolor": GRID})
            results.append(("contacts", f"{left} ↔ {right}: top residue-pair persistence", pair_path, plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})))
    try:
        nodes = pd.read_csv(node_path)
    except (OSError, ValueError, pd.errors.ParserError):
        nodes = pd.DataFrame()
    if {"Chain_Label", "Residue_Label", "WeightedDegree"}.issubset(nodes.columns):
        for chain, table in nodes.groupby("Chain_Label", sort=False):
            table = table.nlargest(30, "WeightedDegree").sort_values("WeightedDegree")
            figure = go.Figure(go.Bar(y=table["Residue_Label"], x=table["WeightedDegree"], orientation="h", marker_color=MAGENTA, hovertemplate="Residue=%{y}<br>Weighted interface degree=%{x:.3f}<extra></extra>"))
            figure.update_layout(title=f"{chain}: interface-residue importance", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 130, "r": 20, "t": 55, "b": 52}, height=max(360, 30 * len(table)), xaxis={"title": "Weighted interface degree", "gridcolor": GRID}, yaxis={"title": "Residue", "gridcolor": GRID})
            results.append(("chemistry", f"{chain}: interface-residue importance", node_path, plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})))
    return results


def categorical_bar_plot(csv_file: Path, title: str, category: str, values: list[str],
                         y_title: str, stacked: bool = False,
                         colors: list[str] | None = None) -> str | None:
    table = _read_table(csv_file)
    if table is None or category not in table or not any(column in table for column in values):
        return None
    table = table.copy().head(80)
    figure = go.Figure()
    for index, column in enumerate(values):
        if column not in table:
            continue
        y = pd.to_numeric(table[column], errors="coerce").fillna(0)
        figure.add_trace(go.Bar(
            x=table[category].astype(str), y=y, name=column.replace("_", " "),
            marker_color=(colors or [CYAN, GREEN, MAGENTA, ORANGE, "#8b5cf6"])[index % len(colors or [CYAN, GREEN, MAGENTA, ORANGE, "#8b5cf6"])],
            hovertemplate=f"{escape(category)}=%{{x}}<br>{escape(column)}=%{{y:.4g}}<extra></extra>",
        ))
    if not figure.data:
        return None
    figure.update_layout(
        title=title, barmode="stack" if stacked else "group", paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 115},
        xaxis={"title": category.replace("_", " "), "gridcolor": GRID, "tickangle": -45},
        yaxis={"title": y_title, "gridcolor": GRID}, height=420, legend={"orientation": "h", "y": 1.13},
    )
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def contact_persistence_plot(events_file: Path, denominator_file: Path,
                             interaction_file: Path | None, protein_name: str) -> str | None:
    """Compute residue-labelled 0–1 persistence from framewise contact events.

    Older pipeline tables can contain local row indices and a micro-fraction.
    The framewise event files retain the original residue labels, so use those
    as the authority for both the x axis and hover text.
    """
    def residue_frames(path: Path) -> tuple[dict[str, set[float]], set[float]]:
        observed: dict[str, set[float]] = {}
        frames: set[float] = set()
        try:
            for chunk in pd.read_csv(path, usecols=["Time_ns", "Residue"], chunksize=100_000):
                chunk["Time_ns"] = pd.to_numeric(chunk["Time_ns"], errors="coerce")
                chunk = chunk.dropna(subset=["Time_ns", "Residue"])
                for residue, times in chunk.groupby("Residue")["Time_ns"]:
                    values = set(times.astype(float).tolist())
                    observed.setdefault(str(residue), set()).update(values)
                    frames.update(values)
        except Exception:
            return {}, set()
        return observed, frames

    observed, _ = residue_frames(events_file)
    _, all_frames = residue_frames(denominator_file)
    if not observed or not all_frames:
        return None
    type_colors = {"H-bond": CYAN, "Ionic": MAGENTA, "Hydrophobic": ORANGE}
    types: dict[str, str] = {}
    if interaction_file and interaction_file.is_file():
        try:
            interactions = pd.read_csv(interaction_file, usecols=["Residue", "Type"])
            for residue, values in interactions.dropna().groupby("Residue")["Type"]:
                counts = values.astype(str).value_counts()
                # Use the chemistry observed most often for this residue.
                # A previous "any H-bond wins" priority made transient
                # H-bonds visually dominate the persistence chart.
                types[str(residue)] = str(counts.index[0]) if not counts.empty else "Contact"
        except Exception:
            pass
    rows = sorted(((residue, len(times) / len(all_frames), types.get(residue, "Contact"))
                   for residue, times in observed.items()), key=lambda item: item[0])
    if not rows:
        return None
    residues, fractions, interaction_types = map(list, zip(*rows))
    colors = [type_colors.get(kind, GREEN) for kind in interaction_types]
    figure = go.Figure(go.Bar(
        x=residues, y=fractions, marker_color=colors,
        customdata=np.column_stack([interaction_types, np.asarray(fractions) * 100]),
        hovertemplate="Residue=%{x}<br>Contact persistence=%{y:.3f} (%{customdata[1]:.1f} of sampled frames)<br>Most frequent interaction chemistry=%{customdata[0]}<extra></extra>",
    ))
    figure.update_layout(
        title=f"{protein_name} ligand-contact persistence", paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"},
        margin={"l": 60, "r": 20, "t": 55, "b": 115}, height=420,
        xaxis={"title": f"{protein_name} residue", "gridcolor": GRID, "tickangle": -45},
        yaxis={"title": "Contact persistence (fraction of sampled frames)", "range": [0, 1], "gridcolor": GRID},
    )
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def distance_carousel_html(analysis: Path) -> str:
    """Render the precomputed residue set as native, hoverable Plotly bars."""
    directory = analysis / "RESIDUE_Contact_Distance"
    residues = [image.stem.removeprefix("dist_") for image in sorted(directory.glob("dist_*.png"))] if directory.is_dir() else []
    source = analysis / "AllContacts_Framewise.csv"
    if not residues or not source.is_file():
        return ""
    wanted = set(residues)
    samples: dict[str, list[float]] = {residue: [] for residue in residues}
    # Exact source data stays downloadable; retain a generous, deterministic
    # browser sample from every residue so publication remains lightweight.
    try:
        for chunk in pd.read_csv(source, usecols=["Residue", "Distance"], chunksize=100_000):
            chunk["Distance"] = pd.to_numeric(chunk["Distance"], errors="coerce")
            for residue, values in chunk.dropna().groupby("Residue")["Distance"]:
                residue = str(residue)
                if residue not in wanted or len(samples[residue]) >= 12_000:
                    continue
                samples[residue].extend(values.iloc[:12_000 - len(samples[residue])].astype(float).tolist())
    except Exception:
        return ""
    slides: list[str] = []
    for index, residue in enumerate(residues):
        values = np.asarray(samples[residue], dtype=float)
        if values.size < 2:
            continue
        counts, edges = np.histogram(values, bins=32)
        centres = (edges[:-1] + edges[1:]) / 2
        shares = counts / counts.sum() * 100
        figure = go.Figure(go.Bar(
            x=centres, y=counts, width=np.diff(edges) * 0.92, marker_color=CYAN,
            customdata=np.column_stack([edges[:-1], edges[1:], shares]),
            hovertemplate=(f"Residue={escape(residue)}<br>Distance interval=%{{customdata[0]:.2f}}–%{{customdata[1]:.2f}} Å"
                           "<br>Observed contact events=%{y}<br>Share of displayed events=%{customdata[2]:.2f}%<extra></extra>"),
        ))
        figure.update_layout(
            title=f"{residue} ligand–residue distance distribution", paper_bgcolor=PANEL, plot_bgcolor=PANEL,
            font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, height=480,
            margin={"l": 65, "r": 25, "t": 60, "b": 55},
            xaxis={"title": "Ligand–residue distance (Å)", "gridcolor": GRID},
            yaxis={"title": "Observed contact events", "gridcolor": GRID},
            annotations=[{"text": f"n={len(values):,} displayed events · mean={values.mean():.2f} Å · median={np.median(values):.2f} Å", "showarrow": False, "x": 0.5, "xref": "paper", "y": 1.12, "yref": "paper", "font": {"color": MUTED}}],
        )
        graph = plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})
        slides.append(f'<div class="distance-slide{" is-active" if not slides else ""}" data-distance-slide>{graph}</div>')
    if not slides:
        return ""
    source_link = source.name
    return f'''<style>.distance-carousel{{border-radius:21px;padding:20px;margin-top:18px;background:linear-gradient(155deg,#10172b,#060912);border:1px solid #23375d;box-shadow:0 18px 50px rgba(0,0,0,.35)}}.distance-carousel-controls{{display:flex;align-items:center;gap:10px;justify-content:center;margin:0 0 12px}}.distance-carousel-controls button{{border:1px solid #285a79;border-radius:10px;background:#0a253a;color:#cfe8ff;padding:9px 12px;font:800 .82rem inherit;cursor:pointer}}.distance-carousel-controls button:hover{{background:#116887;color:#fff;border-color:#35d0ff}}.distance-carousel-controls span{{min-width:150px;text-align:center;color:#7aa6d9;font-size:.88rem}}.distance-carousel-stage{{border:1px solid rgba(53,208,255,.28);border-radius:16px;background:#070b14;overflow:hidden}}.distance-slide{{display:none;margin:0;padding:7px}}.distance-slide.is-active{{display:block}}</style><section class="distance-carousel" aria-label="Per-residue ligand contact-distance distributions">
<div class="section-head"><div><div class="eyebrow">Interactive residue maps</div><h2>Residue interaction distance carousel</h2><p>Hover each cyan bar for its distance interval, observed-event count, and share of the displayed events. The mean and median summarize the selected residue above the plot.</p></div></div>
<div class="distance-carousel-controls"><button type="button" data-distance-prev aria-label="Previous residue map">← Previous</button><span data-distance-status></span><button type="button" data-distance-next aria-label="Next residue map">Next →</button></div>
<div class="distance-carousel-stage">{"".join(slides)}</div><a href="{escape(source_link)}" download style="display:inline-block;color:#35d0ff;font-size:.8rem;font-weight:800;margin:10px 3px 5px">Download exact source CSV for all residue distance distributions</a>{guide_html("contact_frequency")}</section><script>(()=>{{const root=document.currentScript.previousElementSibling,slides=[...root.querySelectorAll('[data-distance-slide]')],status=root.querySelector('[data-distance-status]');let index=0;const render=()=>{{slides.forEach((slide,i)=>slide.classList.toggle('is-active',i===index));status.textContent=`${{index+1}} / ${{slides.length}} · ${{slides[index].querySelector('.plotly-graph-div').layout.title.text}}`;requestAnimationFrame(()=>slides[index].querySelectorAll('.js-plotly-plot').forEach(plot=>Plotly.Plots.resize(plot)))}};root.querySelector('[data-distance-prev]').addEventListener('click',()=>{{index=(index-1+slides.length)%slides.length;render()}});root.querySelector('[data-distance-next]').addEventListener('click',()=>{{index=(index+1)%slides.length;render()}});render()}})();</script>'''


def contact_events_plot(csv_file: Path, title: str, distance_mode: bool = False) -> str | None:
    """Show individual calculated contact events with residue/atom hover context."""
    table = _read_table(csv_file)
    required = {"Time_ns", "Residue", "Distance"}
    if table is None or not required.issubset(table.columns):
        return None
    table = table.copy()
    table["Time_ns"] = pd.to_numeric(table["Time_ns"], errors="coerce")
    table["Distance"] = pd.to_numeric(table["Distance"], errors="coerce")
    table = table.dropna(subset=["Time_ns", "Distance", "Residue"])
    table = _sample_table(table, MAX_EVENT_POINTS)
    if table.empty:
        return None
    y = table["Distance"] if distance_mode else table["Residue"].astype(str)
    custom = np.column_stack([
        table.get("LigAtomName", pd.Series("–", index=table.index)).astype(str),
        table.get("ProtAtomName", pd.Series("–", index=table.index)).astype(str),
        table["Distance"].round(3),
    ])
    figure = go.Figure(go.Scattergl(
        x=table["Time_ns"], y=y, mode="markers", marker={"size": 5, "color": table["Distance"], "colorscale": "Turbo", "showscale": True, "colorbar": {"title": "Distance (Å)"}},
        customdata=custom,
        hovertemplate="Time=%{x:.3f} ns<br>Residue=%{y}<br>Ligand atom=%{customdata[0]}<br>Protein atom=%{customdata[1]}<br>Distance=%{customdata[2]} Å<extra></extra>",
    ))
    figure.update_layout(
        title=title, paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"},
        margin={"l": 85, "r": 35, "t": 55, "b": 52}, height=440,
        xaxis={"title": "Time (ns)", "gridcolor": GRID}, yaxis={"title": "Distance (Å)" if distance_mode else "Protein residue", "gridcolor": GRID},
    )
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def distance_histogram_plot(csv_file: Path) -> str | None:
    table = _read_table(csv_file)
    if table is None or "Distance" not in table:
        return None
    values = pd.to_numeric(table["Distance"], errors="coerce").dropna()
    if values.empty:
        return None
    figure = go.Figure(go.Histogram(x=values, nbinsx=70, marker_color=CYAN,
        hovertemplate="Distance=%{x:.3f} Å<br>Events=%{y}<extra></extra>"))
    figure.update_layout(title="Composite ligand–residue distance distribution", paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 52}, height=380,
        xaxis={"title": "Distance (Å)", "gridcolor": GRID}, yaxis={"title": "Observed events", "gridcolor": GRID})
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def interaction_degree_plot(csv_file: Path) -> str | None:
    table = _read_table(csv_file)
    if table is None or not {"Time_ns", "Residue"}.issubset(table.columns):
        return None
    degree = table.dropna(subset=["Time_ns", "Residue"]).groupby("Time_ns")["Residue"].nunique().reset_index(name="Interaction_degree")
    # Render directly to avoid creating an on-disk derived analysis artifact.
    figure = go.Figure(go.Scatter(x=degree["Time_ns"], y=degree["Interaction_degree"], mode="lines", line={"color": MAGENTA, "width": 2},
        hovertemplate="Time=%{x:.3f} ns<br>Unique contacting residues=%{y}<extra></extra>"))
    figure.update_layout(title="Interaction degree over time", paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 52}, height=380,
        xaxis={"title": "Time (ns)", "gridcolor": GRID}, yaxis={"title": "Unique contacting residues", "gridcolor": GRID})
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def persistence_distance_plot(contacts_file: Path, frequency_file: Path) -> str | None:
    contacts, frequency = _read_table(contacts_file), _read_table(frequency_file)
    if contacts is None or frequency is None or not {"Residue", "Distance"}.issubset(contacts.columns) or not {"Residue", "PercentFrames"}.issubset(frequency.columns):
        return None
    means = contacts.assign(Distance=pd.to_numeric(contacts["Distance"], errors="coerce")).groupby("Residue", as_index=False)["Distance"].mean()
    merged = means.merge(frequency[["Residue", "PercentFrames"]], on="Residue", how="inner").dropna()
    if merged.empty:
        return None
    figure = go.Figure(go.Scatter(x=merged["Distance"], y=merged["PercentFrames"], mode="markers+text", text=merged["Residue"], textposition="top center",
        marker={"size": 10, "color": merged["PercentFrames"], "colorscale": "Viridis", "showscale": True, "colorbar": {"title": "% frames"}},
        hovertemplate="Residue=%{text}<br>Mean distance=%{x:.3f} Å<br>Contact frequency=%{y:.2f}%<extra></extra>"))
    figure.update_layout(title="Contact persistence versus mean distance", paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 65, "r": 25, "t": 55, "b": 52}, height=410,
        xaxis={"title": "Mean ligand–residue distance (Å)", "gridcolor": GRID}, yaxis={"title": "Contact frequency (%)", "gridcolor": GRID})
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def sampled_indices(length: int, maximum: int = MAX_HEATMAP_TIMEPOINTS) -> np.ndarray:
    """Evenly select a matrix axis while retaining its beginning and end."""
    if length <= maximum:
        return np.arange(length, dtype=int)
    return np.linspace(0, length - 1, maximum, dtype=int)


def heatmap_plot(csv_file: Path, kind: str) -> str | None:
    """Render the scientifically meaningful matrix outputs as dedicated heatmaps."""
    try:
        frame = pd.read_csv(csv_file)
    except Exception:
        return None
    if frame.empty:
        return None
    if kind == "dssp":
        values = frame.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        if values.ndim != 2 or not np.isfinite(values).any():
            return None
        frame_idx = sampled_indices(values.shape[0])
        # DSSP rows are trajectory frames, columns are residue indices.
        z = values[frame_idx, :].T
        figure = go.Figure(go.Heatmap(
            z=z, x=frame_idx, y=np.arange(z.shape[0]), zmin=0, zmax=2,
            colorscale=[[0, "#ffd700"], [0.333, "#ffd700"], [0.334, "#1e88e5"], [0.666, "#1e88e5"], [0.667, "#d32f2f"], [1, "#d32f2f"]],
            colorbar={"title": "SSE", "tickvals": [0, 1, 2], "ticktext": ["Coil", "Sheet", "Helix"]},
            hovertemplate="Frame=%{x}<br>Residue=%{y}<br>SSE=%{z}<extra></extra>",
        ))
        title, x_title, y_title = "Secondary-structure evolution", "Sampled trajectory frame", "Residue index"
    elif kind == "contacts":
        if frame.shape[1] < 2:
            return None
        residues = frame.iloc[:, 0].astype(str).to_numpy()
        values = frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        if values.ndim != 2 or not np.isfinite(values).any():
            return None
        time_idx = sampled_indices(values.shape[1])
        # A residue either meets the contact definition in a frame/bin or it
        # does not.  Use a binary matrix and two solid colors rather than a
        # misleading continuous intensity scale.
        z = (values[:, time_idx] > 0).astype(int)
        raw_times = frame.columns[1:]
        x = pd.to_numeric(raw_times[time_idx], errors="coerce").to_numpy(dtype=float)
        figure = go.Figure(go.Heatmap(
            z=z, x=x, y=residues, zmin=0, zmax=1,
            colorscale=[[0, "#081326"], [0.4999, "#081326"], [0.5, "#ff2bd6"], [1, "#ff2bd6"]],
            colorbar={"title": "Contact", "tickmode": "array", "tickvals": [0, 1], "ticktext": ["No", "Yes"]},
            hovertemplate="Time=%{x:.3f} ns<br>Residue=%{y}<br>Contact present=%{z}<extra></extra>",
        ))
        title, x_title, y_title = "Ligand–residue contact presence over time", "Time (ns)", "Protein residue"
    elif kind == "interactions":
        # A residue can make more than one kind of ligand interaction during a
        # time bin.  A single heatmap cell must therefore use an explicit,
        # documented priority rather than implying that categories blend:
        # H-bond (3) > Ionic (2) > Hydrophobic (1) > no interaction (0).
        # This is the same priority contract used by 3A's static E4 heatmap.
        required = {"Time_ns", "Residue", "Type"}
        if not required.issubset(frame.columns):
            return None
        priority = {"Hydrophobic": 1, "Ionic": 2, "H-bond": 3}
        data = frame.loc[:, ["Time_ns", "Residue", "Type"]].copy()
        data["Time_ns"] = pd.to_numeric(data["Time_ns"], errors="coerce")
        data["InteractionCode"] = data["Type"].map(priority)
        data = data.dropna(subset=["Time_ns", "Residue", "InteractionCode"])
        if data.empty:
            return None
        data["Time_bin"] = (data["Time_ns"] / 2.0).round() * 2.0
        matrix = data.pivot_table(
            index="Residue", columns="Time_bin", values="InteractionCode",
            aggfunc="max", fill_value=0,
        )
        if matrix.empty:
            return None
        # Keep the browser response compact without changing the full CSV
        # export.  Two-ns bins plus sampled columns are sufficient for a
        # legible overview; hover still reports the exact displayed bin.
        time_idx = sampled_indices(matrix.shape[1])
        matrix = matrix.iloc[:, time_idx]
        z = matrix.to_numpy(dtype=float)
        labels = np.asarray(["None", "Hydrophobic", "Ionic", "H-bond"], dtype=object)
        customdata = labels[np.clip(z.astype(int), 0, 3)]
        # Deliberately stepped colors: interaction classes are categorical,
        # not a physical quantity to be interpolated as a gradient.
        discrete_scale = [
            [0.0, "#101827"], [0.2499, "#101827"],
            [0.25, ORANGE], [0.4999, ORANGE],
            [0.50, MAGENTA], [0.7499, MAGENTA],
            [0.75, CYAN], [1.0, CYAN],
        ]
        figure = go.Figure(go.Heatmap(
            z=z, x=matrix.columns.to_numpy(dtype=float), y=matrix.index.astype(str),
            zmin=0, zmax=3, colorscale=discrete_scale, customdata=customdata,
            colorbar={
                "title": "Interaction type", "tickmode": "array",
                "tickvals": [0, 1, 2, 3],
                "ticktext": ["None", "Hydrophobic", "Ionic", "H-bond"],
            },
            hovertemplate="Time bin=%{x:.1f} ns<br>Residue=%{y}<br>Type=%{customdata}<extra></extra>",
        ))
        title, x_title, y_title = "Ligand interaction type over time", "Time (ns; 2 ns bins)", "Protein residue"
    else:
        return None
    figure.update_layout(
        title=title, paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 75, "r": 25, "t": 55, "b": 48},
        xaxis={"title": x_title, "gridcolor": GRID, "zerolinecolor": GRID},
        yaxis={"title": y_title, "gridcolor": GRID, "zerolinecolor": GRID}, height=420,
    )
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def detailed_interaction_plot(csv_file: Path) -> str | None:
    """Show every observed interaction type without collapsing co-occurrences."""
    try:
        frame = pd.read_csv(csv_file)
    except Exception:
        return None
    required = {"Time_ns", "Residue", "Type"}
    if not required.issubset(frame.columns):
        return None
    codes = {"Hydrophobic": 1, "Ionic": 2, "H-bond": 3}
    data = frame.loc[:, ["Time_ns", "Residue", "Type"]].copy()
    data["Time_ns"] = pd.to_numeric(data["Time_ns"], errors="coerce")
    data["Code"] = data["Type"].map(codes)
    data = data.dropna(subset=["Time_ns", "Residue", "Code"])
    if data.empty:
        return None
    data["Time_bin"] = (data["Time_ns"] / 2.0).round() * 2.0
    data["Row"] = data["Residue"].astype(str) + " · " + data["Type"].astype(str)
    matrix = data.pivot_table(index="Row", columns="Time_bin", values="Code", aggfunc="max", fill_value=0)
    if matrix.empty:
        return None
    matrix = matrix.iloc[:, sampled_indices(matrix.shape[1])]
    z = matrix.to_numpy(dtype=float)
    colors = [[0.0, "#101827"], [0.2499, "#101827"], [0.25, ORANGE], [0.4999, ORANGE], [0.50, MAGENTA], [0.7499, MAGENTA], [0.75, CYAN], [1.0, CYAN]]
    figure = go.Figure(go.Heatmap(
        z=z, x=matrix.columns.to_numpy(dtype=float), y=matrix.index.astype(str), zmin=0, zmax=3,
        colorscale=colors, colorbar={"title": "Type", "tickmode": "array", "tickvals": [0, 1, 2, 3], "ticktext": ["None", "Hydrophobic", "Ionic", "H-bond"]},
        hovertemplate="Time bin=%{x:.1f} ns<br>Interaction=%{y}<br>Code=%{z}<extra></extra>",
    ))
    figure.update_layout(
        title="Specific ligand interaction events over time", paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 135, "r": 25, "t": 55, "b": 48},
        xaxis={"title": "Time (ns; 2 ns bins)", "gridcolor": GRID, "zerolinecolor": GRID},
        yaxis={"title": "Residue · interaction type", "gridcolor": GRID, "zerolinecolor": GRID}, height=520,
    )
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def dashboard_card(graph: str, name: str, download_href: str = "", download_label: str = "Download CSV", note: str = "") -> str:
    """Wrap each analysis graphic in a short, collapsible scientific guide."""
    link = f'<a href="{escape(download_href)}" download>{escape(download_label)}</a>' if download_href else ""
    note_html = f'<p class="heatmap-note">{escape(note)}</p>' if note else ""
    return f'<article class="card">{graph}{guide_html(graph_guide_key(name))}{note_html}{link}</article>'


def section_html(title: str, eyebrow: str, text: str, cards: list[str]) -> str:
    if not cards:
        return ""
    return f'<section class="gallery-section"><div class="section-head"><div><div class="eyebrow">{escape(eyebrow)}</div><h2>{escape(title)}</h2><p>{escape(text)}</p></div></div><div class="grid">{"".join(cards)}</div></section>'


def tabbed_sections(sections: list[tuple[str, str, str]]) -> str:
    """Build one dashboard surface with tabs, without duplicating plot markup."""
    visible = [(ident, label, body) for ident, label, body in sections if body]
    if not visible:
        return ""
    tabs = "".join(
        f'<button type="button" role="tab" id="tab-{ident}" aria-controls="panel-{ident}" '
        f'aria-selected="{str(index == 0).lower()}" data-dashboard-tab="{ident}">{escape(label)}</button>'
        for index, (ident, label, _) in enumerate(visible)
    )
    panels = "".join(
        f'<section id="panel-{ident}" role="tabpanel" aria-labelledby="tab-{ident}" '
        f'class="dashboard-panel{" is-active" if index == 0 else ""}"{ "" if index == 0 else " hidden"}>{body}</section>'
        for index, (ident, _, body) in enumerate(visible)
    )
    return f'<section id="analysis-tabs" class="dashboard-tabs-shell"><div class="dashboard-tabs" role="tablist" aria-label="Result analyses">{tabs}</div>{panels}</section>'


def build_document(system_name: str, context: PartnerContext, cards: list[str], dashboard: str) -> str:
    """Return the Kyle-hosted, self-contained single-run reading room."""
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{escape(system_name)} · PyMACS results</title><script>{plotlyjs()}</script><script src="/static/3Dmol-min.js"></script><style>
:root{{--bg:#05060a;--panel:#0b1020;--grid:#1c2a4a;--text:#cfe8ff;--muted:#7aa6d9;--cyan:#35d0ff;--green:#18ff8b;--yellow:#ffe66d;--magenta:#ff2bd6;--orange:#ff8a3d}}*{{box-sizing:border-box}}html{{scroll-behavior:smooth}}body{{margin:0;background:radial-gradient(circle at 90% 0,rgba(53,208,255,.15),transparent 30%),radial-gradient(circle at 4% 25%,rgba(255,43,214,.10),transparent 27%),var(--bg);color:var(--text);font-family:Inter,ui-sans-serif,system-ui,-apple-system,sans-serif}}.page{{max-width:1560px;margin:auto;padding:32px 24px 70px}}.hero,.card,.ligand-card,.gallery-section,.structure-card,.dashboard-tabs-shell{{background:linear-gradient(155deg,#10172b,#060912);border:1px solid #23375d;box-shadow:0 18px 50px rgba(0,0,0,.35)}}.hero{{position:relative;overflow:hidden;border-radius:25px;padding:34px;margin-bottom:20px}}.hero:before{{content:"";position:absolute;inset:-48%;background:conic-gradient(from 180deg,transparent,rgba(53,208,255,.16),transparent,rgba(255,43,214,.13),transparent);animation:slow-spin 28s linear infinite;pointer-events:none}}.hero>*{{position:relative;z-index:1}}@keyframes slow-spin{{to{{transform:rotate(360deg)}}}}.eyebrow{{color:var(--green);font-size:.75rem;font-weight:900;letter-spacing:.16em;text-transform:uppercase}}h1{{font-size:clamp(2.5rem,5.2vw,5.2rem);letter-spacing:-.065em;line-height:.92;margin:.35rem 0 1rem}}h2{{letter-spacing:-.025em;margin:.2rem 0 .45rem;font-size:1.5rem}}p{{line-height:1.55}}.sub,.section-head p,.muted,.asset-note{{color:var(--muted);max-width:960px}}.pills,.portal-actions{{display:flex;gap:9px;flex-wrap:wrap;margin-top:18px}}.pill,.portal-action{{border:1px solid rgba(53,208,255,.28);padding:8px 11px;border-radius:999px;background:rgba(53,208,255,.08);font-size:.88rem}}.portal-actions{{margin-top:12px}}.portal-action{{color:var(--text);font-weight:800;text-decoration:none;transition:transform .16s ease,background .16s ease}}.portal-action:hover{{transform:translateY(-1px);background:rgba(53,208,255,.20)}}.portal-action--primary{{background:linear-gradient(90deg,#18b7e6,#1988e6);border-color:transparent;color:#03101b}}details.guide,.overview{{margin:11px 0 2px;border:1px solid rgba(53,208,255,.22);border-radius:13px;background:rgba(5,10,23,.68)}}details.guide summary,.overview summary{{cursor:pointer;padding:10px 12px;color:var(--cyan);font-weight:800;font-size:.88rem}}details.guide p,.overview p{{margin:0;padding:0 12px 11px;color:var(--muted);font-size:.87rem}}.overview{{margin-top:22px;max-width:1020px}}.overview ol{{margin:0;padding:0 28px 15px;color:var(--muted);line-height:1.55}}.gallery-section,.ligand-card,.structure-card{{border-radius:21px;padding:20px;margin-top:18px}}.section-head{{border-left:3px solid var(--green);padding-left:13px;margin-bottom:16px}}.interface-pair{{margin:12px 0;border:1px solid rgba(53,208,255,.24);border-radius:15px;background:rgba(5,10,23,.58);padding:0 14px 14px}}.interface-pair summary{{cursor:pointer;padding:14px 2px 10px;color:var(--cyan);font-weight:900;font-size:1rem}}.interface-pair summary span{{color:var(--muted);font-size:.8rem;font-weight:700;margin-left:8px}}.interface-preview-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px}}.interface-preview{{margin:0;border:1px solid rgba(53,208,255,.25);border-radius:14px;overflow:hidden;background:#05070d}}.interface-preview img{{display:block;width:100%;height:250px;object-fit:contain;background:#fff}}.interface-preview figcaption{{padding:9px 11px;color:var(--muted);font-size:.84rem;font-weight:800}}.dashboard-tabs-shell{{margin-top:18px;border-radius:21px;overflow:hidden}}.dashboard-tabs{{display:flex;gap:7px;overflow:auto;padding:14px 16px 0;border-bottom:1px solid #23375d}}.dashboard-tabs button,.structure-toolbar button{{border:1px solid #285a79;border-radius:10px 10px 0 0;background:#0a253a;color:#cfe8ff;padding:10px 13px;font:800 .82rem inherit;white-space:nowrap;cursor:pointer}}.dashboard-tabs button[aria-selected="true"],.dashboard-tabs button:hover{{background:#116887;color:#fff;border-color:var(--cyan)}}.dashboard-panel{{padding:2px 20px 20px}}.dashboard-panel[hidden]{{display:none}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:16px}}.card{{border-radius:17px;padding:13px;overflow:hidden}}.card .js-plotly-plot{{cursor:zoom-in}}.card a{{display:inline-block;color:var(--cyan);font-size:.8rem;font-weight:800;margin:8px 10px 5px}}.heatmap-note{{margin:7px 10px;color:var(--muted);font-size:.82rem;line-height:1.4}}.ligand-grid{{display:grid;grid-template-columns:minmax(300px,1fr) minmax(300px,1fr);gap:18px;align-items:stretch}}.ligand-2d,.ligand-3d,.complex-3d{{min-height:330px;border:1px solid rgba(255,230,109,.30);border-radius:17px;background:radial-gradient(circle at center,rgba(255,230,109,.08),rgba(5,6,10,.7) 70%);overflow:hidden;position:relative}}.ligand-2d svg{{display:block;width:100%;height:330px;max-width:100%}}.ligand-3d,.complex-3d{{border-color:rgba(53,208,255,.30)}}.complex-3d{{height:min(62vw,650px);min-height:410px}}.structure-toolbar{{display:flex;gap:8px;flex-wrap:wrap;margin:0 0 11px}}.structure-toolbar button{{border-radius:10px;padding:8px 11px}}.structure-toolbar button:hover,.structure-toolbar button.is-active{{background:#116887;color:#fff;border-color:var(--cyan)}}.viewer-label{{color:var(--yellow);font-size:.78rem;font-weight:900;text-transform:uppercase;letter-spacing:.1em;margin:2px 0 8px}}.missing-visual{{padding:24px;color:var(--muted)}}.asset-note{{font-size:.8rem;margin:.7rem 0 0}}code{{display:block;white-space:pre-wrap;overflow-wrap:anywhere;padding:0 12px 12px;color:var(--yellow)}}.graph-modal{{display:none;position:fixed;z-index:9999;inset:0;padding:3vh 3vw;background:rgba(0,0,0,.78);backdrop-filter:blur(7px)}}.graph-modal.open{{display:flex;align-items:center;justify-content:center}}.graph-modal-panel{{position:relative;width:min(1500px,94vw);height:min(900px,92vh);padding:20px 16px 12px;border:1px solid var(--cyan);background:#080d1b;box-shadow:0 0 70px rgba(53,208,255,.18);border-radius:20px}}.graph-modal-close{{position:absolute;right:17px;top:13px;z-index:2;border:1px solid var(--cyan);border-radius:999px;background:#101b33;color:var(--text);font-weight:800;padding:8px 12px;cursor:pointer}}.graph-modal-plot{{width:100%;height:100%}}@media(max-width:700px){{.page{{padding:18px 12px 42px}}.grid,.ligand-grid{{grid-template-columns:1fr}}.hero{{padding:24px 19px}}.dashboard-panel{{padding:2px 12px 12px}}.complex-3d{{height:68vw;min-height:330px}}.graph-modal{{padding:1vh 1vw}}.graph-modal-panel{{width:98vw;height:96vh;padding:14px 4px 4px}}}}
</style></head><body><main class="page"><section class="hero"><div class="eyebrow">PyMACS · single-system results</div><h1>{escape(system_name)}</h1><p class="sub">{escape({"ligand": "Start with global stability, then inspect ligand pose and pocket chemistry.", "peptide": "Start with global stability, then inspect peptide conformation and peptide–protein interface persistence.", "protein_chain": "Start with global and chain stability, then inspect persistent chain interfaces.", "none": "Start with global stability, fold retention, and secondary-structure evolution."}[context.partner_kind])} Every graph is a native Plotly analysis: hover for values, use the legend to isolate traces, and double-click for a large view.</p><div class="pills"><span class="pill">{escape({"ligand": "Ligand / PROTAC", "peptide": "Peptide", "protein_chain": "Chains / interfaces", "none": "Protein / chain"}[context.partner_kind])}: {escape(context.partner_label)}</span><span class="pill">{len(cards)} interactive analyses</span><span class="pill">Kyle-local, fast delivery</span></div><!--PORTAL_ACTIONS--><details class="overview"><summary>ⓘ How to read this result—start here</summary><ol><li><b>Check the fold:</b> protein RMSD, radius of gyration, and secondary-structure evolution should be considered together.</li><li><b>Check the partner:</b> interpret {escape(context.relationship_label)} RMSD/RMSF and contact persistence together—not in isolation.</li><li><b>Check the chemistry:</b> solid-color interaction types show which residues make hydrophobic, ionic, or hydrogen-bond interactions over time.</li><li><b>Do not overclaim:</b> these trajectories test stability of a modeled system; they do not by themselves establish binding affinity or experimental activity.</li></ol></details></section>{dashboard}</main><div class="graph-modal" id="graph-modal" aria-hidden="true"><div class="graph-modal-panel"><button class="graph-modal-close" id="graph-modal-close" type="button">Close ×</button><div class="graph-modal-plot" id="graph-modal-plot"></div></div></div><script>
(() => {{
 const modal=document.getElementById('graph-modal'),destination=document.getElementById('graph-modal-plot');
 const close=()=>{{modal.classList.remove('open');modal.setAttribute('aria-hidden','true');Plotly.purge(destination)}};
 const open=source=>{{const layout=Object.assign({{}},source.layout,{{autosize:true,height:null,margin:{{l:70,r:35,t:65,b:60}}}});Plotly.react(destination,source.data,layout,{{responsive:true,displaylogo:false}});modal.classList.add('open');modal.setAttribute('aria-hidden','false');requestAnimationFrame(()=>Plotly.Plots.resize(destination))}};
 document.querySelectorAll('.js-plotly-plot').forEach(source=>source.addEventListener('dblclick',event=>{{event.preventDefault();event.stopImmediatePropagation();open(source)}},true));
 document.getElementById('graph-modal-close').addEventListener('click',close);modal.addEventListener('click',event=>{{if(event.target===modal)close()}});document.addEventListener('keydown',event=>{{if(event.key==='Escape')close()}});
 document.querySelectorAll('[data-dashboard-tab]').forEach(button=>button.addEventListener('click',()=>{{const target=button.dataset.dashboardTab;document.querySelectorAll('[data-dashboard-tab]').forEach(item=>item.setAttribute('aria-selected',String(item===button)));document.querySelectorAll('.dashboard-panel').forEach(panel=>{{const active=panel.id===`panel-${{target}}`;panel.hidden=!active;panel.classList.toggle('is-active',active);if(active)panel.querySelectorAll('.js-plotly-plot').forEach(plot=>Plotly.Plots.resize(plot))}})}}));
 document.querySelectorAll('.js-ligand-3d').forEach(el=>{{try{{const raw=atob(el.dataset.structureB64),viewer=$3Dmol.createViewer(el,{{backgroundColor:'#05060a',antialias:true}});viewer.addModel(raw,el.dataset.format);viewer.setStyle({{}},{{stick:{{radius:.105,color:'#ffe66d'}},sphere:{{scale:.16,color:'#ffe66d'}}}});viewer.setStyle({{elem:'N'}},{{stick:{{radius:.105,color:'#35d0ff'}},sphere:{{scale:.16,color:'#35d0ff'}}}});viewer.setStyle({{elem:'O'}},{{stick:{{radius:.105,color:'#ff3b5c'}},sphere:{{scale:.16,color:'#ff3b5c'}}}});viewer.setStyle({{elem:'Cl'}},{{stick:{{radius:.105,color:'#18ff8b'}},sphere:{{scale:.20,color:'#18ff8b'}}}});viewer.zoomTo();viewer.render();addEventListener('resize',()=>viewer.resize())}}catch(_){{el.textContent='3D ligand preview unavailable.'}}}});
 document.querySelectorAll('.js-complex-3d').forEach(el=>{{try{{const viewer=$3Dmol.createViewer(el,{{backgroundColor:'#05060a',antialias:true}}),raw=atob(el.dataset.structureB64),ligand=el.dataset.ligand,pdbLigand=el.dataset.pdbLigand||ligand,finalLigand=el.dataset.finalLigandB64?atob(el.dataset.finalLigandB64):'',chains=el.dataset.partnerChains.split(',').filter(Boolean),pdbModel=viewer.addModel(raw,'pdb'),ligandModel=finalLigand?viewer.addModel(finalLigand,'mol2'):null;const ligandStyle={{stick:{{radius:.20,colorscheme:'Jmol'}},sphere:{{scale:.25,colorscheme:'Jmol'}}}},clear=()=>{{pdbModel.setStyle({{}},{{}});if(ligandModel)ligandModel.setStyle({{}},{{}})}},paintLigand=()=>{{if(ligandModel)ligandModel.setStyle({{}},ligandStyle)}};const render=(mode,chain='')=>{{clear();pdbModel.setStyle({{hetflag:false}},{{cartoon:{{color:'spectrum'}}}});if(mode==='complex'){{pdbModel.setStyle({{hetflag:true}},{{stick:{{radius:.16,colorscheme:'Jmol'}}}});if(ligandModel)pdbModel.setStyle({{resn:pdbLigand}},{{}});paintLigand();viewer.zoomTo()}}if(mode==='chain'&&chain){{pdbModel.setStyle({{}},{{cartoon:{{color:'#34495e',opacity:.18}}}});pdbModel.setStyle({{chain}},{{cartoon:{{color:'#35d0ff'}},stick:{{radius:.10}}}});viewer.zoomTo({{chain}})}}if(mode==='partner'&&chains.length){{pdbModel.setStyle({{}},{{cartoon:{{color:'#34495e',opacity:.18}}}});pdbModel.setStyle({{chain:chains[0]}},{{cartoon:{{color:'#ffe66d'}},stick:{{radius:.10}}}});viewer.zoomTo({{chain:chains[0]}})}}if(mode==='interface'&&chains.length>1){{const sel={{chain:chains[0]}};pdbModel.setStyle({{}},{{cartoon:{{color:'#34495e',opacity:.22}}}});pdbModel.setStyle({{chain:chains[0]}},{{cartoon:{{color:'#35d0ff'}}}});pdbModel.setStyle({{within:{{distance:5,sel}}}},{{cartoon:{{color:'#ffe66d'}},stick:{{radius:.10}}}});viewer.zoomTo({{within:{{distance:7,sel}}}})}}if(mode==='pocket'&&pdbLigand){{pdbModel.setStyle({{}},{{cartoon:{{color:'#34495e',opacity:.28}}}});pdbModel.setStyle({{within:{{distance:5,sel:{{resn:pdbLigand}}}}}},{{cartoon:{{color:'#35d0ff'}},stick:{{radius:.10}}}});pdbModel.setStyle({{resn:pdbLigand}},{{}});paintLigand();viewer.zoomTo({{resn:pdbLigand}})}}if(mode==='protein')viewer.zoomTo();viewer.render();el.parentElement.querySelectorAll('[data-structure-view]').forEach(button=>button.classList.toggle('is-active',button.dataset.structureView===mode&&(mode!=='chain'||button.dataset.chainId===chain)))}};el.parentElement.querySelectorAll('[data-structure-view]').forEach(button=>button.addEventListener('click',()=>render(button.dataset.structureView,button.dataset.chainId||'')));render('complex');addEventListener('resize',()=>viewer.resize())}}catch(_){{el.textContent='3D final-structure preview unavailable.'}}}});
}})();
</script></body></html>'''


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a one-run interactive PyMACS result dashboard.")
    parser.add_argument("--root", default=".", help="MD system directory")
    parser.add_argument("--out", default="Analysis_Results/PyMACS_Result_Dashboard.html", help="HTML destination relative to root")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    analysis = root / "Analysis_Results"
    output = root / args.out
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = read_metadata(root)
    system_name = str(metadata.get("system_name") or root.name)
    interface_rin = interface_rin_context(analysis)
    has_interface_rin = bool(interface_rin["protein_detail"] or interface_rin["peptide_detail"])
    setup = {}
    simulation = {}
    try:
        run_config = json.loads((root / "pymacs_run.json").read_text())
        setup = run_config.get("setup", {})
        simulation = run_config.get("simulation", {})
    except (OSError, ValueError, TypeError):
        pass
    context = partner_context(root, analysis, setup, metadata, simulation)
    ligand_name = context.partner_label if context.partner_kind == "ligand" else "No ligand label"
    chain_labels = []
    for entry in str(setup.get("chain_map") or "").split(","):
        chain, separator, label = entry.strip().partition(":")
        if chain and separator and label:
            chain_labels.append(label.strip())
    protein_name = " / ".join(component_display_labels(root, context)) or " / ".join(chain_labels) or "Protein"
    # Full-complex RMSD is meaningful for multiple macromolecular chains or
    # a PROTAC assembly.  A small-molecule ligand is not an extra protein
    # chain, so a one-chain ligand–protein job must not show this card.
    is_multichain_system = context.partner_kind == "protein_chain" or len(chain_labels) > 1 or str(setup.get("mode") or "").casefold() == "protac"
    csv_files = sorted(path for path in analysis.rglob("*.csv") if path.is_file())
    # Every card below is explicitly mapped to a scientific table.  Never
    # revive the legacy arbitrary-column fallback for raw analysis tables.
    gallery_csvs: set[Path] = set()
    colours = [CYAN, GREEN, MAGENTA, ORANGE]
    heatmap_cards = []
    groups: dict[str, list[str]] = {"overview": [], "stability": [], "flexibility": [], "secondary": [], "contacts": [], "chemistry": [], "distances": [], "coupled": [], "structure": []}
    def add(group: str, graph: str | None, title: str, source: Path, label: str = "Download exact source CSV", note: str = "") -> None:
        if graph:
            groups[group].append(dashboard_card(graph, title, source.relative_to(output.parent).as_posix(), label, note))
    # Matrix outputs are promoted to dedicated scientific heatmaps.  They are
    # deliberately sampled along time so the report stays quick to publish and
    # open; the full matrix CSV remains downloadable from the card.
    dssp_file = analysis / "DSSP_encoded_numeric.csv"
    if dssp_file.is_file():
        graph = heatmap_plot(dssp_file, "dssp")
        if graph:
            heatmap_cards.append(dashboard_card(graph, dssp_file.stem, dssp_file.relative_to(output.parent).as_posix(), "Download full SSE matrix CSV"))
    interaction_file = analysis / "InteractionTypes_Framewise.csv"
    if interaction_file.is_file():
        graph = heatmap_plot(interaction_file, "interactions")
        if graph:
            heatmap_cards.append(dashboard_card(
                graph, interaction_file.stem, interaction_file.relative_to(output.parent).as_posix(),
                "Download full interaction-event CSV",
                "Solid colors are categorical: H-bond takes precedence over ionic, then hydrophobic, when more than one interaction occurs in a residue/time bin.",
            ))
        detailed = detailed_interaction_plot(interaction_file)
        if detailed:
            heatmap_cards.append(dashboard_card(
                detailed, "Specific ligand interaction events",
                interaction_file.relative_to(output.parent).as_posix(),
                "Download full interaction-event CSV",
                "Every residue/type pair has its own row, so simultaneous interaction types are visible rather than collapsed into a single priority class.",
            ))
    for contact_file in sorted(analysis.glob("*_contact_map_residue_time.csv")):
        graph = heatmap_plot(contact_file, "contacts")
        if graph:
            heatmap_cards.append(dashboard_card(graph, contact_file.stem, contact_file.relative_to(output.parent).as_posix(), "Download full contact matrix CSV"))

    # The former portal tabs are all rendered here from their scientific CSVs.
    # Each card selects the semantic dimensions of its source rather than
    # joining arbitrary numeric columns into a generic line chart.
    stability_specs = [
        ("Protein_RMSD.csv", "Protein RMSD", "Time_ps", ["Protein_RMSD"], "RMSD (Å)"),
        ("LIG_Complex_RMSD_Overlay.csv", "Protein–ligand RMSD overlay", "Time_ns", ["Protein_RMSD", "Ligand_RMSD"], "RMSD (Å)"),
        ("Radius_of_Gyration_Protein_Ligand_Complex.csv", "Radius of gyration", "Time_ns", ["Rg_Protein_Å", "Rg_Ligand_Å", "Rg_Complex_Å"], "Radius of gyration (Å)"),
        ("native_contact_q.csv", f"{context.relationship_label.title()} Q", "Time_ns", ["Q_native_contact_fraction"], "Contact fraction"),
        ("protein_native_contact_q.csv", "Protein fold/native-contact Q", "Time_ns", ["Q_protein_native_contact_fraction"], "Contact fraction"),
    ]
    if is_multichain_system:
        stability_specs.insert(1, ("FullComplex_RMSD.csv", "Full-complex RMSD", "Time_ps", ["RMSD_Å"], "RMSD (Å)"))
    for file_name, title, x, ys, y_title in stability_specs:
        source = analysis / file_name
        colors = [GREEN, ORANGE, MAGENTA] if "overlay" in title.lower() or "gyration" in title.lower() else [GREEN]
        if source.is_file(): add("stability", multiseries_plot(source, title.replace("Protein", protein_name), x, ys, y_title, colors=colors), title, source)
    protein_rmsd = analysis / "Protein_RMSD.csv"
    if protein_rmsd.is_file():
        add("overview", multiseries_plot(protein_rmsd, f"{protein_name} RMSD", "Time_ps", ["Protein_RMSD"], "RMSD (Å)", colors=[GREEN]), f"{protein_name} RMSD", protein_rmsd)
    if interaction_file.is_file():
        add("overview", heatmap_plot(interaction_file, "interactions"), "Interaction chemistry overview", interaction_file, "Download full interaction-event CSV")
    for source in sorted(analysis.glob("RMSD_*.csv")):
        if source.name not in {"RMSD_BTK.csv"}: continue
        add("stability", multiseries_plot(source, f"{protein_name} chain RMSD", "Time_ps", ["RMSD_Å"], "RMSD (Å)", colors=[GREEN]), f"{protein_name} chain RMSD", source)
    for source in sorted(analysis.glob("*_Ligand_RMSD.csv")):
        partner_rmsd = "Peptide RMSD" if context.partner_kind == "peptide" else "Ligand RMSD"
        add("stability", multiseries_plot(source, partner_rmsd, "Time_ps", ["Ligand_RMSD"], "RMSD (Å)", colors=[ORANGE]), partner_rmsd, source)
    peptide_rmsd = analysis / "Peptide_RMSD.csv"
    if peptide_rmsd.is_file():
        add("stability", multiseries_plot(peptide_rmsd, f"{context.partner_label} peptide RMSD", "Time_ps", ["RMSD_A"], "RMSD (Å)", colors=[ORANGE]), "Peptide RMSD", peptide_rmsd,
            note="Peptide coordinates are aligned on the peptide selection. Interpret this with the interface persistence panel rather than as a binding affinity measure.")

    for source in sorted(analysis.glob("RMSF_*_per_residue.csv")):
        add("flexibility", categorical_bar_plot(source, f"{protein_name} RMSF per residue", "Residue", ["RMSF_Å"], "RMSF (Å)", colors=[GREEN]), f"{protein_name} RMSF per residue", source)
    for source in sorted(analysis.glob("RMSF_*.csv")):
        if source.name.endswith("_per_residue.csv"): continue
        add("flexibility", categorical_bar_plot(source, f"{protein_name} backbone-atom RMSF", "Atom_Index", ["RMSF_Å"], "RMSF (Å)", colors=[GREEN]), f"{protein_name} backbone-atom RMSF", source, note="This is the atom-resolved backbone selection used by the analysis; use the adjacent residue plot for residue-level interpretation.")
    for source in sorted(analysis.glob("*_Ligand_RMSF.csv")):
        partner_rmsf = "Peptide residue RMSF" if context.partner_kind == "peptide" else "Ligand RMSF"
        add("flexibility", categorical_bar_plot(source, partner_rmsf, "Atom", ["RMSF"], "RMSF (Å)", colors=[ORANGE]), partner_rmsf, source)
    peptide_rmsf = analysis / "Peptide_RMSF.csv"
    if peptide_rmsf.is_file():
        add("flexibility", categorical_bar_plot(peptide_rmsf, f"{context.partner_label} peptide residue RMSF", "Residue", ["RMSF_A"], "RMSF (Å)", colors=[ORANGE]), "Peptide residue RMSF", peptide_rmsf,
            note="Higher values locate mobile peptide residues; compare contacting positions with the interface-node map.")

    sse_fractions = analysis / "SSE_Fractions.csv"
    if sse_fractions.is_file(): add("secondary", multiseries_plot(sse_fractions, "SSE fractions", "Time_ns", ["Helix", "Sheet", "Coil"], "Fraction"), "SSE fractions", sse_fractions)
    persistence = analysis / "Residue_SSE_Persistence.csv"
    if persistence.is_file(): add("secondary", categorical_bar_plot(persistence, "Residue SSE persistence", "Residue_Index", ["Helix_Persistence", "Sheet_Persistence", "Coil_Persistence"], "Fraction", stacked=True), "Residue SSE persistence", persistence)
    for source in sorted(analysis.glob("SSE_chain_*_fractions.csv")):
        add("secondary", multiseries_plot(source, source.stem.replace("_", " "), "Time_ns", ["Helix", "Sheet", "Coil"], "Fraction"), source.stem.replace("_", " "), source)

    frequency = analysis / "LIG_contact_frequency_filtered.csv"
    if frequency.is_file(): add("contacts", categorical_bar_plot(frequency, "Ligand-contact frequency", "Residue", ["PercentFrames"], "Frames with contact (%)"), "Ligand-contact frequency", frequency)
    filtered = analysis / "FilteredContacts_Framewise.csv"
    all_contacts = analysis / "AllContacts_Framewise.csv"
    persistence_graph = contact_persistence_plot(filtered, all_contacts, interaction_file if interaction_file.is_file() else None, protein_name) if filtered.is_file() and all_contacts.is_file() else None
    if persistence_graph: add("contacts", persistence_graph, "Ligand-contact persistence", filtered, note="Persistence is recomputed from sampled frames (0–1) using the original residue labels. A bar color now means the most frequently observed classified chemistry for that residue: H-bond blue, ionic pink, hydrophobic orange.")
    if filtered.is_file(): add("contacts", contact_events_plot(filtered, "Specific ligand contact events"), "Specific ligand contact events", filtered, note="Hover a marker for the contacting ligand atom, protein atom, distance, residue, and time.")

    summary = analysis / "InteractionTypes_Summary.csv"
    if summary.is_file(): add("chemistry", categorical_bar_plot(summary, "Interaction-type abundance by residue", "Residue", ["H-bond", "Hydrophobic", "Ionic"], "Observed interaction events", stacked=True), "Interaction-type abundance by residue", summary)
    if interaction_file.is_file():
        add("chemistry", detailed_interaction_plot(interaction_file), "Specific interaction-type events", interaction_file, note="Each residue/type pairing stays separate, so simultaneous chemistry is never collapsed into one line.")
        table = _read_table(interaction_file)
        if table is not None and "Type" in table:
            counts = table["Type"].value_counts().rename_axis("Type").reset_index(name="Events")
            figure = go.Figure(go.Bar(x=counts["Type"], y=counts["Events"], marker_color=[CYAN, MAGENTA, ORANGE][:len(counts)], hovertemplate="Interaction type=%{x}<br>Observed events=%{y}<extra></extra>"))
            figure.update_layout(title="Interaction-type distribution", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 52}, height=380, xaxis={"title": "Interaction type", "gridcolor": GRID}, yaxis={"title": "Observed events", "gridcolor": GRID})
            add("chemistry", plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False}), "Interaction-type distribution", interaction_file)

    if all_contacts.is_file():
        add("distances", distance_histogram_plot(all_contacts), "Composite ligand–residue distance distribution", all_contacts)
        add("distances", contact_events_plot(all_contacts, "All ligand–residue distance events", distance_mode=True), "All ligand–residue distance events", all_contacts, note="Hover provides the exact atom pair, residue, time, and observed separation.")
        add("distances", interaction_degree_plot(all_contacts), "Interaction degree over time", all_contacts)
    if filtered.is_file(): add("distances", contact_events_plot(filtered, "Persistent-contact distance explorer", distance_mode=True), "Persistent-contact distance explorer", filtered)
    if filtered.is_file() and frequency.is_file(): add("distances", persistence_distance_plot(filtered, frequency), "Contact persistence versus mean distance", filtered, note="Contact frequency is joined to each residue's mean observed filtered-contact distance.")

    coupling = analysis / "LIG_SSE_Contact_Coupling_Timeseries.csv"
    if coupling.is_file():
        add("coupled", multiseries_plot(coupling, "SSE versus ligand-contact dynamics", "Time_ns", ["Helix_Fraction", "Sheet_Fraction", "Coil_Fraction", "Raw_Contacts"], "Fraction / contact count"), "SSE versus ligand-contact dynamics", coupling)
        add("coupled", multiseries_plot(coupling, "Smoothed contacts", "Time_ns", ["Raw_Contacts", "Smoothed_Contacts"], "Contact count"), "Smoothed contacts", coupling)
        add("coupled", multiseries_plot(coupling, "Rolling correlation", "Time_ns", ["Rolling_Corr_Contacts_vs_Sheet"], "Correlation"), "Rolling correlation", coupling)
    changepoints = analysis / "LIG_ChangePoint_Summary.csv"
    if changepoints.is_file(): add("coupled", categorical_bar_plot(changepoints, "Change-point segments", "Segment", ["MeanContacts", "MeanHelix", "MeanSheet", "MeanCoil"], "Mean value"), "Change-point segments", changepoints, note="Hover the bars to compare the average contact and secondary-structure values for each detected trajectory segment.")
    for group, title, source, graph in interface_rin_plots(analysis):
        add(group, graph, title, source, note="This panel is isolated to the labelled chain or chain pair; it is not pooled across unrelated interfaces.")
    for index, csv_file in enumerate(csv_files):
        # Keep raw frame-wise tables and alternate selections in the download
        # bundle.  The main gallery is intentionally canonical and therefore
        # never presents several near-identical RMSD/RMSF panels as separate
        # primary conclusions.
        if csv_file not in gallery_csvs:
            continue
        graph = numeric_plot(csv_file, colours[index % len(colours)])
        if graph:
            card = dashboard_card(graph, csv_file.stem, csv_file.relative_to(output.parent).as_posix())
            stem = csv_file.stem.lower()
            if "ligand" in stem or "contact" in stem or "interaction" in stem:
                groups["contacts"].append(card)
            elif "rmsf" in stem:
                groups["flexibility"].append(card)
            elif "sse" in stem or "dssp" in stem:
                groups["secondary"].append(card)
            else:
                groups["overview"].append(card)
    # Matrix outputs are the fastest way to understand whether a ligand keeps
    # the starting pocket chemistry and whether the protein retains its fold,
    # so they deliberately lead the report instead of being buried in exports.
    groups["secondary"] = heatmap_cards[:1] + groups["secondary"]
    groups["contacts"] = heatmap_cards[1:] + groups["contacts"]
    cards = [card for values in groups.values() for card in values]
    if not cards and not has_interface_rin:
        cards.append('<article class="card"><h2>Analysis data is still being written</h2><p>No tabular analysis output was available when this report was generated. Re-run this report after 3A finishes.</p></article>')
    # Publish the chemistry drawing separately as well as embedding it in the
    # self-contained report.  The Kyle results dossier can then show the same
    # fast, recognisable ligand preview before a scientist opens the full HTML.
    ligand_structure, ligand_smiles = ligand_inputs(root, ligand_name)
    published_ligand_svg = ligand_svg(ligand_structure, ligand_smiles)
    if published_ligand_svg:
        (output.parent / "Ligand_2D.svg").write_text(published_ligand_svg, encoding="utf-8")
    overview_title = "First look: interface and fold" if context.is_interface else "First look: fold and binding"
    overview_text = (
        "Read global structure and the first persistent chain interface before opening individual chain-pair views."
        if context.is_interface else ("Read global structure and interaction chemistry before drawing conclusions about the ligand." if context.partner_kind == "ligand" else "Read protein stability, fold retention, and secondary structure together.")
    )
    overview = section_html(overview_title, "Start here", overview_text, groups["overview"])
    # Interface-only systems should never acquire a misleading empty ligand card.
    ligand_panel = ligand_html(root, ligand_name) if context.partner_kind == "ligand" else ""
    network_explorer = interactive_networx_html(analysis)
    interface_evidence = interface_evidence_html(analysis, context, interface_rin)
    structure_panel = final_structure_html(root, context)
    stability = section_html("Stability", "Trajectory stability", f"Compare global protein, complex, {context.relationship_label}, compactness, and native-contact trends.", groups["stability"])
    flexibility = section_html("Flexibility", "Residue and partner mobility", f"Inspect protein and {context.partner_label if context.partner_kind != 'none' else 'chain'} RMSF at atom and residue levels.", groups["flexibility"])
    secondary = section_html("Secondary structure", "Fold evolution", "Read the residue-time secondary-structure maps alongside global SSE fractions and persistence.", groups["secondary"])
    contacts = section_html(context.contact_tab_label, "Binding/interface evidence", f"Use {context.relationship_label} contact frequency, persistence, and event views together.", groups["contacts"])
    chemistry = section_html("Contact chemistry", f"{context.relationship_label.title()} interaction types", "Hydrogen-bond, ionic, and hydrophobic interaction events remain separate and hoverable.", groups["chemistry"])
    distances = section_html("Distances", "Partner-residue separations", f"Explore every calculated {context.relationship_label} separation without collapsing pair context.", groups["distances"]) + distance_carousel_html(analysis)
    coupled = section_html("Coupled dynamics", "Structure and contacts", f"Compare secondary structure with {context.relationship_label} contacts and detected trajectory segments.", groups["coupled"])
    dashboard = tabbed_sections([
        ("overview", "Overview", overview + interface_rin["overview"] + network_explorer),
        ("stability", "Stability", stability),
        ("flexibility", "Flexibility", flexibility),
        ("secondary", "Secondary structure", secondary),
        ("contacts", context.contact_tab_label, ligand_panel + contacts + network_explorer + (interface_rin["peptide_detail"] if context.partner_kind == "peptide" else interface_rin["protein_detail"] if context.partner_kind == "protein_chain" else "")),
        ("chemistry", "Contact chemistry", chemistry),
        ("distances", "Distances", distances),
        ("coupled", "Coupled dynamics", coupled),
        ("interface-evidence", "Interface maps", interface_evidence),
        ("interfaces", "Protein interfaces", "" if context.partner_kind == "protein_chain" else interface_rin["protein_detail"]),
        ("peptide-interfaces", "Peptide interfaces", "" if context.partner_kind == "peptide" else interface_rin["peptide_detail"]),
        ("structure", "Structure", structure_panel),
    ])
    html = build_document(system_name, context, cards, dashboard)
    output.write_text(html, encoding="utf-8")
    manifest = build_results_manifest(root, analysis, output, system_name, ligand_name)
    manifest["partner_context"] = {
        "partner_kind": context.partner_kind,
        "partner_label": context.partner_label,
        "partner_chain_ids": list(context.partner_chain_ids),
        "interface_pairs": [list(pair) for pair in context.interface_pairs],
        "system_class": context.system_class,
        "reference_contact_definition": "first analyzed-frame reference contacts when a Q metric is available",
    }
    (analysis / RESULTS_MANIFEST_NAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    # Retain a stable compatibility entry point for existing external links.
    if heatmap_cards:
        (output.parent / "PyMACS_Heatmap_Explorer.html").write_text(html, encoding="utf-8")
    print(f"Wrote {output}")
    return

    # Historical report body retained below only for source compatibility with
    # older generated artifacts.  The guided gallery above is the active view.
    html = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{escape(system_name)} · PyMACS results</title><script>{plotlyjs()}</script><style>
*{{box-sizing:border-box}} body{{margin:0;background:radial-gradient(circle at 90% 0,rgba(53,208,255,.15),transparent 30%),radial-gradient(circle at 4% 25%,rgba(255,43,214,.10),transparent 25%),{BG};color:{TEXT};font-family:Inter,ui-sans-serif,system-ui,-apple-system,sans-serif}}.page{{max-width:1500px;margin:auto;padding:34px 24px 64px}}.hero,.card{{background:linear-gradient(155deg,#10172b,#060912);border:1px solid #23375d;box-shadow:0 18px 50px rgba(0,0,0,.35)}}.hero{{border-radius:24px;padding:30px;margin-bottom:20px}}.eyebrow{{color:{GREEN};font-size:.75rem;font-weight:900;letter-spacing:.16em}}h1{{font-size:clamp(2.3rem,5vw,4.9rem);letter-spacing:-.06em;line-height:.95;margin:.25rem 0 1rem}}.sub,.meta{{color:{MUTED};line-height:1.55}}.pills{{display:flex;gap:10px;flex-wrap:wrap;margin-top:18px}}.pill{{border:1px solid rgba(53,208,255,.3);padding:8px 11px;border-radius:999px;background:rgba(53,208,255,.08)}}.heatmap-link{{display:block;margin-top:20px;border:1px solid rgba(24,255,139,.72);border-radius:15px;padding:16px 18px;background:linear-gradient(100deg,rgba(24,255,139,.18),rgba(53,208,255,.10));color:{TEXT};text-decoration:none;box-shadow:0 0 30px rgba(24,255,139,.08)}}.heatmap-link strong{{display:block;color:{GREEN};font-size:1.06rem}}.heatmap-link span{{display:block;margin-top:5px;color:{MUTED};font-size:.92rem}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(390px,1fr));gap:16px}}.card{{border-radius:17px;padding:13px;overflow:hidden}}.card .js-plotly-plot{{cursor:zoom-in}}.card a{{display:inline-block;color:{CYAN};font-size:.8rem;font-weight:800;margin:0 10px 8px}}.heatmap-note{{margin:2px 10px 8px;color:{MUTED};font-size:.82rem;line-height:1.4}}.graph-modal{{display:none;position:fixed;z-index:9999;inset:0;padding:3vh 3vw;background:rgba(0,0,0,.78);backdrop-filter:blur(7px)}}.graph-modal.open{{display:flex;align-items:center;justify-content:center}}.graph-modal-panel{{position:relative;width:min(1500px,94vw);height:min(900px,92vh);padding:20px 16px 12px;border:1px solid #36d0ff;background:#080d1b;box-shadow:0 0 70px rgba(53,208,255,.18);border-radius:20px}}.graph-modal-close{{position:absolute;right:17px;top:13px;z-index:2;border:1px solid #35d0ff;border-radius:999px;background:#101b33;color:{TEXT};font-weight:800;padding:8px 12px;cursor:pointer}}.graph-modal-plot{{width:100%;height:100%}}@media(max-width:600px){{.page{{padding:20px 12px 42px}}.grid{{grid-template-columns:1fr}}.graph-modal{{padding:1vh 1vw}}.graph-modal-panel{{width:98vw;height:96vh;padding:14px 4px 4px}}}}
</style></head><body><main class="page"><section class="hero"><div class="eyebrow">PYMACS · SINGLE-SYSTEM INTERACTIVE RESULTS</div><h1>{escape(system_name)}</h1><p class="sub">Interactive analysis plots generated automatically from this run’s local PyMACS CSV outputs. Double-click any graph for a larger reading view; use the comparative dashboard when comparing related simulations.</p><div class="pills"><span class="pill">Ligand / PROTAC: {escape(ligand_name)}</span><span class="pill">{len(cards)} interactive plots</span><span class="pill">Analysis_Results</span></div>{heatmap_link}</section><section class="grid">{''.join(cards)}</section></main><div class="graph-modal" id="graph-modal" aria-hidden="true"><div class="graph-modal-panel"><button class="graph-modal-close" id="graph-modal-close" type="button">Close ×</button><div class="graph-modal-plot" id="graph-modal-plot"></div></div></div><script>
(() => {{
  const modal = document.getElementById('graph-modal');
  const destination = document.getElementById('graph-modal-plot');
  const close = () => {{ modal.classList.remove('open'); modal.setAttribute('aria-hidden', 'true'); Plotly.purge(destination); }};
  const open = (source) => {{
    const layout = Object.assign({{}}, source.layout, {{autosize: true, height: null, margin: {{l: 70, r: 35, t: 65, b: 60}}}});
    Plotly.react(destination, source.data, layout, {{responsive: true, displaylogo: false}});
    modal.classList.add('open'); modal.setAttribute('aria-hidden', 'false');
    requestAnimationFrame(() => Plotly.Plots.resize(destination));
  }};
  document.querySelectorAll('.js-plotly-plot').forEach((source) => source.addEventListener('dblclick', (event) => {{ event.preventDefault(); event.stopImmediatePropagation(); open(source); }}, true));
  document.getElementById('graph-modal-close').addEventListener('click', close);
  modal.addEventListener('click', (event) => {{ if (event.target === modal) close(); }});
  document.addEventListener('keydown', (event) => {{ if (event.key === 'Escape') close(); }});
}})();
</script></body></html>'''
    output.write_text(html, encoding="utf-8")
    if heatmap_cards:
        main_grid = f"<section class=\"grid\">{''.join(cards)}</section>"
        heatmap_grid = f"<section class=\"grid\">{''.join(heatmap_cards)}</section>"
        explorer = html.replace(main_grid, heatmap_grid)
        explorer = explorer.replace(heatmap_link, '<a class="heatmap-link" href="PyMACS_Result_Dashboard.html">← Back to analysis plots</a>')
        explorer = explorer.replace(f'{len(cards)} interactive plots', f'{len(heatmap_cards)} interactive heatmaps')
        (output.parent / "PyMACS_Heatmap_Explorer.html").write_text(explorer, encoding="utf-8")
    print(f"Wrote {output}")


def plotlyjs() -> str:
    # Keep the deliverable self-contained for Kyle/offline viewing.
    from plotly.offline import get_plotlyjs
    return get_plotlyjs()


if __name__ == "__main__":
    main()
