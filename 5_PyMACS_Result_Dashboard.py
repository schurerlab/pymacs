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
import shutil
import subprocess
import textwrap
import zipfile
from dataclasses import dataclass
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot
from plotly.subplots import make_subplots


BG, PANEL, GRID, TEXT, MUTED = "#05060a", "#0b1020", "#1c2a4a", "#cfe8ff", "#7aa6d9"
CYAN, GREEN, MAGENTA, ORANGE = "#35d0ff", "#18ff8b", "#ff2bd6", "#ff8a3d"
# Keep chain identities visually stable across every dossier.  The first two
# colors are the requested Miami orange and neon Miami green; later chains
# remain deliberately high-contrast on the dark report background.
CHAIN_COLORS = (ORANGE, GREEN, MAGENTA, CYAN, "#ffe66d", "#a78bfa")
# Plotly positions a horizontal legend in paper coordinates.  Keeping it below
# the plotting area is the only placement that remains reliable when a card is
# narrow, a title is long, or legend labels wrap onto a second row.
LEGEND_BELOW = {"orientation": "h", "x": 0, "xanchor": "left", "y": -0.22, "yanchor": "top"}
# A browser report should open immediately even when its source CSV contains a
# full 500 ns trajectory at sub-ps sampling.  The CSV download remains exact;
# this only limits the points drawn in the browser.
MAX_RENDER_POINTS = 2000
MAX_HEATMAP_TIMEPOINTS = 250
MAX_EVENT_POINTS = 12000
MAX_TABLE_ROWS = 50000
TRAJECTORY_PREVIEW_FRAMES = 151

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


def component_key_html(root: Path) -> str:
    """Compact dossier-wide key for multi-component color semantics."""
    entries = atomindex_components(root)
    if len(entries) < 2:
        return ""
    chips = []
    for index, entry in enumerate(entries):
        color = CHAIN_COLORS[index % len(CHAIN_COLORS)]
        chain = entry.get("current_chain") or entry.get("source_chain") or "?"
        chips.append(
            f'<span class="component-key-chip" style="--component-color:{color}"><i></i>'
            f'Component {index + 1} · {escape(entry["label"])} <small>chain {escape(chain)}</small></span>'
        )
    return ('<style>.component-key{display:flex;gap:9px;flex-wrap:wrap;margin:0 0 15px}.component-key-chip{display:inline-flex;align-items:center;gap:7px;padding:7px 10px;border:1px solid color-mix(in srgb,var(--component-color) 52%,#285a79);border-radius:999px;background:#091423;color:#dff5ff;font-size:.79rem;font-weight:850}.component-key-chip i{width:9px;height:9px;border-radius:50%;background:var(--component-color);box-shadow:0 0 10px var(--component-color)}.component-key-chip small{color:#8daac7;font-weight:700}</style><aside class="component-key" aria-label="Component color key">' + "".join(chips) + '</aside>')


def legacy_receptor_label(root: Path) -> str:
    """Recover a receptor name for a legacy ligand run without atomIndex.

    Some imported runs retained a trustworthy historical-source path but lost
    their atomIndex and chain map.  The immediate parent of the system in that
    path is normally the receptor/PDB directory (for example ``LEDGF_9g94``).
    Use only that narrowly structured name; generic cache folders such as
    ``historical`` and ``Rank3`` must never become a scientific label.
    """
    candidates: list[str] = []
    try:
        run = json.loads((root / "pymacs_run.json").read_text(encoding="utf-8"))
        historical_source = str(run.get("historical_source") or "").strip()
        if historical_source:
            candidates.append(Path(historical_source).parent.name)
    except (OSError, ValueError, TypeError):
        pass
    candidates.extend(parent.name for parent in list(root.parents)[:6])
    for candidate in candidates:
        match = re.fullmatch(r"([A-Za-z][A-Za-z0-9-]{1,30})_[0-9][A-Za-z0-9]{3}", candidate)
        if match:
            return f"{match.group(1)} (protein)"
    # A final coordinate PDB still retains the real receptor chain even when
    # its project-level name was lost during an older import.  “Receptor
    # (chain A)” is deliberately an honest label: it is more useful than a
    # bare “Protein”, without inventing a target name that the source cannot
    # support.
    for pdb_name in ("Final_Trajectory.pdb", "binding_pocket_only.pdb"):
        try:
            with (root / pdb_name).open(encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith(("ATOM  ", "HETATM")):
                        chain = line[21:22].strip()
                        if chain:
                            return f"Receptor (chain {chain})"
        except OSError:
            continue
    return ""


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


def ligand_evidence_present(analysis: Path) -> bool:
    """Return whether 3B/3C produced a small-molecule evidence bundle.

    This is deliberately based on analysed outputs, rather than the historical
    submission label.  Older published runs occasionally lack ligand metadata;
    legacy runs were sometimes saved with ``mode=protein``
    or without the original ligand metadata even though their ligand analyses
    completed successfully.
    """
    return any((
        (analysis / "InteractionTypes_Framewise.csv").is_file(),
        (analysis / "FilteredContacts_Framewise.csv").is_file(),
        (analysis / "AllContacts_Framewise.csv").is_file(),
        (analysis / "LIG_Complex_RMSD_Overlay.csv").is_file(),
        (analysis / "Radius_of_Gyration_Protein_Ligand_Complex.csv").is_file(),
        any(analysis.glob("*_Ligand_RMSD.csv")),
        any(analysis.glob("*_contact_frequency_filtered.csv")),
    ))


def ligand_label_from_evidence(root: Path, setup: dict, metadata: dict) -> str:
    """Choose the best honest display name when legacy metadata is incomplete."""
    placeholders = {"", "none", "no ligand label", "proteinonly"}
    for value in (
        metadata.get("ligand_name"), metadata.get("ligand_residue"),
        setup.get("ligand"), setup.get("compound_name"),
    ):
        label = str(value or "").strip()
        if label.casefold() not in placeholders:
            return label
    # A parameterized CGenFF file is the next-best authoritative identity.
    structure, _ = ligand_inputs(root, "")
    if structure:
        stem = structure.name.split(".", 1)[0].strip()
        if stem:
            return stem
    # Do not falsely call a small molecule a protein merely because its old
    # record lost a name.  "LIG" is the transparent generic display fallback.
    return "LIG"


def promote_ligand_context_when_evidenced(
    root: Path, analysis: Path, setup: dict, metadata: dict, context: PartnerContext,
) -> tuple[PartnerContext, bool]:
    """Prevent stale protein-only metadata from suppressing ligand reports.

    Peptide classifications stay authoritative: a peptide system is not
    converted merely because it contains a contact table.  For every other
    context, the complete small-molecule analysis bundle wins over an old
    protein/apo label.
    """
    evidence = ligand_evidence_present(analysis)
    if evidence and context.partner_kind not in {"ligand", "peptide"}:
        return (
            PartnerContext(
                "ligand", ligand_label_from_evidence(root, setup, metadata),
                context.partner_chain_ids, context.interface_pairs, "ligand_complex",
            ),
            True,
        )
    return context, evidence


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
    peptide_explorer = analysis / "Peptide_NETWORX" / "Peptide_NETWORX_Explorer.html"
    if peptide_explorer.is_file():
        canonical.append(peptide_explorer.relative_to(root).as_posix())
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
        "peptide_networx": next((path for path in canonical if path.endswith("Peptide_NETWORX_Explorer.html")), None),
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
    "distance_trace": (
        "Closest peptide–protein separation over time",
        "At every sampled frame this trace keeps only the shortest measured gap between the peptide and protein. It is a nearest-neighbour readout, not an average of all distances.",
        "A mostly steady low band means the interface stays close. Brief tall spikes mean the closest measured gap temporarily increased; compare them with the contact timeline before calling a dissociation.",
    ),
    "distance_distribution": (
        "Residue separation distributions",
        "Each violin collects the measured gaps for one protein residue whenever that residue was recorded as contacting the peptide. Wider sections mean that distance occurred more often among those recorded events.",
        "Look first at the middle line and box: smaller typical gaps and a compact shape suggest a consistently close partner. These are observed contact events only, not every trajectory frame.",
    ),
    "distance_ranking": (
        "Persistent-contact distance ranking",
        "This ranks protein residues using their typical measured gap, prioritizing residues that recur in the contact evidence. It condenses many individual contact events into an easy comparison.",
        "Shorter bars mean a smaller typical observed separation. Use the hover values and contact-persistence panels to distinguish a close rare event from a close recurring partner.",
    ),
    "contact_episode": (
        "Uninterrupted observed contact episodes",
        "This measures the longest run of sampled frames in which each protein residue was recorded in contact with the peptide.",
        "Longer bars mean a contact remained present across more consecutive sampled frames. It is a sampled-trajectory continuity measure, not a chemical residence time or binding lifetime.",
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


def guide_html(key: str, prompt: str = "ⓘ ELI5: what am I looking at?") -> str:
    title, meaning, expected = GRAPH_GUIDES.get(key, GRAPH_GUIDES["generic"])
    return (
        f'<details class="guide"><summary>{escape(prompt)}</summary>'
        f'<p><b>What this is:</b> {escape(title)}. {escape(meaning)}</p>'
        f'<p><b>How to read it:</b> {escape(expected)}</p>'
        '<p><b>Keep in mind:</b> One graph is a clue, not a verdict. Check the matching stability, contact, and structure panels before making a biological claim.</p></details>'
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
    """Build a heavy-atom display MOL2 with final MD coordinates and CGenFF bonds.

    GROMACS final PDB snapshots generally omit ``CONECT`` records.  Browsers
    therefore guess small-molecule bonds from distances, which loses aromatic
    bond orders and can make a ligand look chemically wrong.  CGenFF MOL2
    files contain the authoritative atom names and bonds, but their coordinates
    are pre-MD.  This joins those two sources without altering either input.
    The emitted model deliberately excludes explicit hydrogens: they clutter a
    complex viewer and are not useful for the final-structure overview.  Bond
    records (including ``ar`` aromatic bonds) are retained between heavy atoms.

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
        # A CGenFF MOL2 has explicit hydrogens.  Keep the PDB matching strict
        # (so we only accept a fully identified ligand), then make a separate
        # heavy-atom model for 3D display and renumber its bonds.
        display_atoms: list[tuple[int, list[str]]] = []
        old_to_new: dict[int, int] = {}
        for _line_index, fields in atom_records:
            atom_type = fields[5].split(".", 1)[0].upper()
            atom_name = fields[1].upper()
            if atom_type == "H" or atom_name.startswith("H"):
                continue
            old_index = int(fields[0])
            old_to_new[old_index] = len(display_atoms) + 1
            display_atoms.append((old_index, fields))
        if not display_atoms:
            return None, ""

        bond_end = next(
            (index for index in range(bond_start + 1, len(mol2_lines))
             if mol2_lines[index].strip().startswith("@<TRIPOS>")),
            len(mol2_lines),
        )
        display_bonds: list[tuple[int, int, str]] = []
        for index in range(bond_start + 1, bond_end):
            fields = mol2_lines[index].split()
            if not fields:
                continue
            if len(fields) < 4:
                return None, ""
            try:
                source, target = int(fields[1]), int(fields[2])
            except ValueError:
                return None, ""
            if source in old_to_new and target in old_to_new:
                display_bonds.append((old_to_new[source], old_to_new[target], fields[3]))

        molecule_marker = next(
            (index for index, line in enumerate(mol2_lines) if line.strip() == "@<TRIPOS>MOLECULE"),
            None,
        )
        molecule_name = cgenff_mol2.stem
        if molecule_marker is not None:
            for line in mol2_lines[molecule_marker + 1:atom_start - 1]:
                if line.strip():
                    molecule_name = line.strip()
                    break
        rebuilt = [
            "@<TRIPOS>MOLECULE", molecule_name,
            f"{len(display_atoms):5d} {len(display_bonds):5d} 0 0 0",
            "SMALL", "GASTEIGER", "", "@<TRIPOS>ATOM",
        ]
        for _old_index, fields in display_atoms:
            x, y, z = coordinates[fields[1].upper()]
            rebuilt.append("{:<7} {:<8} {:>10.4f} {:>10.4f} {:>10.4f} {}".format(
                old_to_new[int(fields[0])], fields[1], x, y, z, " ".join(fields[5:]),
            ))
        rebuilt.append("@<TRIPOS>BOND")
        rebuilt.extend(
            f"{index:>6} {source:>5} {target:>5} {bond_type}"
            for index, (source, target, bond_type) in enumerate(display_bonds, start=1)
        )
        rebuilt.extend([
            "@<TRIPOS>SUBSTRUCTURE",
            f"     1 {resname:<8} 1 GROUP             0 ****  ****    0 ROOT",
        ])
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


def interactive_networx_html(analysis: Path, *, ligand_expected: bool = False) -> str:
    """Embed the explorer source directly; never navigate to an attachment.

    A ligand dossier must never silently look complete when its atom-resolved
    explorer artifact has not yet been produced.  The status panel keeps that
    distinction visible during recovery/backfill without fabricating a graph.
    """
    explorers = sorted((analysis / "NETWORX").glob("*_interactive.html"))
    if not explorers:
        if not ligand_expected:
            return ""
        events = analysis / "InteractionTypes_Framewise.csv"
        if events.is_file() and events.stat().st_size:
            detail = (
                "Atom-level interaction events are available, but the interactive "
                "NETWORX artifact has not been generated or published yet. Run the "
                "analysis-only 3B NETWORX stage, then regenerate this dashboard."
            )
        else:
            detail = (
                "This run has no atom-level interaction-event table. Regenerate the "
                "completed-trajectory interaction analysis (3A) before running the "
                "analysis-only 3B NETWORX stage."
            )
        return f'''<section class="gallery-section explorer-pending"><div class="section-head"><div><div class="eyebrow">Persistent interaction explorer</div><h2>Interaction explorer pending</h2><p>{detail}</p></div></div></section>'''
    explorer = explorers[0]
    try:
        document = explorer.read_text(encoding="utf-8")
    except OSError:
        return ""
    # srcdoc makes the explorer part of this report.  Its CSS/IDs/scripts are
    # deliberately isolated, so the Overview and Contacts copies can coexist.
    return f'''<section class="gallery-section"><div class="section-head"><div><div class="eyebrow">Persistent interaction explorer</div><h2>Compare retained interactions</h2><p>Choose the displayed occupancy threshold to inspect the chemically classified residue–ligand-atom interactions. The larger embedded canvas supports node rearrangement and SVG export without opening a separate page.</p></div></div><iframe title="Interactive ligand interaction explorer" srcdoc="{escape(document, quote=True)}" sandbox="allow-scripts allow-downloads" style="display:block;width:100%;height:1320px;border:1px solid #285a79;border-radius:16px;background:#05070d" loading="lazy"></iframe></section>'''


def peptide_networx_html(analysis: Path, *, peptide_expected: bool = False) -> str:
    """Embed the dedicated sequence-preserving peptide topology when published."""
    explorer = analysis / "Peptide_NETWORX" / "Peptide_NETWORX_Explorer.html"
    if not explorer.is_file():
        if not peptide_expected:
            return ""
        return '''<section class="gallery-section explorer-pending"><div class="section-head"><div><div class="eyebrow">Peptide topology</div><h2>Peptide NETWORX pending</h2><p>The peptide interface data are available, but the sequence-preserving explorer has not been published. Re-run the 3C interface stage; it now writes this presentation-only view from the completed interface outputs.</p></div></div></section>'''
    try:
        document = explorer.read_text(encoding="utf-8")
    except OSError:
        return ""
    return f'''<section class="gallery-section peptide-networx"><div class="section-head"><div><div class="eyebrow">Peptide interaction topology</div><h2>Full sequence, persistent receptor contacts</h2><p>The peptide backbone is fixed in N → C order; only contact edges and highlighting respond to the occupancy threshold. This schematic preserves sequence context and does not represent a 3D conformation.</p></div></div><iframe title="Peptide NETWORX interaction explorer" srcdoc="{escape(document, quote=True)}" sandbox="allow-scripts allow-downloads" style="display:block;width:100%;height:1180px;border:1px solid #285a79;border-radius:16px;background:#05070d" loading="lazy"></iframe></section>'''


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
        # The multiview explorer owns its state in a select control rather
        # than a global ``pair`` variable.  Set it after its inline scripts
        # have initialized so each dashboard panel opens its own chain pair.
        focused = focused.replace(
            "</html>",
            f"<script>addEventListener('DOMContentLoaded',()=>{{const p=document.querySelector('#pair');if(p){{p.value={json.dumps(pair_id)};p.dispatchEvent(new Event('change'));}}}});</script></html>",
            1,
        )
        if compact:
            focused = focused.replace(".controls{display:flex", "#pairs{display:none}.controls{display:flex", 1)
        return (
            f'<iframe title="{escape(title)}" srcdoc="{escape(focused, quote=True)}" '
            'sandbox="allow-scripts allow-same-origin allow-downloads" style="display:block;width:100%;height:1420px;border:1px solid #285a79;'
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


def interface_svg_snapshots(explorer: Path, threshold: float = .5) -> list[tuple[str, str]]:
    """Create thresholded, chain-coloured SVG previews from explorer data."""
    try:
        text = explorer.read_text(encoding="utf-8")
        match = re.search(r"<script>const DATA=(.*?);\n\(\(\)=>", text, re.DOTALL)
        data = json.loads(match.group(1)) if match else {}
        pair = (data.get("pairs") or [])[0]
        left, right = pair["left"], pair["right"]
        edges = [edge for edge in data.get("edges", []) if edge.get("fraction", 0) >= threshold and
                 {edge.get("a", "").split(":", 1)[0], edge.get("b", "").split(":", 1)[0]} == {left, right}]
        # Component identity comes from the chain prefix, rather than the source
        # table's arbitrary edge direction.
        left_nodes = list(dict.fromkeys(node for edge in edges for node in (edge["a"], edge["b"])
                                        if node.split(":", 1)[0] == left))
        right_nodes = list(dict.fromkeys(node for edge in edges for node in (edge["a"], edge["b"])
                                         if node.split(":", 1)[0] == right))
        if not edges or not left_nodes or not right_nodes:
            return []
        def short(label: str) -> str:
            return label.split(":")[-1]
        def component_label(label: str) -> str:
            return label.split(" ", 1)[0]
        def edge_tooltip(edge: dict) -> str:
            evidence = edge.get("evidence") or []
            distances = [float(row["distance_a"]) for row in evidence if row.get("distance_a") is not None]
            distance = f"\nClosest reported distance: {min(distances):.2f} Å" if distances else ""
            return (f"{short(edge['a'])} ↔ {short(edge['b'])}\n"
                    f"Interaction: {edge.get('type', 'contact')}\n"
                    f"Persistence: {edge.get('fraction', 0):.1%}{distance}\n"
                    f"Atomic interaction observations: {len(evidence)}")
        node_edges = {node: [edge for edge in edges if node in (edge["a"], edge["b"])]
                      for node in left_nodes + right_nodes}
        def node_tooltip(node: str) -> str:
            contacts = node_edges[node]
            strongest = max((float(edge.get("fraction", 0)) for edge in contacts), default=0)
            component = "Component 1 / chain A" if node in left_nodes else "Component 2 / chain B"
            partners: dict[str, dict] = {}
            for edge in contacts:
                partner = edge["b"] if edge["a"] == node else edge["a"]
                entry = partners.setdefault(partner, {"fraction": 0.0, "types": set(), "distances": []})
                entry["fraction"] = max(entry["fraction"], float(edge.get("fraction", 0)))
                entry["types"].add(str(edge.get("type", "contact")))
                entry["distances"].extend(float(row["distance_a"]) for row in (edge.get("evidence") or []) if row.get("distance_a") is not None)
            ordered = sorted(partners.items(), key=lambda item: item[1]["fraction"], reverse=True)
            partner_lines = []
            for partner, entry in ordered[:7]:
                geometry = f" · {min(entry['distances']):.2f} Å" if entry["distances"] else ""
                partner_lines.append(f"  {short(partner)} · {'/'.join(sorted(entry['types']))} · {entry['fraction']:.1%}{geometry}")
            if len(ordered) > 7:
                partner_lines.append(f"  +{len(ordered) - 7} additional persistent partners")
            return (f"{short(node)}\n{component}\n"
                    f"Persistent residue partners: {len(partners)} · strongest {strongest:.1%}\n"
                    "PARTNERS · interaction · persistence · closest distance\n" + "\n".join(partner_lines))
        def y(index: int, count: int) -> int:
            return 128 + round(index * 220 / max(count - 1, 1))
        all_nodes = left_nodes + right_nodes
        def ring(nodes: list[str], cx: int, cy: int, rx: int, ry: int) -> dict[str, tuple[int, int]]:
            return {node: (round(cx + rx * np.cos(-np.pi / 2 + 2 * np.pi * index / max(len(nodes), 1))), round(cy + ry * np.sin(-np.pi / 2 + 2 * np.pi * index / max(len(nodes), 1)))) for index, node in enumerate(nodes)}
        layouts = {
            "Force-directed": ring(all_nodes, 280, 242, 205, 118),
            "Bipartite": {**{node: (92, y(index, len(left_nodes))) for index, node in enumerate(left_nodes)}, **{node: (468, y(index, len(right_nodes))) for index, node in enumerate(right_nodes)}},
            "Concentric": {**ring(left_nodes, 280, 242, 205, 118), **ring(right_nodes, 280, 242, 112, 68)},
            "Circular": ring(all_nodes, 280, 242, 205, 118),
            "Chain-separated": {**{node: (96 + index * 368 / max(len(left_nodes) - 1, 1), 145) for index, node in enumerate(left_nodes)}, **{node: (96 + index * 368 / max(len(right_nodes) - 1, 1), 332) for index, node in enumerate(right_nodes)}},
        }
        snapshots = []
        for layout, positions in layouts.items():
            lines = ''.join(f'<path d="M{x1} {y1} L{x2} {y2}" stroke="#35d0ff" stroke-opacity="{.3 + .7 * edge["fraction"]:.2f}" stroke-width="{1 + 5 * edge["fraction"]:.1f}" stroke-linecap="round" style="cursor:help" data-interface-tooltip-kind="edge" data-interface-tooltip="{escape(edge_tooltip(edge), quote=True)}" aria-label="{escape(edge_tooltip(edge), quote=True)}"/>' for edge in edges for (x1, y1), (x2, y2) in [(positions[edge["a"]], positions[edge["b"]])])
            nodes = ''.join(f'<g style="cursor:help" tabindex="0" data-interface-tooltip-kind="node" data-interface-component="{"one" if node in left_nodes else "two"}" data-interface-tooltip="{escape(node_tooltip(node), quote=True)}" aria-label="{escape(node_tooltip(node), quote=True)}"><circle cx="{x}" cy="{yy}" r="15" fill="{ORANGE if node in left_nodes else GREEN}"/><text x="{x}" y="{yy + 3}" text-anchor="middle" font-size="7" font-weight="800" pointer-events="none">{escape(short(node))}</text></g>' for node, (x, yy) in positions.items())
            svg = f'<svg viewBox="0 0 560 400" role="img" aria-label="{escape(layout)} interface view"><rect data-interface-background="true" width="560" height="400" fill="#07111e"/><path d="M18 101H542" stroke="#35d0ff" stroke-opacity=".28"/><path d="M18 104H188" stroke="{ORANGE}" stroke-opacity=".75"/><path d="M372 104H542" stroke="{GREEN}" stroke-opacity=".75"/><text data-interface-export-title="true" x="18" y="30" fill="#d9ecff" font-size="13" font-weight="800">{escape(layout)} · {escape(left)} ↔ {escape(right)} · ≥{threshold:.0%}</text><circle cx="24" cy="58" r="6" fill="{ORANGE}"/><text x="38" y="62" fill="{ORANGE}" font-size="10" font-weight="800">COMPONENT 1 · {escape(component_label(left))}</text><circle cx="222" cy="58" r="6" fill="{GREEN}"/><text x="236" y="62" fill="{GREEN}" font-size="10" font-weight="800">COMPONENT 2 · {escape(component_label(right))}</text><text data-interface-export-remove="true" x="18" y="84" fill="#8aa9bd" font-size="8">HOVER HUD · residue: partners + strongest persistence · cyan edge: contact chemistry + geometry</text>{lines}{nodes}</svg>'
            snapshots.append((svg, f"{layout} · ≥{threshold:.0%} persistence"))
        return snapshots
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return []


def interface_evidence_html(analysis: Path, context: PartnerContext, rin: dict) -> str:
    """Make interface-network availability obvious; never silently hide a peptide relationship."""
    if not context.is_interface:
        return ""
    interface_dir = analysis / "Interface_RIN"
    explorer = interface_dir / "Interface_RIN_Explorer.html"
    if explorer.is_file():
        threshold_sets: list[str] = []
        available_thresholds: list[str] = []
        for threshold in (.5, .75, .9):
            previews = [f'<figure class="interface-preview">{svg}<figcaption>{escape(label)}</figcaption></figure>'
                        for svg, label in interface_svg_snapshots(explorer, threshold)]
            percent = f"{threshold:.0%}"
            if previews:
                slides = ''.join(f'<article class="interface-slide"{" hidden" if index else ""}>{preview}</article>' for index, preview in enumerate(previews))
                carousel = (
                    f'<div class="interface-carousel"><div class="interface-carousel-track">{slides}</div>'
                    f'<div class="interface-carousel-controls"><button type="button" data-interface-carousel="previous">← Previous</button>'
                    f'<span data-interface-carousel-label>1 / {len(previews)}</span>'
                    '<button type="button" data-interface-carousel="next">Next →</button>'
                    '<button type="button" data-interface-svg-download>Download SVG (transparent)</button></div></div>'
                )
            else:
                carousel = f'<div class="empty-analysis"><p>No residue-pair contact reaches the selected ≥{percent} persistence threshold.</p></div>'
            threshold_sets.append(
                f'<div class="interface-threshold-set" data-interface-threshold-set="{percent}"'
                f'{" hidden" if threshold != .5 else ""}>{carousel}</div>'
            )
            available_thresholds.append(percent)
        controls = ''.join(
            f'<button type="button" data-interface-threshold="{percent}" aria-pressed="{str(percent == "50%").lower()}">≥{percent}</button>'
            for percent in available_thresholds
        )
        static_maps = (f'<div class="interface-threshold-controls" role="group" aria-label="Minimum contact persistence">{controls}</div>'
                       f'{"".join(threshold_sets)}') if threshold_sets else ""
        return (
            '<section class="gallery-section"><div class="section-head"><div><div class="eyebrow">Interface evidence</div>'
            f'<h2>{escape(context.relationship_label.title())} residue maps</h2>'
            '<p>Choose a minimum contact-persistence threshold: 50% shows the broad stable interface, 75% emphasizes stronger contacts, and 90% isolates its core. Component 1 is Miami orange, component 2 is neon green, and cyan edges are shared contacts.</p></div></div>'
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


def representative_snapshots_html(root: Path, analysis: Path) -> str:
    """Publish early/middle/late coordinate snapshots when the final XTC exists.

    This is deliberately optional: it never alters a trajectory or blocks the
    dossier if a historical run retained only its final PDB.
    """
    trajectory, topology = root / "Final_Trajectory.xtc", root / "Final_Trajectory.pdb"
    target = analysis / "Representative_Snapshots"
    names = (("Early", "early"), ("Middle", "middle"), ("Late", "late"))
    if trajectory.is_file() and topology.is_file() and not all((target / f"{slug}.pdb").is_file() for _, slug in names):
        try:
            import MDAnalysis as mda
            universe = mda.Universe(str(topology), str(trajectory))
            count = len(universe.trajectory)
            if count:
                target.mkdir(parents=True, exist_ok=True)
                for (label, slug), frame_index in zip(names, (0, count // 2, count - 1)):
                    universe.trajectory[frame_index]
                    universe.atoms.write(str(target / f"{slug}.pdb"))
        except Exception:
            pass
    available = [(label, target / f"{slug}.pdb") for label, slug in names if (target / f"{slug}.pdb").is_file()]
    if not available:
        return ""
    cards = "".join(
        f'<a class="snapshot-link" href="{escape(path.relative_to(analysis).as_posix())}" download><b>{escape(label)} trajectory frame</b><span>Download PDB snapshot</span></a>'
        for label, path in available
    )
    return f'''<style>.snapshot-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(155px,1fr));gap:10px;margin-top:14px}}.snapshot-link{{display:block;text-decoration:none;border:1px solid #27476a;border-radius:12px;padding:12px;background:#0b1524;color:#cfe8ff}}.snapshot-link:hover{{border-color:#35d0ff;background:#102945}}.snapshot-link b,.snapshot-link span{{display:block}}.snapshot-link span{{color:#7aa6d9;font-size:.82rem;margin-top:4px}}</style><section class="snapshot-grid" aria-label="Representative trajectory snapshots">{cards}</section>'''


def trajectory_playback_html(root: Path, analysis: Path, context: PartnerContext) -> str:
    """Publish a lazy browser-preview trajectory plus its original local assets.

    The full XTC is retained as a downloadable research asset.  Browsers load
    only a centered, solvent-free, evenly sampled multi-model PDB after an
    explicit click, so opening an ordinary dossier never downloads trajectory
    coordinates or allocates WebGL memory for them.
    """
    topology, trajectory = root / "Final_Trajectory.pdb", root / "Final_Trajectory.xtc"
    if not (topology.is_file() and trajectory.is_file()):
        return (
            '<section class="trajectory-card"><div class="eyebrow">Trajectory playback</div>'
            '<h2>Interactive simulation playback is not available yet</h2>'
            '<p>The final topology PDB and XTC must both be published for this run. The ordinary Structure tab remains available.</p></section>'
        )
    asset_dir = analysis / "Trajectory_Playback"
    preview = asset_dir / "centered_trajectory_preview.pdb"
    # A PBC-repaired preview, when made during analysis, is the only
    # trajectory that should be opened by default in an interactive viewer.
    # The raw XTC remains downloadable but can show a molecule on opposite
    # sides of the periodic box (a display artifact, not a dissociation).
    pbc_preview = asset_dir / "PyMACS_centered_pbc_preview.pdb"
    metadata = asset_dir / "preview_metadata.json"
    source_files = ("Final_Trajectory.pdb", "Final_Trajectory.xtc", "binding_pocket_only.pdb", "binding_pocket_only.xtc")
    asset_dir.mkdir(parents=True, exist_ok=True)
    # Keep the exact source files close to the published report.  A hard link
    # avoids duplicate storage on Triton; copy2 is the portable fallback.
    for name in source_files:
        source, target = root / name, asset_dir / name
        if not source.is_file():
            continue
        try:
            if target.exists() and target.stat().st_size == source.stat().st_size:
                continue
            target.unlink(missing_ok=True)
            try:
                os.link(source, target)
            except OSError:
                shutil.copy2(source, target)
        except OSError:
            pass

    visualization_source = pbc_preview if pbc_preview.is_file() else preview
    # A .pse serialises every loaded XTC frame and can easily become several
    # gigabytes.  Instead publish a portable PyMOL *session bundle*: opening
    # this small startup script from the extracted bundle loads the four exact
    # source files, preserves their identities, and starts with the duplicate
    # pocket object hidden.  This gives researchers an immediately usable
    # PyMOL workspace without silently replacing the full trajectories with a
    # down-sampled browser preview.
    pymol_loader = asset_dir / "Open_in_PyMOL.pml"
    pymol_bundle = asset_dir / "PyMACS_PyMOL_full_trajectory_session.zip"
    pymol_session = asset_dir / "PyMACS_centered_trajectory_session.pse"
    available_sources = [name for name in source_files if (asset_dir / name).is_file()]
    if len(available_sources) == len(source_files):
        pymol_loader.write_text(textwrap.dedent("""\
            # PyMACS full-trajectory session
            # Keep this file beside the bundled files, then open it in PyMOL.
            # The PBC-repaired preview is the default visual because raw XTC
            # coordinates can wrap a bound partner across the simulation box.
            # The four exact raw source files are retained in this bundle.
            reinitialize
            load {visualization_name}, pbc_centered_visualization
            load Final_Trajectory.pdb, raw_final_topology
            load binding_pocket_only.pdb, raw_binding_pocket_topology
            hide everything
            show cartoon, pbc_centered_visualization
            color marine, pbc_centered_visualization and chain A
            color tv_yellow, pbc_centered_visualization and chain B
            show sticks, pbc_centered_visualization and chain B
            set stick_radius, 0.16, pbc_centered_visualization and chain B
            set all_states, on
            set cache_frames, 1
            set movie_auto_interpolate, off
            disable raw_final_topology
            disable raw_binding_pocket_topology
            group PyMACS_trajectories, pbc_centered_visualization raw_final_topology raw_binding_pocket_topology
            zoom pbc_centered_visualization
            frame 1
        """.format(visualization_name=visualization_source.name)), encoding="utf-8")
        try:
            bundle_is_current = pymol_bundle.is_file() and pymol_bundle.stat().st_mtime >= max(
                *((asset_dir / name).stat().st_mtime for name in available_sources), visualization_source.stat().st_mtime
            ) and pymol_bundle.stat().st_mtime >= pymol_loader.stat().st_mtime
            if not bundle_is_current:
                temporary_bundle = pymol_bundle.with_suffix(".zip.partial")
                with zipfile.ZipFile(temporary_bundle, "w", compression=zipfile.ZIP_STORED) as archive:
                    archive.write(pymol_loader, pymol_loader.name)
                    archive.write(visualization_source, visualization_source.name)
                    for name in source_files:
                        archive.write(asset_dir / name, name)
                temporary_bundle.replace(pymol_bundle)
        except (OSError, zipfile.BadZipFile) as error:
            print(f"PyMOL session bundle generation skipped: {error}")

    preview_info: dict[str, object] = {}
    needs_preview = not preview.is_file() or preview.stat().st_mtime < max(topology.stat().st_mtime, trajectory.stat().st_mtime)
    if needs_preview:
        temporary = preview.with_suffix(".pdb.partial")
        try:
            import MDAnalysis as mda
            from MDAnalysis.analysis import align

            universe = mda.Universe(str(topology), str(trajectory))
            shown = universe.select_atoms("protein")
            backbone = universe.select_atoms("protein and backbone")
            if not shown.n_atoms or not backbone.n_atoms or not len(universe.trajectory):
                raise ValueError("no protein/peptide atoms available for trajectory preview")
            universe.trajectory[0]
            reference = backbone.positions.copy()
            reference_center = reference.mean(axis=0)
            indices = np.unique(np.linspace(0, len(universe.trajectory) - 1, min(TRAJECTORY_PREVIEW_FRAMES, len(universe.trajectory)), dtype=int))
            times_ns: list[float] = []
            with mda.coordinates.PDB.PDBWriter(str(temporary), multiframe=True) as writer:
                for frame_index in indices:
                    universe.trajectory[int(frame_index)]
                    moving = backbone.positions.copy()
                    moving_center = moving.mean(axis=0)
                    rotation, _ = align.rotation_matrix(moving - moving_center, reference - reference_center)
                    shown.positions = (shown.positions - moving_center) @ rotation.T + reference_center
                    writer.write(shown)
                    times_ns.append(float(universe.trajectory.time) / 1000.0)
            temporary.replace(preview)
            preview_info = {"frames": int(len(indices)), "source_frames": int(len(universe.trajectory)), "times_ns": times_ns, "atoms": int(shown.n_atoms), "centering": "protein-backbone aligned to the first sampled frame"}
            metadata.write_text(json.dumps(preview_info, indent=2), encoding="utf-8")
        except Exception as error:
            print(f"Trajectory preview generation skipped: {error}")
            temporary.unlink(missing_ok=True)
    if not preview_info:
        try:
            preview_info = json.loads(metadata.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            preview_info = {}
    if not preview.is_file():
        return '<section class="trajectory-card"><div class="eyebrow">Trajectory playback</div><h2>Trajectory preview preparation did not complete</h2><p>The original PDB/XTC files are retained, but the browser preview could not be generated on this host.</p></section>'
    # Prefer the GROMACS PBC-repaired coordinates whenever they are present;
    # otherwise retain compatibility with historical dossiers that only have
    # the alignment-only preview.
    browser_preview = pbc_preview if pbc_preview.is_file() else preview
    preview_href = browser_preview.relative_to(analysis).as_posix()
    # 3Dmol uses an absolutely positioned canvas; anchor it to this stage
    # instead of allowing it to escape to the top of the dossier page.
    raw_links = '''<style>.trajectory-stage{position:relative}.trajectory-component-key{display:flex;gap:10px;flex-wrap:wrap;margin:4px 0 12px;color:#cfe8ff;font-size:.82rem;font-weight:800}.trajectory-component-key span{display:inline-flex;align-items:center;gap:6px}.trajectory-component-key i{width:10px;height:10px;border-radius:50%;display:inline-block}</style><div class="trajectory-component-key"><span><i style="background:#ff8a3d"></i>Component 1 · chain A</span><span><i style="background:#18ff8b"></i>Component 2 · chain B</span></div>'''
    if pymol_session.is_file():
        raw_links += (
            f'<a href="{escape(pymol_session.relative_to(analysis).as_posix())}" download>'
            'Download ready-to-open PyMOL session (.pse)</a>'
        )
    if pymol_bundle.is_file():
        raw_links += (
            f'<a href="{escape(pymol_bundle.relative_to(analysis).as_posix())}" download>'
            'Download PyMOL full-session bundle (.zip)</a>'
        )
    raw_links += "".join(
        f'<a href="{escape((asset_dir / name).relative_to(analysis).as_posix())}" download>Download {escape(name)}</a>'
        for name in source_files if (asset_dir / name).is_file()
    )
    peptide_chains = ",".join(context.partner_chain_ids) if context.partner_kind == "peptide" else ""
    component_styles = []
    for index, entry in enumerate(atomindex_components(root)):
        chain = entry.get("current_chain") or entry.get("source_chain")
        if chain:
            component_styles.append({"chain": chain, "label": entry["label"], "color": CHAIN_COLORS[index % len(CHAIN_COLORS)]})
    component_key = "".join(
        f'<span><i style="background:{escape(item["color"])}"></i>{escape(item["label"])} · chain {escape(item["chain"])}</span>'
        for item in component_styles
    )
    playback_note = ("Orange is the first protein component and neon green is the second; the same colors identify these components throughout this dossier."
                     if len(component_styles) > 1 else "Component colors match the chain-resolved plots in this dossier.")
    frame_count = int(preview_info.get("frames", 0))
    source_frame_count = int(preview_info.get("source_frames", frame_count))
    time_values = preview_info.get("times_ns", [])
    last_time = float(time_values[-1]) if isinstance(time_values, list) and time_values else 0.0
    return f'''<style>.trajectory-card{{border-radius:21px;padding:20px;margin-top:18px;background:linear-gradient(155deg,#10172b,#060912);border:1px solid #23375d;box-shadow:0 18px 50px rgba(0,0,0,.35)}}.trajectory-toolbar{{display:flex;align-items:center;gap:9px;flex-wrap:wrap;margin:14px 0 10px}}.trajectory-toolbar button{{border:1px solid #285a79;border-radius:10px;background:#0a253a;color:#dff5ff;padding:9px 12px;font:800 .82rem inherit;cursor:pointer}}.trajectory-toolbar button:hover:not(:disabled){{background:#116887;border-color:#35d0ff}}.trajectory-toolbar button:disabled{{opacity:.46;cursor:not-allowed}}.trajectory-toolbar input{{accent-color:#35d0ff;min-width:min(430px,60vw);flex:1}}.trajectory-stage{{height:min(62vw,650px);min-height:410px;border:1px solid rgba(53,208,255,.30);border-radius:17px;background:#05060a;overflow:hidden}}.trajectory-stage:empty:before{{content:'Load playback to view the centered trajectory';display:grid;place-items:center;height:100%;color:#7aa6d9}}.trajectory-status{{color:#9cb6d0;font-size:.86rem;font-weight:800}}.trajectory-downloads{{display:flex;flex-wrap:wrap;gap:8px;margin-top:13px}}.trajectory-downloads a{{color:#35d0ff;border:1px solid #285a79;border-radius:9px;padding:7px 9px;text-decoration:none;font-size:.78rem;font-weight:800}}</style><section class="trajectory-card" data-trajectory-preview="{escape(preview_href)}" data-peptide-chains="{escape(peptide_chains)}" data-frame-times="{escape(json.dumps(time_values, separators=(',', ':')))}"><div class="section-head"><div><div class="eyebrow">Trajectory playback · on demand</div><h2>Watch the centered simulation, frame by frame</h2><p>This is a solvent-free, {frame_count}-frame preview spanning 0–{last_time:.2f} ns. It is centered and backbone-aligned to the first sampled frame so the complex moves—not the camera. It does not load until you press the button.</p></div></div><div class="trajectory-toolbar"><button type="button" data-trajectory-load>Load interactive playback</button><button type="button" data-trajectory-play disabled>▶ Play</button><button type="button" data-trajectory-reset disabled>Reset view</button><input type="range" min="0" max="{max(0, frame_count - 1)}" value="0" step="1" aria-label="Trajectory frame" data-trajectory-frame disabled><span class="trajectory-status" data-trajectory-status>Not loaded · ordinary dossier load remains fast</span></div><div class="trajectory-stage" data-trajectory-stage></div><details class="guide"><summary>ⓘ ELI5: what am I watching?</summary><p><b>What this is:</b> a deliberately lightweight set of evenly spaced frames from the final MD trajectory. The protein backbone is aligned to frame 1 so ordinary drifting through the simulation box does not look like biological movement.</p><p><b>How to read it:</b> use play for an overview, then pause and drag/zoom to inspect a frame. Components use the dossier palette: chain A is Miami orange and chain B is neon green. Compare apparent motion with the stability and contact plots instead of treating an animation as a binding measurement.</p><p><b>Keep in mind:</b> this is a visual preview, not all {25001:,} raw frames. The exact final and pocket PDB/XTC files are available below.</p></details><div class="trajectory-downloads">{raw_links}</div></section><script>(()=>{{const root=document.currentScript.previousElementSibling,load=root.querySelector('[data-trajectory-load]'),play=root.querySelector('[data-trajectory-play]'),reset=root.querySelector('[data-trajectory-reset]'),slider=root.querySelector('[data-trajectory-frame]'),status=root.querySelector('[data-trajectory-status]'),stage=root.querySelector('[data-trajectory-stage]'),times=JSON.parse(root.dataset.frameTimes||'[]'),componentPalette={{A:'#ff8a3d',B:'#18ff8b',C:'#ff2bd6',D:'#35d0ff',E:'#ffe66d',F:'#a78bfa'}};let viewer,model,frame=0,timer=null;const label=()=>{{let time=Number(times[frame]??0);status.textContent=`Frame ${{frame+1}} / ${{Number(slider.max)+1}} · ${{time.toFixed(2)}} ns`}};const show=next=>{{frame=Math.max(0,Math.min(Number(slider.max),Number(next)));slider.value=frame;viewer.setFrame(frame);viewer.render();label()}};const stop=()=>{{if(timer){{clearInterval(timer);timer=null}}play.textContent='▶ Play'}};load.addEventListener('click',async()=>{{if(!window.$3Dmol){{status.textContent='Local 3D viewer library is unavailable.';return}}load.disabled=true;status.textContent='Loading centered trajectory preview…';try{{const response=await fetch(root.dataset.trajectoryPreview);if(!response.ok)throw new Error('preview fetch failed');const pdb=await response.text();viewer=$3Dmol.createViewer(stage,{{backgroundColor:'#05060a',antialias:true}});model=viewer.addModelsAsFrames(pdb,'pdb');model.setStyle({{}},{{cartoon:{{color:'#496985'}}}});Object.entries(componentPalette).forEach(([chain,color])=>model.setStyle({{chain}},{{cartoon:{{color}}}}));root.dataset.peptideChains.split(',').filter(Boolean).forEach(chain=>model.setStyle({{chain}},{{cartoon:{{color:'#ffe66d'}},stick:{{radius:.16,color:'#ffdf85'}}}}));viewer.zoomTo();viewer.render();slider.disabled=play.disabled=reset.disabled=false;show(0);addEventListener('resize',()=>viewer&&viewer.resize())}}catch(error){{status.textContent='Could not load the trajectory preview. The exact PDB/XTC files remain downloadable.';load.disabled=false}}}});slider.addEventListener('input',()=>{{if(viewer){{stop();show(slider.value)}}}});play.addEventListener('click',()=>{{if(!viewer)return;if(timer){{stop();return}}play.textContent='❚❚ Pause';timer=setInterval(()=>{{let next=frame>=Number(slider.max)?0:frame+1;show(next)}},115)}});reset.addEventListener('click',()=>{{if(viewer){{viewer.zoomTo();viewer.render()}}}})}})();</script>'''


def final_structure_html(root: Path, context: PartnerContext) -> str:
    """Embed the final complex or apo PDB in the report's 3Dmol viewer."""
    candidates = (root / "Final_Trajectory.pdb", root / "binding_pocket_only.pdb")
    structure = next((path for path in candidates if path.is_file()), None)
    # Historical runs often retained the submitted complex but not a final
    # trajectory snapshot.  It remains a useful, explicitly labelled starting
    # structure for visual inspection; never imply that it is an MD endpoint.
    fallback = False
    if structure is None:
        structure = next((path for path in sorted(root.glob("*.pdb"))
                          if path.name not in {"protein.pdb", "protein_pdb2gmx_ready.pdb",
                                               "protein_pdbfixer_ready.pdb", "selected_input_source.pdb"}), None)
        fallback = structure is not None
    if structure is None:
        return (
            '<section class="structure-card"><div class="eyebrow">Final structure</div>'
            '<h2>3D complex / apo viewer</h2><p class="muted">No final PDB was published for this run.</p></section>'
        )
    try:
        payload = base64.b64encode(structure.read_bytes()).decode("ascii")
    except OSError:
        return ""
    metadata_ligand = "".join(ch for ch in context.partner_label if ch.isalnum())[:4] if context.partner_kind == "ligand" else ""
    ligand = metadata_ligand
    final_ligand_payload = ""
    pdb_ligand_resname = ligand
    # A number of older dossiers do not preserve the ligand classification in
    # their run metadata.  A CGenFF MOL2 plus a strict final-PDB atom-name match
    # is stronger evidence than that missing label, so it is a safe fallback.
    cgenff_mol2, _ = ligand_inputs(root, metadata_ligand)
    if cgenff_mol2:
        # Keep the final MD pose from the PDB while recovering the explicit
        # CGenFF bond graph.  If the metadata only says "Ligand", the strict
        # atom-name match below finds the actual 3-letter PDB residue (LIG is
        # the common fallback) without guessing a chemical identity.
        final_ligand, matched_resname = mol2_with_final_pdb_coordinates(structure, cgenff_mol2, metadata_ligand)
        if final_ligand:
            final_ligand_payload = base64.b64encode(final_ligand.encode("utf-8")).decode("ascii")
            pdb_ligand_resname = matched_resname
            ligand = matched_resname
    chain_buttons = "".join(f'<button type="button" data-structure-view="chain" data-chain-id="{escape(chain_id)}">Chain {escape(chain_id)}</button>' for chain_id in context.partner_chain_ids)
    if context.partner_kind == "ligand" or final_ligand_payload:
        focused_button = '<button type="button" data-structure-view="pocket">Ligand + pocket</button>'
    elif context.partner_kind == "peptide":
        focused_button = '<button type="button" data-structure-view="partner">Peptide</button><button type="button" data-structure-view="interface">Peptide interface</button>'
    elif context.partner_kind == "protein_chain":
        focused_button = '<button type="button" data-structure-view="interface">Selected interface</button>'
    else:
        focused_button = ''
    protein_label = "Protein/receptor" if context.partner_kind == "peptide" else "Protein"
    bond_note = " Ligand bonds use the matching CGenFF topology." if final_ligand_payload else ""
    snapshots = representative_snapshots_html(root, root / "Analysis_Results")
    source_note = ("submitted/input complex; no final trajectory PDB was published" if fallback
                   else "final published coordinate snapshot")
    return f'''<section class="structure-card"><div class="section-head"><div><div class="eyebrow">Structure viewer</div><h2>3D complex presentation</h2><p>Drag to rotate and scroll to zoom. This viewer is independent from the contact explorer and controls its own peptide, protein, and side-chain display.{bond_note}</p></div></div><div class="structure-toolbar" role="group" aria-label="Structure focus"><button type="button" data-structure-view="complex">Full complex</button><button type="button" data-structure-view="protein">{protein_label}</button>{chain_buttons}{focused_button}</div><div class="structure-control-grid"><label>Protein color <select data-structure-control="protein-color"><option value="slate">Slate blue</option><option value="spectrum">Spectrum</option><option value="chain">By chain</option><option value="mint">Mint</option><option value="violet">Violet</option></select></label><label>Peptide color <select data-structure-control="peptide-color"><option value="gold">Gold</option><option value="magenta">Magenta</option><option value="cyan">Cyan</option><option value="element">By element</option></select></label><label>Peptide side chains <select data-structure-control="peptide-sidechains"><option value="stick">Sticks</option><option value="line">Lines</option><option value="off">Hidden</option></select></label><label>Protein side chains <select data-structure-control="protein-sidechains"><option value="off">Hidden</option><option value="line">Lines</option><option value="stick">Sticks</option></select></label></div><div class="complex-3d js-complex-3d" data-ligand="{escape(ligand)}" data-pdb-ligand="{escape(pdb_ligand_resname)}" data-final-ligand-b64="{final_ligand_payload}" data-partner-kind="{escape(context.partner_kind)}" data-partner-chains="{escape(','.join(context.partner_chain_ids))}" data-structure-b64="{payload}"></div><p class="asset-note">Source: {escape(structure.name)} · {source_note}.</p>{snapshots}</section>'''


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
            font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 55, "r": 20, "t": 55, "b": 120},
            xaxis={"title": "Time (ns)", "gridcolor": GRID, "zerolinecolor": GRID},
            yaxis={"title": "Fraction of residues", "range": [0, 1], "gridcolor": GRID, "zerolinecolor": GRID},
            hovermode="x unified", height=360, legend=LEGEND_BELOW,
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
    is_picoseconds = x_column.lower() == "time_ps"
    if is_picoseconds:
        # Published dashboards always use nanoseconds; legacy CSVs may still
        # retain a ps column for provenance.
        table[x_column] = table[x_column] / 1000.0
        x_label = "Time (ns)"
    else:
        x_label = x_column.replace("_", " ")
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
            hovertemplate=f"{escape(x_label)}=%{{x:.4g}}<br>{escape(column)}=%{{y:.4g}}<extra></extra>",
        ))
    if not figure.data:
        return None
    figure.update_layout(
        title=title, paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"},
        margin={"l": 60, "r": 20, "t": 55, "b": 120}, hovermode="x unified", height=380,
        xaxis={"title": x_label, "gridcolor": GRID},
        yaxis={"title": y_title, "gridcolor": GRID},
        legend=LEGEND_BELOW,
    )
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def _chain_label_from_metric(source: Path, prefix: str, suffix: str = ".csv") -> str:
    """Recover PyMACS' component label without pretending it is a PDB chain ID."""
    name = source.name
    return name[len(prefix):-len(suffix)] if name.startswith(prefix) and name.endswith(suffix) else source.stem


def chain_rmsd_plot(sources: list[Path], *, title: str = "Chain-resolved RMSD",
                    colors: list[str] | None = None) -> str | None:
    """Overlay one explicitly aligned anchor trace per polymer component."""
    figure = go.Figure()
    palette = colors or CHAIN_COLORS
    for source in sources:
        table = _read_table(source)
        if table is None or "RMSD_Å" not in table.columns:
            continue
        time_column = "Time_ns" if "Time_ns" in table.columns else "Time_ps" if "Time_ps" in table.columns else None
        if time_column is None:
            continue
        table = table.copy()
        table[time_column] = pd.to_numeric(table[time_column], errors="coerce")
        table["RMSD_Å"] = pd.to_numeric(table["RMSD_Å"], errors="coerce")
        table = _sample_table(table.dropna(subset=[time_column, "RMSD_Å"]).sort_values(time_column))
        if len(table) < 2:
            continue
        label = _chain_label_from_metric(source, "RMSD_")
        selection = str(table.get("Atom_Selection", pd.Series(["per-chain atoms"])).iloc[0])
        visible_index = len(figure.data)
        figure.add_trace(go.Scatter(
            x=table[time_column] if time_column == "Time_ns" else table[time_column] / 1000.0, y=table["RMSD_Å"], mode="lines", name=f"{label} ({selection})",
            line={"color": palette[visible_index % len(palette)], "width": 2.5},
            hovertemplate=f"Chain/component={escape(label)}<br>Time=%{{x:.3f}} ns<br>RMSD=%{{y:.3f}} Å<br>Selection={escape(selection)}<extra></extra>",
        ))
    if not figure.data:
        return None
    figure.update_layout(title=title, paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 120},
        hovermode="x unified", height=400, xaxis={"title": "Time (ns)", "gridcolor": GRID},
        yaxis={"title": "RMSD (Å)", "gridcolor": GRID}, legend=LEGEND_BELOW)
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def chain_rmsf_plot(sources: list[Path], *, title: str = "Chain-resolved residue RMSF",
                    colors: list[str] | None = None) -> str | None:
    """Plot chain-local positions while exposing the source residue numbers on hover.

    Different chains often have independent numbering (or non-1 starts).  A
    shared absolute x-axis therefore invents a visual relationship and causes
    misleading connecting lines.  Each trace begins at its own local position;
    hover always reports the deposited/source residue identifier.
    """
    figure = go.Figure()
    palette = colors or CHAIN_COLORS
    for source in sources:
        table = _read_table(source)
        if table is None or not {"Residue", "RMSF_Å"}.issubset(table.columns):
            continue
        table = table.copy()
        table["RMSF_Å"] = pd.to_numeric(table["RMSF_Å"], errors="coerce")
        table = table.dropna(subset=["RMSF_Å"]).reset_index(drop=True)
        if len(table) < 2:
            continue
        label = _chain_label_from_metric(source, "RMSF_", "_per_residue.csv")
        selection = str(table.get("Anchor_Atom", pd.Series(["residue anchor"])).iloc[0])
        original_residue = table.get("Original_Residue", table["Residue"]).astype(str)
        visible_index = len(figure.data)
        figure.add_trace(go.Scatter(
            x=np.arange(1, len(table) + 1), y=table["RMSF_Å"], mode="lines+markers", name=f"{label} ({selection})",
            customdata=np.column_stack([original_residue]),
            line={"color": palette[visible_index % len(palette)], "width": 2.5}, marker={"size": 4},
            hovertemplate=f"Chain/component={escape(label)}<br>Chain-local position=%{{x}}<br>Original residue=%{{customdata[0]}}<br>RMSF=%{{y:.3f}} Å<br>Anchor={escape(selection)}<extra></extra>",
        ))
    if not figure.data:
        return None
    figure.update_layout(title=title, paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 120},
        hovermode="closest", height=420, xaxis={"title": "Chain-local residue position (see hover for original residue)", "gridcolor": GRID},
        yaxis={"title": "RMSF (Å)", "gridcolor": GRID}, legend=LEGEND_BELOW)
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def _interface_time_data(analysis: Path) -> tuple[pd.DataFrame, Path]:
    """Load the canonical pair-level interface timeline without inventing zeros."""
    source = analysis / "Interface_RIN" / "interface_contact_timeseries.csv"
    try:
        table = pd.read_csv(source)
    except (OSError, ValueError, pd.errors.ParserError):
        return pd.DataFrame(), source
    required = {"Time_ns", "ContactCount", "MinDistance_A"}
    if not required.issubset(table.columns):
        return pd.DataFrame(), source
    table = table.copy()
    for column in required:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table = table.dropna(subset=list(required)).sort_values("Time_ns")
    return table, source


def _coupled_figure(figure: go.Figure, title: str, *, height: int = 720) -> str:
    """Apply one spacious, legible style to the coupled-dynamics panels."""
    figure.update_layout(
        title=title, paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"},
        margin={"l": 72, "r": 68, "t": 78, "b": 82}, height=height,
        legend=LEGEND_BELOW, hovermode="x unified",
    )
    figure.update_xaxes(gridcolor=GRID, title="Time (ns)")
    figure.update_yaxes(gridcolor=GRID)
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def coupled_engagement_plot(analysis: Path) -> str | None:
    """Show interface contact count, closest gap, and native-fold retention together."""
    interface, _ = _interface_time_data(analysis)
    q_source = analysis / "protein_native_contact_q.csv"
    if interface.empty or not q_source.is_file():
        return None
    q_data = _read_table(q_source)
    if q_data is None or not {"Time_ns", "Q_protein_native_contact_fraction"}.issubset(q_data.columns):
        return None
    q_data = q_data.copy()
    q_data["Time_ns"] = pd.to_numeric(q_data["Time_ns"], errors="coerce")
    q_data["Q_protein_native_contact_fraction"] = pd.to_numeric(q_data["Q_protein_native_contact_fraction"], errors="coerce")
    q_data = q_data.dropna().sort_values("Time_ns")
    interface, q_data = _sample_table(interface), _sample_table(q_data)
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=.08,
                        subplot_titles=("Interface contacts", "Closest interface separation", "Protein native-contact retention (Q)"))
    fig.add_trace(go.Scattergl(x=interface["Time_ns"], y=interface["ContactCount"], mode="lines", name="Interface contacts", line={"color": CYAN, "width": 2}, hovertemplate="Time=%{x:.3f} ns<br>Residue-pair contacts=%{y:.0f}<extra></extra>"), row=1, col=1)
    fig.add_trace(go.Scattergl(x=interface["Time_ns"], y=interface["MinDistance_A"], mode="lines", name="Closest gap", line={"color": CYAN, "width": 2}, hovertemplate="Time=%{x:.3f} ns<br>Closest interface gap=%{y:.3f} Å<extra></extra>"), row=2, col=1)
    fig.add_trace(go.Scattergl(x=q_data["Time_ns"], y=q_data["Q_protein_native_contact_fraction"], mode="lines", name="Native-contact Q", line={"color": CYAN, "width": 2}, hovertemplate="Time=%{x:.3f} ns<br>Native-contact Q=%{y:.3f}<extra></extra>"), row=3, col=1)
    fig.update_yaxes(title="Contacts", row=1, col=1)
    fig.update_yaxes(title="Gap (Å)", row=2, col=1)
    fig.update_yaxes(title="Q fraction", range=[0, 1.05], row=3, col=1)
    return _coupled_figure(fig, "Shared interface engagement and fold retention", height=790)


def coupled_motion_binding_plot(analysis: Path) -> str | None:
    """Relate shared interface contacts to each component's compactness and motion."""
    interface, _ = _interface_time_data(analysis)
    if interface.empty:
        return None
    interface = interface[["Time_ns", "ContactCount"]].sort_values("Time_ns")
    traces: list[tuple[str, pd.DataFrame, str, str]] = []
    for index, source in enumerate(sorted(analysis.glob("Radius_of_Gyration_*.csv"))):
        if source.name == "Radius_of_Gyration_Protein.csv":
            continue
        table = _read_table(source)
        if table is None or not {"Time_ns", "Rg_Component_Å"}.issubset(table.columns):
            continue
        label = source.stem.removeprefix("Radius_of_Gyration_")
        table = table[["Time_ns", "Rg_Component_Å"]].copy()
        table["Time_ns"] = pd.to_numeric(table["Time_ns"], errors="coerce")
        table["Rg_Component_Å"] = pd.to_numeric(table["Rg_Component_Å"], errors="coerce")
        traces.append((label, table.dropna().sort_values("Time_ns"), "Rg_Component_Å", CHAIN_COLORS[index % len(CHAIN_COLORS)]))
    rmsd_sources = sorted(source for source in analysis.glob("RMSD_*.csv") if source.name != "RMSD_Protein.csv")
    for index, source in enumerate(rmsd_sources):
        table = _read_table(source)
        if table is None or not {"Time_ns", "RMSD_Å"}.issubset(table.columns):
            continue
        label = source.stem.removeprefix("RMSD_")
        table = table[["Time_ns", "RMSD_Å"]].copy()
        table["Time_ns"] = pd.to_numeric(table["Time_ns"], errors="coerce")
        table["RMSD_Å"] = pd.to_numeric(table["RMSD_Å"], errors="coerce")
        traces.append((label, table.dropna().sort_values("Time_ns"), "RMSD_Å", CHAIN_COLORS[index % len(CHAIN_COLORS)]))
    if not traces:
        return None
    fig = make_subplots(rows=2, cols=1, vertical_spacing=.12,
                        subplot_titles=("Component compactness versus interface contacts", "Component RMSD versus interface contacts"))
    for label, table, metric, color in traces:
        row = 1 if metric == "Rg_Component_Å" else 2
        merged = pd.merge_asof(interface, table, on="Time_ns", direction="nearest")
        merged = _sample_table(merged.dropna(subset=[metric]))
        if merged.empty:
            continue
        fig.add_trace(go.Scattergl(x=merged["ContactCount"], y=merged[metric], mode="markers", name=f"{label} {'Rg' if row == 1 else 'RMSD'}", marker={"color": color, "size": 5, "opacity": .42}, hovertemplate=f"Component={escape(label)}<br>Interface contacts=%{{x:.0f}}<br>{'Rg' if row == 1 else 'RMSD'}=%{{y:.3f}} Å<extra></extra>"), row=row, col=1)
    fig.update_xaxes(title="Interface residue-pair contacts", row=1, col=1)
    fig.update_xaxes(title="Interface residue-pair contacts", row=2, col=1)
    fig.update_yaxes(title="Radius of gyration (Å)", row=1, col=1)
    fig.update_yaxes(title="RMSD (Å)", row=2, col=1)
    fig.update_layout(title="Component motion versus interface engagement", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 72, "r": 32, "t": 78, "b": 120}, height=700, legend=LEGEND_BELOW)
    fig.update_xaxes(gridcolor=GRID); fig.update_yaxes(gridcolor=GRID)
    return plot(fig, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def coupled_chemistry_plot(analysis: Path) -> str | None:
    """Bin validated atom-level chemistry into a readable time-resolved signal."""
    source = analysis / "Interface_RIN" / "ValidatedInterfaceInteractionEvents.csv"
    if not source.is_file():
        return None
    pieces = []
    try:
        for chunk in pd.read_csv(source, usecols=["Frame", "Time_ns", "InteractionClass"], chunksize=200_000):
            chunk = chunk[chunk["InteractionClass"].isin(["H-bond", "Ionic", "Hydrophobic"])].copy()
            if not chunk.empty:
                chunk["Time_ns"] = pd.to_numeric(chunk["Time_ns"], errors="coerce")
                chunk["Frame"] = pd.to_numeric(chunk["Frame"], errors="coerce")
                pieces.append(chunk.dropna())
    except (OSError, ValueError, pd.errors.ParserError):
        return None
    if not pieces:
        return None
    events = pd.concat(pieces, ignore_index=True)
    events["Time bin (ns)"] = np.floor(events["Time_ns"] / 2.0) * 2.0
    counts = events.groupby(["Time bin (ns)", "InteractionClass"])["Frame"].nunique().unstack(fill_value=0).sort_index()
    if counts.empty:
        return None
    fig = go.Figure()
    for kind, color in (("H-bond", CYAN), ("Ionic", MAGENTA), ("Hydrophobic", ORANGE)):
        if kind in counts:
            fig.add_trace(go.Scatter(x=counts.index, y=counts[kind], stackgroup="chemistry", mode="lines", name=kind, line={"color": color, "width": 1.6}, hovertemplate=f"Time bin=%{{x:.1f}} ns<br>{kind} frames=%{{y:.0f}}<extra></extra>"))
    fig.update_layout(title="Validated interface chemistry through time", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 68, "r": 30, "t": 65, "b": 110}, height=450, legend=LEGEND_BELOW, hovermode="x unified", xaxis={"title": "Time (ns; 2 ns bins)", "gridcolor": GRID}, yaxis={"title": "Frames with validated chemistry", "gridcolor": GRID})
    return plot(fig, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def coupled_rewiring_plot(analysis: Path) -> str | None:
    """Quantify whether the same residue-pair network survives from frame to frame."""
    source = analysis / "Interface_RIN" / "interface_residue_pair_framewise.csv.gz"
    if not source.is_file():
        return None
    frames: dict[int, tuple[float, set[str]]] = {}
    try:
        for chunk in pd.read_csv(source, compression="gzip", usecols=["Frame", "Time_ns", "Pair_ID"], chunksize=200_000):
            chunk["Frame"] = pd.to_numeric(chunk["Frame"], errors="coerce")
            chunk["Time_ns"] = pd.to_numeric(chunk["Time_ns"], errors="coerce")
            for frame, block in chunk.dropna().groupby("Frame", sort=False):
                frame_id = int(frame)
                prior = frames.get(frame_id)
                pairs = set(block["Pair_ID"].astype(str))
                frames[frame_id] = (float(block["Time_ns"].iloc[0]), pairs if prior is None else prior[1] | pairs)
    except (OSError, ValueError, pd.errors.ParserError):
        return None
    if len(frames) < 2:
        return None
    rows = []
    previous: set[str] | None = None
    for _, (time, pairs) in sorted(frames.items()):
        shared = len(pairs & previous) if previous is not None else 0
        union = len(pairs | previous) if previous is not None else 0
        rows.append((time, len(pairs), shared / union if union else np.nan, len(pairs - (previous or set())), len((previous or set()) - pairs)))
        previous = pairs
    table = pd.DataFrame(rows, columns=["Time_ns", "Active pairs", "Pair retention", "New pairs", "Lost pairs"])
    table = _sample_table(table)
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=.11, subplot_titles=("Active residue-pair contacts", "Frame-to-frame interface rewiring"))
    fig.add_trace(go.Scattergl(x=table["Time_ns"], y=table["Active pairs"], mode="lines", name="Active pairs", line={"color": CYAN, "width": 2}, hovertemplate="Time=%{x:.3f} ns<br>Active residue pairs=%{y:.0f}<extra></extra>"), row=1, col=1)
    fig.add_trace(go.Scattergl(x=table["Time_ns"], y=table["Pair retention"], mode="lines", name="Pair retention", line={"color": CYAN, "width": 2}, hovertemplate="Time=%{x:.3f} ns<br>Shared pair retention=%{y:.1%}<extra></extra>"), row=2, col=1)
    fig.add_trace(go.Scattergl(x=table["Time_ns"], y=table["New pairs"], mode="lines", name="New pairs", line={"color": GREEN, "width": 1.4}, hovertemplate="Time=%{x:.3f} ns<br>New pairs=%{y:.0f}<extra></extra>"), row=2, col=1)
    fig.add_trace(go.Scattergl(x=table["Time_ns"], y=table["Lost pairs"], mode="lines", name="Lost pairs", line={"color": ORANGE, "width": 1.4}, hovertemplate="Time=%{x:.3f} ns<br>Lost pairs=%{y:.0f}<extra></extra>"), row=2, col=1)
    fig.update_yaxes(title="Pairs", row=1, col=1); fig.update_yaxes(title="Retention / pair count", row=2, col=1)
    return _coupled_figure(fig, "Interface contact rewiring", height=690)


def coupled_sse_binding_plot(sse_source: Path, interface: pd.DataFrame, component_label: str, color: str) -> str | None:
    """Give one component a local SSE lane while retaining the shared contact reference."""
    table = _read_table(sse_source)
    required = {"Time_ns", "Helix", "Sheet", "Coil"}
    if table is None or interface.empty or not required.issubset(table.columns):
        return None
    table = table[list(required)].copy()
    for column in required:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table = _sample_table(table.dropna().sort_values("Time_ns"))
    interface = _sample_table(interface[["Time_ns", "ContactCount"]].dropna().sort_values("Time_ns"))
    if len(table) < 2:
        return None
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    for column, line_color in (("Helix", CYAN), ("Sheet", GREEN), ("Coil", MAGENTA)):
        fig.add_trace(go.Scattergl(x=table["Time_ns"], y=table[column], mode="lines", name=column, line={"color": line_color, "width": 1.8}, hovertemplate=f"Component={escape(component_label)}<br>Time=%{{x:.3f}} ns<br>{column} fraction=%{{y:.3f}}<extra></extra>"), secondary_y=False)
    fig.add_trace(go.Scattergl(x=interface["Time_ns"], y=interface["ContactCount"], mode="lines", name="Interface contacts", line={"color": color, "width": 2.5, "dash": "dot"}, hovertemplate="Time=%{x:.3f} ns<br>Interface contacts=%{y:.0f}<extra></extra>"), secondary_y=True)
    fig.update_layout(title=f"{component_label}: secondary structure versus interface engagement", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 68, "r": 68, "t": 65, "b": 110}, height=450, legend=LEGEND_BELOW, hovermode="x unified")
    fig.update_xaxes(title="Time (ns)", gridcolor=GRID); fig.update_yaxes(title="SSE fraction", range=[0, 1.05], gridcolor=GRID, secondary_y=False); fig.update_yaxes(title="Interface contacts", gridcolor=GRID, secondary_y=True)
    return plot(fig, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


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
        yaxis={"title": y_title, "gridcolor": GRID}, height=420, legend=LEGEND_BELOW,
    )
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def chain_sse_persistence_plot(csv_file: Path, chain_label: str) -> str | None:
    """Render one chain-local SSE persistence profile with source IDs on hover."""
    table = _read_table(csv_file)
    required = {"Chain_Local_Position", "Original_Residue", "Helix_Persistence", "Sheet_Persistence", "Coil_Persistence"}
    if table is None or not required.issubset(table.columns):
        return None
    table = table.copy().head(250)
    position = pd.to_numeric(table["Chain_Local_Position"], errors="coerce")
    table = table.loc[position.notna()].copy()
    if table.empty:
        return None
    figure = go.Figure()
    for column, color, label in (
        ("Coil_Persistence", MAGENTA, "Coil persistence"),
        ("Sheet_Persistence", GREEN, "Sheet persistence"),
        ("Helix_Persistence", CYAN, "Helix persistence"),
    ):
        values = pd.to_numeric(table[column], errors="coerce").fillna(0)
        figure.add_trace(go.Bar(
            x=table["Chain_Local_Position"], y=values, name=label, marker_color=color,
            customdata=np.column_stack([table["Original_Residue"].astype(str)]),
            hovertemplate=f"Chain={escape(chain_label)}<br>Chain-local position=%{{x}}<br>Original residue=%{{customdata[0]}}<br>{label}=%{{y:.3f}}<extra></extra>",
        ))
    figure.update_layout(
        title=f"Chain {chain_label} residue SSE persistence", barmode="stack",
        paper_bgcolor=PANEL, plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"},
        margin={"l": 60, "r": 20, "t": 55, "b": 120}, height=430,
        xaxis={"title": "Chain-local residue position (see hover for original residue)", "gridcolor": GRID},
        yaxis={"title": "Fraction", "range": [0, 1], "gridcolor": GRID}, legend=LEGEND_BELOW,
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


def interface_scorecard_html(analysis: Path) -> str:
    """Plain-language, five-row answer to: what holds this interface together?"""
    source = analysis / "Interface_RIN" / "interface_residue_pair_contacts.csv"
    table = _read_table(source)
    needed = {"ChainA_Label", "ChainB_Label", "ResidueA", "ResidueB", "ContactFraction", "MeanMinDistance_A"}
    if table is None or not needed.issubset(table.columns):
        return ""
    table = table.copy()
    peptide_a = table["ChainA_Label"].astype(str).str.contains(r"\(peptide,", case=False, regex=True)
    peptide_b = table["ChainB_Label"].astype(str).str.contains(r"\(peptide,", case=False, regex=True)
    table = table[peptide_a | peptide_b].copy()
    if table.empty:
        return ""
    table["Peptide position"] = pd.Series(np.where(peptide_a, table["ResidueA"], table["ResidueB"]), index=table.index).astype(str).str.rsplit(":", n=1).str[-1]
    table["Protein partner"] = pd.Series(np.where(peptide_a, table["ResidueB"], table["ResidueA"]), index=table.index).astype(str).str.rsplit(":", n=1).str[-1]
    # Interface source tables can retain waters/ions as peptide-chain records.
    # They are useful provenance, but never meaningful peptide positions.
    canonical_peptide = r"^(?:ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL)-?\d+$"
    table = table[table["Peptide position"].str.match(canonical_peptide, case=False, na=False)].copy()
    table["ContactFraction"] = pd.to_numeric(table["ContactFraction"], errors="coerce")
    table["MeanMinDistance_A"] = pd.to_numeric(table["MeanMinDistance_A"], errors="coerce")
    table = table.dropna(subset=["ContactFraction", "MeanMinDistance_A"]).sort_values("ContactFraction", ascending=False).head(5)
    if table.empty:
        return ""
    rows = "".join(
        f'<tr><td>{escape(row["Peptide position"])}</td><td>{escape(row["Protein partner"])}</td>'
        f'<td>{row["ContactFraction"] * 100:.1f}%</td><td>{row["MeanMinDistance_A"]:.2f} Å</td></tr>'
        for _, row in table.iterrows()
    )
    metric_note = ('<details class="guide"><summary>ⓘ Why can this table and the colored network disagree?</summary>'
                   '<p><b>This table uses geometric proximity:</b> it counts a pair whenever any peptide/protein heavy atoms are within the 4 Å interface cutoff.</p>'
                   '<p><b>The interactive network is stricter:</b> it shows a specific classified interaction (hydrogen bond, ionic, hydrophobic, or other) only when that interaction passes its current occupancy filter. A pair can stay physically close in every frame but form a particular hydrogen bond in fewer than 50% of frames.</p>'
                   '<p>So a row here is an interface-nearness clue. The network answers the narrower question: which <em>specific interaction type</em> remains persistent at the selected threshold?</p></details>')
    return f'''<style>.interface-scorecard,.representative-frames{{margin:0 0 18px;padding:18px;border:1px solid rgba(24,255,139,.32);border-radius:17px;background:linear-gradient(145deg,#0d2130,#08101c)}}.interface-scorecard h2,.representative-frames h2{{margin:.2rem 0 .35rem}}.interface-scorecard p,.representative-frames p{{color:#9cb6d0;margin:.25rem 0 .8rem;max-width:900px}}.interface-table{{overflow:auto}}.interface-table table{{width:100%;border-collapse:collapse;font-size:.9rem}}.interface-table th{{color:#7aa6d9;text-align:left;text-transform:uppercase;font-size:.7rem;letter-spacing:.06em}}.interface-table td,.interface-table th{{padding:9px;border-bottom:1px solid #203653;white-space:nowrap}}.interface-table td:first-child{{color:#ffe66d;font-weight:900}}</style><section class="interface-scorecard"><div><div class="eyebrow">Start here · geometric interface anchors</div><h2>Which residue pairs stay physically near each other?</h2><p>These are the five peptide–protein pairs that were within the 4 Å heavy-atom interface cutoff most often. This is a proximity summary, not the same filter used by the colored interaction network below.</p></div><div class="interface-table"><table><thead><tr><th>Peptide</th><th>Protein partner</th><th>Geometrically near</th><th>Typical closest gap</th></tr></thead><tbody>{rows}</tbody></table></div>{metric_note}{guide_html("contact_frequency")}</section>'''


def peptide_importance_plot(analysis: Path) -> str | None:
    """Rank peptide residues by the number and strength of retained partners."""
    source = analysis / "Interface_RIN" / "interface_node_summary.csv"
    table = _read_table(source)
    needed = {"Chain_Label", "Residue_Label", "Degree", "WeightedDegree", "MaxContactFraction"}
    if table is None or not needed.issubset(table.columns):
        return None
    table = table[table["Chain_Label"].astype(str).str.contains(r"\(peptide,", case=False, regex=True)].copy()
    # Keep the explanatory ranking aligned with the visible peptide sequence:
    # omit non-polymer entries such as terminal SOL/water records.
    canonical_peptide = r":(?:ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL)-?\d+$"
    table = table[table["Residue_Label"].astype(str).str.contains(canonical_peptide, case=False, regex=True)].copy()
    if table.empty:
        return None
    table["Degree"] = pd.to_numeric(table["Degree"], errors="coerce").fillna(0)
    table["WeightedDegree"] = pd.to_numeric(table["WeightedDegree"], errors="coerce").fillna(0)
    table["MaxContactFraction"] = pd.to_numeric(table["MaxContactFraction"], errors="coerce").fillna(0)
    table["Position"] = table["Residue_Label"].astype(str).str.rsplit(":", n=1).str[-1]
    table = table.sort_values("WeightedDegree", ascending=True)
    figure = go.Figure(go.Bar(
        y=table["Position"], x=table["WeightedDegree"], orientation="h", marker_color=MAGENTA,
        customdata=np.column_stack([table["Degree"], table["MaxContactFraction"] * 100]),
        hovertemplate="Peptide residue=%{y}<br>Combined contact strength=%{x:.2f}<br>Different protein partners=%{customdata[0]:.0f}<br>Strongest single partner=%{customdata[1]:.1f}% of frames<extra></extra>",
    ))
    figure.update_layout(title="Which peptide positions matter most to this interface?", paper_bgcolor=PANEL, plot_bgcolor=PANEL,
                         font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, height=max(390, 27 * len(table)),
                         margin={"l": 95, "r": 25, "t": 65, "b": 55},
                         xaxis={"title": "Combined contact strength (higher = more or stronger recurring partners)", "gridcolor": GRID},
                         yaxis={"title": "Peptide residue", "gridcolor": GRID})
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def peptide_contact_regime_plot(analysis: Path) -> str | None:
    """Show which peptide positions have any recorded receptor contact over time."""
    source = analysis / "Interface_RIN" / "interface_residue_pair_presence.npz"
    if not source.is_file():
        return None
    try:
        # This is a locally generated PyMACS archive; label arrays are stored
        # as NumPy object arrays, so loading them requires pickle support.
        payload = np.load(source, allow_pickle=True)
        matrix = np.asarray(payload["contact_matrix"], dtype=np.uint8)
        left = [str(value) for value in np.asarray(payload["residue_a_labels"]).tolist()]
        right = [str(value) for value in np.asarray(payload["residue_b_labels"]).tolist()]
        times = np.asarray(payload["times_ns"], dtype=float)
    except (OSError, KeyError, ValueError):
        return None
    if matrix.ndim != 2 or not matrix.size or matrix.shape[1] != len(times):
        return None
    peptide_pattern = re.compile(r"\(peptide, chain [^)]+\):([A-Z]{3})(-?\d+)$", re.I)
    rows: dict[str, list[int]] = {}
    for index, (a, b) in enumerate(zip(left, right)):
        peptide = a if peptide_pattern.search(a) else b if peptide_pattern.search(b) else ""
        match = peptide_pattern.search(peptide)
        if not match:
            continue
        label = f"{match.group(1).upper()}{match.group(2)}"
        if match.group(1).upper() not in {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}:
            continue
        rows.setdefault(label, []).append(index)
    if not rows:
        return None
    ordered = sorted(rows, key=lambda label: int(re.search(r"-?\d+$", label).group()))
    values = np.vstack([(matrix[indices].max(axis=0) > 0).astype(float) for indices in (rows[label] for label in ordered)])
    sample = sampled_indices(values.shape[1])
    figure = go.Figure(go.Heatmap(z=values[:, sample], x=times[sample], y=ordered, zmin=0, zmax=1,
                                  colorscale=[[0, "#101827"], [.49, "#101827"], [.5, CYAN], [1, GREEN]],
                                  colorbar={"title": "Any contact", "tickmode": "array", "tickvals": [0, 1], "ticktext": ["No", "Yes"]},
                                  hovertemplate="Time=%{x:.3f} ns<br>Peptide residue=%{y}<br>Any recorded protein contact=%{z}<extra></extra>"))
    figure.update_layout(title="When does each peptide position touch the protein?", paper_bgcolor=PANEL, plot_bgcolor=PANEL,
                         font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, height=max(390, 26 * len(ordered)),
                         margin={"l": 90, "r": 25, "t": 65, "b": 55},
                         xaxis={"title": "Time (ns)", "gridcolor": GRID}, yaxis={"title": "Peptide residue", "gridcolor": GRID})
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def representative_interface_frames_html(analysis: Path) -> str:
    """Name easy-to-check trajectory moments from published interface data."""
    source = analysis / "Interface_RIN" / "interface_contact_timeseries.csv"
    table = _read_table(source)
    needed = {"Time_ns", "Frame", "ContactCount", "MinDistance_A"}
    if table is None or not needed.issubset(table.columns):
        return ""
    table = table.copy()
    for column in ("Time_ns", "Frame", "ContactCount", "MinDistance_A"):
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table = table.dropna(subset=list(needed))
    if table.empty:
        return ""
    selected = [("Most-connected moment", table.loc[table["ContactCount"].idxmax()]),
                ("Typical moment", table.iloc[(table["ContactCount"] - table["ContactCount"].median()).abs().argsort().iloc[0]]),
                ("Closest moment", table.loc[table["MinDistance_A"].idxmin()])]
    cards = "".join(f'<div class="frame-card"><b>{escape(name)}</b><span>Frame {int(row["Frame"]):,} · {row["Time_ns"]:.2f} ns</span><small>{int(row["ContactCount"])} residue-pair contacts · closest gap {row["MinDistance_A"]:.2f} Å</small></div>' for name, row in selected)
    return f'''<style>.frame-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px}}.frame-card{{border:1px solid #285a79;border-radius:12px;background:#07111e;padding:12px}}.frame-card b,.frame-card span,.frame-card small{{display:block}}.frame-card b{{color:#35d0ff}}.frame-card span{{margin-top:5px;font-weight:800}}.frame-card small{{margin-top:5px;color:#9cb6d0;line-height:1.35}}</style><section class="representative-frames"><div class="eyebrow">Three easy checkpoints</div><h2>Which moments should I inspect?</h2><p>These are real frames from the analyzed trajectory. Use their time or frame number when checking the Structure and interface panels.</p><div class="frame-grid">{cards}</div>{guide_html("contact_frequency")}</section>'''


def distance_carousel_html(analysis: Path) -> str:
    """A compact, evidence-first gallery for interpreting interface distances.

    Every panel is derived from recorded contact events or the explicitly
    published interface time series.  Missing frames are never converted to
    zero-distance observations, which prevents solvent/placeholder artifacts
    from becoming visually persuasive spokes or traces.
    """
    events_path = analysis / "AllContacts_Framewise.csv"
    interface_events_path = analysis / "Interface_RIN" / "interface_residue_pair_framewise.csv.gz"
    if events_path.is_file():
        try:
            events = pd.read_csv(events_path, usecols=lambda column: column in {"Time_ns", "Residue", "Distance"})
        except (OSError, ValueError, pd.errors.ParserError):
            return ""
    elif interface_events_path.is_file():
        # Protein complexes and peptide complexes both persist exact pairwise
        # interface distances in the 3C framewise table.  Retain peptide-only
        # labels for peptide dossiers, but retain labelled residue *pairs* for
        # protein–protein systems instead of leaving their Distance tab blank.
        try:
            raw = pd.read_csv(interface_events_path, compression="gzip", usecols=lambda column: column in {
                "Time_ns", "ChainA_Label", "ChainB_Label", "ResidueA", "ResidueB", "MinDistance_A"
            })
        except (OSError, ValueError, pd.errors.ParserError):
            return ""
        needed = {"Time_ns", "ChainA_Label", "ChainB_Label", "ResidueA", "ResidueB", "MinDistance_A"}
        if not needed.issubset(raw.columns):
            return ""
        peptide_a = raw["ChainA_Label"].astype(str).str.contains(r"\(peptide,", case=False, regex=True)
        peptide_b = raw["ChainB_Label"].astype(str).str.contains(r"\(peptide,", case=False, regex=True)
        amino = r":(?:ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL)-?\d+$"
        valid_a = raw["ResidueA"].astype(str).str.contains(amino, case=False, regex=True)
        valid_b = raw["ResidueB"].astype(str).str.contains(amino, case=False, regex=True)
        is_peptide_interface = bool((peptide_a | peptide_b).any())
        if is_peptide_interface:
            raw = raw[(peptide_a & valid_a) | (peptide_b & valid_b)].copy()
            events = pd.DataFrame({
                "Time_ns": raw["Time_ns"],
                "Residue": np.where(peptide_a.loc[raw.index], raw["ResidueB"], raw["ResidueA"]),
                "Distance": raw["MinDistance_A"],
            })
        else:
            raw = raw[valid_a & valid_b].copy()
            events = pd.DataFrame({
                "Time_ns": raw["Time_ns"],
                "Residue": raw["ResidueA"].astype(str) + " ↔ " + raw["ResidueB"].astype(str),
                "Distance": raw["MinDistance_A"],
            })
        if events.empty:
            return ""
        events_path = interface_events_path
    else:
        return ""
    required = {"Time_ns", "Residue", "Distance"}
    if not required.issubset(events.columns):
        return ""
    events["Time_ns"] = pd.to_numeric(events["Time_ns"], errors="coerce")
    events["Distance"] = pd.to_numeric(events["Distance"], errors="coerce")
    events = events.dropna(subset=["Time_ns", "Residue", "Distance"])
    events = events[(events["Distance"] > 0) & np.isfinite(events["Distance"])]
    if events.empty:
        return ""

    def graph_div(figure: go.Figure) -> str:
        figure.update_layout(
            paper_bgcolor=PANEL, plot_bgcolor=PANEL,
            font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"},
            margin={"l": 68, "r": 26, "t": 72, "b": 58}, height=455,
        )
        return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})

    slides: list[tuple[str, str, str]] = []
    timeseries = analysis / "Interface_RIN" / "interface_contact_timeseries.csv"
    try:
        trace = pd.read_csv(timeseries)
    except (OSError, ValueError, pd.errors.ParserError):
        trace = pd.DataFrame()
    if {"Time_ns", "MinDistance_A"}.issubset(trace.columns):
        trace = trace.copy()
        trace["Time_ns"] = pd.to_numeric(trace["Time_ns"], errors="coerce")
        trace["MinDistance_A"] = pd.to_numeric(trace["MinDistance_A"], errors="coerce")
        trace = trace.dropna(subset=["Time_ns", "MinDistance_A"])
        trace = trace[(trace["MinDistance_A"] > 0) & np.isfinite(trace["MinDistance_A"])]
    if trace.empty:
        trace = events.groupby("Time_ns", as_index=False)["Distance"].min().rename(columns={"Distance": "MinDistance_A"})
        trace_label = "minimum observed contact-event separation"
    else:
        trace = trace.sort_values("Time_ns")
        trace_label = "closest interface separation"
    if len(trace) > MAX_RENDER_POINTS:
        trace = trace.iloc[np.linspace(0, len(trace) - 1, MAX_RENDER_POINTS, dtype=int)]
    median = float(trace["MinDistance_A"].median())
    q05, q95 = trace["MinDistance_A"].quantile([.05, .95]).astype(float)
    figure = go.Figure()
    figure.add_trace(go.Scattergl(x=trace["Time_ns"], y=trace["MinDistance_A"], mode="lines", name=trace_label,
                                  line={"color": ORANGE, "width": 1.8},
                                  hovertemplate="Time=%{x:.3f} ns<br>Closest separation=%{y:.3f} Å<extra></extra>"))
    figure.add_hline(y=median, line_dash="dash", line_color=CYAN, annotation_text=f"median {median:.2f} Å")
    figure.update_layout(title="1 · Closest interface separation over time",
                         xaxis={"title": "Time (ns)", "gridcolor": GRID},
                         yaxis={"title": "Closest separation (Å)", "gridcolor": GRID},
                         annotations=[{"text": f"Observed range (5th–95th percentile): {q05:.2f}–{q95:.2f} Å · no missing frame is plotted as zero", "showarrow": False, "x": .5, "xref": "paper", "y": 1.16, "yref": "paper", "font": {"color": MUTED}}])
    slides.append(("Closest separation", graph_div(figure), "distance_trace"))

    per_residue = events.groupby("Residue")["Distance"].agg(["count", "median", "mean"]).sort_values("count", ascending=False).head(12)
    selected = events[events["Residue"].isin(per_residue.index)].copy()
    figure = go.Figure()
    for index, residue in enumerate(per_residue.index):
        values = selected.loc[selected["Residue"] == residue, "Distance"].to_numpy(dtype=float)
        if len(values) > MAX_EVENT_POINTS // max(1, len(per_residue)):
            values = values[np.linspace(0, len(values) - 1, MAX_EVENT_POINTS // len(per_residue), dtype=int)]
        figure.add_trace(go.Violin(x=values, y=[str(residue)] * len(values), orientation="h", name=str(residue), showlegend=False,
                                   box_visible=True, meanline_visible=True, points=False,
                                   line_color=[CYAN, GREEN, MAGENTA, ORANGE, "#a78bfa"][index % 5],
                                   hovertemplate=f"Residue={escape(str(residue))}<br>Observed separation=%{{x:.3f}} Å<extra></extra>"))
    figure.update_layout(title="2 · Separation distributions for the 12 most-observed contacts",
                        xaxis={"title": "Observed contact separation (Å)", "gridcolor": GRID},
                        yaxis={"title": "Residue or residue pair", "gridcolor": GRID},
                         annotations=[{"text": "Violin width reflects recorded contact-event samples only; boxes show median and quartiles.", "showarrow": False, "x": .5, "xref": "paper", "y": 1.16, "yref": "paper", "font": {"color": MUTED}}])
    slides.append(("Residue distributions", graph_div(figure), "distance_distribution"))

    frequency = next(iter(sorted(analysis.glob("*_contact_frequency_filtered.csv"))), None)
    rank = per_residue.reset_index().rename(columns={"index": "Residue"})
    rank["MedianDistance_A"] = rank["median"]
    rank["ObservedEvents"] = rank["count"]
    rank = rank[["Residue", "MedianDistance_A", "ObservedEvents"]]
    if frequency is not None:
        try:
            occupancy = pd.read_csv(frequency, usecols=["Residue", "PercentFrames"])
            occupancy["PercentFrames"] = pd.to_numeric(occupancy["PercentFrames"], errors="coerce")
            rank = rank.merge(occupancy, on="Residue", how="left")
        except (OSError, ValueError, pd.errors.ParserError):
            pass
    if "PercentFrames" not in rank:
        rank["PercentFrames"] = np.nan
    rank = rank.sort_values(["PercentFrames", "ObservedEvents"], ascending=False, na_position="last").head(15).sort_values("MedianDistance_A")
    hover = np.column_stack([rank["ObservedEvents"], rank["PercentFrames"].fillna(-1)])
    figure = go.Figure(go.Bar(y=rank["Residue"].astype(str), x=rank["MedianDistance_A"], orientation="h", marker_color=GREEN,
                              customdata=hover,
                              hovertemplate="Contact=%{y}<br>Median separation=%{x:.3f} Å<br>Observed contact events=%{customdata[0]:.0f}<br>Geometric occupancy=%{customdata[1]:.1f}%<extra></extra>"))
    figure.update_layout(title="3 · Persistent-contact distance ranking",
                        xaxis={"title": "Median observed separation (Å)", "gridcolor": GRID},
                        yaxis={"title": "Residue or residue pair", "gridcolor": GRID},
                         annotations=[{"text": "Ranked by published geometric-contact occupancy when available; otherwise by observed event count.", "showarrow": False, "x": .5, "xref": "paper", "y": 1.16, "yref": "paper", "font": {"color": MUTED}}])
    slides.append(("Persistent-contact ranking", graph_div(figure), "distance_ranking"))

    frame_times = np.sort(events["Time_ns"].unique())
    default_step = float(np.median(np.diff(frame_times))) if len(frame_times) > 1 else 0.0
    episodes: list[tuple[str, float]] = []
    for residue, table in events.groupby("Residue"):
        times = np.sort(table["Time_ns"].unique())
        if not len(times):
            continue
        start = previous = float(times[0])
        tolerance = default_step * 1.5 + 1e-9
        for current in times[1:]:
            current = float(current)
            if default_step and current - previous <= tolerance:
                previous = current
                continue
            episodes.append((str(residue), previous - start + default_step))
            start = previous = current
        episodes.append((str(residue), previous - start + default_step))
    lifetime = pd.DataFrame(episodes, columns=["Residue", "Duration_ns"])
    if not lifetime.empty:
        lifetime = lifetime[lifetime["Duration_ns"] >= 0].sort_values("Duration_ns", ascending=False).head(30)
        figure = go.Figure(go.Bar(y=lifetime["Residue"], x=lifetime["Duration_ns"], orientation="h", marker_color=MAGENTA,
                                  hovertemplate="Contact=%{y}<br>Continuous observed contact episode=%{x:.3f} ns<extra></extra>"))
        figure.update_layout(title="4 · Longest uninterrupted observed contact episodes",
                             xaxis={"title": "Episode duration at sampled-frame resolution (ns)", "gridcolor": GRID},
                             yaxis={"title": "Residue or residue pair", "gridcolor": GRID},
                             annotations=[{"text": "Episodes break when a residue is absent from the recorded contact-event frame; this is not an inferred binding lifetime.", "showarrow": False, "x": .5, "xref": "paper", "y": 1.16, "yref": "paper", "font": {"color": MUTED}}])
        slides.append(("Contact episodes", graph_div(figure), "contact_episode"))
    if not slides:
        return ""
    rendered = "".join(
        f'<div class="distance-slide{" is-active" if index == 0 else ""}" data-distance-slide data-title="{escape(title)}">'
        f'{content}{guide_html(guide, "ⓘ Explain this graph")}</div>'
        for index, (title, content, guide) in enumerate(slides)
    )
    return f'''<style>.distance-carousel{{border-radius:21px;padding:20px;margin-top:18px;background:linear-gradient(155deg,#10172b,#060912);border:1px solid #23375d;box-shadow:0 18px 50px rgba(0,0,0,.35)}}.distance-carousel-controls{{display:flex;align-items:center;gap:10px;justify-content:center;margin:0 0 12px;flex-wrap:wrap}}.distance-carousel-controls button{{border:1px solid #285a79;border-radius:10px;background:#0a253a;color:#cfe8ff;padding:9px 12px;font:800 .82rem inherit;cursor:pointer}}.distance-carousel-controls button:hover{{background:#116887;color:#fff;border-color:#35d0ff}}.distance-carousel-controls span{{min-width:205px;text-align:center;color:#7aa6d9;font-size:.88rem;font-weight:700}}.distance-carousel-stage{{border:1px solid rgba(53,208,255,.28);border-radius:16px;background:#070b14;overflow:hidden}}.distance-slide{{display:none;margin:0;padding:7px}}.distance-slide.is-active{{display:block}}.distance-evidence{{margin:12px 3px 3px;color:#a9c1dc;font-size:.84rem;line-height:1.45}}</style><section class="distance-carousel" aria-label="Interface distance analysis">
<div class="section-head"><div><div class="eyebrow">Distance evidence · observed events only</div><h2>Interface distance analysis</h2><p>Four complementary views make separation, persistence, residue specificity, and contact continuity readable without treating a missing contact as a zero-distance measurement.</p></div></div>
<div class="distance-carousel-controls"><button type="button" data-distance-prev aria-label="Previous distance view">← Previous</button><span data-distance-status></span><button type="button" data-distance-next aria-label="Next distance view">Next →</button></div>
<div class="distance-carousel-stage">{rendered}</div><p class="distance-evidence">All summaries are derived from <code>{escape(events_path.name)}</code>{" and " + escape(timeseries.name) if timeseries.is_file() else ""}. Use the exact source files in the report downloads for full-resolution values.</p></section><script>(()=>{{const root=document.currentScript.previousElementSibling,slides=[...root.querySelectorAll('[data-distance-slide]')],status=root.querySelector('[data-distance-status]');let index=0;const render=()=>{{slides.forEach((slide,i)=>slide.classList.toggle('is-active',i===index));status.textContent=`${{index+1}} / ${{slides.length}} · ${{slides[index].dataset.title}}`;requestAnimationFrame(()=>slides[index].querySelectorAll('.js-plotly-plot').forEach(plot=>Plotly.Plots.resize(plot)))}};root.querySelector('[data-distance-prev]').addEventListener('click',()=>{{index=(index-1+slides.length)%slides.length;render()}});root.querySelector('[data-distance-next]').addEventListener('click',()=>{{index=(index+1)%slides.length;render()}});render()}})();</script>'''


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


def _event_contact_sets(events_file: Path) -> tuple[pd.DataFrame, dict[float, set[str]]]:
    """Read a bounded contact-event table and index its residue states by time."""
    table = _read_table(events_file)
    if table is None or not {"Time_ns", "Residue"}.issubset(table.columns):
        return pd.DataFrame(), {}
    table = table.copy()
    table["Time_ns"] = pd.to_numeric(table["Time_ns"], errors="coerce")
    table = table.dropna(subset=["Time_ns", "Residue"])
    if table.empty:
        return pd.DataFrame(), {}
    states = {
        float(time): set(rows["Residue"].astype(str))
        for time, rows in table.groupby("Time_ns", sort=True)
    }
    return table, states


def interaction_explorer_summary(analysis: Path) -> tuple[list[dict], float | None]:
    """Return published typed edges and the strongest observed typed occupancy.

    The interactive explorer intentionally serializes only edges that reach its
    lowest display threshold (10%).  Its empty-state message retains the
    strongest *observed* edge below that threshold, which is essential context
    for a scorecard: a geometric contact can be persistent while no individual
    residue--atom--chemistry assignment is persistent enough to display.
    """
    explorer = next(iter(sorted((analysis / "NETWORX").glob("*_interactive.html"))), None)
    if not explorer or not explorer.is_file():
        return [], None
    try:
        source = explorer.read_text(encoding="utf-8")
        match = re.search(r"const EDGES=(\[.*?\]);\s*const ", source, flags=re.DOTALL)
        edges = json.loads(match.group(1)) if match else []
        values = [float(edge.get("fraction", 0.0)) for edge in edges]
        maximum = max(values, default=0.0)
        # This literal is injected by 3B for sparse systems.  It is the
        # occupancy before 3B applies the 10% explorer display threshold.
        sparse = re.search(r"Strongest observed edge:\s*([0-9.]+)%", source)
        if sparse:
            maximum = max(maximum, float(sparse.group(1)) / 100.0)
        return edges, maximum
    except (OSError, ValueError, AttributeError):
        return [], None


def ligand_scorecard_html(analysis: Path, protein_name: str) -> str:
    """A plain-language binding summary derived only from published analyses."""
    metrics: list[tuple[str, str, str]] = []
    protein_rmsd = _read_table(analysis / "Protein_RMSD.csv")
    if protein_rmsd is not None and "Protein_RMSD" in protein_rmsd:
        values = pd.to_numeric(protein_rmsd["Protein_RMSD"], errors="coerce").dropna()
        if not values.empty:
            metrics.append(("Protein drift", f"{values.median():.2f} Å", "median backbone RMSD"))
    overlay = next(iter(sorted(analysis.glob("*_Complex_RMSD_Overlay.csv"))), None)
    if overlay and overlay.is_file():
        table = _read_table(overlay)
        if table is not None and "Ligand_RMSD" in table:
            values = pd.to_numeric(table["Ligand_RMSD"], errors="coerce").dropna()
            if not values.empty:
                metrics.append(("Ligand motion", f"{values.median():.2f} Å", "median ligand RMSD"))
    frequency = next(iter(sorted(analysis.glob("*_contact_frequency_filtered.csv"))), None)
    if frequency and frequency.is_file():
        table = _read_table(frequency)
        if table is not None and {"Residue", "PercentFrames"}.issubset(table.columns):
            values = pd.to_numeric(table["PercentFrames"], errors="coerce")
            if values.notna().any():
                leader = table.loc[values.idxmax(), "Residue"]
                metrics.append(("Most continuously proximal", f"{leader}", f"{values.max():.1f}% geometric-contact occupancy"))
                metrics.append(("Persistent proximal residues", str(int((values >= 50).sum())), "at least 50% geometric-contact occupancy"))
    # Reuse the resolved edge records published by 3B's interaction explorer.
    # They already apply its hydrogen/virtual-site remapping and exclude
    # unresolved atom labels, so this value must exactly match the largest
    # percentage a reader can see in the explorer.
    resolved_edges, strongest_typed = interaction_explorer_summary(analysis)
    if resolved_edges:
        strongest = max(resolved_edges, key=lambda edge: float(edge.get("fraction", 0.0)))
        occupancy = 100.0 * float(strongest.get("fraction", 0.0))
        metrics.append(("Strongest typed interaction", str(strongest.get("residue", "—")), f"{occupancy:.1f}% {str(strongest.get('type', 'typed interaction')).lower()} · ligand atom {str(strongest.get('atom', '—'))}"))
    elif strongest_typed is not None:
        metrics.append(("Strongest typed interaction", f"{strongest_typed * 100:.1f}%", "no exact residue–atom/type pair reaches the 10% explorer threshold"))
    all_contacts = analysis / "AllContacts_Framewise.csv"
    table = _read_table(all_contacts)
    if table is not None and "Distance" in table:
        distances = pd.to_numeric(table["Distance"], errors="coerce").dropna()
        if not distances.empty:
            metrics.append(("Closest contact", f"{distances.min():.2f} Å", "minimum observed atom-pair distance"))
    q_file = analysis / "native_contact_q.csv"
    q_table = _read_table(q_file)
    if q_table is not None and "Q_native_contact_fraction" in q_table:
        values = pd.to_numeric(q_table["Q_native_contact_fraction"], errors="coerce").dropna()
        if not values.empty:
            metrics.append(("Ligand native Q", f"{values.median():.2f}", "median native-contact fraction"))
    if not metrics:
        return ""
    cards = "".join(
        f'<div class="ligand-score"><strong>{escape(label)}</strong><b>{escape(value)}</b><span>{escape(detail)}</span></div>'
        for label, value, detail in metrics[:7]
    )
    return f'''<style>.ligand-scorecard{{display:grid;grid-template-columns:repeat(auto-fit,minmax(155px,1fr));gap:10px}}.ligand-score{{border:1px solid #27476a;border-radius:14px;padding:13px;background:linear-gradient(145deg,#0d1b2c,#080c14)}}.ligand-score strong{{display:block;color:#7aa6d9;font-size:.74rem;letter-spacing:.06em;text-transform:uppercase}}.ligand-score b{{display:block;color:#35d0ff;font-size:1.45rem;margin:5px 0;line-height:1.1}}.ligand-score span{{display:block;color:#c0d1e5;font-size:.8rem;line-height:1.35}}.scorecard-info{{grid-column:1/-1;border:1px solid rgba(53,208,255,.28);border-radius:11px;background:#091424;color:#b9cce3;padding:0 10px;font-size:.82rem;line-height:1.45}}.scorecard-info summary{{cursor:pointer;color:#35d0ff;font-weight:800;padding:9px 0;list-style:none}}.scorecard-info summary::-webkit-details-marker{{display:none}}.scorecard-info p{{margin:0 0 10px}}</style><section class="ligand-scorecard"><div style="grid-column:1/-1"><div class="eyebrow">Binding scorecard</div><p style="margin:.2rem 0 .8rem">A first-pass summary for {escape(protein_name)}. These are descriptive trajectory metrics, not affinity estimates.</p></div>{cards}<details class="scorecard-info"><summary>ⓘ Why the two interaction occupancies can differ</summary><p><b>Geometric-contact occupancy</b> means a residue is within the contact-distance cutoff of any ligand atom. <b>Typed interaction occupancy</b> is stricter: the same residue–ligand-atom pair must be classified as the same hydrogen-bond, ionic, hydrophobic, or water-bridge interaction in that sampled frame. A residue can therefore be continuously proximal without any one typed pair reaching the explorer's 10% display threshold. The typed-interaction card is taken from the same resolved interaction records displayed in the explorer below.</p></details></section>'''


def contact_state_timeline_plot(events_file: Path) -> str | None:
    table, states = _event_contact_sets(events_file)
    if table.empty or len(states) < 2:
        return None
    times = np.asarray(sorted(states), dtype=float)
    counts = np.asarray([len(states[time]) for time in times], dtype=float)
    if len(times) > MAX_RENDER_POINTS:
        keep = np.linspace(0, len(times) - 1, MAX_RENDER_POINTS, dtype=int)
        times, counts = times[keep], counts[keep]
    low, high = np.quantile(counts, [1 / 3, 2 / 3])
    state = np.where(counts <= low, 0, np.where(counts >= high, 2, 1))
    names = np.asarray(["low-contact", "moderate-contact", "high-contact"])[state]
    figure = go.Figure(go.Scatter(
        x=times, y=counts, mode="lines", line={"color": CYAN, "width": 2},
        fill="tozeroy", fillcolor="rgba(53,208,255,.14)", customdata=names,
        hovertemplate="Time=%{x:.3f} ns<br>Unique contacting residues=%{y}<br>State=%{customdata}<extra></extra>",
    ))
    figure.update_layout(title="Ligand contact-state timeline", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 52}, height=360, xaxis={"title": "Time (ns)", "gridcolor": GRID}, yaxis={"title": "Unique contacting residues", "gridcolor": GRID})
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def rmsd_distribution_plot(overlay_file: Path) -> str | None:
    table = _read_table(overlay_file)
    if table is None:
        return None
    series = [(column, color) for column, color in (("Protein_RMSD", GREEN), ("Ligand_RMSD", ORANGE)) if column in table]
    if not series:
        return None
    figure = go.Figure()
    for column, color in series:
        values = pd.to_numeric(table[column], errors="coerce").dropna()
        if values.empty:
            continue
        figure.add_trace(go.Histogram(x=values, name=column.replace("_", " "), opacity=.62, marker_color=color, nbinsx=45))
    if not figure.data:
        return None
    figure.update_layout(title="RMSD distribution across the trajectory", barmode="overlay", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 120}, height=380, xaxis={"title": "RMSD (Å)", "gridcolor": GRID}, yaxis={"title": "Sampled frames", "gridcolor": GRID}, legend=LEGEND_BELOW)
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def rmsd_rg_scatter_plot(overlay_file: Path, rg_file: Path) -> str | None:
    rmsd, rg = _read_table(overlay_file), _read_table(rg_file)
    if rmsd is None or rg is None or not {"Time_ns", "Protein_RMSD"}.issubset(rmsd.columns) or not {"Time_ns", "Rg_Protein_Å"}.issubset(rg.columns):
        return None
    left = rmsd.loc[:, ["Time_ns", "Protein_RMSD"]].copy()
    right = rg.loc[:, ["Time_ns", "Rg_Protein_Å"]].copy()
    for frame in (left, right): frame["Time_ns"] = pd.to_numeric(frame["Time_ns"], errors="coerce")
    merged = pd.merge_asof(left.dropna().sort_values("Time_ns"), right.dropna().sort_values("Time_ns"), on="Time_ns", direction="nearest")
    merged = _sample_table(merged.dropna())
    if len(merged) < 3:
        return None
    figure = go.Figure(go.Scatter(x=merged["Protein_RMSD"], y=merged["Rg_Protein_Å"], mode="markers", marker={"size": 6, "color": merged["Time_ns"], "colorscale": "Viridis", "showscale": True, "colorbar": {"title": "Time (ns)"}}, text=merged["Time_ns"], hovertemplate="Protein RMSD=%{x:.3f} Å<br>Protein Rg=%{y:.3f} Å<br>Time=%{text:.3f} ns<extra></extra>"))
    figure.update_layout(title="Protein RMSD versus compactness", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 65, "r": 25, "t": 55, "b": 52}, height=400, xaxis={"title": "Protein RMSD (Å)", "gridcolor": GRID}, yaxis={"title": "Protein radius of gyration (Å)", "gridcolor": GRID})
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def rmsf_persistence_scatter_plot(rmsf_file: Path, frequency_file: Path, interaction_file: Path | None) -> str | None:
    rmsf, frequency = _read_table(rmsf_file), _read_table(frequency_file)
    if rmsf is None or frequency is None or not {"Residue", "RMSF_Å"}.issubset(rmsf.columns) or not {"Residue", "PercentFrames"}.issubset(frequency.columns):
        return None
    # RMSF files from older runs encode residue number only (e.g. ``54``),
    # whereas contact tables use chemistry-aware labels (e.g. ``GLU54``).
    # Join on the parsed residue index while preserving the informative contact
    # label for plotting and hover text.
    rmsf = rmsf[["Residue", "RMSF_Å"]].copy()
    frequency = frequency[["Residue", "PercentFrames"]].copy()
    rmsf["Residue_Index"] = pd.to_numeric(rmsf["Residue"].astype(str).str.extract(r"(\d+)", expand=False), errors="coerce")
    frequency["Residue_Index"] = pd.to_numeric(frequency["Residue"].astype(str).str.extract(r"(\d+)", expand=False), errors="coerce")
    frequency = frequency.dropna(subset=["Residue_Index"]).drop_duplicates("Residue_Index")
    merged = rmsf.merge(frequency[["Residue_Index", "Residue", "PercentFrames"]], on="Residue_Index", how="inner", suffixes=("_RMSF", ""))
    for column in ("RMSF_Å", "PercentFrames"): merged[column] = pd.to_numeric(merged[column], errors="coerce")
    merged = merged.dropna()
    if merged.empty:
        return None
    chemistry = pd.Series("Contact", index=merged.index)
    if interaction_file and interaction_file.is_file():
        interactions = _read_table(interaction_file)
        if interactions is not None and {"Residue", "Type"}.issubset(interactions.columns):
            dominant = interactions.groupby("Residue")["Type"].agg(lambda values: values.astype(str).value_counts().index[0])
            chemistry = merged["Residue"].map(dominant).fillna("Contact")
    palette = {"H-bond": CYAN, "Ionic": MAGENTA, "Hydrophobic": ORANGE, "Contact": GREEN}
    figure = go.Figure()
    for kind in chemistry.unique():
        selected = chemistry == kind
        figure.add_trace(go.Scatter(x=merged.loc[selected, "RMSF_Å"], y=merged.loc[selected, "PercentFrames"], mode="markers+text", text=merged.loc[selected, "Residue"], textposition="top center", name=str(kind), marker={"size": 10, "color": palette.get(str(kind), GREEN)}, hovertemplate="Residue=%{text}<br>RMSF=%{x:.3f} Å<br>Contact occupancy=%{y:.1f}%<extra></extra>"))
    figure.update_layout(title="Residue flexibility versus ligand-contact persistence", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 65, "r": 25, "t": 55, "b": 120}, height=420, xaxis={"title": "Residue RMSF (Å)", "gridcolor": GRID}, yaxis={"title": "Contact occupancy (%)", "gridcolor": GRID}, legend=LEGEND_BELOW)
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def contact_lifetime_plot(events_file: Path) -> str | None:
    _table, states = _event_contact_sets(events_file)
    if len(states) < 3:
        return None
    times = np.asarray(sorted(states), dtype=float)
    step = float(np.median(np.diff(times))) if len(times) > 1 else 0.0
    if step <= 0:
        return None
    runs: dict[str, list[float]] = {}
    active: dict[str, tuple[float, float]] = {}
    for time in times:
        current = states[float(time)]
        for residue in list(active):
            if residue not in current:
                start, previous = active.pop(residue)
                runs.setdefault(residue, []).append(max(step, previous - start + step))
        for residue in current:
            active.setdefault(residue, (float(time), float(time)))
            start, _ = active[residue]
            active[residue] = (start, float(time))
    for residue, (start, previous) in active.items():
        runs.setdefault(residue, []).append(max(step, previous - start + step))
    rows = sorted(((residue, float(np.median(values)), float(max(values))) for residue, values in runs.items()), key=lambda row: row[1], reverse=True)[:25]
    if not rows:
        return None
    residues, medians, maxima = map(list, zip(*rows))
    figure = go.Figure([go.Bar(y=residues, x=medians, orientation="h", name="Median uninterrupted lifetime", marker_color=CYAN), go.Bar(y=residues, x=maxima, orientation="h", name="Longest uninterrupted lifetime", marker_color=GREEN, opacity=.58)])
    figure.update_layout(title="Ligand-contact uninterrupted lifetimes", barmode="overlay", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 105, "r": 20, "t": 55, "b": 120}, height=max(380, 24 * len(rows)), xaxis={"title": "Lifetime (ns; sampled-frame resolution)", "gridcolor": GRID}, yaxis={"title": "Residue", "gridcolor": GRID}, legend=LEGEND_BELOW)
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def contact_turnover_plot(events_file: Path) -> str | None:
    _table, states = _event_contact_sets(events_file)
    if len(states) < 2:
        return None
    times = sorted(states)
    previous: set[str] = set()
    rows = []
    for time in times:
        current = states[time]
        rows.append((time, len(current), len(current - previous), len(previous - current)))
        previous = current
    frame = pd.DataFrame(rows, columns=["Time_ns", "Active", "Gained", "Lost"])
    frame = _sample_table(frame)
    figure = go.Figure()
    for column, color in (("Active", CYAN), ("Gained", GREEN), ("Lost", MAGENTA)):
        figure.add_trace(go.Scatter(x=frame["Time_ns"], y=frame[column], mode="lines", name=column, line={"color": color, "width": 2}, stackgroup="turnover" if column != "Active" else None))
    figure.update_layout(title="Ligand-contact turnover over time", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 120}, height=380, xaxis={"title": "Time (ns)", "gridcolor": GRID}, yaxis={"title": "Residue contacts", "gridcolor": GRID}, legend=LEGEND_BELOW)
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def contact_fingerprint_html(events_file: Path, frequency_file: Path, interaction_file: Path | None) -> str:
    events, frequency = _read_table(events_file), _read_table(frequency_file)
    if events is None or frequency is None or not {"Residue", "Distance"}.issubset(events.columns) or not {"Residue", "PercentFrames"}.issubset(frequency.columns):
        return ""
    events = events.copy(); events["Distance"] = pd.to_numeric(events["Distance"], errors="coerce")
    distance = events.groupby("Residue")["Distance"].agg(["mean", "min"])
    if "LigAtomName" in events:
        atoms = events.groupby("Residue")["LigAtomName"].agg(
            lambda values: ", ".join(sorted(set(values.dropna().astype(str)))[:3]) or "—"
        )
    else:
        atoms = pd.Series(dtype=str)
    table = frequency[["Residue", "PercentFrames"]].copy().merge(distance, left_on="Residue", right_index=True, how="left").merge(atoms.rename("Ligand atoms"), left_on="Residue", right_index=True, how="left")
    table["Dominant chemistry"] = "Contact"
    if interaction_file and interaction_file.is_file():
        interactions = _read_table(interaction_file)
        if interactions is not None and {"Residue", "Type"}.issubset(interactions.columns):
            dominant = interactions.groupby("Residue")["Type"].agg(lambda values: values.astype(str).value_counts().index[0])
            table["Dominant chemistry"] = table["Residue"].map(dominant).fillna("Contact")
    table = table.sort_values("PercentFrames", ascending=False).head(30)
    if table.empty:
        return ""
    rows = "".join(f"<tr><td>{escape(str(row['Residue']))}</td><td>{float(row['PercentFrames']):.1f}%</td><td>{escape(str(row['Dominant chemistry']))}</td><td>{float(row['mean']):.2f} Å</td><td>{float(row['min']):.2f} Å</td><td>{escape(str(row['Ligand atoms']))}</td></tr>" for _, row in table.iterrows())
    return f'''<style>.fingerprint{{overflow:auto}}.fingerprint table{{width:100%;border-collapse:collapse;font-size:.86rem}}.fingerprint th{{color:#7aa6d9;text-align:left;text-transform:uppercase;font-size:.72rem;letter-spacing:.06em}}.fingerprint td,.fingerprint th{{padding:9px;border-bottom:1px solid #203653;white-space:nowrap}}.fingerprint td:first-child{{color:#35d0ff;font-weight:800}}</style><div class="fingerprint"><table><thead><tr><th>Residue</th><th>Occupancy</th><th>Dominant chemistry</th><th>Mean distance</th><th>Minimum distance</th><th>Ligand atoms</th></tr></thead><tbody>{rows}</tbody></table></div>'''


def chemistry_composition_plot(interaction_file: Path) -> str | None:
    table = _read_table(interaction_file)
    if table is None or not {"Time_ns", "Type"}.issubset(table.columns):
        return None
    table = table.copy(); table["Time_ns"] = pd.to_numeric(table["Time_ns"], errors="coerce"); table = table.dropna(subset=["Time_ns", "Type"])
    if table.empty:
        return None
    bins = np.linspace(float(table["Time_ns"].min()), float(table["Time_ns"].max()), min(120, max(20, int(np.sqrt(len(table))))) + 1)
    table["bin"] = pd.cut(table["Time_ns"], bins=bins, include_lowest=True)
    grouped = table.groupby(["bin", "Type"], observed=True).size().unstack(fill_value=0)
    x = np.asarray([interval.mid for interval in grouped.index], dtype=float)
    figure = go.Figure()
    colors = {"H-bond": CYAN, "Ionic": MAGENTA, "Hydrophobic": ORANGE, "Water bridge": GREEN}
    for kind in grouped.columns:
        figure.add_trace(go.Scatter(x=x, y=grouped[kind], mode="lines", stackgroup="chemistry", name=str(kind), line={"width": .8, "color": colors.get(str(kind), GREEN)}, hovertemplate="Time=%{x:.3f} ns<br>Events=%{y}<extra>" + escape(str(kind)) + "</extra>"))
    figure.update_layout(title="Interaction chemistry composition over time", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 60, "r": 20, "t": 55, "b": 120}, height=390, xaxis={"title": "Time (ns)", "gridcolor": GRID}, yaxis={"title": "Observed interaction events per time bin", "gridcolor": GRID}, legend=LEGEND_BELOW)
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def residue_chemistry_heatmap(interaction_file: Path) -> str | None:
    table = _read_table(interaction_file)
    if table is None or not {"Residue", "Type"}.issubset(table.columns):
        return None
    total = table.groupby("Residue").size().rename("Total")
    matrix = table.groupby(["Residue", "Type"]).size().unstack(fill_value=0).div(total, axis=0)
    matrix = matrix.loc[matrix.max(axis=1).sort_values(ascending=False).head(40).index]
    if matrix.empty:
        return None
    figure = go.Figure(go.Heatmap(z=matrix.to_numpy(), x=matrix.columns.astype(str), y=matrix.index.astype(str), colorscale="Turbo", zmin=0, zmax=1, colorbar={"title": "Fraction of events"}, hovertemplate="Residue=%{y}<br>Interaction=%{x}<br>Fraction=%{z:.1%}<extra></extra>"))
    figure.update_layout(title="Residue × interaction-chemistry occupancy", paper_bgcolor=PANEL, plot_bgcolor=PANEL, font={"color": TEXT, "family": "Inter, ui-sans-serif, system-ui"}, margin={"l": 100, "r": 25, "t": 55, "b": 52}, height=max(400, 19 * len(matrix)), xaxis={"title": "Interaction chemistry", "gridcolor": GRID}, yaxis={"title": "Protein residue", "gridcolor": GRID})
    return plot(figure, include_plotlyjs=False, output_type="div", config={"responsive": True, "displaylogo": False})


def sampled_indices(length: int, maximum: int = MAX_HEATMAP_TIMEPOINTS) -> np.ndarray:
    """Evenly select a matrix axis while retaining its beginning and end."""
    if length <= maximum:
        return np.arange(length, dtype=int)
    return np.linspace(0, length - 1, maximum, dtype=int)


def heatmap_plot(csv_file: Path, kind: str, title_override: str | None = None) -> str | None:
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
    if title_override:
        title = title_override
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


def component_lanes_section_html(root: Path, title: str, eyebrow: str, text: str,
                                 shared_cards: list[str], component_cards: dict[int, list[str]]) -> str:
    """Render a hierarchy for complex dossiers: shared first, then each chain.

    The cyan lane contains genuine whole-complex or interface observables.  A
    component plot is never silently mixed into that lane, even if the visual
    grid happens to reflow on a smaller display.
    """
    if not shared_cards and not any(component_cards.values()):
        return ""
    entries = atomindex_components(root)
    lanes = []
    if shared_cards:
        lanes.append(
            '<section class="component-analysis-lane component-analysis-lane--shared">'
            '<header><span class="component-analysis-kicker">Complex / shared</span>'
            '<h3>Whole-complex and interface analysis</h3>'
            '<p>These observables are calculated across the complete assembly or its defined interface. Cyan marks them as shared rather than chain-specific.</p>'
            '</header><div class="grid">' + "".join(shared_cards) + '</div></section>'
        )
    for component_index in sorted(component_cards):
        cards = component_cards[component_index]
        if not cards:
            continue
        index = component_index - 1
        entry = entries[index] if index < len(entries) else {}
        label = str(entry.get("label") or f"component {component_index}")
        chain = str(entry.get("current_chain") or entry.get("source_chain") or "?")
        accent = CHAIN_COLORS[index % len(CHAIN_COLORS)]
        lanes.append(
            f'<section class="component-analysis-lane" style="--component-accent:{accent}"><header>'
            f'<span class="component-analysis-kicker">Component {component_index}</span>'
            f'<h3>{escape(label)} <small>chain {escape(chain)}</small></h3>'
            '<p>Chain-local calculation: positions and values belong to this component only and never continue into another chain.</p>'
            '</header><div class="grid">' + "".join(cards) + '</div></section>'
        )
    return ('<style>.component-analysis-lanes{display:grid;gap:19px}.component-analysis-lane{--component-accent:#35d0ff;border:1px solid #285a79;border-left:5px solid var(--component-accent);border-radius:18px;padding:16px;background:linear-gradient(135deg,#0b1a28,#080b14 48%)}.component-analysis-lane--shared{--component-accent:#35d0ff}.component-analysis-lane header{padding:0 2px 13px;margin-bottom:15px;border-bottom:1px solid #23375d}.component-analysis-lane h3{margin:3px 0;color:#cfe8ff;font-size:1.25rem}.component-analysis-lane h3 small{color:#8daac7;font-size:.75em;font-weight:700}.component-analysis-lane p{margin:0;color:#7aa6d9;font-size:.88rem;max-width:820px}.component-analysis-kicker{display:block;color:var(--component-accent);font-size:.72rem;font-weight:900;letter-spacing:.14em;text-transform:uppercase}</style><section class="gallery-section component-lanes-section"><div class="section-head"><div><div class="eyebrow">'
            + escape(eyebrow) + '</div><h2>' + escape(title) + '</h2><p>' + escape(text) + '</p></div></div><div class="component-analysis-lanes">'
            + "".join(lanes) + '</div></section>')


def chain_secondary_section_html(cards_by_chain: dict[str, list[str]]) -> str:
    """Give each complex component a discrete, chain-local SSE reading lane."""
    blocks = []
    for index, (chain, cards) in enumerate(cards_by_chain.items()):
        if not cards:
            continue
        accent = CHAIN_COLORS[index % len(CHAIN_COLORS)]
        blocks.append(
            f'<section class="chain-dynamics-group" style="--chain-accent:{accent}"><header>'
            f'<span class="chain-dynamics-kicker">Component {index + 1}</span><h3>Chain {escape(chain)}</h3>'
            '<p>Read this component independently: residue positions, SSE evolution, and persistence never continue into another chain.</p>'
            '</header><div class="chain-dynamics-grid">' + "".join(cards) + "</div></section>"
        )
    if not blocks:
        return ""
    return ('<style>.chain-dynamics-stack{display:grid;gap:19px}.chain-dynamics-group{border:1px solid #285a79;border-left:5px solid var(--chain-accent);border-radius:18px;padding:16px;background:linear-gradient(135deg,#0b1a28,#080b14 48%)}.chain-dynamics-group header{padding:0 2px 13px;border-bottom:1px solid #23375d}.chain-dynamics-group h3{margin:3px 0;color:#cfe8ff;font-size:1.35rem}.chain-dynamics-group p{margin:0;color:#7aa6d9;font-size:.88rem;max-width:780px}.chain-dynamics-kicker{display:block;color:var(--chain-accent);font-size:.72rem;font-weight:900;letter-spacing:.14em;text-transform:uppercase}.chain-dynamics-grid{display:grid;grid-template-columns:repeat(3,minmax(290px,1fr));gap:15px;margin-top:15px;overflow-x:auto}.chain-dynamics-grid .card{min-width:290px}@media(max-width:700px){.chain-dynamics-grid{grid-template-columns:1fr}.chain-dynamics-grid .card{min-width:0}}</style><section class="gallery-section chain-secondary-section"><div class="section-head"><div>'
            '<div class="eyebrow">Fold evolution · chain-resolved</div><h2>Secondary structure by chain</h2>'
            '<p>Each component has its own reading lane. Compare chains only after reading their local maps and persistence profiles.</p>'
            '</div></div><div class="chain-dynamics-stack">' + "".join(blocks) + "</div></section>")


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
:root{{--bg:#05060a;--panel:#0b1020;--grid:#1c2a4a;--text:#cfe8ff;--muted:#7aa6d9;--cyan:#35d0ff;--green:#18ff8b;--yellow:#ffe66d;--magenta:#ff2bd6;--orange:#ff8a3d}}*{{box-sizing:border-box}}html{{scroll-behavior:smooth}}body{{margin:0;background:radial-gradient(circle at 90% 0,rgba(53,208,255,.15),transparent 30%),radial-gradient(circle at 4% 25%,rgba(255,43,214,.10),transparent 27%),var(--bg);color:var(--text);font-family:Inter,ui-sans-serif,system-ui,-apple-system,sans-serif}}.page{{max-width:1560px;margin:auto;padding:32px 24px 70px}}.hero,.card,.ligand-card,.gallery-section,.structure-card,.dashboard-tabs-shell{{background:linear-gradient(155deg,#10172b,#060912);border:1px solid #23375d;box-shadow:0 18px 50px rgba(0,0,0,.35)}}.hero{{position:relative;overflow:hidden;border-radius:25px;padding:34px;margin-bottom:20px}}.hero:before{{content:"";position:absolute;inset:-48%;background:conic-gradient(from 180deg,transparent,rgba(53,208,255,.16),transparent,rgba(255,43,214,.13),transparent);animation:slow-spin 28s linear infinite;pointer-events:none}}.hero>*{{position:relative;z-index:1}}@keyframes slow-spin{{to{{transform:rotate(360deg)}}}}.eyebrow{{color:var(--green);font-size:.75rem;font-weight:900;letter-spacing:.16em;text-transform:uppercase}}h1{{font-size:clamp(2.5rem,5.2vw,5.2rem);letter-spacing:-.065em;line-height:.92;margin:.35rem 0 1rem}}h2{{letter-spacing:-.025em;margin:.2rem 0 .45rem;font-size:1.5rem}}p{{line-height:1.55}}.sub,.section-head p,.muted,.asset-note{{color:var(--muted);max-width:960px}}.pills,.portal-actions{{display:flex;gap:9px;flex-wrap:wrap;margin-top:18px}}.pill,.portal-action{{border:1px solid rgba(53,208,255,.28);padding:8px 11px;border-radius:999px;background:rgba(53,208,255,.08);font-size:.88rem}}.portal-actions{{margin-top:12px}}.portal-action{{color:var(--text);font-weight:800;text-decoration:none;transition:transform .16s ease,background .16s ease}}.portal-action:hover{{transform:translateY(-1px);background:rgba(53,208,255,.20)}}.portal-action--primary{{background:linear-gradient(90deg,#18b7e6,#1988e6);border-color:transparent;color:#03101b}}details.guide,.overview{{margin:11px 0 2px;border:1px solid rgba(53,208,255,.22);border-radius:13px;background:rgba(5,10,23,.68)}}details.guide summary,.overview summary{{cursor:pointer;padding:10px 12px;color:var(--cyan);font-weight:800;font-size:.88rem}}details.guide p,.overview p{{margin:0;padding:0 12px 11px;color:var(--muted);font-size:.87rem}}.overview{{margin-top:22px;max-width:1020px}}.overview ol{{margin:0;padding:0 28px 15px;color:var(--muted);line-height:1.55}}.gallery-section,.ligand-card,.structure-card{{border-radius:21px;padding:20px;margin-top:18px}}.section-head{{border-left:3px solid var(--green);padding-left:13px;margin-bottom:16px}}.interface-pair{{margin:12px 0;border:1px solid rgba(53,208,255,.24);border-radius:15px;background:rgba(5,10,23,.58);padding:0 14px 14px}}.interface-pair summary{{cursor:pointer;padding:14px 2px 10px;color:var(--cyan);font-weight:900;font-size:1rem}}.interface-pair summary span{{color:var(--muted);font-size:.8rem;font-weight:700;margin-left:8px}}.interface-preview-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px}}.interface-preview{{margin:0;border:1px solid rgba(53,208,255,.25);border-radius:14px;overflow:hidden;background:#05070d}}.interface-preview img{{display:block;width:100%;height:250px;object-fit:contain;background:#fff}}.interface-preview svg{{display:block;width:100%;height:auto;min-height:330px}}.interface-preview figcaption{{padding:9px 11px;color:var(--muted);font-size:.84rem;font-weight:800}}.dashboard-tabs-shell{{margin-top:18px;border-radius:21px;overflow:hidden}}.dashboard-tabs{{display:flex;gap:7px;overflow:auto;padding:14px 16px 0;border-bottom:1px solid #23375d}}.dashboard-tabs button,.structure-toolbar button{{border:1px solid #285a79;border-radius:10px 10px 0 0;background:#0a253a;color:#cfe8ff;padding:10px 13px;font:800 .82rem inherit;white-space:nowrap;cursor:pointer}}.dashboard-tabs button[aria-selected="true"],.dashboard-tabs button:hover{{background:#116887;color:#fff;border-color:var(--cyan)}}.dashboard-panel{{padding:2px 20px 20px}}.dashboard-panel[hidden]{{display:none}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:16px}}.card{{border-radius:17px;padding:13px;overflow:hidden}}.card .js-plotly-plot{{cursor:zoom-in}}.card a{{display:inline-block;color:var(--cyan);font-size:.8rem;font-weight:800;margin:8px 10px 5px}}.heatmap-note{{margin:7px 10px;color:var(--muted);font-size:.82rem;line-height:1.4}}.ligand-grid{{display:grid;grid-template-columns:minmax(300px,1fr) minmax(300px,1fr);gap:18px;align-items:stretch}}.ligand-2d,.ligand-3d,.complex-3d{{min-height:330px;border:1px solid rgba(255,230,109,.30);border-radius:17px;background:radial-gradient(circle at center,rgba(255,230,109,.08),rgba(5,6,10,.7) 70%);overflow:hidden;position:relative}}.ligand-2d svg{{display:block;width:100%;height:330px;max-width:100%}}.ligand-3d,.complex-3d{{border-color:rgba(53,208,255,.30)}}.complex-3d{{height:min(62vw,650px);min-height:410px}}.structure-toolbar{{display:flex;gap:8px;flex-wrap:wrap;margin:0 0 11px}}.structure-toolbar button{{border-radius:10px;padding:8px 11px}}.structure-toolbar button:hover,.structure-toolbar button.is-active{{background:#116887;color:#fff;border-color:var(--cyan)}}.viewer-label{{color:var(--yellow);font-size:.78rem;font-weight:900;text-transform:uppercase;letter-spacing:.1em;margin:2px 0 8px}}.missing-visual{{padding:24px;color:var(--muted)}}.asset-note{{font-size:.8rem;margin:.7rem 0 0}}code{{display:block;white-space:pre-wrap;overflow-wrap:anywhere;padding:0 12px 12px;color:var(--yellow)}}.graph-modal{{display:none;position:fixed;z-index:9999;inset:0;padding:3vh 3vw;background:rgba(0,0,0,.78);backdrop-filter:blur(7px)}}.graph-modal.open{{display:flex;align-items:center;justify-content:center}}.graph-modal-panel{{position:relative;width:min(1500px,94vw);height:min(900px,92vh);padding:20px 16px 12px;border:1px solid var(--cyan);background:#080d1b;box-shadow:0 0 70px rgba(53,208,255,.18);border-radius:20px}}.graph-modal-close{{position:absolute;right:17px;top:13px;z-index:2;border:1px solid var(--cyan);border-radius:999px;background:#101b33;color:var(--text);font-weight:800;padding:8px 12px;cursor:pointer}}.graph-modal-plot{{width:100%;height:100%}}@media(max-width:700px){{.page{{padding:18px 12px 42px}}.grid,.ligand-grid{{grid-template-columns:1fr}}.hero{{padding:24px 19px}}.dashboard-panel{{padding:2px 12px 12px}}.complex-3d{{height:68vw;min-height:330px}}.graph-modal{{padding:1vh 1vw}}.graph-modal-panel{{width:98vw;height:96vh;padding:14px 4px 4px}}}}
</style></head><body><main class="page"><section class="hero"><div class="eyebrow">PyMACS · single-system results</div><h1>{escape(system_name)}</h1><p class="sub">{escape({"ligand": "Start with global stability, then inspect ligand pose and pocket chemistry.", "peptide": "Start with global stability, then inspect peptide conformation and peptide–protein interface persistence.", "protein_chain": "Start with global and chain stability, then inspect persistent chain interfaces.", "none": "Start with global stability, fold retention, and secondary-structure evolution."}[context.partner_kind])} Every graph is a native Plotly analysis: hover for values, use the legend to isolate traces, and double-click for a large view.</p><div class="pills"><span class="pill">{escape({"ligand": "Ligand / PROTAC", "peptide": "Peptide", "protein_chain": "Chains / interfaces", "none": "Protein / chain"}[context.partner_kind])}: {escape(context.partner_label)}</span><span class="pill">{len(cards)} interactive analyses</span><span class="pill">Kyle-local, fast delivery</span></div><!--PORTAL_ACTIONS--><details class="overview"><summary>ⓘ How to read this result—start here</summary><ol><li><b>Check the fold:</b> protein RMSD, radius of gyration, and secondary-structure evolution should be considered together.</li><li><b>Check the partner:</b> interpret {escape(context.relationship_label)} RMSD/RMSF and contact persistence together—not in isolation.</li><li><b>Check the chemistry:</b> solid-color interaction types show which residues make hydrophobic, ionic, or hydrogen-bond interactions over time.</li><li><b>Do not overclaim:</b> these trajectories test stability of a modeled system; they do not by themselves establish binding affinity or experimental activity.</li></ol></details></section>{dashboard}</main><div class="graph-modal" id="graph-modal" aria-hidden="true"><div class="graph-modal-panel"><button class="graph-modal-close" id="graph-modal-close" type="button">Close ×</button><div class="graph-modal-plot" id="graph-modal-plot"></div></div></div><script>
(() => {{
 const modal=document.getElementById('graph-modal'),destination=document.getElementById('graph-modal-plot');
 const close=()=>{{modal.classList.remove('open');modal.setAttribute('aria-hidden','true');Plotly.purge(destination)}};
 const open=source=>{{const layout=Object.assign({{}},source.layout,{{autosize:true,height:null,margin:{{l:70,r:35,t:65,b:60}}}});Plotly.react(destination,source.data,layout,{{responsive:true,displaylogo:false}});modal.classList.add('open');modal.setAttribute('aria-hidden','false');requestAnimationFrame(()=>Plotly.Plots.resize(destination))}};
 document.querySelectorAll('.js-plotly-plot').forEach(source=>source.addEventListener('dblclick',event=>{{event.preventDefault();event.stopImmediatePropagation();open(source)}},true));
 document.getElementById('graph-modal-close').addEventListener('click',close);modal.addEventListener('click',event=>{{if(event.target===modal)close()}});document.addEventListener('keydown',event=>{{if(event.key==='Escape')close()}});
 document.querySelectorAll('[data-dashboard-tab]').forEach(button=>button.addEventListener('click',()=>{{const target=button.dataset.dashboardTab;document.querySelectorAll('[data-dashboard-tab]').forEach(item=>item.setAttribute('aria-selected',String(item===button)));document.querySelectorAll('.dashboard-panel').forEach(panel=>{{const active=panel.id===`panel-${{target}}`;panel.hidden=!active;panel.classList.toggle('is-active',active);if(active)panel.querySelectorAll('.js-plotly-plot').forEach(plot=>Plotly.Plots.resize(plot))}})}}));
 document.querySelectorAll('.js-ligand-3d').forEach(el=>{{try{{const raw=atob(el.dataset.structureB64),viewer=$3Dmol.createViewer(el,{{backgroundColor:'#05060a',antialias:true}});viewer.addModel(raw,el.dataset.format);viewer.setStyle({{}},{{stick:{{radius:.105,color:'#ffe66d'}},sphere:{{scale:.16,color:'#ffe66d'}}}});viewer.setStyle({{elem:'N'}},{{stick:{{radius:.105,color:'#35d0ff'}},sphere:{{scale:.16,color:'#35d0ff'}}}});viewer.setStyle({{elem:'O'}},{{stick:{{radius:.105,color:'#ff3b5c'}},sphere:{{scale:.16,color:'#ff3b5c'}}}});viewer.setStyle({{elem:'Cl'}},{{stick:{{radius:.105,color:'#18ff8b'}},sphere:{{scale:.20,color:'#18ff8b'}}}});viewer.zoomTo();viewer.render();addEventListener('resize',()=>viewer.resize())}}catch(_){{el.textContent='3D ligand preview unavailable.'}}}});
 document.querySelectorAll('.js-complex-3d').forEach(el=>{{try{{const viewer=$3Dmol.createViewer(el,{{backgroundColor:'#05060a',antialias:true}}),raw=atob(el.dataset.structureB64),ligand=el.dataset.ligand,pdbLigand=el.dataset.pdbLigand||ligand,finalLigand=el.dataset.finalLigandB64?atob(el.dataset.finalLigandB64):'',chains=el.dataset.partnerChains.split(',').filter(Boolean),pdbModel=viewer.addModel(raw,'pdb'),ligandModel=finalLigand?viewer.addModel(finalLigand,'mol2'):null;const ligandStyle={{stick:{{radius:.20,colorscheme:'Jmol',singleBonds:false,aromaticStyle:'solid'}},sphere:{{scale:.25,colorscheme:'Jmol'}}}},clear=()=>{{pdbModel.setStyle({{}},{{}});if(ligandModel)ligandModel.setStyle({{}},{{}})}},paintLigand=()=>{{if(ligandModel)ligandModel.setStyle({{}},ligandStyle)}};const render=(mode,chain='')=>{{clear();pdbModel.setStyle({{hetflag:false}},{{cartoon:{{color:'spectrum'}}}});if(mode==='complex'){{pdbModel.setStyle({{hetflag:true}},{{stick:{{radius:.16,colorscheme:'Jmol'}}}});if(ligandModel)pdbModel.setStyle({{resn:pdbLigand}},{{}});paintLigand();viewer.zoomTo()}}if(mode==='chain'&&chain){{pdbModel.setStyle({{}},{{cartoon:{{color:'#34495e',opacity:.18}}}});pdbModel.setStyle({{chain}},{{cartoon:{{color:'#35d0ff'}},stick:{{radius:.10}}}});viewer.zoomTo({{chain}})}}if(mode==='partner'&&chains.length){{pdbModel.setStyle({{}},{{cartoon:{{color:'#34495e',opacity:.18}}}});pdbModel.setStyle({{chain:chains[0]}},{{cartoon:{{color:'#ffe66d'}},stick:{{radius:.10}}}});viewer.zoomTo({{chain:chains[0]}})}}if(mode==='interface'&&chains.length>1){{const sel={{chain:chains[0]}};pdbModel.setStyle({{}},{{cartoon:{{color:'#34495e',opacity:.22}}}});pdbModel.setStyle({{chain:chains[0]}},{{cartoon:{{color:'#35d0ff'}}}});pdbModel.setStyle({{within:{{distance:5,sel}}}},{{cartoon:{{color:'#ffe66d'}},stick:{{radius:.10}}}});viewer.zoomTo({{within:{{distance:7,sel}}}})}}if(mode==='pocket'&&pdbLigand){{pdbModel.setStyle({{}},{{cartoon:{{color:'#34495e',opacity:.28}}}});pdbModel.setStyle({{within:{{distance:5,sel:{{resn:pdbLigand}}}}}},{{cartoon:{{color:'#35d0ff'}},stick:{{radius:.10}}}});pdbModel.setStyle({{resn:pdbLigand}},{{}});paintLigand();viewer.zoomTo({{resn:pdbLigand}})}}if(mode==='protein')viewer.zoomTo();viewer.render();el.parentElement.querySelectorAll('[data-structure-view]').forEach(button=>button.classList.toggle('is-active',button.dataset.structureView===mode&&(mode!=='chain'||button.dataset.chainId===chain)))}};el.parentElement.querySelectorAll('[data-structure-view]').forEach(button=>button.addEventListener('click',()=>render(button.dataset.structureView,button.dataset.chainId||'')));render('complex');addEventListener('resize',()=>viewer.resize())}}catch(_){{el.textContent='3D final-structure preview unavailable.'}}}});
}})();
</script><script>
/* Independent presentation viewer for the dashboard Structure tab. */
document.querySelectorAll('.js-complex-3d').forEach(el=>{{
  try {{
    el.replaceChildren();
    const panel=el.closest('.structure-card'), raw=atob(el.dataset.structureB64), viewer=$3Dmol.createViewer(el,{{backgroundColor:'#05060a',antialias:true}}), model=viewer.addModel(raw,'pdb'), chains=el.dataset.partnerChains.split(',').filter(Boolean), proteinColor={{slate:'#496985',mint:'#4fd1a5',violet:'#a78bfa'}}, peptideColor={{gold:'#ffe66d',magenta:'#ff70d2',cyan:'#35d0ff'}}, state={{mode:'complex',chain:'',protein:'slate',peptide:'gold',pepSide:'stick',proSide:'off'}};
    const peptideSel=()=>chains.length?{{chain:chains[0]}}:{{}};
    const proteinSel=()=>chains.length?{{chain:chains[0],invert:true}}:{{hetflag:false}};
    const residueTone=atom=>{{const name=String(atom.resn||'').toUpperCase();if(['ASP','GLU'].includes(name))return '#ff5a78';if(['LYS','ARG','HIS'].includes(name))return '#5bc0ff';if(['SER','THR','ASN','GLN','TYR','CYS'].includes(name))return '#c084fc';if(['ALA','VAL','ILE','LEU','MET','PHE','TRP','PRO','GLY'].includes(name))return '#ffcf5b';return '#e6edf7'}};
    const sideStyle=(kind,color)=>{{const paint=color==='residue'?{{colorfunc:residueTone}}:{{color}};return kind==='stick'?{{stick:{{radius:.16,...paint}}}}:kind==='line'?{{line:{{linewidth:1.3,...paint}}}}:{{}}}};
    const render=()=>{{
      model.setStyle({{}},{{cartoon:{{color:state.protein==='spectrum'?'spectrum':proteinColor[state.protein]||undefined,colorscheme:state.protein==='chain'?'chain':undefined}}}});
      if(chains.length) {{ model.setStyle(peptideSel(),{{cartoon:{{color:state.peptide==='residue'?undefined:peptideColor[state.peptide],colorfunc:state.peptide==='residue'?residueTone:undefined}},...sideStyle(state.pepSide,state.peptide)}}); }}
      if(state.proSide!=='off') model.setStyle(proteinSel(),sideStyle(state.proSide,'#8fb7d8'));
      if(state.mode==='chain'&&state.chain) {{model.setStyle({{}},{{cartoon:{{color:'#24364c',opacity:.25}}}});model.setStyle({{chain:state.chain}},{{cartoon:{{color:'#35d0ff'}},...sideStyle(state.proSide,'#35d0ff')}});viewer.zoomTo({{chain:state.chain}})}}
      else if(state.mode==='partner'&&chains.length) viewer.zoomTo(peptideSel());
      else if(state.mode==='interface'&&chains.length) viewer.zoomTo({{within:{{distance:7,sel:peptideSel()}}}});
      else viewer.zoomTo();
      viewer.render(); panel.querySelectorAll('[data-structure-view]').forEach(b=>b.classList.toggle('is-active',b.dataset.structureView===state.mode&&(state.mode!=='chain'||b.dataset.chainId===state.chain)));
    }};
    const residueOption=panel.querySelector('[data-structure-control="peptide-color"] option[value="element"]');if(residueOption){{residueOption.value='residue';residueOption.textContent='By residue type'}};
    panel.querySelectorAll('[data-structure-view]').forEach(b=>b.addEventListener('click',()=>{{state.mode=b.dataset.structureView;state.chain=b.dataset.chainId||'';render()}}));
    panel.querySelectorAll('[data-structure-control]').forEach(control=>control.addEventListener('change',()=>{{const key=control.dataset.structureControl;if(key==='protein-color')state.protein=control.value;if(key==='peptide-color')state.peptide=control.value;if(key==='peptide-sidechains')state.pepSide=control.value;if(key==='protein-sidechains')state.proSide=control.value;render()}}));
    render(); addEventListener('resize',()=>viewer.resize());
  }} catch (_) {{ el.textContent='3D structure preview unavailable.'; }}
}});
document.querySelectorAll('.structure-card').forEach(card=>{{
  const controls=card.querySelector('.structure-control-grid'); if(!controls)return;
  Object.assign(controls.style,{{display:'grid',gridTemplateColumns:'repeat(auto-fit,minmax(175px,1fr))',gap:'10px',margin:'0 0 13px',padding:'12px',border:'1px solid rgba(53,208,255,.24)',borderRadius:'14px',background:'linear-gradient(135deg,#09192b,#10122b)'}});
  controls.querySelectorAll('label').forEach(label=>{{Object.assign(label.style,{{display:'grid',gap:'6px',color:'#9bbce0',fontSize:'.76rem',fontWeight:'800'}});const select=label.querySelector('select');Object.assign(select.style,{{appearance:'none',padding:'9px 32px 9px 11px',borderRadius:'9px',border:'1px solid #285a79',background:'#07111e',color:'#dff5ff',fontWeight:'800',cursor:'pointer'}})}});
}});
document.querySelectorAll('.interface-carousel').forEach(carousel=>{{
  const slides=[...carousel.querySelectorAll('.interface-slide')], label=carousel.querySelector('[data-interface-carousel-label]'); let index=0;
  Object.assign(carousel.style,{{border:'1px solid rgba(53,208,255,.3)',borderRadius:'16px',overflow:'hidden',background:'#07111e'}});
  carousel.querySelectorAll('.interface-slide').forEach(slide=>slide.style.padding='14px');
  carousel.querySelector('.interface-carousel-controls').style.cssText='display:flex;justify-content:center;align-items:center;gap:14px;padding:12px;border-top:1px solid #23375d;color:#b7d5ed;font-weight:800';
  carousel.querySelectorAll('[data-interface-carousel]').forEach(button=>{{button.style.cssText='border:1px solid #285a79;border-radius:9px;background:#0a253a;color:#dff5ff;padding:8px 12px;font:800 .82rem inherit;cursor:pointer';button.addEventListener('click',()=>{{index=(index+(button.dataset.interfaceCarousel==='next'?1:-1)+slides.length)%slides.length;slides.forEach((slide,i)=>slide.hidden=i!==index);label.textContent=`${{index+1}} / ${{slides.length}}`}})}});
  const download=carousel.querySelector('[data-interface-svg-download]'); if(download){{
    download.style.cssText='border:1px solid #18ff8b;border-radius:9px;background:#08291d;color:#d9fff0;padding:8px 12px;font:800 .82rem inherit;cursor:pointer';
    download.addEventListener('click',()=>{{
      const svg=carousel.querySelector('.interface-slide:not([hidden]) svg'); if(!svg)return;
      const copy=svg.cloneNode(true);
      copy.querySelectorAll('[data-interface-background],[data-interface-export-remove]').forEach(node=>node.remove());
      copy.querySelectorAll('[data-interface-export-title]').forEach(node=>node.setAttribute('fill','#10233b'));
      copy.setAttribute('xmlns','http://www.w3.org/2000/svg');
      const blob=new Blob([new XMLSerializer().serializeToString(copy)],{{type:'image/svg+xml;charset=utf-8'}});
      const link=document.createElement('a'); link.href=URL.createObjectURL(blob); link.download=`pymacs-interface-${{carousel.closest('[data-interface-threshold-set]')?.dataset.interfaceThresholdSet?.replace('%','pct')||'map'}}-${{index+1}}.svg`; link.click(); URL.revokeObjectURL(link.href);
    }});
  }}
}});
document.querySelectorAll('.interface-threshold-controls').forEach(controls=>{{
  Object.assign(controls.style,{{display:'flex',gap:'8px',flexWrap:'wrap',margin:'0 0 12px'}});
  const section=controls.closest('.gallery-section'); const buttons=[...controls.querySelectorAll('[data-interface-threshold]')];
  const palette={{'50%':{{color:'#ff6adf',glow:'rgba(255,106,223,.34)',panel:'#32112d'}},'75%':{{color:'#35d0ff',glow:'rgba(53,208,255,.34)',panel:'#0b2940'}},'90%':{{color:'#18ff8b',glow:'rgba(24,255,139,.34)',panel:'#08291d'}}}};
  const activate=threshold=>{{
    section.querySelectorAll('[data-interface-threshold-set]').forEach(deck=>deck.hidden=deck.dataset.interfaceThresholdSet!==threshold);
    buttons.forEach(button=>{{const active=button.dataset.interfaceThreshold===threshold,theme=palette[button.dataset.interfaceThreshold];button.setAttribute('aria-pressed',active);Object.assign(button.style,{{border:`1px solid ${{theme.color}}`,borderRadius:'999px',background:active?theme.panel:'#0a1928',boxShadow:active?`0 0 16px ${{theme.glow}}`:'none',color:active?theme.color:'#9bbce0',padding:'8px 13px',font:'900 .82rem inherit',letterSpacing:'.03em',cursor:'pointer',transition:'background .16s,box-shadow .16s,color .16s'}})}});
  }};
  buttons.forEach(button=>button.addEventListener('click',()=>activate(button.dataset.interfaceThreshold)));
  activate(buttons.find(button=>button.getAttribute('aria-pressed')==='true')?.dataset.interfaceThreshold||buttons[0]?.dataset.interfaceThreshold);
}});
(() => {{
  const hud=document.createElement('div'); hud.id='interface-map-hud'; hud.setAttribute('role','status'); document.body.append(hud);
  Object.assign(hud.style,{{position:'fixed',zIndex:'10001',maxWidth:'390px',maxHeight:'270px',overflowY:'auto',padding:'11px 13px',border:'1px solid #35d0ff',borderRadius:'10px',background:'linear-gradient(135deg,rgba(4,12,27,.97),rgba(8,34,52,.97))',boxShadow:'0 0 26px rgba(53,208,255,.34),inset 0 0 18px rgba(53,208,255,.08)',color:'#dff5ff',font:'700 12px/1.35 ui-monospace,SFMono-Regular,Menlo,monospace',whiteSpace:'pre-line',pointerEvents:'none',opacity:'0',transform:'translateY(4px)',transition:'opacity 55ms linear,transform 55ms linear'}});
  const themes={{one:{{accent:'#ff8a3d',glow:'rgba(255,138,61,.38)',background:'linear-gradient(135deg,rgba(40,18,7,.98),rgba(70,31,10,.98))'}},two:{{accent:'#18ff8b',glow:'rgba(24,255,139,.36)',background:'linear-gradient(135deg,rgba(4,34,20,.98),rgba(7,60,32,.98))'}},edge:{{accent:'#35d0ff',glow:'rgba(53,208,255,.34)',background:'linear-gradient(135deg,rgba(4,12,27,.98),rgba(8,34,52,.98))'}}}};
  const paint=target=>{{const theme=target.dataset.interfaceTooltipKind==='edge'?themes.edge:(themes[target.dataset.interfaceComponent]||themes.edge);hud.style.borderColor=theme.accent;hud.style.boxShadow=`0 0 26px ${{theme.glow}},inset 0 0 18px ${{theme.glow}}`;hud.style.background=theme.background}};
  const place=event=>{{const x=Math.min(innerWidth-405,event.clientX+16),y=Math.min(innerHeight-285,event.clientY+16);hud.style.left=`${{Math.max(10,x)}}px`;hud.style.top=`${{Math.max(10,y)}}px`}};
  document.querySelectorAll('[data-interface-tooltip]').forEach(target=>{{
    const show=event=>{{paint(target);hud.textContent=target.dataset.interfaceTooltip;place(event);hud.style.opacity='1';hud.style.transform='translateY(0)'}};
    target.addEventListener('pointerenter',show);target.addEventListener('pointermove',place);target.addEventListener('pointerleave',()=>{{hud.style.opacity='0';hud.style.transform='translateY(4px)'}});
    target.addEventListener('focus',()=>{{paint(target);hud.textContent=target.dataset.interfaceTooltip;hud.style.left='18px';hud.style.top='18px';hud.style.opacity='1';hud.style.transform='translateY(0)'}});target.addEventListener('blur',()=>{{hud.style.opacity='0';hud.style.transform='translateY(4px)'}});
  }});
}})();
</script><script>
/* Forward chain-aware explorer selections between independently sandboxed panels. */
addEventListener('message',event=>{{const data=event.data||{{}};if(data.type!=='pymacs-rin-selection'||!data.residue)return;document.querySelectorAll('iframe[title*="interface explorer"]').forEach(frame=>{{if(frame.contentWindow!==event.source)frame.contentWindow?.postMessage(data,'*')}})}});
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
    context, has_ligand_evidence = promote_ligand_context_when_evidenced(
        root, analysis, setup, metadata, context,
    )
    ligand_name = context.partner_label if context.partner_kind == "ligand" else "No ligand label"
    chain_labels = []
    for entry in str(setup.get("chain_map") or "").split(","):
        chain, separator, label = entry.strip().partition(":")
        if chain and separator and label:
            chain_labels.append(label.strip())
    # A complex must never present a pooled residue axis as if chains were a
    # continuous sequence.  Metric files are a useful fallback for historical
    # runs whose metadata did not preserve a chain map.
    is_multichain_system = (
        context.partner_kind == "protein_chain"
        or len(chain_labels) > 1
        or len(list(analysis.glob("RMSD_*.csv"))) > 1
    )
    protein_name = (
        " / ".join(component_display_labels(root, context))
        or " / ".join(chain_labels)
        or legacy_receptor_label(root)
        or "Protein"
    )
    # Full-complex RMSD is meaningful for multiple macromolecular chains or
    # a PROTAC assembly.  A small-molecule ligand is not an extra protein
    # chain, so a one-chain ligand–protein job must not show this card.
    is_multichain_system = is_multichain_system or str(setup.get("mode") or "").casefold() == "protac"
    csv_files = sorted(path for path in analysis.rglob("*.csv") if path.is_file())
    # The context resolver above promotes legacy records when their completed
    # small-molecule evidence bundle proves that they are ligand systems.
    has_ligand_evidence = has_ligand_evidence or context.partner_kind == "ligand"
    # Every card below is explicitly mapped to a scientific table.  Never
    # revive the legacy arbitrary-column fallback for raw analysis tables.
    gallery_csvs: set[Path] = set()
    colours = [CYAN, GREEN, MAGENTA, ORANGE]
    heatmap_cards = []
    secondary_heatmap_cards = []
    chain_secondary_cards: dict[str, list[str]] = {}
    groups: dict[str, list[str]] = {"overview": [], "stability": [], "flexibility": [], "secondary": [], "contacts": [], "chemistry": [], "distances": [], "coupled": [], "structure": []}
    component_cards: dict[str, dict[int, list[str]]] = {group: {} for group in groups}
    def add(group: str, graph: str | None, title: str, source: Path, label: str = "Download exact source CSV", note: str = "", component: int | None = None) -> None:
        if graph:
            card = dashboard_card(graph, title, source.relative_to(output.parent).as_posix(), label, note)
            if is_multichain_system and component is not None:
                component_cards[group].setdefault(component, []).append(card)
            else:
                groups[group].append(card)
    # Matrix outputs are promoted to dedicated scientific heatmaps.  They are
    # deliberately sampled along time so the report stays quick to publish and
    # open; the full matrix CSV remains downloadable from the card.
    dssp_file = analysis / "DSSP_encoded_numeric.csv"
    if dssp_file.is_file() and not is_multichain_system:
        graph = heatmap_plot(dssp_file, "dssp")
        if graph:
            secondary_heatmap_cards.append(dashboard_card(graph, dssp_file.stem, dssp_file.relative_to(output.parent).as_posix(), "Download full SSE matrix CSV"))
    # Protein-complex chains have independent residue numbering.  Publish a
    # heatmap for each chain rather than letting a single pooled matrix imply
    # that the final residue of one component neighbours the first of another.
    for source in sorted(analysis.glob("DSSP_chain_*.csv")):
        chain_label = source.stem.removeprefix("DSSP_chain_")
        graph = heatmap_plot(source, "dssp", f"Chain {chain_label} secondary-structure evolution")
        if graph:
            card = dashboard_card(
                graph, f"Chain {chain_label} SSE evolution",
                source.relative_to(output.parent).as_posix(), "Download chain SSE matrix CSV",
                "This map is chain-local: residue positions do not continue into another protein chain.",
            )
            if is_multichain_system:
                chain_secondary_cards.setdefault(chain_label, []).append(card)
            else:
                secondary_heatmap_cards.append(card)
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
    # Ligand reports start with a compact decision aid, then retain every
    # underlying raw CSV as a download from its dedicated visual card.
    if has_ligand_evidence:
        scorecard = ligand_scorecard_html(analysis, protein_name)
        if scorecard:
            source = analysis / "Protein_RMSD.csv"
            groups["overview"].append(dashboard_card(scorecard, "Binding scorecard", source.relative_to(output.parent).as_posix() if source.is_file() else "", note="Read this as a trajectory summary, then use the stability, contacts, chemistry, and distance tabs for the supporting evidence."))
    if context.partner_kind == "peptide":
        interface_scorecard = interface_scorecard_html(analysis)
        if interface_scorecard:
            groups["overview"].append(interface_scorecard)
        representative_frames = representative_interface_frames_html(analysis)
        if representative_frames:
            groups["overview"].append(representative_frames)
        interface_source = analysis / "Interface_RIN" / "interface_node_summary.csv"
        importance = peptide_importance_plot(analysis)
        if importance and interface_source.is_file():
            add("contacts", importance, "Peptide residue importance ranking", interface_source,
                note="ELI5: taller bars mean this peptide position repeatedly reaches more protein partners or reaches them more persistently. This ranks simulated interface roles; it is not an experimental mutation result.")
        regime_source = analysis / "Interface_RIN" / "interface_residue_pair_presence.npz"
        regime = peptide_contact_regime_plot(analysis)
        if regime and regime_source.is_file():
            add("contacts", regime, "Peptide contact timeline", regime_source,
                label="Download exact contact-presence matrix",
                note="ELI5: each bright cell means that peptide position touched at least one protein residue in that sampled frame. Dark means no recorded contact in that frame, not a zero-distance measurement.")
    for contact_file in sorted(analysis.glob("*_contact_map_residue_time.csv")):
        graph = heatmap_plot(contact_file, "contacts")
        if graph:
            heatmap_cards.append(dashboard_card(graph, contact_file.stem, contact_file.relative_to(output.parent).as_posix(), "Download full contact matrix CSV"))

    # The former portal tabs are all rendered here from their scientific CSVs.
    # Each card selects the semantic dimensions of its source rather than
    # joining arbitrary numeric columns into a generic line chart.
    stability_specs = [
        ("Protein_RMSD.csv", "Whole-complex RMSD", "Time_ps", ["Protein_RMSD"], "RMSD (Å)"),
        ("Radius_of_Gyration_Protein.csv", "Whole-complex radius of gyration", "Time_ns", ["Rg_Protein_Å"], "Radius of gyration (Å)"),
        ("LIG_Complex_RMSD_Overlay.csv", "Protein–ligand RMSD overlay", "Time_ns", ["Protein_RMSD", "Ligand_RMSD"], "RMSD (Å)"),
        ("Radius_of_Gyration_Protein_Ligand_Complex.csv", "Radius of gyration", "Time_ns", ["Rg_Protein_Å", "Rg_Ligand_Å", "Rg_Complex_Å"], "Radius of gyration (Å)"),
        ("native_contact_q.csv", f"{context.relationship_label.title()} Q", "Time_ns", ["Q_native_contact_fraction"], "Contact fraction"),
        ("protein_native_contact_q.csv", "Protein fold/native-contact Q", "Time_ns", ["Q_protein_native_contact_fraction"], "Contact fraction"),
    ]
    if is_multichain_system:
        stability_specs.insert(1, ("FullComplex_RMSD.csv", "Full-complex RMSD", "Time_ps", ["RMSD_Å"], "RMSD (Å)"))
    for file_name, title, x, ys, y_title in stability_specs:
        source = analysis / file_name
        # Whole-complex observables are neutral cyan: reserve Miami orange and
        # neon green exclusively for the individually resolved components.
        colors = [CYAN, GREEN, MAGENTA] if "overlay" in title.lower() or "gyration" in title.lower() else [CYAN]
        if source.is_file(): add("stability", multiseries_plot(source, title.replace("Protein", protein_name), x, ys, y_title, colors=colors), title, source)
    if is_multichain_system:
        component_rg_sources = sorted(
            source for source in analysis.glob("Radius_of_Gyration_*.csv")
            if source.name not in {"Radius_of_Gyration_Protein.csv", "Radius_of_Gyration_Protein_Ligand_Complex.csv"}
        )
        for component_index, source in enumerate(component_rg_sources):
            label = source.stem.removeprefix("Radius_of_Gyration_")
            add("stability", multiseries_plot(source, f"Component {component_index + 1}: {label} radius of gyration", "Time_ns", ["Rg_Component_Å"], "Radius of gyration (Å)", colors=[CHAIN_COLORS[component_index % len(CHAIN_COLORS)]]),
                f"Component {component_index + 1}: {label} radius of gyration", source,
                note="All atoms in this component are used for Rg. Compare each independent component with the whole-complex Rg rather than treating it as a backbone RMSD measure.", component=component_index + 1)
    protein_rmsd = analysis / "Protein_RMSD.csv"
    if protein_rmsd.is_file():
        add("overview", multiseries_plot(protein_rmsd, f"{protein_name} RMSD", "Time_ps", ["Protein_RMSD"], "RMSD (Å)", colors=[GREEN]), f"{protein_name} RMSD", protein_rmsd)
    if interaction_file.is_file():
        add("overview", heatmap_plot(interaction_file, "interactions"), "Interaction chemistry overview", interaction_file, "Download full interaction-event CSV")
    chain_rmsd_sources = sorted(
        source for source in analysis.glob("RMSD_*.csv")
        if source.name not in {"RMSD_Protein.csv"}
    )
    if chain_rmsd_sources:
        add("stability", chain_rmsd_plot(chain_rmsd_sources), "Chain-resolved RMSD", chain_rmsd_sources[0],
            note="Each trace is measured after alignment on that component's declared anchor selection (Cα for proteins). Colors identify components: first = Miami orange, second = neon Miami green.")
        add("overview", chain_rmsd_plot(chain_rmsd_sources), "Chain-resolved RMSD", chain_rmsd_sources[0],
            note="Start with this: each color is one component, fitted and measured on its own Cα anchors.")
    for source in sorted(analysis.glob("*_Ligand_RMSD.csv")):
        partner_rmsd = "Peptide RMSD" if context.partner_kind == "peptide" else "Ligand RMSD"
        add("stability", multiseries_plot(source, partner_rmsd, "Time_ps", ["Ligand_RMSD"], "RMSD (Å)", colors=[ORANGE]), partner_rmsd, source)
    peptide_rmsd = analysis / "Peptide_RMSD.csv"
    if peptide_rmsd.is_file():
        add("stability", multiseries_plot(peptide_rmsd, f"{context.partner_label} peptide RMSD", "Time_ps", ["RMSD_A"], "RMSD (Å)", colors=[ORANGE]), "Peptide RMSD", peptide_rmsd,
            note="Peptide coordinates are aligned on the peptide selection. Interpret this with the interface persistence panel rather than as a binding affinity measure.")

    if has_ligand_evidence:
        overlay_source = next(iter(sorted(analysis.glob("*_Complex_RMSD_Overlay.csv"))), None)
        rg_source = analysis / "Radius_of_Gyration_Protein_Ligand_Complex.csv"
        if overlay_source:
            add("stability", rmsd_distribution_plot(overlay_source), "RMSD distribution", overlay_source, note="The distribution shows how broadly the ligand and protein sampled structural deviation across the entire trajectory.")
            if rg_source.is_file():
                add("stability", rmsd_rg_scatter_plot(overlay_source, rg_source), "Protein RMSD versus compactness", rg_source, note="Points are sampled frames colored by time; coherent clusters indicate recurring structural regimes.")

    # RMSF_Protein_per_residue.csv is a legacy pooled trace.  It must never
    # claim a component palette slot or be presented beside chain-local data.
    chain_rmsf_sources = sorted(
        source for source in analysis.glob("RMSF_*_per_residue.csv")
        if source.name != "RMSF_Protein_per_residue.csv"
    )
    if chain_rmsf_sources:
        for component_index, source in enumerate(chain_rmsf_sources):
            label = _chain_label_from_metric(source, "RMSF_", "_per_residue.csv")
            add("flexibility", chain_rmsf_plot([source], title=f"Component {component_index + 1}: {label} residue RMSF", colors=[CHAIN_COLORS[component_index % len(CHAIN_COLORS)]]),
                f"Component {component_index + 1}: {label} residue RMSF", source,
                note="This component is plotted independently with its own chain-local residue positions. Hover reports the original residue ID and anchor atom used for RMSF.", component=component_index + 1)
        add("overview", chain_rmsf_plot(chain_rmsf_sources), "Chain-resolved residue RMSF", chain_rmsf_sources[0],
            note="Comparison only: each color is an independent component beginning at local position 1; hover preserves the source residue identifier.")
    for source in sorted(analysis.glob("RMSF_*.csv")):
        if source.name.endswith("_per_residue.csv"): continue
        add("flexibility", categorical_bar_plot(source, f"{protein_name} backbone-atom RMSF", "Atom_Index", ["RMSF_Å"], "RMSF (Å)", colors=[GREEN]), f"{protein_name} backbone-atom RMSF", source, note="This is the atom-resolved backbone selection used by the analysis; use the adjacent residue plot for residue-level interpretation.")
    for source in sorted(analysis.glob("*_Ligand_RMSF.csv")):
        partner_rmsf = "Peptide residue RMSF" if context.partner_kind == "peptide" else "Ligand RMSF"
        add("flexibility", categorical_bar_plot(source, partner_rmsf, "Atom", ["RMSF"], "RMSF (Å)", colors=[ORANGE]), partner_rmsf, source)
    peptide_rmsf = analysis / "Peptide_RMSF.csv"
    if peptide_rmsf.is_file():
        add("flexibility", categorical_bar_plot(peptide_rmsf, f"{context.partner_label} peptide residue RMSF", "Residue", ["RMSF_A"], "RMSF (Å)", colors=[ORANGE]), "Peptide residue RMSF", peptide_rmsf,
            note="Higher values locate mobile peptide residues; compare contacting positions with the sequence-preserving peptide map. Terminal residues often fluctuate more because they have fewer structural constraints, so interpret an isolated terminal spike in sequence context.")

    frequency_for_ligand = next(iter(sorted(analysis.glob("*_contact_frequency_filtered.csv"))), None)
    if has_ligand_evidence and frequency_for_ligand:
        rmsf_residue = next(iter(sorted(analysis.glob("RMSF_*_per_residue.csv"))), None)
        if rmsf_residue:
            add("flexibility", rmsf_persistence_scatter_plot(rmsf_residue, frequency_for_ligand, interaction_file if interaction_file.is_file() else None), "Residue flexibility versus ligand-contact persistence", frequency_for_ligand, note="Stable anchors tend toward high contact occupancy with lower residue RMSF; this is descriptive, not an affinity score.")

    sse_fractions = analysis / "SSE_Fractions.csv"
    if sse_fractions.is_file() and not is_multichain_system:
        add("secondary", multiseries_plot(sse_fractions, "SSE fractions", "Time_ns", ["Helix", "Sheet", "Coil"], "Fraction"), "SSE fractions", sse_fractions)
    persistence = analysis / "Residue_SSE_Persistence.csv"
    if persistence.is_file() and not is_multichain_system:
        add("secondary", categorical_bar_plot(persistence, "Residue SSE persistence", "Residue_Index", ["Helix_Persistence", "Sheet_Persistence", "Coil_Persistence"], "Fraction", stacked=True), "Residue SSE persistence", persistence)
    for source in sorted(analysis.glob("SSE_chain_*_fractions.csv")):
        chain_label = source.stem.removeprefix("SSE_chain_").removesuffix("_fractions")
        graph = multiseries_plot(source, f"Chain {chain_label} SSE fractions", "Time_ns", ["Helix", "Sheet", "Coil"], "Fraction")
        if is_multichain_system and graph:
            chain_secondary_cards.setdefault(chain_label, []).append(dashboard_card(graph, f"Chain {chain_label} SSE fractions", source.relative_to(output.parent).as_posix()))
        else:
            add("secondary", graph, source.stem.replace("_", " "), source)
    for source in sorted(analysis.glob("Residue_SSE_Persistence_chain_*.csv")):
        chain_label = source.stem.removeprefix("Residue_SSE_Persistence_chain_")
        graph = chain_sse_persistence_plot(source, chain_label)
        if is_multichain_system and graph:
            chain_secondary_cards.setdefault(chain_label, []).append(dashboard_card(
                graph, f"Chain {chain_label} residue SSE persistence", source.relative_to(output.parent).as_posix(),
                note="Bars are chain-local. Hover reports the deposited residue identifier, so no numbering is merged across chains."))
        else:
            add("secondary", graph, f"Chain {chain_label} residue SSE persistence", source)

    frequency = frequency_for_ligand or (analysis / "LIG_contact_frequency_filtered.csv")
    if frequency.is_file(): add("contacts", categorical_bar_plot(frequency, "Ligand-contact frequency", "Residue", ["PercentFrames"], "Frames with contact (%)"), "Ligand-contact frequency", frequency)
    filtered = analysis / "FilteredContacts_Framewise.csv"
    all_contacts = analysis / "AllContacts_Framewise.csv"
    persistence_graph = contact_persistence_plot(filtered, all_contacts, interaction_file if interaction_file.is_file() else None, protein_name) if filtered.is_file() and all_contacts.is_file() else None
    if persistence_graph: add("contacts", persistence_graph, "Ligand-contact persistence", filtered, note="Persistence is recomputed from sampled frames (0–1) using the original residue labels. A bar color now means the most frequently observed classified chemistry for that residue: H-bond blue, ionic pink, hydrophobic orange.")
    if filtered.is_file(): add("contacts", contact_events_plot(filtered, "Specific ligand contact events"), "Specific ligand contact events", filtered, note="Hover a marker for the contacting ligand atom, protein atom, distance, residue, and time.")
    if has_ligand_evidence and filtered.is_file():
        add("overview", contact_state_timeline_plot(filtered), "Ligand contact-state timeline", filtered, note="The line counts unique contacting residues at each sampled time. Low, moderate, and high states are relative to this trajectory's observed range.")
        add("contacts", contact_lifetime_plot(filtered), "Ligand-contact uninterrupted lifetimes", filtered, note="Consecutive observed frames are grouped into uninterrupted contact episodes at the sampled-frame resolution.")
        add("contacts", contact_turnover_plot(filtered), "Ligand-contact turnover", filtered, note="Shows active, newly gained, and lost contacting residues per sampled frame.")
        if frequency_for_ligand:
            fingerprint = contact_fingerprint_html(filtered, frequency_for_ligand, interaction_file if interaction_file.is_file() else None)
            if fingerprint:
                groups["contacts"].append(dashboard_card(fingerprint, "Residue contact fingerprint", filtered.relative_to(output.parent).as_posix(), note="Ranked residue summary: occupancy, dominant chemistry, distances, and the ligand atoms observed in the contact events."))

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
        if has_ligand_evidence:
            add("chemistry", chemistry_composition_plot(interaction_file), "Interaction chemistry composition over time", interaction_file, note="Counts are binned event observations, so use the persistent interaction explorer for exact residue–atom occupancy.")
            add("chemistry", residue_chemistry_heatmap(interaction_file), "Residue × interaction chemistry occupancy", interaction_file, note="Each row is normalized within that residue's observed interaction events.")

    # Distance evidence is intentionally consolidated in the four-view
    # carousel below.  Keeping one evidence-first reading path prevents the
    # same contact table from becoming a noisy wall of near-duplicate plots.

    coupling = analysis / "LIG_SSE_Contact_Coupling_Timeseries.csv"
    if coupling.is_file():
        add("coupled", multiseries_plot(coupling, "SSE versus ligand-contact dynamics", "Time_ns", ["Helix_Fraction", "Sheet_Fraction", "Coil_Fraction", "Raw_Contacts"], "Fraction / contact count"), "SSE versus ligand-contact dynamics", coupling)
        add("coupled", multiseries_plot(coupling, "Smoothed contacts", "Time_ns", ["Raw_Contacts", "Smoothed_Contacts"], "Contact count"), "Smoothed contacts", coupling)
        add("coupled", multiseries_plot(coupling, "Rolling correlation", "Time_ns", ["Rolling_Corr_Contacts_vs_Sheet"], "Correlation"), "Rolling correlation", coupling)
    changepoints = analysis / "LIG_ChangePoint_Summary.csv"
    if changepoints.is_file(): add("coupled", categorical_bar_plot(changepoints, "Change-point segments", "Segment", ["MeanContacts", "MeanHelix", "MeanSheet", "MeanCoil"], "Mean value"), "Change-point segments", changepoints, note="Hover the bars to compare the average contact and secondary-structure values for each detected trajectory segment.")
    # Protein complexes preserve a different, richer coupling bundle than
    # peptide/ligand jobs.  These graphs deliberately reuse saved analysis
    # products rather than requiring an MD rerun: interface geometry, native
    # fold retention, component Rg/RMSD, validated chemistry, and exact
    # residue-pair framewise presence.
    if is_multichain_system:
        interface_time, interface_time_source = _interface_time_data(analysis)
        q_source = analysis / "protein_native_contact_q.csv"
        graph = coupled_engagement_plot(analysis)
        if graph and q_source.is_file():
            add("coupled", graph, "Shared interface engagement and fold retention", interface_time_source,
                note="ELI5: read this as the complex's handshake monitor. Contacts show how many residue pairs touch, the gap shows how close the proteins get, and Q shows how much of the original protein fold remains. These are shared cyan measurements, not a single-chain result.")
        graph = coupled_motion_binding_plot(analysis)
        if graph:
            add("coupled", graph, "Component motion versus interface engagement", interface_time_source,
                note="ELI5: each dot is one sampled moment. If a component's orange or green cloud shifts as the cyan interface-contact count changes, its shape or motion is coupled to binding; a diffuse cloud means there is no simple one-to-one relationship.")
        chemistry_source = analysis / "Interface_RIN" / "ValidatedInterfaceInteractionEvents.csv"
        graph = coupled_chemistry_plot(analysis)
        if graph and chemistry_source.is_file():
            add("coupled", graph, "Validated interface chemistry through time", chemistry_source,
                note="ELI5: this shows which molecular glue is present over time. Hydrogen bonds, ionic interactions, and hydrophobic contacts are counted from validated atom-level events in 2 ns windows; taller colored areas mean that kind of glue was found in more frames.")
        rewiring_source = analysis / "Interface_RIN" / "interface_residue_pair_framewise.csv.gz"
        graph = coupled_rewiring_plot(analysis)
        if graph and rewiring_source.is_file():
            add("coupled", graph, "Interface contact rewiring", rewiring_source,
                label="Download exact framewise residue-pair contacts",
                note="ELI5: a stable interface keeps many of the same residue pairs from one frame to the next. Retention near 1 means the contact network is being reused; spikes in new or lost pairs mean the interface is rearranging.")
        if not interface_time.empty:
            for component_index, source in enumerate(sorted(analysis.glob("SSE_chain_*_fractions.csv")), start=1):
                label = source.stem.removeprefix("SSE_chain_").removesuffix("_fractions")
                graph = coupled_sse_binding_plot(source, interface_time, f"Component {component_index}: {label}", CHAIN_COLORS[(component_index - 1) % len(CHAIN_COLORS)])
                add("coupled", graph, f"Component {component_index}: {label} SSE versus interface engagement", source,
                    note="ELI5: this component's helix, sheet, and coil fractions are compared with the shared interface-contact line. Look for sustained changes that happen together, rather than a single noisy frame.", component=component_index)
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
    groups["secondary"] = secondary_heatmap_cards + groups["secondary"]
    groups["contacts"] = heatmap_cards + groups["contacts"]
    cards = ([card for values in groups.values() for card in values]
             + [card for lanes in component_cards.values() for values in lanes.values() for card in values]
             + [card for values in chain_secondary_cards.values() for card in values])
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
        if context.is_interface else ("Read global structure and interaction chemistry before drawing conclusions about the ligand." if has_ligand_evidence else "Read protein stability, fold retention, and secondary structure together.")
    )
    overview = (component_lanes_section_html(root, overview_title, "Start here", overview_text,
                                              groups["overview"], component_cards["overview"])
                if is_multichain_system else section_html(overview_title, "Start here", overview_text, groups["overview"]))
    # Interface-only systems should never acquire a misleading empty ligand card.
    ligand_panel = ligand_html(root, ligand_name) if context.partner_kind == "ligand" else ""
    network_explorer = interactive_networx_html(analysis, ligand_expected=has_ligand_evidence)
    peptide_topology = peptide_networx_html(analysis, peptide_expected=context.partner_kind == "peptide")
    interface_evidence = interface_evidence_html(analysis, context, interface_rin)
    structure_panel = final_structure_html(root, context)
    trajectory_panel = trajectory_playback_html(root, analysis, context)
    component_key = component_key_html(root) if is_multichain_system else ""
    stability = component_key + (component_lanes_section_html(root, "Stability", "Trajectory stability", f"Compare global protein, complex, {context.relationship_label}, compactness, and native-contact trends. Cyan denotes whole-complex observables; orange and neon green denote the individual components.", groups["stability"], component_cards["stability"])
                               if is_multichain_system else section_html("Stability", "Trajectory stability", f"Compare global protein, complex, {context.relationship_label}, compactness, and native-contact trends.", groups["stability"]))
    flexibility = component_key + (component_lanes_section_html(root, "Flexibility", "Residue and partner mobility", f"Inspect protein and {context.partner_label if context.partner_kind != 'none' else 'chain'} RMSF at atom and residue levels.", groups["flexibility"], component_cards["flexibility"])
                                 if is_multichain_system else section_html("Flexibility", "Residue and partner mobility", f"Inspect protein and {context.partner_label if context.partner_kind != 'none' else 'chain'} RMSF at atom and residue levels.", groups["flexibility"]))
    secondary = component_key + ((chain_secondary_section_html(chain_secondary_cards) if is_multichain_system else "") or section_html("Secondary structure", "Fold evolution", "Read the residue-time secondary-structure maps alongside global SSE fractions and persistence.", groups["secondary"]))
    contacts = component_key + (component_lanes_section_html(root, context.contact_tab_label, "Binding/interface evidence", f"Use {context.relationship_label} contact frequency, persistence, and event views together.", groups["contacts"], component_cards["contacts"])
                              if is_multichain_system else section_html(context.contact_tab_label, "Binding/interface evidence", f"Use {context.relationship_label} contact frequency, persistence, and event views together.", groups["contacts"]))
    chemistry = component_key + (component_lanes_section_html(root, "Contact chemistry", f"{context.relationship_label.title()} interaction types", "Hydrogen-bond, ionic, and hydrophobic interaction events remain separate and hoverable.", groups["chemistry"], component_cards["chemistry"])
                               if is_multichain_system else section_html("Contact chemistry", f"{context.relationship_label.title()} interaction types", "Hydrogen-bond, ionic, and hydrophobic interaction events remain separate and hoverable.", groups["chemistry"]))
    distances = component_key + (component_lanes_section_html(root, "Distances", "Partner-residue separations", f"Explore every calculated {context.relationship_label} separation without collapsing pair context.", groups["distances"], component_cards["distances"])
                                 if is_multichain_system else section_html("Distances", "Partner-residue separations", f"Explore every calculated {context.relationship_label} separation without collapsing pair context.", groups["distances"])) + distance_carousel_html(analysis)
    coupled = component_key + (component_lanes_section_html(root, "Coupled dynamics", "Structure and contacts", f"Compare secondary structure with {context.relationship_label} contacts and detected trajectory segments.", groups["coupled"], component_cards["coupled"])
                               if is_multichain_system else section_html("Coupled dynamics", "Structure and contacts", f"Compare secondary structure with {context.relationship_label} contacts and detected trajectory segments.", groups["coupled"]))
    dashboard = tabbed_sections([
        ("overview", "Overview", component_key + overview + peptide_topology + interface_rin["overview"] + network_explorer),
        ("stability", "Stability", stability),
        ("flexibility", "Flexibility", flexibility),
        ("secondary", "Secondary structure", secondary),
        ("contacts", context.contact_tab_label, ligand_panel + contacts + network_explorer + peptide_topology + (interface_rin["peptide_detail"] if context.partner_kind == "peptide" else interface_rin["protein_detail"] if context.partner_kind == "protein_chain" else "")),
        ("chemistry", "Contact chemistry", chemistry),
        ("distances", "Distances", distances),
        ("coupled", "Coupled dynamics", coupled),
        ("interface-evidence", "Interface maps", interface_evidence),
        ("interfaces", "Protein interfaces", "" if context.partner_kind == "protein_chain" else interface_rin["protein_detail"]),
        ("peptide-interfaces", "Peptide interfaces", "" if context.partner_kind == "peptide" else interface_rin["peptide_detail"]),
        ("structure", "Structure", structure_panel),
        ("trajectory", "Trajectory playback", trajectory_panel),
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
