#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MD ANALYSIS FIGUREBOOK BUILDER
Joseph-Michael Edition — Rg AND/OR SAFE VERSION

✔ Handles ligand-bound and ligandless/protein-centric reports
✔ Figure 8 uses intelligent AND/OR logic:
    - ligand mode: prefer Radius_of_Gyration_Overlay.png, fallback to Radius_of_Gyration_Protein.png
    - protein mode: prefer Radius_of_Gyration_Protein.png, fallback to Radius_of_Gyration_Overlay.png
✔ All text wraps
✔ No overflow
✔ Title page guaranteed content
✔ Publication-ready PDF
✔ Skips ligand-only figures automatically for ligandless runs
"""

import os
import fnmatch
import tempfile
import subprocess
import re
import csv
from typing import List, Optional, Tuple
from pathlib import Path

from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from pymacs_component_utils import load_system_registry

from reportlab.platypus import (
    Paragraph, PageBreak, Image as RLImage,
    BaseDocTemplate, Frame, PageTemplate
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.units import inch
from reportlab.graphics import renderPDF
from svglib.svglib import svg2rlg


TITLE_PAGE_DESCRIPTION_LIGAND = """
Molecular Dynamics Analysis of {LIG} in Complex with {Protein}

Trajectory Generated with GROMACS and the Customized PyMACS Analysis Framework

This report presents a comprehensive molecular dynamics (MD) analysis of a ligand-bound protein system, integrating structural, dynamic, and interaction-level metrics into a unified, publication-ready dataset. The simulation was performed using GROMACS, and all post-processing was executed using our in-house PyMACS analysis suite—a customized Python-based framework purpose-built for high-resolution drug binding trajectory evaluation.

The PyMACS pipeline couples curated structural bioinformatics algorithms with automated figure generation to quantify global stability, chain-specific motions, ligand flexibility, secondary-structure evolution, and the geometry and persistence of residue-level interactions. Together, these metrics allow detailed characterization of how the ligand engages the protein interface, how the interface adapts over time, and which structural features contribute most strongly to sustained binding.

Across this report, each figure provides a distinct analytical lens:

Global and chain-resolved RMSD/RMSF reveal the stability, convergence, and dynamic hotspots of the complex.

Ligand RMSD/RMSF track binding-mode stability and internal flexibility.

Secondary-structure analyses (DSSP) highlight folding transitions, helix/sheet integrity, and ligand-associated structural rearrangements.

Time-resolved contact maps, frequency ranking, and interaction-type distributions capture the chemical and temporal landscape of ligand recognition.

Distance distributions, persistence metrics, and interaction timelines decode the geometric and kinetic properties of individual contacts.

Interaction networks provide structural and quantitative visualization of the binding interface, highlighting key residues that anchor or modulate ligand engagement.

By integrating these orthogonal analyses, this report delivers a mechanistic portrait of the ligand–protein system: how stable the complex remains over time, how the ligand adapts within its pocket, which residues dominate affinity, and which structural rearrangements accompany sustained binding.

This document is intended as a rigorous, standalone reference for comparing drug designs, validating docking hypotheses, guiding mutagenesis or scaffold optimization, and supporting manuscript-quality mechanistic interpretation.
"""

TITLE_PAGE_DESCRIPTION_PROTEIN = """
Molecular Dynamics Analysis of {Protein}

Trajectory Generated with GROMACS and the Customized PyMACS Analysis Framework

This report presents a comprehensive molecular dynamics (MD) analysis of a protein-centric simulation, integrating structural and dynamical metrics into a unified, publication-ready dataset. The simulation was performed using GROMACS, and all post-processing was executed using our in-house PyMACS analysis suite—a customized Python-based framework designed for high-resolution trajectory evaluation of apo proteins, protein-only systems, and peptide–peptide or protein–protein assemblies where ligand-specific analyses are not applicable.

The PyMACS pipeline quantifies global structural stability, chain-specific motions, residue-level flexibility, compactness, and secondary-structure evolution over time. Together, these metrics provide a detailed view of how the protein architecture behaves throughout the simulation, which regions remain structurally anchored, and which regions undergo adaptive motions or conformational rearrangement.

Across this report, each figure provides a distinct analytical lens:

Global and chain-resolved RMSD/RMSF reveal stability, convergence, and dynamic hotspots.

Radius of gyration reports on overall compactness and breathing motions of the protein assembly.

Secondary-structure analyses (DSSP) highlight helix/sheet integrity, unfolding events, and structural transitions over time.

Residue-level persistence plots identify which regions remain structurally stable and which residues participate in flexible or adaptive motions.

By integrating these orthogonal analyses, this report delivers a mechanistic portrait of the simulated protein system: how stable the structure remains, which domains or chains move most strongly, whether compactness is preserved, and which secondary-structure elements remain persistent or undergo reorganization.

This document is intended as a rigorous, standalone reference for evaluating protein stability, comparing apo and bound states, interpreting domain motions, supporting manuscript-quality mechanistic discussion, and guiding subsequent modeling or experimental follow-up.
"""

TITLE_PAGE_DESCRIPTION_BIOLOGICAL = """
Molecular Dynamics Analysis of {Protein}

Trajectory Generated with GROMACS and the Customized PyMACS Analysis Framework

This report presents a comprehensive molecular dynamics (MD) analysis of a biological polymer system, integrating protein and nucleic-acid structural metrics into a unified, publication-ready dataset. The simulation was performed using GROMACS, and all post-processing was executed using our in-house PyMACS analysis suite—a customized Python-based framework designed for high-resolution trajectory evaluation of protein:RNA, protein:DNA, RNA-only, DNA-only, and mixed biomolecular assemblies where ligand-specific analyses are not applicable.

The PyMACS pipeline quantifies global structural stability, chain-specific motions, residue-level flexibility, compactness, secondary-structure evolution, and chain-to-chain biomolecular contacts over time. Together, these metrics provide a detailed view of how the polymer assembly behaves throughout the simulation, which chains remain structurally anchored, which segments undergo adaptive motions, and how strongly the biomolecular partners remain associated during the trajectory.

Across this report, each figure provides a distinct analytical lens:

Global and chain-resolved RMSD/RMSF reveal stability, convergence, and dynamic hotspots across protein and nucleic-acid chains.

Radius of gyration reports on overall compactness and breathing motions of the biomolecular assembly.

Secondary-structure analyses (DSSP) capture folding transitions, helix/sheet integrity, and structural persistence for the protein components.

Chain–chain biomolecule interaction summaries and distance traces report how strongly the polymer partners remain associated, helping identify persistent interfaces versus transient contacts.

By integrating these orthogonal analyses, this report delivers a mechanistic portrait of the simulated biological system: how stable the assembly remains, which chains or domains move most strongly, whether compactness is preserved, and whether the biomolecular interface remains persistent or reorganizes over time.

This document is intended as a rigorous, standalone reference for evaluating biological macromolecular assemblies, comparing apo and bound-like states, interpreting chain-level motions, and supporting manuscript-quality mechanistic discussion.
"""

PROTAC_MANIFEST_CSV = os.path.join("Analysis_Results", "PROTAC", "QC", "protac_figure_manifest.csv")


GENERIC_CATEGORY_NOTES = {
    "protein_rmsd": "This figure summarizes RMSD-based structural drift for the protein components over the trajectory and helps assess equilibration, stability, and large-scale conformational change.",
    "protein_rmsf": "This figure reports residue-level flexibility and highlights dynamic hotspots, mobile loops, and rigid structural cores across the protein components.",
    "ligand_rmsd": "This figure tracks the ligand / PROTAC positional stability over time relative to its starting conformation.",
    "ligand_rmsf": "This figure reports ligand / PROTAC internal flexibility at the atom level across the simulation.",
    "radius_of_gyration": "This figure reports radius of gyration as a compactness metric for the full complex and component proteins over time.",
    "contacts_timeseries": "This figure summarizes time-resolved ligand contact persistence across the protein components and helps compare which partner side engages the PROTAC most strongly.",
    "contact_heatmap": "This heatmap summarizes which residues contact the ligand over time and is ordered to highlight persistent interface hotspots.",
    "interaction_network": "This network view highlights the most persistent residue contacts to the ligand, emphasizing the strongest interaction hotspots while suppressing low-frequency clutter.",
    "interaction_composition": "This figure breaks residue-contact events into broad interaction classes to compare how each protein component contributes hydrophobic, ionic, aromatic, and hydrophilic contacts.",
    "pca": "This PCA projection summarizes dominant collective protein motions sampled during the trajectory.",
}


def safe_input(prompt: str, default: str = "") -> str:
    try:
        value = input(prompt).strip()
        return value if value else default
    except EOFError:
        return default


def prepare_image_for_pdf(img_path: str, max_width_px: int = 1600) -> Optional[str]:
    if not os.path.exists(img_path):
        return None

    try:
        img = Image.open(img_path)
        w, h = img.size

        if w > max_width_px:
            scale = max_width_px / float(w)
            img = img.resize((int(w * scale), int(h * scale)), Image.LANCZOS)

        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".jpg")
        img.convert("RGB").save(tmp.name, "JPEG", quality=92)
        return tmp.name
    except Exception as exc:
        print(f"⚠️ Failed to prepare image {img_path}: {exc}")
        return None


def detect_report_mode() -> str:
    print("\nSelect report type:")
    print("  1. Ligand-bound / small-molecule / peptide-bound")
    print("  2. Ligandless / protein-centric / peptide:peptide / protein:protein")
    print("  3. Biological system / protein:RNA / protein:DNA / RNA:DNA")
    choice = safe_input("Enter choice [1/2/3, default 2]: ", "2")
    if choice == "1":
        return "ligand"
    if choice == "3":
        return "biological"
    return "protein"


def detect_ligand_mol2() -> str:
    mol2 = sorted([f for f in os.listdir(".") if f.endswith(".cgenff.mol2")])
    if not mol2:
        return safe_input("Enter ligand code: ").upper()

    if len(mol2) == 1:
        code = mol2[0].replace(".cgenff.mol2", "")
        if safe_input(f"Detected {code}. Use? [Y/n]: ", "y").lower() in ("", "y", "yes"):
            return code.upper()

    for i, f in enumerate(mol2, 1):
        print(f"{i}) {f}")
    sel = safe_input("Pick # or ENTER to type manually: ", "")
    if sel.isdigit() and 1 <= int(sel) <= len(mol2):
        return mol2[int(sel) - 1].replace(".cgenff.mol2", "").upper()
    return safe_input("Enter ligand code: ").upper()


def chain_label_from_index(idx: int) -> str:
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if idx < len(letters):
        return letters[idx]
    return f"{letters[(idx // 26) - 1]}{letters[idx % 26]}"


def parse_atomindex() -> Tuple[List[str], List[str]]:
    if not os.path.exists("atomIndex.txt"):
        return ["Protein"], ["A"]

    proteins = []
    chains = []

    with open("atomIndex.txt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            proteins.append(line.split()[0])
            chains.append(chain_label_from_index(len(chains)))

    if not proteins:
        return ["Protein"], ["A"]

    return proteins, chains


def load_caption_file_text() -> Optional[str]:
    for candidate in ("4_MDfigs.txt",):
        if os.path.exists(candidate):
            with open(candidate) as f:
                return f.read()
    return None


def load_graph_notes_map(path: str = "4_GraphNotes.txt") -> dict:
    if not os.path.exists(path):
        return {}

    text = Path(path).read_text(encoding="utf-8")
    blocks = [block.strip() for block in re.split(r"\n\s*\n+", text) if block.strip()]
    notes = {}
    for block in blocks:
        lines = [line.strip() for line in block.splitlines() if line.strip()]
        if len(lines) < 2:
            continue
        key = lines[0]
        notes[key] = "\n".join(lines[1:])
    return notes


def normalize_bool(value) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def detect_protac_manifest() -> Optional[str]:
    return PROTAC_MANIFEST_CSV if os.path.exists(PROTAC_MANIFEST_CSV) else None


def load_protac_manifest_rows(manifest_csv: str) -> List[dict]:
    rows = []
    with open(manifest_csv, "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            row["include_in_pdf"] = normalize_bool(row.get("include_in_pdf", "true"))
            try:
                row["sort_order"] = int(float(row.get("sort_order", 0)))
            except ValueError:
                row["sort_order"] = 0
            rows.append(row)
    rows.sort(key=lambda row: (row["sort_order"], row.get("figure_title", "")))
    return rows


def note_for_manifest_row(row: dict, notes_map: dict) -> str:
    description_key = row.get("description_key", "").strip()
    figure_category = row.get("figure_category", "").strip()
    notes_hint = row.get("notes_hint", "").strip()

    if description_key and description_key in notes_map:
        return notes_map[description_key]
    if figure_category and figure_category in notes_map:
        return notes_map[figure_category]
    base = GENERIC_CATEGORY_NOTES.get(figure_category, "This figure documents one aspect of the trajectory analysis generated by the PyMACS workflow.")
    if notes_hint:
        return f"{base}\n\nHint: {notes_hint}"
    return base


def resolve_possible_files(target: str, base: str = "Analysis_Results") -> List[str]:
    matches = []
    if not os.path.exists(base):
        return matches

    for root, _, files in os.walk(base):
        for f in files:
            if f == target or fnmatch.fnmatch(f, target):
                matches.append(os.path.join(root, f))

    return sorted(set(matches))


def find_first_existing(candidates: List[str], base: str = "Analysis_Results") -> Tuple[Optional[str], Optional[str]]:
    for candidate in candidates:
        direct = os.path.join(base, candidate)
        if os.path.exists(direct):
            return candidate, direct
        matches = resolve_possible_files(candidate, base=base)
        if matches:
            return candidate, matches[0]
    return None, None


def generate_rdkit_hero_image_as_svg(mol2_file: str, width: int = 900, height: int = 900, use_color: bool = True) -> str:
    result = subprocess.run(
        ["obabel", mol2_file, "-osmi", "-xk"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )

    smiles = result.stdout.strip().split()[0]
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES: {smiles}")

    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    Chem.rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.bondLineWidth = 2
    options.fixedBondLength = 25
    options.dotsPerAngstrom = 100
    if not use_color:
        options.useBWAtomPalette()

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    out = tempfile.NamedTemporaryFile(delete=False, suffix=".svg").name
    with open(out, "w") as f:
        f.write(drawer.GetDrawingText())

    return out


def add_header(canvas, doc):
    canvas.saveState()
    if doc.page == 1:
        header_text = "PyMACS Molecular Dynamics Simulation Report"
        canvas.setFont("Helvetica-Bold", 14)
    else:
        header_text = f"PyMACS Molecular Dynamics Simulation Report — Page {doc.page}"
        canvas.setFont("Helvetica-Bold", 10)
    canvas.drawCentredString(doc.pagesize[0] / 2, doc.pagesize[1] - 30, header_text)
    canvas.restoreState()


def add_ligand_to_pages(canvas, doc, small_ligand_svg: Optional[str] = None, y_offset: int = 20):
    if not small_ligand_svg:
        return
    if doc.page % 2 == 1:
        canvas.saveState()
        drawing = svg2rlg(small_ligand_svg)
        if drawing is not None:
            sx = (0.8 * inch) / max(drawing.width, 1)
            sy = (0.8 * inch) / max(drawing.height, 1)
            drawing.scale(sx, sy)
            renderPDF.draw(drawing, canvas, 40, y_offset)
        canvas.restoreState()


def build_title_page(flow, styles, mode: str, pretty_ligand: Optional[str], protein_title: str):
    if mode == "ligand":
        title_text = f"MD ANALYSIS REPORT — {pretty_ligand} Bound to {protein_title}"
        description = TITLE_PAGE_DESCRIPTION_LIGAND.format(LIG=pretty_ligand, Protein=protein_title)
    elif mode == "biological":
        title_text = f"MD ANALYSIS REPORT — Biological System: {protein_title}"
        description = TITLE_PAGE_DESCRIPTION_BIOLOGICAL.format(Protein=protein_title)
    else:
        title_text = f"MD ANALYSIS REPORT — {protein_title}"
        description = TITLE_PAGE_DESCRIPTION_PROTEIN.format(Protein=protein_title)

    flow.append(Paragraph(title_text, styles["TitleBig"]))
    flow.append(Paragraph(description.replace("\n", "<br/>"), styles["Body"]))
    flow.append(PageBreak())


def get_distance_histogram_caption(i: int) -> Tuple[str, str]:
    suffix = chr(97 + i)
    title = f"Figure 30{suffix}. Interaction Distance Distribution"
    body = (
        "Narrow peaks suggest stable geometry; broad peaks indicate dynamic or transient contacts. "
        "These histograms help distinguish structurally anchored contacts from interactions that sample "
        "multiple conformational states during the simulation."
    )
    return title, body


def get_rg_metadata_for_mode(mode: str):
    overlay_caption = (
        "Figure 8. Global compactness dynamics of the protein–ligand complex measured by radius of gyration.",
        "Time evolution of the radius of gyration (Rg) for the protein backbone, bound ligand, and full complex over the course of the molecular dynamics simulation. Protein Rg reports on global folding stability and breathing motions, while ligand Rg captures intramolecular collapse or extension during binding. The complex Rg reflects overall assembly compactness. Together, these metrics provide a complementary view of conformational stability and structural organization beyond RMSD-based analyses."
    )
    protein_caption = (
        "Figure 8. Global compactness dynamics of the protein measured by radius of gyration.",
        "Time evolution of the radius of gyration (Rg) for the protein over the course of the molecular dynamics simulation. Rg reports on the global compactness of the structure and can reveal folding stability, breathing motions, opening/closing behavior, or expansion of flexible domains. Stable Rg suggests maintenance of the overall fold, while systematic increases or broad fluctuations may indicate conformational loosening or sampling of more open states."
    )

    if mode == "ligand":
        return ["Radius_of_Gyration_Overlay.png", "Radius_of_Gyration_Protein.png"], {
            "Radius_of_Gyration_Overlay.png": overlay_caption,
            "Radius_of_Gyration_Protein.png": protein_caption,
        }
    return ["Radius_of_Gyration_Protein.png", "Radius_of_Gyration_Overlay.png"], {
        "Radius_of_Gyration_Protein.png": protein_caption,
        "Radius_of_Gyration_Overlay.png": (
            "Figure 8. Global compactness dynamics measured by radius of gyration.",
            "This figure tracks the radius of gyration (Rg) of the simulated system over time. In protein-centric runs, the protein trace provides the most relevant measure of global compactness and structural breathing. Stable values suggest preserved folding and limited expansion or collapse, whereas sustained shifts may indicate compaction, loosening, or domain-level rearrangement."
        ),
    }


def load_main_captions_from_file(ligand: str, chains: List[str], proteins: List[str], mode: str):
    text = load_caption_file_text()
    if not text:
        return []

    blocks = text.split("---")
    out = []
    rg_seen = False

    for b in blocks:
        lines = [x.strip() for x in b.splitlines() if x.strip()]
        if len(lines) < 2:
            continue

        fname = lines[0].replace("{LIG}", ligand)
        title = lines[1].replace("{LIG}", ligand)
        caption = "<br/>".join(lines[2:]).replace("{LIG}", ligand)

        if fname in {"Radius_of_Gyration_Overlay.png", "Radius_of_Gyration_Protein.png"}:
            if not rg_seen:
                candidates, caption_map = get_rg_metadata_for_mode(mode)
                out.append({
                    "type": "alternatives",
                    "candidates": candidates,
                    "caption_map": caption_map,
                })
                rg_seen = True
            continue

        if "{CHAIN}" in fname:
            for c in chains:
                out.append({
                    "type": "single",
                    "filename": fname.replace("{CHAIN}", str(c)),
                    "title": title.replace("{CHAIN}", str(c)),
                    "caption": caption.replace("{CHAIN}", str(c)),
                })
            continue

        if "{PROTEIN_NAME}" in fname:
            for p in proteins:
                out.append({
                    "type": "single",
                    "filename": fname.replace("{PROTEIN_NAME}", p),
                    "title": title.replace("{PROTEIN_NAME}", p),
                    "caption": caption.replace("{PROTEIN_NAME}", p),
                })
            continue

        out.append({
            "type": "single",
            "filename": fname,
            "title": title,
            "caption": caption,
        })

    if not rg_seen:
        candidates, caption_map = get_rg_metadata_for_mode(mode)
        out.append({
            "type": "alternatives",
            "candidates": candidates,
            "caption_map": caption_map,
        })

    return out


def load_protein_centric_captions(chains: List[str], proteins: List[str]):
    template_entries = [
        ("Protein_RMSD.png", "Figure 1. Global Protein RMSD Over Time",
         "This figure presents the backbone root-mean-square deviation (RMSD) of the protein throughout the molecular dynamics simulation. RMSD reflects the average positional shift of backbone atoms compared to the initial reference structure. A low, steady RMSD indicates structural stability and equilibration, while progressive increases may signify conformational drift, domain movement, or gradual unfolding. Sharp spikes often correspond to flexible loops or transient motions. This plot serves as a primary metric for verifying simulation stability and assessing major structural transitions."),
        ("RMSD_{PROTEIN_NAME}.png", "Figure 2. Structural Stability of Individual Protein Chains",
         "This plot isolates RMSD for each protein chain independently, allowing detailed comparison of their structural behavior. In multichain systems, chains may experience different environments, constraints, or degrees of flexibility. A stable chain shows low RMSD, while a rising or fluctuating RMSD may indicate hinge motions, inter-domain rearrangements, or partial unfolding. This chain-level decomposition reveals asymmetry in protein flexibility and helps interpret localized versus collective motions."),
        ("RMSF_{PROTEIN_NAME}.png", "Figure 3. Atom-Level Flexibility of Individual Protein Chains",
         "This figure shows the root-mean-square fluctuation (RMSF) of individual atoms within the selected protein chain. RMSF quantifies how strongly each atom fluctuates around its time-averaged position throughout the simulation. Higher values often correspond to flexible loops, termini, or solvent-exposed regions, whereas lower values indicate rigid structural cores. This view complements RMSD by resolving local rather than global dynamics."),
        ("RMSF_{PROTEIN_NAME}_per_residue_contacts.png", "Figure 4. Residue-Level Flexibility Profile",
         "This figure shows the root-mean-square fluctuation (RMSF) of individual protein residues for the analyzed chain, calculated from Cα atoms over the molecular dynamics trajectory. RMSF quantifies the extent of positional fluctuation around each residue’s average structure, with higher values typically corresponding to flexible loops, termini, or regulatory regions, and lower values indicating more rigid secondary-structure elements such as α-helices or β-sheets. Breaks in the plotted line generally reflect discontinuities in residue numbering or unresolved structural segments rather than true interruptions in dynamics. This plot highlights local flexibility hotspots and helps identify candidate hinge regions, exposed loops, or mobile termini."),
        ("FullComplex_RMSD.png", "Figure 5. Global RMSD of the Full Simulated Protein System",
         "This plot reports RMSD for the entire simulated system, including all protein chains present in the trajectory. It captures global changes such as domain motions, multimer rearrangements, or cooperative structural shifts. Stable full-system RMSD indicates overall convergence and equilibration, whereas significant increases may reflect large-scale breathing motions or long-timescale adaptations."),
        ("DSSP_Heatmap.png", "Figure 7. Secondary Structure Evolution Over Time (DSSP Heatmap)",
         "This heatmap visualizes how secondary-structure elements—α-helices, β-sheets, and coil regions—change throughout the simulation. Each row corresponds to a residue and each column to a timepoint. Continuous colored bands indicate stable structural motifs, while transitions or disruptions signal unfolding events, local melting, or structural rearrangements. DSSP analysis is particularly valuable for identifying whether the protein retains its folded architecture or undergoes structural adaptation."),
        ("SSE_Fractions.png", "Figure 8. Global Secondary Structure Composition",
         "This plot shows the fraction of residues adopting helix, sheet, or coil conformations across the trajectory. Stability in overall composition suggests preserved structural integrity, while shifts may indicate progressive conformational change, tightening, or partial disordering. This figure provides a concise global summary of the protein’s structural state."),
        ("DSSP_chain_{CHAIN}_heatmap.png", "Figure 9. Chain-Resolved Secondary Structure Dynamics",
         "This figure maps secondary-structure changes over time for an individual protein chain. It clarifies whether structural transitions observed in the global analysis originate from a specific chain or are distributed across the complex. Chain-specific DSSP patterns help identify asymmetric responses, cooperative structural shifts, or localized unfolding."),
        ("SSE_chain_{CHAIN}_fractions.png", "Figure 10. Chain-Resolved Secondary Structure Composition",
         "This figure summarizes the fraction of helix, sheet, and coil states over time for an individual chain. Comparing these profiles across chains helps distinguish localized structural changes from system-wide trends and highlights whether one chain is especially dynamic or structurally stable."),
        ("Residue_Helix_Persistence.png", "Figure 11a. Residue-Level Helix Persistence",
         "This figure provides a residue-by-residue analysis of how consistently each residue remains in a helical state throughout the simulation. Persistence refers to the fraction of simulation frames in which a residue maintains a particular structural assignment. High helical persistence indicates stable α-helical segments or rigid structural cores."),
        ("Residue_Sheet_Persistence.png", "Figure 11b. Residue-Level Sheet Persistence",
         "This figure provides a residue-by-residue analysis of how consistently each residue remains in a β-sheet state throughout the simulation. High sheet persistence indicates residues that form stable β-structure elements, whereas low persistence suggests structural plasticity or transitions into alternative conformational states."),
        ("Residue_Flexibility.png", "Figure 11c. Residue-Level Flexibility Index",
         "This figure shows a residue-level flexibility index derived from secondary-structure persistence. Residues with higher values spend more time in flexible or coil-like states and may correspond to hinges, loops, termini, or adaptive regions that participate in conformational rearrangement."),
    ]

    out = []
    for fname, title, caption in template_entries:
        if "{CHAIN}" in fname:
            for c in chains:
                out.append({
                    "type": "single",
                    "filename": fname.replace("{CHAIN}", str(c)),
                    "title": title.replace("{CHAIN}", str(c)),
                    "caption": caption.replace("{CHAIN}", str(c)),
                })
        elif "{PROTEIN_NAME}" in fname:
            for p in proteins:
                out.append({
                    "type": "single",
                    "filename": fname.replace("{PROTEIN_NAME}", p),
                    "title": title.replace("{PROTEIN_NAME}", p),
                    "caption": caption.replace("{PROTEIN_NAME}", p),
                })
        else:
            out.append({
                "type": "single",
                "filename": fname,
                "title": title,
                "caption": caption,
            })

    candidates, caption_map = get_rg_metadata_for_mode("protein")
    out.insert(5, {
        "type": "alternatives",
        "candidates": candidates,
        "caption_map": caption_map,
    })

    return out


def load_biological_system_captions(chains: List[str], proteins: List[str]):
    out = load_protein_centric_captions(chains, proteins)

    biomolecule_entries = [
        (
            "Biomolecule_Chain_Interaction_Heatmap.png",
            "Figure 12. Biomolecule Chain–Chain Contact Fractions",
            "This heatmap summarizes the fraction of simulation frames in which each biomolecular chain pair remained in contact under the configured distance cutoff. Higher values indicate persistent interfaces between protein and nucleic-acid partners, whereas lower values suggest transient association or limited chain–chain engagement.",
        ),
        (
            "Biomolecule_Chain_Interaction_Distances.png",
            "Figure 13. Biomolecule Chain–Chain Minimum Distance Over Time",
            "This figure tracks the minimum anchor-atom distance between each biomolecular chain pair over time. Stable low distances indicate sustained association of the polymer partners, while increases or large fluctuations suggest interface breathing, partial dissociation, or conformational rearrangement.",
        ),
    ]

    for filename, title, caption in biomolecule_entries:
        out.append({
            "type": "single",
            "filename": filename,
            "title": title,
            "caption": caption,
        })

    return out


def add_figure_and_caption(flow, styles, image_path: str, title: str, caption: str):
    img = prepare_image_for_pdf(image_path)
    if not img:
        print(f"⚠️ Unable to process figure: {image_path}")
        return False

    flow.append(RLImage(img, width=9.2 * inch, height=6.0 * inch))
    flow.append(PageBreak())
    flow.append(Paragraph(title, styles["CaptionTitle"]))
    flow.append(Paragraph(caption.replace("\n", "<br/>"), styles["Caption"]))
    flow.append(PageBreak())
    return True


def build_protac_pdf_from_manifest(manifest_csv: str):
    rows = load_protac_manifest_rows(manifest_csv)
    included_rows = [row for row in rows if row.get("include_in_pdf")]
    notes_map = load_graph_notes_map()

    ligand_code = included_rows[0].get("ligand_code", "PROTAC") if included_rows else "PROTAC"
    protein_components = []
    for row in included_rows:
        protein_components.extend([x.strip() for x in row.get("protein_components", "").split(",") if x.strip()])
    protein_components = list(dict.fromkeys(protein_components))
    protein_title = " & ".join(protein_components) if protein_components else "PROTAC Complex"

    output_name = safe_input(
        "Output PDF name [default PROTAC_MD_ANALYSIS_FIGUREBOOK.pdf]: ",
        "PROTAC_MD_ANALYSIS_FIGUREBOOK.pdf",
    )
    if not output_name.lower().endswith(".pdf"):
        output_name += ".pdf"

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="TitleBig", fontSize=22, spaceAfter=24, alignment=1, leading=28))
    styles.add(ParagraphStyle(name="Body", fontSize=11, leading=14, spaceAfter=12))
    styles.add(ParagraphStyle(name="CaptionTitle", fontSize=13, leading=16, spaceAfter=8, fontName="Helvetica-Bold"))
    styles.add(ParagraphStyle(name="Caption", fontSize=10, leading=13, spaceBefore=4))

    doc = BaseDocTemplate(
        output_name,
        pagesize=landscape(letter),
        rightMargin=50,
        leftMargin=50,
        topMargin=50,
        bottomMargin=16,
    )

    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height - 40, id="normal")
    doc.addPageTemplates([
        PageTemplate(id="all_pages", frames=frame, onPage=lambda canvas, d: add_header(canvas, d))
    ])

    flow = []
    build_title_page(flow, styles, "ligand", ligand_code, protein_title)

    any_figures_added = False
    for row in included_rows:
        figure_path = row.get("figure_path", "")
        if not figure_path or not os.path.exists(figure_path):
            print(f"⚠️ Manifest figure missing on disk: {figure_path}")
            continue

        caption = note_for_manifest_row(row, notes_map)
        title = row.get("figure_title", Path(figure_path).name)
        added = add_figure_and_caption(flow, styles, figure_path, title, caption)
        any_figures_added = any_figures_added or added

    if not any_figures_added:
        flow.append(
            Paragraph(
                "No manifest-listed PROTAC figures were found on disk. The report was created with a title page only.",
                styles["Body"],
            )
        )

    doc.build(flow)
    print(f"📘 {output_name} COMPLETE (PROTAC manifest mode)")


def build_pdf():
    mode = detect_report_mode()
    system_registry = load_system_registry(".") or {}

    ligand = None
    pretty_ligand = None
    small_ligand_svg = None

    proteins, chains = parse_atomindex()
    protein_title = " & ".join(proteins)

    if mode == "protein" and (system_registry.get("has_rna") or system_registry.get("has_dna")):
        print(
            "🧬 Detected a nucleic-acid-containing system from pymacs_system_registry.json. "
            "Consider selecting report type 3 for the biomolecule-aware report."
        )

    if mode == "ligand":
        ligand = detect_ligand_mol2()
        pretty_ligand = safe_input(f"Pretty ligand name for {ligand} [default {ligand}]: ", ligand).strip() or ligand

        mol2_file = f"{ligand}.cgenff.mol2"
        if os.path.exists(mol2_file):
            try:
                small_ligand_svg = generate_rdkit_hero_image_as_svg(mol2_file, width=300, height=300, use_color=True)
            except Exception as exc:
                print(f"⚠️ Unable to generate ligand corner image: {exc}")
                small_ligand_svg = None

        captions = load_main_captions_from_file(ligand, chains, proteins, mode="ligand")
    elif mode == "biological":
        captions = load_biological_system_captions(chains, proteins)
    else:
        captions = load_protein_centric_captions(chains, proteins)

    output_name = safe_input("Output PDF name [default MD_ANALYSIS_FIGUREBOOK.pdf]: ", "MD_ANALYSIS_FIGUREBOOK.pdf")
    if not output_name.lower().endswith(".pdf"):
        output_name += ".pdf"

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="TitleBig", fontSize=22, spaceAfter=24, alignment=1, leading=28))
    styles.add(ParagraphStyle(name="Body", fontSize=11, leading=14, spaceAfter=12))
    styles.add(ParagraphStyle(name="CaptionTitle", fontSize=13, leading=16, spaceAfter=8, fontName="Helvetica-Bold"))
    styles.add(ParagraphStyle(name="Caption", fontSize=10, leading=13, spaceBefore=4))

    doc = BaseDocTemplate(
        output_name,
        pagesize=landscape(letter),
        rightMargin=50,
        leftMargin=50,
        topMargin=50,
        bottomMargin=16,
    )

    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height - 40, id="normal")
    doc.addPageTemplates([
        PageTemplate(
            id="all_pages",
            frames=frame,
            onPage=lambda canvas, d: (
                add_header(canvas, d),
                add_ligand_to_pages(canvas, d, small_ligand_svg, y_offset=20),
            ),
        )
    ])

    flow = []
    build_title_page(flow, styles, mode, pretty_ligand, protein_title)

    any_figures_added = False

    for entry in captions:
        if entry["type"] == "single":
            filename = entry["filename"]
            direct = os.path.join("Analysis_Results", filename)
            files = [direct] if os.path.exists(direct) else resolve_possible_files(filename)

            if not files:
                print(f"⚠️ Figure not found: {filename}")
                continue

            for f in files:
                added = add_figure_and_caption(flow, styles, f, entry["title"], entry["caption"])
                any_figures_added = any_figures_added or added

        elif entry["type"] == "alternatives":
            matched_name, matched_path = find_first_existing(entry["candidates"], base="Analysis_Results")
            if not matched_path:
                print(f"⚠️ Figure not found (alternatives tried): {', '.join(entry['candidates'])}")
                continue

            title, caption = entry["caption_map"][matched_name]
            added = add_figure_and_caption(flow, styles, matched_path, title, caption)
            any_figures_added = any_figures_added or added

    if mode == "ligand":
        dist_dir = os.path.join("Analysis_Results", "RESIDUE_Contact_Distance")
        if os.path.exists(dist_dir):
            dists = sorted(
                f for f in os.listdir(dist_dir)
                if f.lower().endswith((".png", ".jpg", ".jpeg")) and f.startswith("dist_")
            )
            for i, dist_file in enumerate(dists):
                img = os.path.join(dist_dir, dist_file)
                title, caption = get_distance_histogram_caption(i)
                added = add_figure_and_caption(flow, styles, img, title, caption)
                any_figures_added = any_figures_added or added
        else:
            print(f"ℹ️ Distance histogram directory not found: {dist_dir}")

    if not any_figures_added:
        flow.append(Paragraph(
            "No matching figures were found in Analysis_Results. The report was created with a title page only.",
            styles["Body"],
        ))

    doc.build(flow)
    print(f"📘 {output_name} COMPLETE")


if __name__ == "__main__":
    manifest_csv = detect_protac_manifest()
    if manifest_csv:
        build_protac_pdf_from_manifest(manifest_csv)
    else:
        build_pdf()
