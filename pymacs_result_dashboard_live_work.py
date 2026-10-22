#!/usr/bin/env python3
"""Build a self-contained interactive PyMACS result dossier for one MD run.

This intentionally complements (rather than replaces) 00_Pymacs_X_analysis.py.
The latter is a comparative apo/ligand dashboard; this program is the automatic
per-system report that follows a normal 3A/PROTAC analysis job.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import subprocess
import textwrap
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
        "Q is the fraction of protein–ligand heavy-atom contacts from the first analyzed frame that remain within the stated cutoff.",
        "Higher, sustained Q means the starting modeled interface is retained. Q evaluates persistence relative to that model; it does not establish experimental binding affinity.",
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


def read_metadata(root: Path) -> dict:
    try:
        return json.loads((root / "pymacs_run.json").read_text()).get("metadata", {})
    except (OSError, ValueError, TypeError):
        return {}


def report_csv_files(csv_files: list[Path], setup: dict) -> list[Path]:
    """Return the non-redundant CSVs appropriate for the single-run gallery.

    The analysis directory intentionally keeps both global, per-chain, per-atom,
    and overlay tables.  They are all useful scientific exports, but several
    become the *same visible plot* in a single-chain system.  Keep the raw
    files intact and curate only the automatic gallery.
    """
    names = {path.name.lower() for path in csv_files}
    chain_map = str(setup.get("chain_map") or "")
    protein_chains = [item for item in chain_map.split(",") if item.strip()]
    one_protein_chain = len(protein_chains) == 1
    selected: list[Path] = []
    for csv_file in csv_files:
        stem = csv_file.stem.lower()
        # This three-column table contains the same protein trace as
        # Protein_RMSD.csv plus a ligand trace that is already published in
        # the dedicated ligand-RMSD card.  numeric_plot would otherwise pick
        # only its first numeric signal and duplicate Protein_RMSD.
        if stem.endswith("_complex_rmsd_overlay"):
            continue
        # Prefer residue RMSF over atom-index RMSF when both are present.
        # The per-atom export remains downloadable in the data bundle.
        if stem.startswith("rmsf_") and not stem.endswith("_per_residue"):
            residue_variant = f"{csv_file.stem}_per_residue.csv".lower()
            if residue_variant in names:
                continue
        if one_protein_chain:
            # In a one-chain protein, the per-chain RMSD/fraction tables are
            # redundant with global Protein_RMSD and SSE_Fractions.  Keep the
            # global summaries; preserve all originals on disk.
            if stem.startswith("rmsd_"):
                continue
            if stem.startswith("sse_chain_") and stem.endswith("_fractions"):
                continue
        selected.append(csv_file)
    return selected


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


def build_document(system_name: str, ligand_name: str, cards: list[str], galleries: str) -> str:
    """Return the Kyle-hosted, self-contained single-run reading room."""
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{escape(system_name)} · PyMACS results</title><script>{plotlyjs()}</script><script src="/static/3Dmol-min.js"></script><style>
:root{{--bg:#05060a;--panel:#0b1020;--grid:#1c2a4a;--text:#cfe8ff;--muted:#7aa6d9;--cyan:#35d0ff;--green:#18ff8b;--yellow:#ffe66d;--magenta:#ff2bd6;--orange:#ff8a3d}}*{{box-sizing:border-box}}html{{scroll-behavior:smooth}}body{{margin:0;background:radial-gradient(circle at 90% 0,rgba(53,208,255,.15),transparent 30%),radial-gradient(circle at 4% 25%,rgba(255,43,214,.10),transparent 27%),var(--bg);color:var(--text);font-family:Inter,ui-sans-serif,system-ui,-apple-system,sans-serif}}.page{{max-width:1560px;margin:auto;padding:32px 24px 70px}}.hero,.card,.ligand-card,.gallery-section{{background:linear-gradient(155deg,#10172b,#060912);border:1px solid #23375d;box-shadow:0 18px 50px rgba(0,0,0,.35)}}.hero{{position:relative;overflow:hidden;border-radius:25px;padding:34px;margin-bottom:20px}}.hero:before{{content:"";position:absolute;inset:-48%;background:conic-gradient(from 180deg,transparent,rgba(53,208,255,.16),transparent,rgba(255,43,214,.13),transparent);animation:slow-spin 28s linear infinite;pointer-events:none}}.hero>*{{position:relative;z-index:1}}@keyframes slow-spin{{to{{transform:rotate(360deg)}}}}.eyebrow{{color:var(--green);font-size:.75rem;font-weight:900;letter-spacing:.16em;text-transform:uppercase}}h1{{font-size:clamp(2.5rem,5.2vw,5.2rem);letter-spacing:-.065em;line-height:.92;margin:.35rem 0 1rem}}h2{{letter-spacing:-.025em;margin:.2rem 0 .45rem;font-size:1.5rem}}p{{line-height:1.55}}.sub,.section-head p,.muted,.asset-note{{color:var(--muted);max-width:960px}}.pills{{display:flex;gap:9px;flex-wrap:wrap;margin-top:18px}}.pill{{border:1px solid rgba(53,208,255,.28);padding:8px 11px;border-radius:999px;background:rgba(53,208,255,.08);font-size:.88rem}}details.guide,.overview{{margin:11px 0 2px;border:1px solid rgba(53,208,255,.22);border-radius:13px;background:rgba(5,10,23,.68)}}details.guide summary,.overview summary{{cursor:pointer;padding:10px 12px;color:var(--cyan);font-weight:800;font-size:.88rem}}details.guide p,.overview p{{margin:0;padding:0 12px 11px;color:var(--muted);font-size:.87rem}}.overview{{margin-top:22px;max-width:1020px}}.overview ol{{margin:0;padding:0 28px 15px;color:var(--muted);line-height:1.55}}.gallery-section,.ligand-card{{border-radius:21px;padding:20px;margin-top:18px}}.section-head{{border-left:3px solid var(--green);padding-left:13px;margin-bottom:16px}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:16px}}.card{{border-radius:17px;padding:13px;overflow:hidden}}.card .js-plotly-plot{{cursor:zoom-in}}.card a{{display:inline-block;color:var(--cyan);font-size:.8rem;font-weight:800;margin:8px 10px 5px}}.heatmap-note{{margin:7px 10px;color:var(--muted);font-size:.82rem;line-height:1.4}}.ligand-grid{{display:grid;grid-template-columns:minmax(300px,1fr) minmax(300px,1fr);gap:18px;align-items:stretch}}.ligand-2d,.ligand-3d{{min-height:330px;border:1px solid rgba(255,230,109,.30);border-radius:17px;background:radial-gradient(circle at center,rgba(255,230,109,.08),rgba(5,6,10,.7) 70%);overflow:hidden;position:relative}}.ligand-2d svg{{display:block;width:100%;height:330px;max-width:100%}}.ligand-3d{{border-color:rgba(53,208,255,.30)}}.viewer-label{{color:var(--yellow);font-size:.78rem;font-weight:900;text-transform:uppercase;letter-spacing:.1em;margin:2px 0 8px}}.missing-visual{{padding:24px;color:var(--muted)}}.asset-note{{font-size:.8rem;margin:.7rem 0 0}}code{{display:block;white-space:pre-wrap;overflow-wrap:anywhere;padding:0 12px 12px;color:var(--yellow)}}.graph-modal{{display:none;position:fixed;z-index:9999;inset:0;padding:3vh 3vw;background:rgba(0,0,0,.78);backdrop-filter:blur(7px)}}.graph-modal.open{{display:flex;align-items:center;justify-content:center}}.graph-modal-panel{{position:relative;width:min(1500px,94vw);height:min(900px,92vh);padding:20px 16px 12px;border:1px solid var(--cyan);background:#080d1b;box-shadow:0 0 70px rgba(53,208,255,.18);border-radius:20px}}.graph-modal-close{{position:absolute;right:17px;top:13px;z-index:2;border:1px solid var(--cyan);border-radius:999px;background:#101b33;color:var(--text);font-weight:800;padding:8px 12px;cursor:pointer}}.graph-modal-plot{{width:100%;height:100%}}@media(max-width:700px){{.page{{padding:18px 12px 42px}}.grid,.ligand-grid{{grid-template-columns:1fr}}.hero{{padding:24px 19px}}.graph-modal{{padding:1vh 1vw}}.graph-modal-panel{{width:98vw;height:96vh;padding:14px 4px 4px}}}}
</style></head><body><main class="page"><section class="hero"><div class="eyebrow">PyMACS · single-system results</div><h1>{escape(system_name)}</h1><p class="sub">A guided reading room for this simulation. Start with global stability, then inspect ligand pose and pocket chemistry. The raw data remain downloadable, but every published result is visible here first.</p><div class="pills"><span class="pill">Ligand / PROTAC: {escape(ligand_name)}</span><span class="pill">{len(cards)} interactive analyses</span><span class="pill">Kyle-local, fast delivery</span></div><details class="overview"><summary>ⓘ How to read this result—start here</summary><ol><li><b>Check the fold:</b> protein RMSD, radius of gyration, and secondary-structure evolution should be considered together.</li><li><b>Check the pose:</b> interpret ligand RMSD/RMSF with persistent contacts—not in isolation.</li><li><b>Check the chemistry:</b> solid-color interaction types show which residues make hydrophobic, ionic, or hydrogen-bond interactions over time.</li><li><b>Do not overclaim:</b> these trajectories test stability of a modeled system; they do not by themselves establish binding affinity or experimental activity.</li></ol></details></section>{galleries}</main><div class="graph-modal" id="graph-modal" aria-hidden="true"><div class="graph-modal-panel"><button class="graph-modal-close" id="graph-modal-close" type="button">Close ×</button><div class="graph-modal-plot" id="graph-modal-plot"></div></div></div><script>
(() => {{
 const modal=document.getElementById('graph-modal'),destination=document.getElementById('graph-modal-plot');
 const close=()=>{{modal.classList.remove('open');modal.setAttribute('aria-hidden','true');Plotly.purge(destination)}};
 const open=source=>{{const layout=Object.assign({{}},source.layout,{{autosize:true,height:null,margin:{{l:70,r:35,t:65,b:60}}}});Plotly.react(destination,source.data,layout,{{responsive:true,displaylogo:false}});modal.classList.add('open');modal.setAttribute('aria-hidden','false');requestAnimationFrame(()=>Plotly.Plots.resize(destination))}};
 document.querySelectorAll('.js-plotly-plot').forEach(source=>source.addEventListener('dblclick',event=>{{event.preventDefault();event.stopImmediatePropagation();open(source)}},true));
 document.getElementById('graph-modal-close').addEventListener('click',close);modal.addEventListener('click',event=>{{if(event.target===modal)close()}});document.addEventListener('keydown',event=>{{if(event.key==='Escape')close()}});
 document.querySelectorAll('.js-ligand-3d').forEach(el=>{{try{{const raw=atob(el.dataset.structureB64),viewer=$3Dmol.createViewer(el,{{backgroundColor:'#05060a',antialias:true}});viewer.addModel(raw,el.dataset.format);viewer.setStyle({{}},{{stick:{{radius:.105,color:'#ffe66d'}},sphere:{{scale:.16,color:'#ffe66d'}}}});viewer.setStyle({{elem:'N'}},{{stick:{{radius:.105,color:'#35d0ff'}},sphere:{{scale:.16,color:'#35d0ff'}}}});viewer.setStyle({{elem:'O'}},{{stick:{{radius:.105,color:'#ff3b5c'}},sphere:{{scale:.16,color:'#ff3b5c'}}}});viewer.setStyle({{elem:'Cl'}},{{stick:{{radius:.105,color:'#18ff8b'}},sphere:{{scale:.20,color:'#18ff8b'}}}});viewer.zoomTo();viewer.render();addEventListener('resize',()=>viewer.resize())}}catch(_){{el.textContent='3D ligand preview unavailable.'}}}});
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
    setup = {}
    try:
        setup = json.loads((root / "pymacs_run.json").read_text()).get("setup", {})
    except (OSError, ValueError, TypeError):
        pass
    ligand_name = str(metadata.get("ligand_name") or metadata.get("ligand_residue") or setup.get("ligand") or "No ligand label")
    csv_files = sorted(path for path in analysis.rglob("*.csv") if path.is_file())
    colours = [CYAN, GREEN, MAGENTA, ORANGE]
    heatmap_cards = []
    groups: dict[str, list[str]] = {"stability": [], "ligand": [], "protein": [], "other": []}
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
    for index, csv_file in enumerate(report_csv_files(csv_files, setup)):
        graph = numeric_plot(csv_file, colours[index % len(colours)])
        if graph:
            card = dashboard_card(graph, csv_file.stem, csv_file.relative_to(output.parent).as_posix())
            stem = csv_file.stem.lower()
            if "ligand" in stem or "contact" in stem or "interaction" in stem:
                groups["ligand"].append(card)
            elif "rmsd" in stem or "rmsf" in stem or "gyration" in stem or "sse" in stem or "dssp" in stem:
                groups["stability"].append(card)
            elif "protein" in stem or "residue" in stem:
                groups["protein"].append(card)
            else:
                groups["other"].append(card)
    # Matrix outputs are the fastest way to understand whether a ligand keeps
    # the starting pocket chemistry and whether the protein retains its fold,
    # so they deliberately lead the report instead of being buried in exports.
    groups["stability"] = heatmap_cards[:1] + groups["stability"]
    groups["ligand"] = heatmap_cards[1:] + groups["ligand"]
    cards = [card for values in groups.values() for card in values]
    if not cards:
        cards.append('<article class="card"><h2>Analysis data is still being written</h2><p>No tabular analysis output was available when this report was generated. Re-run this report after 3A finishes.</p></article>')
    # Publish the chemistry drawing separately as well as embedding it in the
    # self-contained report.  The Kyle results dossier can then show the same
    # fast, recognisable ligand preview before a scientist opens the full HTML.
    ligand_structure, ligand_smiles = ligand_inputs(root, ligand_name)
    published_ligand_svg = ligand_svg(ligand_structure, ligand_smiles)
    if published_ligand_svg:
        (output.parent / "Ligand_2D.svg").write_text(published_ligand_svg, encoding="utf-8")
    ligand_panel = ligand_html(root, ligand_name)
    galleries = (
        section_html("First look: fold and stability", "Start here", "Read the global structure before drawing conclusions about the ligand.", groups["stability"])
        + ligand_panel
        + section_html("Ligand pose and pocket chemistry", "Binding-site evidence", "Use these together: a pose can move while preserving key interactions, and contacts alone do not establish affinity.", groups["ligand"])
        + section_html("Residue-level protein behavior", "Protein detail", "Inspect local flexibility after checking global stability and the ligand interaction fingerprint.", groups["protein"])
        + section_html("Additional analysis", "Supporting outputs", "These outputs add context and should be interpreted together with the main stability and interaction panels.", groups["other"])
    )
    html = build_document(system_name, ligand_name, cards, galleries)
    output.write_text(html, encoding="utf-8")
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
