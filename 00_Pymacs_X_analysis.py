#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Fast PyMACS comparative dashboard builder.

Purpose
-------
Build the most important apo-vs-ligand and ligand-vs-ligand figures/tables from
existing PyMACS / 3A_analysis / 3B-NETWORX outputs without loading huge
framewise contact files by default.

Normal fast command:
    python 00_Pymacs_X_analysis.py --root . --family HUMAN --open
    python 00_Pymacs_X_analysis.py --root . --family RAT --open
    python 00_Pymacs_X_analysis.py --root . --family MOUSE --open

    python 00_Pymacs_X_analysis_FAST.py --root . --family ####ANY SPECIES####### --open

Optional deep fallback, only if you need to summarize framewise contacts:
    python 00_Pymacs_X_analysis_FAST.py --root . --family HUMAN --use-framewise --open

Expected layout:
    HUMAN_APO/Analysis_Results/
    HUMAN_DRUG1/Analysis_Results/
    HUMAN_DRUG2/Analysis_Results/
    (Any Species Name works; just be consistent.)

Core outputs:
    - Ligand structure cards from .cgenff.mol2/.mol2/.sdf files
    - Protein RMSD overlay
    - Radius of gyration overlay
    - Protein RMSF overlay
    - Apo-relative protein-state deltas
    - Binding-site/contact-residue RMSF delta vs apo
    - Ligand RMSD overlay
    - Contact residue summary heatmap, if summary exists
    - Interaction type composition, if summary exists
    - Ligand-vs-ligand comparison table
    - Ligand stability flags
    - CSV exports
    - Polished HTML dashboard with favicon/header, hover explainers, and export links
"""

from __future__ import annotations

import argparse
import base64
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import webbrowser
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.offline import get_plotlyjs

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors
    RDKIT_AVAILABLE = True
except Exception:
    RDKIT_AVAILABLE = False


TRON_BG = "#05060a"
TRON_PANEL = "#0b1020"
TRON_GRID = "#1c2a4a"
TRON_TEXT = "#cfe8ff"
TRON_MUTED = "#7aa6d9"
TRON_CYAN = "#35d0ff"
TRON_BLUE = "#2b6cff"
TRON_MAGENTA = "#ff2bd6"
TRON_PURPLE = "#8a5cff"
TRON_GREEN = "#18ff8b"
TRON_YELLOW = "#ffe66d"
TRON_ORANGE = "#ff8a3d"
TRON_RED = "#ff3b5c"
TRON_WHITE = "#f5fbff"

# Ligand SVG atom colors for dark-mode visibility.
# RDKit's default atom palette can disappear on the dark dashboard background,
# so we force a high-contrast palette for common medicinal-chemistry atoms.
LIGAND_CARBON_YELLOW = "#ffe66d"
LIGAND_HYDROGEN_SLATE = "#9fb7d5"
LIGAND_NITROGEN_CYAN = "#35d0ff"
LIGAND_OXYGEN_RED = "#ff2b2b"
LIGAND_FLUORINE_MINT = "#7cff9b"
LIGAND_CHLORINE_GREEN = "#00ff66"
LIGAND_BROMINE_ORANGE = "#ff8a3d"
LIGAND_IODINE_PURPLE = "#d36bff"
LIGAND_PHOSPHORUS_ORANGE = "#ffb347"
LIGAND_SULFUR_GOLD = "#ffd84d"

# RDKit commonly uses these defaults; they can be hard to see on dark UI.
RDKIT_DEFAULT_N_BLUE = "#3050F8"
RDKIT_DEFAULT_O_RED = "#FF0D0D"
RDKIT_DEFAULT_HALOGEN_GREEN = "#00CC00"
RDKIT_DEFAULT_BR_BROWN = "#A62929"
RDKIT_DEFAULT_I_PURPLE = "#940094"
RDKIT_DEFAULT_P_ORANGE = "#FF8000"
RDKIT_DEFAULT_S_GOLD = "#CCCC00"

APO_COLOR = TRON_WHITE
COLOR_CYCLE = [
    TRON_CYAN,
    TRON_GREEN,
    TRON_YELLOW,
    TRON_BLUE,
    TRON_ORANGE,
    TRON_PURPLE,
    "#00f5ff",
    "#00ffa2",
    "#ffd166",
    "#5dade2",
]
CROSS_SPECIES_TRACE_COLORS = [
    "#35D0FF",
    "#FF2BD6",
    "#18FF8B",
    "#8A5CFF",
    "#FF8A3D",
    "#2B6CFF",
    "#FF3B5C",
    "#00FFA2",
    "#D36BFF",
    "#FFE66D",
    "#B6FF00",
    "#FF6EC7",
    "#00B7FF",
    "#C77DFF",
    "#7CFFCB",
    "#FF5E00",
]
TYPE_COLORS = {
    "H-bond": TRON_YELLOW,
    "Hydrophobic": TRON_BLUE,
    "Ionic": TRON_GREEN,
    "Water bridge": "#4fb4ff",
    "Polar": "#b6f7ff",
    "Van der Waals": TRON_PURPLE,
    "Other": "#9aa8c7",
}
PLOTLY_CONFIG = {"responsive": True, "displaylogo": False, "scrollZoom": True}

PYMACS_LOGO_URL = "https://github.com/schurerlab/Pymacs/blob/main/images/Pymacs_Logo.png?raw=true"
PYMACS_LOGO_OUTPUT_REL = Path("images") / "favicon.png"


AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HSD", "HSE", "HSP", "HIP", "HID", "HIE",
}


def phase(msg: str) -> None:
    print("\n" + "=" * 92, flush=True)
    print(f"🧩 {msg}", flush=True)
    print("=" * 92, flush=True)


def info(msg: str) -> None:
    print(f"   • {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"⚠️  {msg}", flush=True)


def die(msg: str, code: int = 1) -> None:
    print(f"❌ {msg}", flush=True)
    raise SystemExit(code)


def html_escape(x: Any) -> str:
    import html
    return html.escape("" if x is None else str(x))


def slugify(s: str) -> str:
    return re.sub(r"[^a-zA-Z0-9._-]+", "_", str(s)).strip("_") or "item"


def nice(x: Any, nd: int = 3) -> str:
    try:
        if x is None:
            return ""
        fx = float(x)
        if math.isnan(fx) or math.isinf(fx):
            return ""
        return f"{fx:.{nd}f}"
    except Exception:
        return str(x)


def rel_href(from_file: Path, to_file: Path) -> str:
    return os.path.relpath(to_file, start=from_file.parent).replace("\\", "/")


def first_existing(paths: Iterable[Path]) -> Optional[Path]:
    for p in paths:
        if p and p.exists() and p.is_file():
            return p
    return None


def file_size_mb(path: Optional[Path]) -> Optional[float]:
    if path is None or not path.exists():
        return None
    try:
        return path.stat().st_size / (1024 * 1024)
    except Exception:
        return None


def safe_read_csv(path: Optional[Path], **kwargs) -> Optional[pd.DataFrame]:
    if path is None or not path.exists() or not path.is_file():
        return None
    try:
        return pd.read_csv(path, **kwargs)
    except Exception as e:
        warn(f"Could not read CSV {path}: {e}")
        return None


def read_header(path: Optional[Path]) -> List[str]:
    if path is None or not path.exists():
        return []
    try:
        return list(pd.read_csv(path, nrows=0).columns)
    except Exception:
        return []


def find_first_col_from_names(names: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    if not names:
        return None
    exact = {c: c for c in names}
    for cand in candidates:
        if cand in exact:
            return exact[cand]
    lower = {c.lower(): c for c in names}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    for cand in candidates:
        lc = cand.lower().replace("_", "").replace(" ", "").replace("-", "")
        for col in names:
            cc = col.lower().replace("_", "").replace(" ", "").replace("-", "")
            if lc == cc or lc in cc or cc in lc:
                return col
    return None


def find_first_col(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    if df is None or df.empty:
        return None
    return find_first_col_from_names(list(df.columns), candidates)


def detect_time_ns(df: pd.DataFrame) -> Optional[pd.Series]:
    if df is None or df.empty:
        return None
    for col in ["Time_ns", "time_ns", "t_ns", "Time (ns)", "Time(ns)", "ns"]:
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce")
    for col in ["Time_ps", "time_ps", "t_ps", "Time (ps)", "Time(ps)", "ps"]:
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce") / 1000.0
    for col in ["Time", "time"]:
        if col in df.columns:
            vals = pd.to_numeric(df[col], errors="coerce")
            if vals.dropna().empty:
                continue
            return vals / 1000.0 if vals.dropna().max() > 1000 else vals
    for col in ["Frame", "frame", "Frame_ID", "frame_id", "Frames", "frames"]:
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce")
    return None


def numeric_mean(series: Optional[pd.Series]) -> Optional[float]:
    if series is None:
        return None
    vals = pd.to_numeric(series, errors="coerce").dropna()
    return None if vals.empty else float(vals.mean())


def numeric_std(series: Optional[pd.Series]) -> Optional[float]:
    if series is None:
        return None
    vals = pd.to_numeric(series, errors="coerce").dropna()
    return None if vals.empty else float(vals.std())


def numeric_last(series: Optional[pd.Series]) -> Optional[float]:
    if series is None:
        return None
    vals = pd.to_numeric(series, errors="coerce").dropna()
    return None if vals.empty else float(vals.iloc[-1])


def numeric_last_half_mean(series: Optional[pd.Series]) -> Optional[float]:
    if series is None:
        return None
    vals = pd.to_numeric(series, errors="coerce").dropna()
    if vals.empty:
        return None
    return float(vals.iloc[len(vals)//2:].mean())


def linear_slope(x: Optional[pd.Series], y: Optional[pd.Series]) -> Optional[float]:
    if x is None or y is None:
        return None
    tmp = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if len(tmp) < 3:
        return None
    try:
        slope, _ = np.polyfit(tmp["x"].values.astype(float), tmp["y"].values.astype(float), 1)
        return float(slope)
    except Exception:
        return None


def safe_delta(bound: Any, apo: Any) -> Optional[float]:
    try:
        if bound is None or apo is None:
            return None
        return float(bound) - float(apo)
    except Exception:
        return None


def cosine_similarity_matrix(mat: pd.DataFrame) -> pd.DataFrame:
    if mat is None or mat.empty:
        return pd.DataFrame()
    arr = mat.fillna(0.0).values.astype(float)
    norms = np.linalg.norm(arr, axis=1)
    out = np.zeros((arr.shape[0], arr.shape[0]), dtype=float)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[0]):
            denom = norms[i] * norms[j]
            out[i, j] = float(np.dot(arr[i], arr[j]) / denom) if denom else 0.0
    return pd.DataFrame(out, index=mat.index, columns=mat.index)


def normalize_interaction_type(value: Any) -> str:
    s = str(value).strip() if value is not None else ""
    sl = s.lower()
    if any(k in sl for k in ["hbond", "hydrogen", "h-bond", "hb"]):
        return "H-bond"
    if any(k in sl for k in ["hydroph", "apolar", "nonpolar", "alkyl", "pi-alkyl"]):
        return "Hydrophobic"
    if any(k in sl for k in ["ionic", "salt", "electro", "charge", "cation", "anion"]):
        return "Ionic"
    if any(k in sl for k in ["water", "bridge"]):
        return "Water bridge"
    if "polar" in sl:
        return "Polar"
    if any(k in sl for k in ["vdw", "van der", "proximal"]):
        return "Van der Waals"
    return s if s else "Other"


def is_excluded_interaction_type(value: Any) -> bool:
    s = str(value).strip().lower() if value is not None else ""
    return s in {"total", "all", "sum", "overall", "grand_total", "grand total"}


def build_run_palette(runs: List["RunData"], apo_name: Optional[str] = None) -> Dict[str, str]:
    palette: Dict[str, str] = {}
    idx = 0
    for r in runs:
        if apo_name and r.name == apo_name:
            palette[r.name] = APO_COLOR
        else:
            palette[r.name] = COLOR_CYCLE[idx % len(COLOR_CYCLE)]
            idx += 1
    return palette


def build_cross_species_palette(runs: List["RunData"]) -> Dict[str, str]:
    ordered = sorted(
        runs,
        key=lambda r: (
            str(r.family),
            0 if r.is_apo else 1,
            str(r.ligand_display or r.ligand_code or r.state or r.name),
            str(r.name),
        ),
    )
    palette: Dict[str, str] = {}
    for i, run in enumerate(ordered):
        palette[run.name] = CROSS_SPECIES_TRACE_COLORS[i % len(CROSS_SPECIES_TRACE_COLORS)]
    return palette


def make_vertical_line_trace(xs: List[float], y0: float, y1: float, color: str, name: str) -> go.Scatter:
    xvals: List[Optional[float]] = []
    yvals: List[Optional[float]] = []
    for x in xs:
        xvals.extend([x, x, None])
        yvals.extend([y0, y1, None])
    return go.Scatter(
        x=xvals,
        y=yvals,
        mode="lines",
        name=name,
        line=dict(color=color, width=1, dash="dot"),
        opacity=0.75,
        visible="legendonly",
        hoverinfo="skip",
    )


def cross_species_legend_layout() -> Dict[str, Any]:
    return dict(
        bgcolor="rgba(0,0,0,0)",
        borderwidth=0,
        font=dict(color=TRON_MUTED),
        groupclick="toggleitem",
        itemclick="toggle",
        itemdoubleclick="toggleothers",
    )


def parse_residue_sort_tuple_from_row(row: pd.Series) -> Tuple[str, int, str]:
    chain = str(row.get("Chain") or "")
    try:
        resnum = int(float(row.get("ResNum"))) if pd.notna(row.get("ResNum")) else 10**9
    except Exception:
        resnum = 10**9
    residue = str(row.get("Residue") or "")
    return (chain, resnum, residue)


def parse_residue_label(label: Any) -> Dict[str, Any]:
    raw = "" if label is None else str(label).strip()
    raw = raw.replace("_", ":")
    raw = re.sub(r"\s+", " ", raw)
    chain = ""
    resname = ""
    resnum: Optional[int] = None
    m = re.search(r"(?:(?P<chain>[A-Za-z0-9]+)[:/])?\s*(?P<resname>[A-Za-z]{3,4})\s*[-:]?\s*(?P<num>-?\d+)", raw)
    if m:
        chain = m.group("chain") or ""
        resname = m.group("resname").upper()
        try:
            resnum = int(m.group("num"))
        except Exception:
            resnum = None
    else:
        nums = re.findall(r"-?\d+", raw)
        if nums:
            try:
                resnum = int(nums[-1])
            except Exception:
                resnum = None
        for token in re.findall(r"[A-Za-z]{3,4}", raw):
            if token.upper() in AA3:
                resname = token.upper()
                break
        if ":" in raw:
            first = raw.split(":", 1)[0].strip()
            if first and len(first) <= 3:
                chain = first
    key = raw.upper() if resnum is None else (f"{chain}:{resnum}" if chain else str(resnum))
    label2 = raw if resnum is None else f"{chain + ':' if chain else ''}{resname if resname else 'RES'}{resnum}"
    return {"raw": raw, "chain": chain, "resname": resname, "resnum": resnum, "key": key, "label": label2}


def residue_sort_key(label: Any) -> Tuple[str, int, str]:
    p = parse_residue_label(label)
    return (p["chain"] or "", p["resnum"] if p["resnum"] is not None else 10**9, str(label))


def residue_key_from_values(residue_value: Any = None, chain_value: Any = None, resname_value: Any = None, resnum_value: Any = None) -> Tuple[str, str, str, Optional[int]]:
    chain = "" if chain_value is None or pd.isna(chain_value) else str(chain_value).strip()
    resname = "" if resname_value is None or pd.isna(resname_value) else str(resname_value).strip().upper()
    resnum: Optional[int] = None
    if resnum_value is not None and not pd.isna(resnum_value):
        try:
            resnum = int(float(resnum_value))
        except Exception:
            resnum = None
    if residue_value is not None and not pd.isna(residue_value):
        p = parse_residue_label(residue_value)
        chain = chain or p["chain"]
        resname = resname or p["resname"]
        resnum = resnum if resnum is not None else p["resnum"]
    if resnum is None:
        raw = "" if residue_value is None or pd.isna(residue_value) else str(residue_value)
        return raw.upper(), raw, chain, resnum
    key = f"{chain}:{resnum}" if chain else str(resnum)
    label = f"{chain + ':' if chain else ''}{resname if resname else 'RES'}{resnum}"
    return key, label, chain, resnum


@dataclass
class LigandVisual:
    structure_file: Optional[Path] = None
    svg: Optional[str] = None
    formula: Optional[str] = None
    mw: Optional[float] = None
    smiles: Optional[str] = None
    note: Optional[str] = None
    structure_text_b64: Optional[str] = None
    structure_format: Optional[str] = None
    depiction_file: Optional[Path] = None


@dataclass
class RunData:
    name: str
    run_dir: Path
    analysis_dir: Path
    family: str
    state: str
    is_apo: bool
    ligand_code: Optional[str] = None
    ligand_display: Optional[str] = None
    protein_rmsd: Optional[pd.DataFrame] = None
    protein_rmsf: Optional[pd.DataFrame] = None
    radius_gyration: Optional[pd.DataFrame] = None
    ligand_rmsd: Optional[pd.DataFrame] = None
    ligand_rmsf: Optional[pd.DataFrame] = None
    bound_complex_rmsd: Optional[pd.DataFrame] = None
    residue_frequency: Optional[pd.DataFrame] = None
    type_fraction: Optional[pd.DataFrame] = None
    residue_hover_summary: Optional[pd.DataFrame] = None
    residue_interaction_summary: Optional[pd.DataFrame] = None
    ligand_visual: LigandVisual = field(default_factory=LigandVisual)
    files: Dict[str, Optional[str]] = field(default_factory=dict)
    metrics: Dict[str, Any] = field(default_factory=dict)


def infer_family_and_state(run_name: str) -> Tuple[str, str]:
    parts = [p for p in run_name.split("_") if p]
    if not parts:
        return run_name, run_name
    return parts[0], ("_".join(parts[1:]) if len(parts) > 1 else run_name)


def find_run_dirs(root: Path) -> List[Path]:
    direct = [p for p in sorted(root.iterdir()) if p.is_dir() and (p / "Analysis_Results").is_dir()]
    if direct:
        return direct
    found = [ar.parent for ar in sorted(root.rglob("Analysis_Results")) if ar.is_dir()]
    seen, uniq = set(), []
    for p in found:
        rp = str(p.resolve())
        if rp not in seen:
            uniq.append(p)
            seen.add(rp)
    return uniq


def glob_many(base: Path, patterns: Sequence[str]) -> List[Path]:
    out: List[Path] = []
    if not base.exists():
        return out
    for pat in patterns:
        out.extend(sorted(base.glob(pat)))
    seen, uniq = set(), []
    for p in out:
        rp = str(p.resolve())
        if p.exists() and p.is_file() and rp not in seen:
            uniq.append(p)
            seen.add(rp)
    return uniq


def find_ligand_rmsd_path(analysis_dir: Path) -> Optional[Path]:
    return first_existing(glob_many(analysis_dir, ["*_Ligand_RMSD.csv", "*Ligand*RMSD*.csv", "*ligand*rmsd*.csv"]))


def extract_ligand_code_from_filename(path: Optional[Path]) -> Optional[str]:
    if path is None:
        return None
    stem = path.stem
    for pat in [r"(.+?)_Ligand_RMSD$", r"(.+?)_ligand_rmsd$", r"(.+?)_RMSD$", r"(.+?)_rmsd$"]:
        m = re.match(pat, stem, flags=re.IGNORECASE)
        if m:
            return m.group(1)
    return None


def find_ligand_rmsf_path(analysis_dir: Path, ligand_code: Optional[str]) -> Optional[Path]:
    candidates: List[Path] = []
    if ligand_code:
        candidates.extend(glob_many(analysis_dir, [f"{ligand_code}_rmsf.csv", f"{ligand_code}_RMSF.csv", f"{ligand_code}*rmsf*.csv", f"{ligand_code}*RMSF*.csv"]))
    candidates.extend(glob_many(analysis_dir, ["*_ligand_rmsf.csv", "*Ligand_RMSF.csv", "*ligand*RMSF*.csv", "*ligand*rmsf*.csv"]))
    filtered = [p for p in candidates if "protein" not in p.name.lower() and "backbone" not in p.name.lower() and not p.name.lower().startswith("rmsf_")]
    return first_existing(filtered)


def find_protein_rmsf_paths(analysis_dir: Path, ligand_code: Optional[str]) -> List[Path]:
    exact = glob_many(analysis_dir, ["Protein_RMSF.csv", "protein_rmsf.csv", "Backbone_RMSF.csv", "backbone_rmsf.csv", "rmsf_per_residue.csv", "RMSF_Protein.csv", "GlobalProtein_RMSF.csv"])
    if exact:
        return exact
    candidates = glob_many(analysis_dir, ["RMSF_*.csv", "*Protein*RMSF*.csv", "*protein*rmsf*.csv", "*Backbone*RMSF*.csv", "*backbone*rmsf*.csv"])
    lig = ligand_code.lower() if ligand_code else ""
    return [p for p in candidates if not (lig and lig in p.name.lower()) and "ligand" not in p.name.lower() and "complex" not in p.name.lower()]


def find_residue_summary_path(analysis_dir: Path) -> Optional[Path]:
    patterns = [
        "Residue_Contact_Frequency.csv", "residue_contact_frequency.csv", "ResidueContactFrequency.csv",
        "Contact_Frequency.csv", "contact_frequency.csv", "ContactFrequency.csv", "ResidueContacts_Summary.csv",
        "ResidueContact_Summary.csv", "ContactFrequency_Summary.csv", "FilteredContacts_Summary.csv",
        "*residue*contact*frequency*.csv", "*Residue*Contact*Frequency*.csv", "*contact*frequency*.csv",
        "*Contact*Frequency*.csv", "*residue*summary*.csv", "*Residue*Summary*.csv",
    ]
    candidates = glob_many(analysis_dir, patterns)
    net = analysis_dir / "NETWORX"
    if net.exists():
        candidates.extend(glob_many(net, patterns + ["*summary*.csv"]))
    return first_existing([p for p in candidates if "framewise" not in p.name.lower()])


def find_type_summary_path(analysis_dir: Path) -> Optional[Path]:
    patterns = ["InteractionTypes_Summary.csv", "interaction_types_summary.csv", "Interaction_Type_Summary.csv", "interaction_type_summary.csv", "InteractionTypes.csv", "interaction_types.csv", "*Interaction*Summary*.csv", "*interaction*summary*.csv"]
    candidates = glob_many(analysis_dir, patterns)
    net = analysis_dir / "NETWORX"
    if net.exists():
        candidates.extend(glob_many(net, patterns))
    return first_existing([p for p in candidates if "framewise" not in p.name.lower()])


def find_framewise_interaction_path(analysis_dir: Path) -> Optional[Path]:
    patterns = ["FilteredContacts_Framewise.csv", "InteractionTypes_Framewise.csv", "Contacts_Framewise.csv", "ResidueContacts_Framewise.csv", "AllContacts_Framewise.csv", "*contact*framewise*.csv", "*Contact*Framewise*.csv", "*interaction*framewise*.csv", "*Interaction*Framewise*.csv"]
    candidates = glob_many(analysis_dir, patterns)
    net = analysis_dir / "NETWORX"
    if net.exists():
        candidates.extend(glob_many(net, patterns))
    return first_existing(candidates)


def find_contact_framewise_path(analysis_dir: Path) -> Optional[Path]:
    patterns = ["FilteredContacts_Framewise.csv", "Contacts_Framewise.csv", "ResidueContacts_Framewise.csv", "AllContacts_Framewise.csv", "*Filtered*Framewise*.csv", "*Contact*Framewise*.csv", "*contact*framewise*.csv"]
    candidates = glob_many(analysis_dir, patterns)
    net = analysis_dir / "NETWORX"
    if net.exists():
        candidates.extend(glob_many(net, patterns))
    filtered = [p for p in candidates if "interactiontypes" not in p.name.lower()]
    return first_existing(filtered)


def find_interaction_types_framewise_path(analysis_dir: Path) -> Optional[Path]:
    patterns = ["InteractionTypes_Framewise.csv", "*InteractionTypes*Framewise*.csv", "*interactiontypes*framewise*.csv", "*interaction*framewise*.csv"]
    candidates = glob_many(analysis_dir, patterns)
    net = analysis_dir / "NETWORX"
    if net.exists():
        candidates.extend(glob_many(net, patterns))
    return first_existing(candidates)


# =============================================================================
# Ligand structure discovery and 2D rendering
# =============================================================================

def obabel_mol2_to_sdf_temp(mol2_path: Path) -> Optional[Path]:
    """Fallback conversion for MOL2 files that RDKit cannot parse directly."""
    obabel = shutil.which("obabel")
    if not obabel:
        return None

    tmp_dir = Path(tempfile.mkdtemp(prefix="pymacs_lig2d_"))
    out_sdf = tmp_dir / f"{mol2_path.stem}.sdf"
    try:
        result = subprocess.run(
            [obabel, str(mol2_path), "-O", str(out_sdf)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode == 0 and out_sdf.exists() and out_sdf.stat().st_size > 0:
            return out_sdf
    except Exception:
        return None
    return None


def find_ligand_structure_file(run: RunData) -> Optional[Path]:
    """
    Find the ligand structure file without touching heavy trajectory/contact data.

    Priority:
      1. <ligand>.cgenff.mol2 in the run folder
      2. any *.cgenff.mol2 in the run folder
      3. <ligand>.mol2 / any *.mol2
      4. SDF/MOL/SMI files in common ligand subfolders
    """
    if run.is_apo:
        return None

    ligand_tokens = []
    for token in [run.ligand_code, run.ligand_display, run.state]:
        if token and str(token).upper() != "APO":
            ligand_tokens.append(str(token))

    search_dirs = [
        run.run_dir,
        run.analysis_dir,
        run.run_dir / "SDF",
        run.run_dir / "Ligands",
        run.run_dir / "ligands",
        run.run_dir / "structures",
        run.run_dir / "Structures",
        run.analysis_dir / "NETWORX",
    ]

    priority_patterns: List[str] = []
    for token in ligand_tokens:
        priority_patterns.extend([
            f"{token}.cgenff.mol2",
            f"{token.upper()}.cgenff.mol2",
            f"{token.lower()}.cgenff.mol2",
            f"{token}*.cgenff.mol2",
            f"*{token}*.cgenff.mol2",
            f"{token}.mol2",
            f"{token.upper()}.mol2",
            f"{token.lower()}.mol2",
            f"*{token}*.mol2",
            f"*{token}*.sdf",
            f"*{token}*.mol",
            f"*{token}*.smi",
            f"*{token}*.smiles",
        ])

    fallback_patterns = [
        "*.cgenff.mol2",
        "*.mol2",
        "*.sdf",
        "*.mol",
        "*.smi",
        "*.smiles",
        "*.pdb",
    ]

    for patterns in [priority_patterns, fallback_patterns]:
        for directory in search_dirs:
            if not directory.exists():
                continue
            for pattern in patterns:
                hits = sorted(directory.glob(pattern))
                if hits:
                    return hits[0]
    return None


def read_first_smiles(path: Path) -> Optional[str]:
    try:
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts:
                return parts[0]
    except Exception:
        return None
    return None


def find_matching_perceived_smiles_file(structure_path: Optional[Path]) -> Optional[Path]:
    """
    Given a 3D ligand structure like UNL.cgenff.mol2, prefer UNL.perceived.smi for
    the 2D depiction. This gives a cleaner canonical 2D graph while still using the
    original MOL2 coordinates for the 3D viewer.
    """
    if structure_path is None:
        return None

    parent = structure_path.parent
    name = structure_path.name
    candidates: List[Path] = []

    # Best case: <PREFIX>.cgenff.mol2 -> <PREFIX>.perceived.smi
    if name.lower().endswith('.cgenff.mol2'):
        prefix = name[:-len('.cgenff.mol2')]
        candidates.extend([
            parent / f'{prefix}.perceived.smi',
            parent / f'{prefix}.smi',
            parent / f'{prefix}.smiles',
            parent / f'{prefix}.perceived.smiles',
        ])
    else:
        stem = structure_path.stem
        # Also normalize things like LIG.mol2 -> LIG.perceived.smi
        candidates.extend([
            parent / f'{stem}.perceived.smi',
            parent / f'{stem}.smi',
            parent / f'{stem}.smiles',
            parent / f'{stem}.perceived.smiles',
        ])

    for cand in candidates:
        if cand.exists():
            return cand

    # Broad fallback: any *.perceived.smi in the same folder.
    perceived_hits = sorted(parent.glob('*.perceived.smi'))
    if perceived_hits:
        if len(perceived_hits) == 1:
            return perceived_hits[0]
        # try matching by uppercase token overlap
        structure_tokens = set(re.split(r'[^A-Za-z0-9]+', structure_path.stem.upper()))
        scored = []
        for hit in perceived_hits:
            hit_tokens = set(re.split(r'[^A-Za-z0-9]+', hit.stem.upper()))
            scored.append((len(structure_tokens & hit_tokens), hit.name, hit))
        scored.sort(reverse=True)
        return scored[0][2]

    return None


def rdkit_mol_from_file(path: Path):
    if not RDKIT_AVAILABLE:
        return None

    suffix = path.suffix.lower()
    mol = None

    try:
        if suffix == ".sdf":
            suppl = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
            mol = next((m for m in suppl if m is not None), None)

        elif suffix == ".mol":
            mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=True)

        elif suffix == ".mol2":
            mol = Chem.MolFromMol2File(str(path), removeHs=False, sanitize=True)
            if mol is None:
                mol = Chem.MolFromMol2File(str(path), removeHs=False, sanitize=False)
                if mol is not None:
                    try:
                        Chem.SanitizeMol(
                            mol,
                            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
                        )
                    except Exception:
                        pass

        elif suffix in {".smi", ".smiles"}:
            for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
                line = line.strip()
                if line:
                    mol = Chem.MolFromSmiles(line.split()[0])
                    break

        elif suffix == ".pdb":
            mol = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=True)

    except Exception:
        mol = None

    return mol



def clean_ligand_for_display(mol):
    """
    Display-only cleanup for the 2D card:
      - remove nonpolar C-H hydrogens to reduce clutter
      - keep hydrogens attached to heteroatoms
      - leave the original molecule untouched for formula/MW/3D payload
    """
    if mol is None:
        return None
    try:
        rw = Chem.RWMol(Chem.Mol(mol))
        remove_idxs = []
        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            nbrs = list(atom.GetNeighbors())
            if not nbrs or nbrs[0].GetAtomicNum() == 6:
                remove_idxs.append(atom.GetIdx())
        for idx in sorted(remove_idxs, reverse=True):
            rw.RemoveAtom(idx)
        out = rw.GetMol()
        try:
            Chem.SanitizeMol(out, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        except Exception:
            pass
        return out
    except Exception:
        try:
            return Chem.RemoveHs(mol)
        except Exception:
            return Chem.Mol(mol)


def atom_label_color(atomic_num: int) -> str:
    return {
        1:  LIGAND_HYDROGEN_SLATE,
        6:  LIGAND_CARBON_YELLOW,
        7:  LIGAND_NITROGEN_CYAN,
        8:  LIGAND_OXYGEN_RED,
        9:  LIGAND_FLUORINE_MINT,
        15: LIGAND_PHOSPHORUS_ORANGE,
        16: LIGAND_SULFUR_GOLD,
        17: LIGAND_CHLORINE_GREEN,
        35: LIGAND_BROMINE_ORANGE,
        53: LIGAND_IODINE_PURPLE,
    }.get(int(atomic_num), TRON_TEXT)


def atom_label_text(atom) -> str:
    sym = atom.GetSymbol()
    if sym == 'C':
        return ''
    charge = atom.GetFormalCharge()
    if charge > 0:
        sym += '+' if charge == 1 else f'{charge}+'
    elif charge < 0:
        sym += '-' if charge == -1 else f'{abs(charge)}-'
    return sym


def svg_line(x1, y1, x2, y2, color, width=2.55, extra='') -> str:
    return (
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
        f'stroke="{color}" stroke-width="{width:.2f}" stroke-linecap="round" {extra}/>'
    )


def svg_bond_lines(x1, y1, x2, y2, order, color, is_aromatic=False) -> str:
    """
    Manual 2D bond drawing based on RDKit coordinates. This avoids RDKit's
    terminal-label clipping problem where atoms like Cl can look detached,
    without adding ugly highlight overlays to C-N/C-O bonds.
    """
    dx = x2 - x1
    dy = y2 - y1
    length = math.hypot(dx, dy) or 1.0
    ox = -dy / length * 3.5
    oy = dx / length * 3.5
    if is_aromatic:
        return (
            svg_line(x1, y1, x2, y2, color, width=2.25)
            + svg_line(x1 + ox, y1 + oy, x2 + ox, y2 + oy, color, width=1.10, extra='stroke-dasharray="4 4" opacity="0.72"')
        )
    if order >= 3:
        return (
            svg_line(x1, y1, x2, y2, color, width=2.0)
            + svg_line(x1 + ox, y1 + oy, x2 + ox, y2 + oy, color, width=1.45)
            + svg_line(x1 - ox, y1 - oy, x2 - ox, y2 - oy, color, width=1.45)
        )
    if order >= 2:
        return (
            svg_line(x1 + ox * 0.72, y1 + oy * 0.72, x2 + ox * 0.72, y2 + oy * 0.72, color, width=1.80)
            + svg_line(x1 - ox * 0.72, y1 - oy * 0.72, x2 - ox * 0.72, y2 - oy * 0.72, color, width=1.80)
        )
    return svg_line(x1, y1, x2, y2, color, width=2.45)


def recolor_rdkit_svg(svg: str) -> str:
    """Final SVG color cleanup for RDKit's common default atom colors."""
    replacements = {
        "#000000": LIGAND_CARBON_YELLOW,
        "#000": LIGAND_CARBON_YELLOW,
        "rgb(0,0,0)": LIGAND_CARBON_YELLOW,
        "rgb(0, 0, 0)": LIGAND_CARBON_YELLOW,
        "#3050F8": LIGAND_NITROGEN_CYAN,
        "#3050f8": LIGAND_NITROGEN_CYAN,
        "rgb(48,80,248)": LIGAND_NITROGEN_CYAN,
        "rgb(48, 80, 248)": LIGAND_NITROGEN_CYAN,
        "#0000FF": LIGAND_NITROGEN_CYAN,
        "#0000ff": LIGAND_NITROGEN_CYAN,
        "#FF0D0D": LIGAND_OXYGEN_RED,
        "#ff0d0d": LIGAND_OXYGEN_RED,
        "#FF0000": LIGAND_OXYGEN_RED,
        "#ff0000": LIGAND_OXYGEN_RED,
        "#00CC00": LIGAND_CHLORINE_GREEN,
        "#00cc00": LIGAND_CHLORINE_GREEN,
        "#1FF01F": LIGAND_CHLORINE_GREEN,
        "#1ff01f": LIGAND_CHLORINE_GREEN,
        "#90E050": LIGAND_FLUORINE_MINT,
        "#90e050": LIGAND_FLUORINE_MINT,
        "#A62929": LIGAND_BROMINE_ORANGE,
        "#a62929": LIGAND_BROMINE_ORANGE,
        "#940094": LIGAND_IODINE_PURPLE,
        "#FF8000": LIGAND_PHOSPHORUS_ORANGE,
        "#ff8000": LIGAND_PHOSPHORUS_ORANGE,
        "#CCCC00": LIGAND_SULFUR_GOLD,
        "#cccc00": LIGAND_SULFUR_GOLD,
        "#C6C600": LIGAND_SULFUR_GOLD,
        "#c6c600": LIGAND_SULFUR_GOLD,
    }
    for old_color, new_color in replacements.items():
        svg = svg.replace(f"stroke:{old_color}", f"stroke:{new_color}")
        svg = svg.replace(f"fill:{old_color}", f"fill:{new_color}")
        svg = svg.replace(f"stroke: {old_color}", f"stroke: {new_color}")
        svg = svg.replace(f"fill: {old_color}", f"fill: {new_color}")
        svg = svg.replace(f'stroke="{old_color}"', f'stroke="{new_color}"')
        svg = svg.replace(f'fill="{old_color}"', f'fill="{new_color}"')
        svg = svg.replace(f"stroke='{old_color}'", f"stroke='{new_color}'")
        svg = svg.replace(f"fill='{old_color}'", f"fill='{new_color}'")

    # RDKit sometimes emits style strings with arbitrary whitespace.
    regex_replacements = [
        (r"#000000|rgb\(0,\s*0,\s*0\)", LIGAND_CARBON_YELLOW),
        (r"#3050f8|#0000ff|rgb\(48,\s*80,\s*248\)|rgb\(0,\s*0,\s*255\)", LIGAND_NITROGEN_CYAN),
        (r"#ff0d0d|#ff0000|rgb\(255,\s*13,\s*13\)|rgb\(255,\s*0,\s*0\)", LIGAND_OXYGEN_RED),
        (r"#90e050|rgb\(144,\s*224,\s*80\)", LIGAND_FLUORINE_MINT),
        (r"#1ff01f|#00cc00|rgb\(31,\s*240,\s*31\)|rgb\(0,\s*204,\s*0\)", LIGAND_CHLORINE_GREEN),
        (r"#a62929|rgb\(166,\s*41,\s*41\)", LIGAND_BROMINE_ORANGE),
        (r"#940094|rgb\(148,\s*0,\s*148\)", LIGAND_IODINE_PURPLE),
        (r"#ff8000|rgb\(255,\s*128,\s*0\)", LIGAND_PHOSPHORUS_ORANGE),
        (r"#cccc00|#c6c600|rgb\(204,\s*204,\s*0\)|rgb\(198,\s*198,\s*0\)", LIGAND_SULFUR_GOLD),
    ]
    for pattern, new_color in regex_replacements:
        svg = re.sub(rf"stroke:\s*({pattern})", f"stroke:{new_color}", svg, flags=re.IGNORECASE)
        svg = re.sub(rf"fill:\s*({pattern})", f"fill:{new_color}", svg, flags=re.IGNORECASE)
        svg = re.sub(rf"stroke=['\"]({pattern})['\"]", f"stroke='{new_color}'", svg, flags=re.IGNORECASE)
        svg = re.sub(rf"fill=['\"]({pattern})['\"]", f"fill='{new_color}'", svg, flags=re.IGNORECASE)
    return svg


def render_clean_ligand_svg_manual(mol, width: int = 620, height: int = 360) -> Optional[str]:
    """
    Clean dark-mode 2D ligand SVG using RDKit's native renderer.

    This intentionally no longer manually draws the ring bonds. RDKit's native
    drawer gives cleaner aromatic rings and better double-bond geometry. We still
    remove nonpolar C-H hydrogens, enlarge labels/bonds, and force the atom
    palette for dark-mode visibility.
    """
    if mol is None:
        return None

    display_mol = clean_ligand_for_display(mol)
    if display_mol is None or display_mol.GetNumAtoms() == 0:
        return None

    try:
        # Avoid kekulization failures on CGenFF MOL2 aromatic systems.
        from rdkit.Chem.Draw import rdMolDraw2D
        display_mol = rdMolDraw2D.PrepareMolForDrawing(
            display_mol,
            kekulize=False,
            addChiralHs=False,
            wedgeBonds=True,
        )
    except Exception:
        try:
            AllChem.Compute2DCoords(display_mol)
        except Exception:
            try:
                Chem.rdDepictor.Compute2DCoords(display_mol)
            except Exception:
                return None

    drawer = Draw.MolDraw2DSVG(width, height)
    try:
        opts = drawer.drawOptions()
        opts.clearBackground = False
        opts.bondLineWidth = 2.6
        opts.fixedBondLength = 34
        opts.padding = 0.08
        opts.multipleBondOffset = 0.18
        opts.fixedFontSize = 18
        opts.minFontSize = 14
        opts.maxFontSize = 24
        opts.legendFontSize = 0
        # Keep atom labels legible on the dark card.
        try:
            opts.updateAtomPalette({
                1:  (0.62, 0.72, 0.84),  # H
                6:  (1.00, 0.90, 0.20),  # C / bond skeleton
                7:  (0.21, 0.82, 1.00),  # N
                8:  (1.00, 0.17, 0.17),  # O
                9:  (0.49, 1.00, 0.61),  # F
                15: (1.00, 0.70, 0.28),  # P
                16: (1.00, 0.85, 0.30),  # S
                17: (0.00, 1.00, 0.40),  # Cl
                35: (1.00, 0.54, 0.24),  # Br
                53: (0.83, 0.42, 1.00),  # I
            })
        except Exception:
            pass
    except Exception:
        pass

    try:
        drawer.DrawMolecule(display_mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        if svg.startswith("<?xml"):
            svg = svg.split("?>", 1)[1]
        svg = recolor_rdkit_svg(svg)
        # Make sure it scales correctly inside our card.
        svg = svg.replace('<svg ', '<svg class="pymacs-ligand-svg" ', 1) if 'class="pymacs-ligand-svg"' not in svg[:200] else svg
        return svg
    except Exception:
        return None


def mol_to_3d_payload(mol, structure_path: Optional[Path]) -> Tuple[Optional[str], Optional[str]]:
    """Return base64-encoded structure text and a 3Dmol.js format label."""
    payload_text = None
    payload_format = None
    if mol is not None:
        try:
            m3d = Chem.Mol(mol)
            if m3d.GetNumConformers() == 0:
                try:
                    m3d = Chem.AddHs(m3d, addCoords=True)
                    AllChem.EmbedMolecule(m3d, randomSeed=0xC0FFEE)
                    AllChem.MMFFOptimizeMolecule(m3d, maxIters=200)
                except Exception:
                    pass
            payload_text = Chem.MolToMolBlock(m3d)
            payload_format = 'sdf'
        except Exception:
            payload_text = None
            payload_format = None
    if not payload_text and structure_path is not None and structure_path.exists():
        try:
            payload_text = structure_path.read_text(encoding='utf-8', errors='ignore')
            suffix = structure_path.suffix.lower().lstrip('.')
            payload_format = {'mol2': 'mol2', 'sdf': 'sdf', 'mol': 'sdf', 'pdb': 'pdb'}.get(suffix, suffix or 'sdf')
        except Exception:
            payload_text = None
            payload_format = None
    if not payload_text:
        return None, None
    return base64.b64encode(payload_text.encode('utf-8', errors='ignore')).decode('ascii'), payload_format


def render_ligand_svg(structure_path: Optional[Path]) -> LigandVisual:
    """
    Build ligand visuals using two different sources on purpose:
      - 3D viewer: from the original structure file (preferably .cgenff.mol2)
      - 2D depiction: from the matching .perceived.smi when available

    This avoids the awkward 2D geometry that can happen when trying to derive a
    publication-style 2D depiction directly from a 3D MOL2 coordinate file.
    """
    lv = LigandVisual(structure_file=structure_path)
    if structure_path is None:
        lv.note = 'No ligand structure file found.'
        return lv

    # Always prepare the 3D payload from the original structure path.
    if not RDKIT_AVAILABLE:
        lv.note = f'RDKit not available. Found structure file: {structure_path.name}'
        lv.structure_text_b64, lv.structure_format = mol_to_3d_payload(None, structure_path)
        return lv

    structure_mol = rdkit_mol_from_file(structure_path)
    if structure_mol is None and structure_path.suffix.lower() == '.mol2':
        sdf_fallback = obabel_mol2_to_sdf_temp(structure_path)
        if sdf_fallback is not None:
            structure_mol = rdkit_mol_from_file(sdf_fallback)

    # Prefer a matching .perceived.smi for the 2D picture.
    depiction_path = find_matching_perceived_smiles_file(structure_path)
    lv.depiction_file = depiction_path

    depiction_mol = None
    if depiction_path is not None:
        smi = read_first_smiles(depiction_path)
        if smi:
            try:
                depiction_mol = Chem.MolFromSmiles(smi)
            except Exception:
                depiction_mol = None

    # Fallback to the parsed structure if no perceived smiles was available.
    if depiction_mol is None and structure_mol is not None:
        depiction_mol = Chem.Mol(structure_mol)

    if depiction_mol is None and structure_mol is None:
        lv.note = f'Failed to parse ligand structure: {structure_path.name}'
        lv.structure_text_b64, lv.structure_format = mol_to_3d_payload(None, structure_path)
        return lv

    try:
        lv.svg = render_clean_ligand_svg_manual(depiction_mol)
        if not lv.svg:
            raise ValueError('2D renderer returned empty output')
    except Exception as e:
        try:
            mol2d = clean_ligand_for_display(depiction_mol)
            from rdkit.Chem.Draw import rdMolDraw2D
            mol2d = rdMolDraw2D.PrepareMolForDrawing(
                mol2d,
                kekulize=False,
                addChiralHs=False,
                wedgeBonds=True,
            )
            drawer = Draw.MolDraw2DSVG(620, 360)
            opts = drawer.drawOptions()
            opts.clearBackground = False
            opts.bondLineWidth = 2.6
            opts.fixedBondLength = 34
            opts.padding = 0.08
            opts.multipleBondOffset = 0.18
            opts.fixedFontSize = 18
            opts.minFontSize = 14
            opts.maxFontSize = 24
            try:
                opts.updateAtomPalette({
                    1:  (0.62, 0.72, 0.84),
                    6:  (1.00, 0.90, 0.20),
                    7:  (0.21, 0.82, 1.00),
                    8:  (1.00, 0.17, 0.17),
                    9:  (0.49, 1.00, 0.61),
                    15: (1.00, 0.70, 0.28),
                    16: (1.00, 0.85, 0.30),
                    17: (0.00, 1.00, 0.40),
                    35: (1.00, 0.54, 0.24),
                    53: (0.83, 0.42, 1.00),
                })
            except Exception:
                pass
            drawer.DrawMolecule(mol2d)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            if svg.startswith('<?xml'):
                svg = svg.split('?>', 1)[1]
            lv.svg = recolor_rdkit_svg(svg)
            lv.note = f'Used fallback 2D drawing after renderer issue: {e}'
        except Exception as e2:
            lv.note = f'Failed to draw 2D ligand for {structure_path.name}: {e2}'

    # Prefer chemistry annotations from the original structure if possible.
    props_mol = structure_mol if structure_mol is not None else depiction_mol
    try:
        lv.formula = Chem.rdMolDescriptors.CalcMolFormula(props_mol)
    except Exception:
        pass
    try:
        lv.mw = float(Descriptors.MolWt(props_mol))
    except Exception:
        pass
    try:
        lv.smiles = Chem.MolToSmiles(Chem.RemoveHs(depiction_mol))
    except Exception:
        try:
            lv.smiles = Chem.MolToSmiles(Chem.RemoveHs(props_mol))
        except Exception:
            pass

    lv.structure_text_b64, lv.structure_format = mol_to_3d_payload(structure_mol, structure_path)

    if depiction_path is not None:
        lv.note = (lv.note + ' | ' if lv.note else '') + f'2D from {depiction_path.name}; 3D from {structure_path.name}'
    else:
        lv.note = (lv.note + ' | ' if lv.note else '') + f'2D and 3D both sourced from {structure_path.name}'

    return lv


def load_one_run(run_dir: Path, use_framewise: bool = False, max_framewise_mb: float = 500.0, make_ligand_visual: bool = True, contact_hover_framewise: bool = False, max_hover_framewise_mb: float = 250.0) -> RunData:
    analysis_dir = run_dir / "Analysis_Results"
    name = run_dir.name
    family, state = infer_family_and_state(name)
    is_apo = "APO" in state.upper() or name.upper().endswith("_APO")
    info(f"Loading {name}...")
    ligand_rmsd_path = find_ligand_rmsd_path(analysis_dir)
    ligand_code = extract_ligand_code_from_filename(ligand_rmsd_path) or (None if is_apo else state)
    protein_rmsd_path = first_existing(glob_many(analysis_dir, ["Protein_RMSD.csv", "protein_rmsd.csv", "rmsd_over_time.csv", "RMSD_over_time.csv", "Backbone_RMSD.csv", "GlobalProtein_RMSD.csv", "*Protein*RMSD*.csv", "*protein*rmsd*.csv"]))
    rg_path = first_existing(glob_many(analysis_dir, ["Radius_of_Gyration_Protein_Ligand_Complex.csv", "Radius_of_Gyration_Protein.csv", "Radius_of_Gyration.csv", "Rg.csv", "*Gyration*.csv", "*gyration*.csv", "*Rg*.csv"]))
    bound_complex_path = first_existing(glob_many(analysis_dir, ["*_BoundComplex_RMSD.csv", "*BoundComplex*RMSD*.csv", "*bound*complex*rmsd*.csv"]))
    ligand_rmsf_path = find_ligand_rmsf_path(analysis_dir, ligand_code)
    protein_rmsf_paths = find_protein_rmsf_paths(analysis_dir, ligand_code)
    residue_summary_path = find_residue_summary_path(analysis_dir)
    type_summary_path = find_type_summary_path(analysis_dir)
    framewise_path = find_framewise_interaction_path(analysis_dir)
    contact_framewise_path = find_contact_framewise_path(analysis_dir)
    interaction_types_framewise_path = find_interaction_types_framewise_path(analysis_dir)
    run = RunData(name=name, run_dir=run_dir, analysis_dir=analysis_dir, family=family, state=state, is_apo=is_apo, ligand_code=ligand_code, ligand_display="APO" if is_apo else (ligand_code or state))
    run.files = {
        "protein_rmsd_csv": str(protein_rmsd_path) if protein_rmsd_path else None,
        "protein_rmsf_csvs": ";".join(str(p) for p in protein_rmsf_paths) if protein_rmsf_paths else None,
        "radius_gyration_csv": str(rg_path) if rg_path else None,
        "ligand_rmsd_csv": str(ligand_rmsd_path) if ligand_rmsd_path else None,
        "ligand_rmsf_csv": str(ligand_rmsf_path) if ligand_rmsf_path else None,
        "bound_complex_rmsd_csv": str(bound_complex_path) if bound_complex_path else None,
        "residue_summary_csv": str(residue_summary_path) if residue_summary_path else None,
        "interaction_type_summary_csv": str(type_summary_path) if type_summary_path else None,
        "framewise_interaction_csv_skipped_by_default": str(framewise_path) if framewise_path else None,
        "contact_hover_framewise_csv": str(contact_framewise_path) if contact_framewise_path else None,
        "interaction_type_framewise_csv": str(interaction_types_framewise_path) if interaction_types_framewise_path else None,
    }
    if make_ligand_visual and not run.is_apo:
        ligand_structure_file = find_ligand_structure_file(run)
        run.files["ligand_structure_file"] = str(ligand_structure_file) if ligand_structure_file else None
        run.ligand_visual = render_ligand_svg(ligand_structure_file)
    else:
        run.files["ligand_structure_file"] = None

    run.protein_rmsd = load_time_series_csv(protein_rmsd_path, ["Protein_RMSD", "protein_rmsd", "Backbone_RMSD", "RMSD", "RMSD_A", "RMSD_Å"], "protein_rmsd_A")
    run.radius_gyration = load_time_series_csv(rg_path, ["Rg", "Radius_of_Gyration", "Rg_Protein", "Rg_Protein_Ligand_Complex", "Rg_A", "RadiusOfGyration"], "rg_A")
    run.ligand_rmsd = load_time_series_csv(ligand_rmsd_path, ["Ligand_RMSD", "ligand_rmsd", "RMSD", "RMSD_A", "RMSD_Å"], "ligand_rmsd_A")
    run.bound_complex_rmsd = load_time_series_csv(bound_complex_path, ["BoundComplex_RMSD", "Complex_RMSD", "RMSD", "RMSD_A", "RMSD_Å"], "bound_complex_rmsd_A")
    run.protein_rmsf = load_combined_protein_rmsf(protein_rmsf_paths)
    run.ligand_rmsf = load_rmsf_csv(ligand_rmsf_path, source_label=ligand_rmsf_path.stem) if ligand_rmsf_path else None
    run.residue_frequency = load_residue_summary(residue_summary_path)
    run.type_fraction = load_type_summary(type_summary_path)
    run.residue_hover_summary = None
    run.residue_interaction_summary = None
    if contact_hover_framewise:
        run.residue_hover_summary = load_residue_hover_summary_from_framewise(contact_framewise_path, max_mb=max_hover_framewise_mb)
    if run.residue_hover_summary is None and run.residue_frequency is not None:
        run.residue_hover_summary = standardize_hover_summary_from_residue_frequency(run, source_label=residue_summary_path.name if residue_summary_path else "residue summary")
    run.residue_interaction_summary = load_residue_interaction_summary(type_summary_path, framewise_path=interaction_types_framewise_path, max_mb=max_hover_framewise_mb)
    interaction_meta = {"interaction_rows": None, "unique_contact_residues": int(run.residue_frequency["ResidueKey"].nunique()) if run.residue_frequency is not None else 0, "interaction_frame_count": None, "interaction_load_first_half": None, "interaction_load_last_half": None, "interaction_load_ratio_last_over_first": None, "used_framewise": False, "framewise_file_size_mb": file_size_mb(framewise_path)}
    if use_framewise and (run.residue_frequency is None or run.type_fraction is None):
        size = file_size_mb(framewise_path)
        if framewise_path and size is not None and size > max_framewise_mb:
            warn(f"Skipping framewise parse for {name}: {size:.1f} MB > --max-framewise-mb {max_framewise_mb}")
        elif framewise_path:
            residue_df, type_df, meta = summarize_framewise_chunked(framewise_path)
            if run.residue_frequency is None:
                run.residue_frequency = residue_df
            if run.type_fraction is None:
                run.type_fraction = type_df
            interaction_meta.update(meta)
            interaction_meta["used_framewise"] = True
    run.metrics = build_basic_metrics(run, interaction_meta)
    info(f"{name}: apo={run.is_apo}, protein_rmsd={'Y' if run.protein_rmsd is not None else 'N'}, protein_rmsf={'Y' if run.protein_rmsf is not None else 'N'}, rg={'Y' if run.radius_gyration is not None else 'N'}, ligand_rmsd={'Y' if run.ligand_rmsd is not None else 'N'}, contacts={'Y' if run.residue_frequency is not None else 'N'}")
    return run


def build_basic_metrics(run: RunData, interaction_meta: Dict[str, Any]) -> Dict[str, Any]:
    pr, rg, lr, bc = run.protein_rmsd, run.radius_gyration, run.ligand_rmsd, run.bound_complex_rmsd
    metrics = {
        "family": run.family, "state": run.state, "is_apo": run.is_apo, "ligand_display": run.ligand_display,
        "protein_rmsd_mean_A": numeric_mean(pr["protein_rmsd_A"]) if pr is not None else None,
        "protein_rmsd_std_A": numeric_std(pr["protein_rmsd_A"]) if pr is not None else None,
        "protein_rmsd_final_A": numeric_last(pr["protein_rmsd_A"]) if pr is not None else None,
        "protein_rmsd_last_half_mean_A": numeric_last_half_mean(pr["protein_rmsd_A"]) if pr is not None else None,
        "protein_rmsd_slope_A_per_ns": linear_slope(pr["time_ns"], pr["protein_rmsd_A"]) if pr is not None else None,
        "protein_rmsf_mean_A": numeric_mean(run.protein_rmsf["RMSF_A"]) if run.protein_rmsf is not None else None,
        "protein_rmsf_std_A": numeric_std(run.protein_rmsf["RMSF_A"]) if run.protein_rmsf is not None else None,
        "rg_mean_A": numeric_mean(rg["rg_A"]) if rg is not None else None,
        "rg_std_A": numeric_std(rg["rg_A"]) if rg is not None else None,
        "rg_final_A": numeric_last(rg["rg_A"]) if rg is not None else None,
        "rg_last_half_mean_A": numeric_last_half_mean(rg["rg_A"]) if rg is not None else None,
        "rg_slope_A_per_ns": linear_slope(rg["time_ns"], rg["rg_A"]) if rg is not None else None,
        "ligand_rmsd_mean_A": numeric_mean(lr["ligand_rmsd_A"]) if lr is not None else None,
        "ligand_rmsd_std_A": numeric_std(lr["ligand_rmsd_A"]) if lr is not None else None,
        "ligand_rmsd_final_A": numeric_last(lr["ligand_rmsd_A"]) if lr is not None else None,
        "ligand_rmsd_last_half_mean_A": numeric_last_half_mean(lr["ligand_rmsd_A"]) if lr is not None else None,
        "ligand_rmsd_slope_A_per_ns": linear_slope(lr["time_ns"], lr["ligand_rmsd_A"]) if lr is not None else None,
        "ligand_rmsf_mean_A": numeric_mean(run.ligand_rmsf["RMSF_A"]) if run.ligand_rmsf is not None else None,
        "bound_complex_rmsd_mean_A": numeric_mean(bc["bound_complex_rmsd_A"]) if bc is not None else None,
        "has_protein_rmsd": pr is not None, "has_protein_rmsf": run.protein_rmsf is not None, "has_rg": rg is not None,
        "has_ligand_rmsd": lr is not None, "has_ligand_rmsf": run.ligand_rmsf is not None,
        "has_residue_contacts": run.residue_frequency is not None, "has_interaction_types": run.type_fraction is not None,
        "has_ligand_visual": run.ligand_visual.svg is not None,
    }
    metrics.update(interaction_meta)
    return metrics


def choose_apo(runs: List[RunData], apo_name: Optional[str]) -> Optional[RunData]:
    if apo_name:
        for r in runs:
            if r.name == apo_name:
                return r
        warn(f"Explicit apo run not found: {apo_name}")
    apos = [r for r in runs if r.is_apo]
    if not apos:
        return None
    exact = [r for r in apos if r.state.upper() == "APO"]
    return exact[0] if exact else apos[0]


def build_summary_table(runs: List[RunData]) -> pd.DataFrame:
    df = pd.DataFrame([{"Run": r.name, **r.metrics} for r in runs])
    preferred = ["Run", "family", "state", "is_apo", "ligand_display", "protein_rmsd_mean_A", "protein_rmsd_last_half_mean_A", "protein_rmsd_final_A", "protein_rmsd_slope_A_per_ns", "protein_rmsf_mean_A", "protein_rmsf_std_A", "rg_mean_A", "rg_last_half_mean_A", "rg_final_A", "rg_slope_A_per_ns", "ligand_rmsd_mean_A", "ligand_rmsd_last_half_mean_A", "ligand_rmsd_final_A", "ligand_rmsd_slope_A_per_ns", "ligand_rmsf_mean_A", "unique_contact_residues", "interaction_rows", "used_framewise", "framewise_file_size_mb", "has_protein_rmsd", "has_protein_rmsf", "has_rg", "has_ligand_rmsd", "has_residue_contacts", "has_interaction_types", "has_ligand_visual"]
    cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
    return df[cols] if not df.empty else df


def build_availability_table(runs: List[RunData]) -> pd.DataFrame:
    return pd.DataFrame([{
        "Run": r.name, "State": r.state, "Ligand": r.ligand_display,
        "Protein RMSD": "Yes" if r.protein_rmsd is not None else "No",
        "Protein RMSF": "Yes" if r.protein_rmsf is not None else "No",
        "Rg": "Yes" if r.radius_gyration is not None else "No",
        "Ligand RMSD": "Yes" if r.ligand_rmsd is not None else "No",
        "Ligand RMSF": "Yes" if r.ligand_rmsf is not None else "No",
        "Residue contact summary": "Yes" if r.residue_frequency is not None else "No",
        "Interaction type summary": "Yes" if r.type_fraction is not None else "No",
        "2D ligand": "Yes" if r.ligand_visual.svg is not None else "No",
        "Ligand file": r.ligand_visual.structure_file.name if r.ligand_visual.structure_file else "",
        "Framewise used": "Yes" if r.metrics.get("used_framewise") else "No",
        "Framewise MB": nice(r.metrics.get("framewise_file_size_mb"), 1),
    } for r in runs])


def cross_label(run: RunData) -> str:
    ligand = run.ligand_display or run.ligand_code or run.state or run.name
    return f"{run.family} | {ligand}"


def dominant_interaction_type_for_run(run: RunData) -> Optional[str]:
    if run.type_fraction is None or run.type_fraction.empty or "InteractionType" not in run.type_fraction.columns:
        return None
    work = run.type_fraction.copy()
    work = work[~work["InteractionType"].apply(is_excluded_interaction_type)]
    if work.empty:
        return None
    frac = pd.to_numeric(work["fraction"], errors="coerce") if "fraction" in work.columns else pd.Series(np.nan, index=work.index)
    if frac.dropna().empty and "count" in work.columns:
        frac = pd.to_numeric(work["count"], errors="coerce")
    if frac.dropna().empty:
        return str(work.iloc[0]["InteractionType"])
    return str(work.loc[frac.fillna(-np.inf).idxmax(), "InteractionType"])


def dataframe_to_long_time(df: Optional[pd.DataFrame], value_candidates: Sequence[str], run: RunData, metric_name: str) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()

    time = detect_time_ns(df)
    value_col = find_first_col(df, value_candidates)
    if value_col is None:
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        numeric_cols = [c for c in numeric_cols if str(c).lower() not in {"time", "time_ns", "time_ps", "frame"}]
        if not numeric_cols:
            return pd.DataFrame()
        value_col = numeric_cols[0]

    fallback_time = np.arange(len(df), dtype=float)
    out = pd.DataFrame({
        "species": run.family,
        "run_name": run.name,
        "state": run.state,
        "ligand_code": run.ligand_code,
        "ligand_display": run.ligand_display or run.ligand_code or run.state,
        "is_apo": run.is_apo,
        "label": cross_label(run),
        "metric": metric_name,
        "time_ns": pd.to_numeric(time, errors="coerce") if time is not None else fallback_time,
        "value": pd.to_numeric(df[value_col], errors="coerce"),
    })
    return out.dropna(subset=["value"]).reset_index(drop=True)


def rmsf_to_long(df: Optional[pd.DataFrame], run: RunData, metric_name: str = "protein_rmsf") -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()

    res_col = find_first_col(df, ["Residue", "residue", "ResNum", "resnum", "Residue_Index", "index"])
    val_col = find_first_col(df, ["RMSF", "rmsf", "RMSF_A", "rmsf_A", "Value", "value"])
    if val_col is None:
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        if not numeric_cols:
            return pd.DataFrame()
        val_col = numeric_cols[-1]

    residue_raw = df[res_col] if res_col else pd.Series(np.arange(1, len(df) + 1), index=df.index)
    residue_num = pd.to_numeric(df["ResNum"], errors="coerce") if "ResNum" in df.columns else residue_raw.apply(lambda x: parse_residue_label(x).get("resnum"))
    residue_key = df["ResidueKey"].astype(str) if "ResidueKey" in df.columns else residue_raw.apply(lambda x: parse_residue_label(x).get("key"))
    chain_vals = df["Chain"].astype(str) if "Chain" in df.columns else residue_raw.apply(lambda x: parse_residue_label(x).get("chain") or "")
    residue_index = np.arange(1, len(df) + 1)
    normalized_position = np.linspace(0.0, 1.0, num=len(df)) if len(df) > 1 else np.array([0.0])

    out = pd.DataFrame({
        "species": run.family,
        "run_name": run.name,
        "state": run.state,
        "ligand_code": run.ligand_code,
        "ligand_display": run.ligand_display or run.ligand_code or run.state,
        "is_apo": run.is_apo,
        "label": cross_label(run),
        "metric": metric_name,
        "residue_raw": residue_raw.astype(str),
        "residue_num": pd.to_numeric(residue_num, errors="coerce"),
        "residue_key": residue_key.astype(str),
        "chain": chain_vals.astype(str),
        "residue_index": residue_index,
        "normalized_position": normalized_position,
        "rmsf_A": pd.to_numeric(df[val_col], errors="coerce"),
    })
    return out.dropna(subset=["rmsf_A"]).reset_index(drop=True)


def contact_residue_x_positions_for_run(run: RunData) -> List[float]:
    if run.is_apo or run.residue_frequency is None or run.residue_frequency.empty:
        return []
    tmp = run.residue_frequency.copy()
    xs = pd.to_numeric(tmp.get("ResNum"), errors="coerce") if "ResNum" in tmp.columns else pd.Series(np.nan, index=tmp.index)
    if xs.isna().all() and "Residue" in tmp.columns:
        xs = pd.to_numeric(tmp["Residue"].map(lambda v: parse_residue_label(v).get("resnum")), errors="coerce")
    return sorted({float(v) for v in xs.dropna().tolist()})


def build_contact_hover_records_for_run(run: RunData, rmsf_sub: pd.DataFrame, x_mode: str = "residue") -> pd.DataFrame:
    if run.is_apo or x_mode != "residue" or run.residue_hover_summary is None or run.residue_hover_summary.empty:
        return pd.DataFrame()
    hover = run.residue_hover_summary.copy()
    if "ResidueKey" not in hover.columns:
        return pd.DataFrame()
    if run.residue_interaction_summary is not None and not run.residue_interaction_summary.empty and "ResidueKey" in run.residue_interaction_summary.columns:
        inter = run.residue_interaction_summary.copy()
        join_cols = [c for c in inter.columns if c != "ResidueKey" and c not in hover.columns]
        hover = hover.merge(inter[["ResidueKey"] + join_cols], on="ResidueKey", how="left")
    else:
        for col in ["dominant_interaction_type", "interaction_type_list", "hbond_fraction", "hydrophobic_fraction", "ionic_fraction", "polar_fraction", "water_bridge_fraction", "vdw_fraction"]:
            if col not in hover.columns:
                hover[col] = np.nan

    rmsf_work = rmsf_sub.copy()
    if "residue_key" not in rmsf_work.columns:
        rmsf_work["residue_key"] = rmsf_work["residue_raw"].astype(str).map(lambda v: parse_residue_label(v).get("key"))
    if "chain" not in rmsf_work.columns:
        rmsf_work["chain"] = rmsf_work["residue_raw"].astype(str).map(lambda v: parse_residue_label(v).get("chain") or "")
    rmsf_work["residue_num"] = pd.to_numeric(rmsf_work.get("residue_num"), errors="coerce")
    rmsf_work["residue_index"] = pd.to_numeric(rmsf_work.get("residue_index"), errors="coerce")
    rmsf_work["rmsf_A"] = pd.to_numeric(rmsf_work.get("rmsf_A"), errors="coerce")
    rmsf_work = rmsf_work.dropna(subset=["rmsf_A"]).copy()

    ymax = float(rmsf_work["rmsf_A"].max()) if not rmsf_work.empty else 1.0
    rows = []
    for _, row in hover.iterrows():
        residue_key = str(row.get("ResidueKey") or "")
        residue_label = str(row.get("residue_label") or row.get("Residue") or residue_key)
        chain = "" if pd.isna(row.get("Chain")) else str(row.get("Chain") or "")
        resnum = pd.to_numeric(pd.Series([row.get("ResNum")]), errors="coerce").iloc[0]
        exact = rmsf_work[rmsf_work["residue_key"].astype(str) == residue_key].copy()
        if exact.empty and pd.notna(resnum):
            exact = rmsf_work[pd.to_numeric(rmsf_work["residue_num"], errors="coerce") == float(resnum)].copy()
        if exact.empty:
            exact = rmsf_work[rmsf_work["residue_raw"].astype(str) == residue_label].copy()
        x_val = np.nan
        y_val = np.nan
        if not exact.empty:
            x_num = pd.to_numeric(exact["residue_num"], errors="coerce")
            if x_num.notna().any():
                x_val = float(x_num.dropna().iloc[0])
            else:
                x_idx = pd.to_numeric(exact["residue_index"], errors="coerce")
                if x_idx.notna().any():
                    x_val = float(x_idx.dropna().iloc[0])
            rmsf_vals = pd.to_numeric(exact["rmsf_A"], errors="coerce").dropna()
            y_val = float(rmsf_vals.iloc[0]) if not rmsf_vals.empty else np.nan
        if pd.isna(x_val):
            parsed = parse_residue_label(residue_label).get("resnum")
            if parsed is not None:
                x_val = float(parsed)
        if pd.isna(x_val):
            continue
        rows.append({
            "x": float(x_val),
            "y": y_val,
            "y_top": ymax,
            "label": cross_label(run),
            "species": run.family,
            "run_name": run.name,
            "ligand_display": run.ligand_display or run.ligand_code or run.state,
            "residue_label": residue_label,
            "residue_key": residue_key,
            "chain": chain,
            "resnum": float(resnum) if pd.notna(resnum) else np.nan,
            "contact_fraction": pd.to_numeric(pd.Series([row.get("contact_fraction")]), errors="coerce").iloc[0],
            "contact_percent": pd.to_numeric(pd.Series([row.get("contact_percent")]), errors="coerce").iloc[0],
            "contact_count": pd.to_numeric(pd.Series([row.get("contact_count")]), errors="coerce").iloc[0],
            "contact_frames": pd.to_numeric(pd.Series([row.get("contact_frames")]), errors="coerce").iloc[0],
            "total_frames": pd.to_numeric(pd.Series([row.get("total_frames")]), errors="coerce").iloc[0],
            "dominant_interaction_type": row.get("dominant_interaction_type"),
            "interaction_type_list": row.get("interaction_type_list"),
            "hbond_fraction": pd.to_numeric(pd.Series([row.get("hbond_fraction")]), errors="coerce").iloc[0],
            "hydrophobic_fraction": pd.to_numeric(pd.Series([row.get("hydrophobic_fraction")]), errors="coerce").iloc[0],
            "ionic_fraction": pd.to_numeric(pd.Series([row.get("ionic_fraction")]), errors="coerce").iloc[0],
            "polar_fraction": pd.to_numeric(pd.Series([row.get("polar_fraction")]), errors="coerce").iloc[0],
            "water_bridge_fraction": pd.to_numeric(pd.Series([row.get("water_bridge_fraction")]), errors="coerce").iloc[0],
            "vdw_fraction": pd.to_numeric(pd.Series([row.get("vdw_fraction")]), errors="coerce").iloc[0],
            "mean_distance_A": pd.to_numeric(pd.Series([row.get("mean_distance_A")]), errors="coerce").iloc[0],
            "min_distance_A": pd.to_numeric(pd.Series([row.get("min_distance_A")]), errors="coerce").iloc[0],
            "top_ligand_atom": row.get("top_ligand_atom"),
            "top_protein_atom": row.get("top_protein_atom"),
            "nearest_rmsf_A": y_val,
            "source": row.get("source"),
        })
    out = pd.DataFrame(rows)
    return out.sort_values(["contact_fraction", "residue_label"], ascending=[False, True]).reset_index(drop=True) if not out.empty else out


def _hover_clean(value: Any, kind: str = "text") -> str:
    if value is None:
        return "n/a"
    try:
        if pd.isna(value):
            return "n/a"
    except Exception:
        pass
    if kind == "percent":
        try:
            return f"{float(value):.1f}%"
        except Exception:
            return "n/a"
    if kind == "fraction":
        try:
            return f"{float(value):.3f}"
        except Exception:
            return "n/a"
    if kind == "count":
        try:
            return str(int(round(float(value))))
        except Exception:
            return "n/a"
    if kind == "float":
        try:
            return f"{float(value):.3f}"
        except Exception:
            return "n/a"
    return str(value) if str(value).strip() else "n/a"


def make_contact_residue_hover_trace(records: pd.DataFrame, color: str, name: str, ymax: float) -> go.Scatter:
    xvals: List[Optional[float]] = []
    yvals: List[Optional[float]] = []
    custom_rows: List[List[Any]] = []
    if records is None or records.empty:
        return go.Scatter()
    for _, row in records.iterrows():
        meta = [
            _hover_clean(row.get("label")),
            _hover_clean(row.get("species")),
            _hover_clean(row.get("run_name")),
            _hover_clean(row.get("ligand_display")),
            _hover_clean(row.get("residue_label")),
            _hover_clean(row.get("chain")),
            _hover_clean(row.get("resnum"), kind="count"),
            _hover_clean(row.get("contact_percent"), kind="percent"),
            _hover_clean(row.get("contact_fraction"), kind="fraction"),
            _hover_clean(row.get("contact_frames"), kind="count"),
            _hover_clean(row.get("total_frames"), kind="count"),
            _hover_clean(row.get("contact_count"), kind="count"),
            _hover_clean(row.get("dominant_interaction_type")),
            _hover_clean(row.get("interaction_type_list")),
            _hover_clean(row.get("hbond_fraction"), kind="fraction"),
            _hover_clean(row.get("hydrophobic_fraction"), kind="fraction"),
            _hover_clean(row.get("ionic_fraction"), kind="fraction"),
            _hover_clean(row.get("mean_distance_A"), kind="float"),
            _hover_clean(row.get("min_distance_A"), kind="float"),
            _hover_clean(row.get("top_ligand_atom")),
            _hover_clean(row.get("top_protein_atom")),
            _hover_clean(row.get("nearest_rmsf_A"), kind="float"),
            _hover_clean(row.get("source")),
        ]
        x = float(row["x"])
        xvals.extend([x, x, None])
        yvals.extend([0.0, ymax, None])
        custom_rows.extend([meta, meta, ["n/a"] * len(meta)])
    return go.Scatter(
        x=xvals,
        y=yvals,
        customdata=np.array(custom_rows, dtype=object),
        mode="lines+markers",
        name=name,
        line=dict(color=color, width=1, dash="dot"),
        marker=dict(size=7, symbol="diamond", color=color),
        visible="legendonly",
        opacity=0.75,
        hovertemplate="%{customdata[0]}<br>Species: %{customdata[1]}<br>Run: %{customdata[2]}<br>Ligand: %{customdata[3]}<br>Residue: %{customdata[4]}<br>Chain: %{customdata[5]}<br>Residue number: %{customdata[6]}<br>Contact fraction: %{customdata[7]}<br>Raw contact fraction: %{customdata[8]}<br>Contact frames: %{customdata[9]} / %{customdata[10]}<br>Contact count: %{customdata[11]}<br>Dominant type: %{customdata[12]}<br>Interaction types: %{customdata[13]}<br>H-bond fraction: %{customdata[14]}<br>Hydrophobic fraction: %{customdata[15]}<br>Ionic fraction: %{customdata[16]}<br>Mean distance: %{customdata[17]} A<br>Min distance: %{customdata[18]} A<br>Top ligand atom: %{customdata[19]}<br>Top protein atom: %{customdata[20]}<br>RMSF at/near residue: %{customdata[21]} A<br>Source: %{customdata[22]}<extra></extra>",
    )


def build_contact_hover_export(runs: List[RunData]) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for run in runs:
        if run.is_apo or run.protein_rmsf is None or run.protein_rmsf.empty:
            continue
        rmsf_sub = rmsf_to_long(run.protein_rmsf, run, metric_name="protein_rmsf")
        rec = build_contact_hover_records_for_run(run, rmsf_sub, x_mode="residue")
        if rec is not None and not rec.empty:
            frames.append(rec)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    cols = [
        "species", "run_name", "ligand_display", "residue_label", "residue_key", "chain", "resnum",
        "contact_fraction", "contact_percent", "contact_count", "contact_frames", "total_frames",
        "dominant_interaction_type", "interaction_type_list", "hbond_fraction", "hydrophobic_fraction", "ionic_fraction",
        "mean_distance_A", "min_distance_A", "top_ligand_atom", "top_protein_atom", "source",
    ]
    for col in cols:
        if col not in out.columns:
            out[col] = np.nan
    return out[cols].copy()


def build_cross_species_availability(summary: pd.DataFrame) -> pd.DataFrame:
    if summary is None or summary.empty:
        return pd.DataFrame()
    cols = [
        "species", "run_name", "state", "ligand_code", "ligand_display", "is_apo",
        "available_protein_rmsd", "available_protein_rmsf", "available_rg",
        "available_ligand_rmsd", "available_ligand_rmsf", "available_contacts",
        "used_framewise", "framewise_file_size_mb",
    ]
    return summary[[c for c in cols if c in summary.columns]].copy()


def build_cross_species_index(summary: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "species", "run_name", "state", "ligand_code", "ligand_display", "is_apo",
        "protein_rmsd_mean_A", "protein_rmsd_last_A", "protein_rmsd_last_half_mean_A",
        "protein_rmsf_mean_A", "rg_mean_A", "rg_last_half_mean_A",
        "ligand_rmsd_mean_A", "ligand_rmsd_last_half_mean_A",
        "ligand_rmsf_mean_A",
        "contact_residue_count", "dominant_interaction_type",
        "available_protein_rmsd", "available_protein_rmsf",
        "available_rg", "available_ligand_rmsd",
        "available_ligand_rmsf", "available_contacts",
    ]
    out = summary.copy() if summary is not None else pd.DataFrame()
    for col in cols:
        if col not in out.columns:
            out[col] = None
    return out[cols].copy()


def build_cross_species_ligand_table(runs: List[RunData]) -> pd.DataFrame:
    base = build_pairwise_ligand_table(runs)
    if base is None or base.empty:
        return pd.DataFrame()
    lookup = {r.name: r for r in runs}
    rows = []
    for _, row in base.iterrows():
        a = lookup.get(str(row["Ligand_A_Run"]))
        b = lookup.get(str(row["Ligand_B_Run"]))
        rows.append({
            "Species_A": a.family if a else None,
            "State_A": a.state if a else None,
            "Ligand_A_Label": cross_label(a) if a else row.get("Ligand_A"),
            "Ligand_A_Run": row.get("Ligand_A_Run"),
            "Ligand_A": row.get("Ligand_A"),
            "Species_B": b.family if b else None,
            "State_B": b.state if b else None,
            "Ligand_B_Label": cross_label(b) if b else row.get("Ligand_B"),
            "Ligand_B_Run": row.get("Ligand_B_Run"),
            "Ligand_B": row.get("Ligand_B"),
            "CrossSpeciesPair": bool(a and b and a.family != b.family),
            "ContactResidue_Jaccard": row.get("ContactResidue_Jaccard"),
            "SharedContactResidues_N": row.get("SharedContactResidues_N"),
            "UniqueContactResidues_A_N": row.get("UniqueContactResidues_A_N"),
            "UniqueContactResidues_B_N": row.get("UniqueContactResidues_B_N"),
            "InteractionType_CosineSimilarity": row.get("InteractionType_CosineSimilarity"),
            "Delta_ProteinRMSDMean_A_minus_B": row.get("Delta_ProteinRMSDMean_A_minus_B"),
            "Delta_RgMean_A_minus_B": row.get("Delta_RgMean_A_minus_B"),
            "Delta_LigandRMSDMean_A_minus_B": row.get("Delta_LigandRMSDMean_A_minus_B"),
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["CrossSpeciesPair", "ContactResidue_Jaccard", "InteractionType_CosineSimilarity"], ascending=[False, False, False]).reset_index(drop=True)


def build_cross_metric_summary(summary: pd.DataFrame) -> pd.DataFrame:
    if summary is None or summary.empty:
        return pd.DataFrame()
    metric_cols = [
        "protein_rmsd_last_half_mean_A",
        "protein_rmsf_mean_A",
        "rg_last_half_mean_A",
        "ligand_rmsd_last_half_mean_A",
        "contact_residue_count",
    ]
    id_cols = [c for c in ["species", "run_name", "state", "ligand_code", "ligand_display", "is_apo", "label"] if c in summary.columns]
    present_metrics = [c for c in metric_cols if c in summary.columns]
    if not present_metrics:
        return pd.DataFrame()
    out = summary[id_cols + present_metrics].melt(id_vars=id_cols, value_vars=present_metrics, var_name="metric", value_name="value")
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    return out.dropna(subset=["value"]).reset_index(drop=True)


def build_cross_species_stability(runs: List[RunData]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    apo_delta_frames: List[pd.DataFrame] = []
    binding_summary_frames: List[pd.DataFrame] = []
    for family in sorted({r.family for r in runs}):
        fam_runs = [r for r in runs if r.family == family]
        apo = choose_apo(fam_runs, None)
        if apo is None:
            continue
        apo_delta = build_apo_delta_summary(fam_runs, apo)
        if apo_delta is not None and not apo_delta.empty:
            apo_delta_frames.append(apo_delta.assign(species=family))
        binding_delta = build_binding_site_rmsf_delta(fam_runs, apo)
        binding_summary = build_binding_site_summary(binding_delta)
        if binding_summary is not None and not binding_summary.empty:
            binding_summary_frames.append(binding_summary.assign(species=family))

    apo_delta_all = pd.concat(apo_delta_frames, ignore_index=True) if apo_delta_frames else pd.DataFrame()
    binding_summary_all = pd.concat(binding_summary_frames, ignore_index=True) if binding_summary_frames else pd.DataFrame()
    stability = build_stability_flags(runs, apo_delta_all, binding_summary_all)
    if stability is None or stability.empty:
        return pd.DataFrame(), apo_delta_all, binding_summary_all

    lookup = {r.name: r for r in runs}
    stability.insert(0, "species", stability["Run"].map(lambda name: lookup.get(str(name)).family if lookup.get(str(name)) else None))
    stability.insert(1, "state", stability["Run"].map(lambda name: lookup.get(str(name)).state if lookup.get(str(name)) else None))
    stability.insert(2, "ligand_code", stability["Run"].map(lambda name: lookup.get(str(name)).ligand_code if lookup.get(str(name)) else None))
    stability.insert(3, "label", stability["Run"].map(lambda name: cross_label(lookup.get(str(name))) if lookup.get(str(name)) else str(name)))
    return stability, apo_delta_all, binding_summary_all


def build_cross_species_components(runs: List[RunData], top_residues: int = 40, min_contact_fraction: float = 0.0) -> Dict[str, Any]:
    ordered_runs = sorted(runs, key=lambda r: (r.family, not r.is_apo, r.state, r.name))
    summary_rows: List[Dict[str, Any]] = []
    rmsd_long: List[pd.DataFrame] = []
    rg_long: List[pd.DataFrame] = []
    ligand_rmsd_long: List[pd.DataFrame] = []
    rmsf_long: List[pd.DataFrame] = []
    ligand_rmsf_long: List[pd.DataFrame] = []

    for run in ordered_runs:
        rmsd_long.append(dataframe_to_long_time(run.protein_rmsd, ["RMSD_A", "Protein_RMSD_A", "Backbone_RMSD_A", "RMSD", "rmsd"], run, "protein_rmsd"))
        rg_long.append(dataframe_to_long_time(run.radius_gyration, ["Rg_A", "Radius_Gyration_A", "Radius of gyration", "Rg", "rg"], run, "radius_gyration"))
        ligand_rmsd_long.append(dataframe_to_long_time(run.ligand_rmsd, ["Ligand_RMSD_A", "RMSD_A", "RMSD", "rmsd"], run, "ligand_rmsd"))
        rmsf_long.append(rmsf_to_long(run.protein_rmsf, run, metric_name="protein_rmsf"))
        if run.ligand_rmsf is not None and not run.ligand_rmsf.empty:
            ligand_rmsf_long.append(rmsf_to_long(run.ligand_rmsf, run, metric_name="ligand_rmsf"))

        summary_rows.append({
            "species": run.family,
            "family": run.family,
            "run_name": run.name,
            "state": run.state,
            "ligand_code": run.ligand_code,
            "ligand_display": run.ligand_display or run.ligand_code or run.state,
            "is_apo": run.is_apo,
            "label": cross_label(run),
            "protein_rmsd_mean_A": run.metrics.get("protein_rmsd_mean_A"),
            "protein_rmsd_last_A": run.metrics.get("protein_rmsd_final_A"),
            "protein_rmsd_last_half_mean_A": run.metrics.get("protein_rmsd_last_half_mean_A"),
            "protein_rmsf_mean_A": run.metrics.get("protein_rmsf_mean_A"),
            "rg_mean_A": run.metrics.get("rg_mean_A"),
            "rg_last_half_mean_A": run.metrics.get("rg_last_half_mean_A"),
            "ligand_rmsd_mean_A": run.metrics.get("ligand_rmsd_mean_A"),
            "ligand_rmsd_last_half_mean_A": run.metrics.get("ligand_rmsd_last_half_mean_A"),
            "ligand_rmsf_mean_A": run.metrics.get("ligand_rmsf_mean_A"),
            "contact_residue_count": run.metrics.get("unique_contact_residues"),
            "dominant_interaction_type": dominant_interaction_type_for_run(run),
            "available_protein_rmsd": bool(run.protein_rmsd is not None and not run.protein_rmsd.empty),
            "available_protein_rmsf": bool(run.protein_rmsf is not None and not run.protein_rmsf.empty),
            "available_rg": bool(run.radius_gyration is not None and not run.radius_gyration.empty),
            "available_ligand_rmsd": bool(run.ligand_rmsd is not None and not run.ligand_rmsd.empty),
            "available_ligand_rmsf": bool(run.ligand_rmsf is not None and not run.ligand_rmsf.empty),
            "available_contacts": bool(run.residue_frequency is not None and not run.residue_frequency.empty),
            "used_framewise": bool(run.metrics.get("used_framewise")),
            "framewise_file_size_mb": run.metrics.get("framewise_file_size_mb"),
        })

    summary = pd.DataFrame(summary_rows)
    availability = build_cross_species_availability(summary)
    cross_index = build_cross_species_index(summary)
    cross_metric_summary = build_cross_metric_summary(summary)
    contact_mat = build_contact_matrix(ordered_runs, topn=top_residues, min_contact_fraction=min_contact_fraction)
    stability, apo_delta_all, binding_summary_all = build_cross_species_stability(ordered_runs)
    pairwise = build_cross_species_ligand_table(ordered_runs)
    contact_hover_long = build_contact_hover_export(ordered_runs)
    ligand_runs = [r for r in ordered_runs if not r.is_apo]
    cross_type_mat = build_interaction_type_matrix(ligand_runs)
    if cross_type_mat is not None and not cross_type_mat.empty:
        ligand_lookup = {r.name: r for r in ligand_runs}
        cross_type_mat.index = [cross_label(ligand_lookup[idx]) if idx in ligand_lookup else idx for idx in cross_type_mat.index]
    comps: Dict[str, Any] = {
        "summary": summary,
        "availability": availability,
        "cross_index": cross_index,
        "cross_metric_summary": cross_metric_summary,
        "rmsd_long": pd.concat(rmsd_long, ignore_index=True) if any(not df.empty for df in rmsd_long) else pd.DataFrame(),
        "rmsf_long": pd.concat(rmsf_long, ignore_index=True) if any(not df.empty for df in rmsf_long) else pd.DataFrame(),
        "rg_long": pd.concat(rg_long, ignore_index=True) if any(not df.empty for df in rg_long) else pd.DataFrame(),
        "ligand_rmsd_long": pd.concat(ligand_rmsd_long, ignore_index=True) if any(not df.empty for df in ligand_rmsd_long) else pd.DataFrame(),
        "ligand_rmsf_long": pd.concat(ligand_rmsf_long, ignore_index=True) if any(not df.empty for df in ligand_rmsf_long) else pd.DataFrame(),
        "contact_mat": contact_mat,
        "stability": stability,
        "pairwise": pairwise,
        "apo_delta_all": apo_delta_all,
        "binding_summary_all": binding_summary_all,
        "cross_interaction_type_matrix": cross_type_mat,
        "contact_hover_long": contact_hover_long,
    }
    comps["figs"] = make_cross_species_overlay_figures(ordered_runs, comps)
    return comps


def classify_protein_state_change(row: Dict[str, Any]) -> str:
    rg_delta = row.get("delta_rg_last_half_mean_A") or row.get("delta_rg_mean_A")
    rmsd_delta = row.get("delta_protein_rmsd_last_half_mean_A") or row.get("delta_protein_rmsd_mean_A")
    rmsf_delta = row.get("delta_protein_rmsf_mean_A")
    tags = []
    try:
        if rg_delta is not None and not pd.isna(rg_delta):
            tags.append("more compact" if float(rg_delta) <= -0.20 else "more expanded" if float(rg_delta) >= 0.20 else "similar compactness")
    except Exception:
        pass
    try:
        if rmsd_delta is not None and not pd.isna(rmsd_delta):
            tags.append("lower global RMSD" if float(rmsd_delta) <= -0.20 else "higher global RMSD" if float(rmsd_delta) >= 0.20 else "similar global RMSD")
    except Exception:
        pass
    try:
        if rmsf_delta is not None and not pd.isna(rmsf_delta):
            if float(rmsf_delta) <= -0.10:
                tags.append("globally less flexible")
            elif float(rmsf_delta) >= 0.10:
                tags.append("globally more flexible")
    except Exception:
        pass
    return "; ".join(tags) if tags else "insufficient apo-relative data"


def build_apo_delta_summary(runs: List[RunData], apo: Optional[RunData]) -> pd.DataFrame:
    if apo is None:
        return pd.DataFrame()
    metric_cols = ["protein_rmsd_mean_A", "protein_rmsd_last_half_mean_A", "protein_rmsd_final_A", "protein_rmsf_mean_A", "rg_mean_A", "rg_last_half_mean_A", "rg_final_A", "protein_rmsd_slope_A_per_ns", "rg_slope_A_per_ns"]
    rows = []
    for r in runs:
        if r.name == apo.name:
            continue
        row = {"Run": r.name, "Ligand": r.ligand_display, "ApoRun": apo.name}
        for m in metric_cols:
            a, b = apo.metrics.get(m), r.metrics.get(m)
            row[f"{m}_apo"] = a
            row[f"{m}_bound"] = b
            row[f"delta_{m}"] = safe_delta(b, a)
        row["protein_state_call"] = classify_protein_state_change(row)
        rows.append(row)
    return pd.DataFrame(rows)


def build_rmsf_delta_long(runs: List[RunData], apo: Optional[RunData]) -> pd.DataFrame:
    if apo is None or apo.protein_rmsf is None:
        return pd.DataFrame()
    apo_map = dict(zip(apo.protein_rmsf["ResidueKey"].astype(str), apo.protein_rmsf["RMSF_A"].astype(float)))
    rows = []
    for r in runs:
        if r.name == apo.name or r.protein_rmsf is None:
            continue
        for _, rr in r.protein_rmsf.iterrows():
            key = str(rr["ResidueKey"])
            if key not in apo_map:
                continue
            apo_val, bound_val = float(apo_map[key]), float(rr["RMSF_A"])
            rows.append({"Run": r.name, "Ligand": r.ligand_display, "Residue": rr["Residue"], "ResidueKey": key, "Apo_RMSF_A": apo_val, "Bound_RMSF_A": bound_val, "Delta_RMSF_A": bound_val - apo_val, "Direction": "More flexible" if bound_val > apo_val else "More rigid", "Abs_Delta_RMSF_A": abs(bound_val - apo_val)})
    out = pd.DataFrame(rows)
    return out.sort_values(["Run", "Abs_Delta_RMSF_A"], ascending=[True, False]).reset_index(drop=True) if not out.empty else out


def build_binding_site_rmsf_delta(runs: List[RunData], apo: Optional[RunData]) -> pd.DataFrame:
    if apo is None or apo.protein_rmsf is None:
        return pd.DataFrame()
    apo_map = dict(zip(apo.protein_rmsf["ResidueKey"].astype(str), apo.protein_rmsf["RMSF_A"].astype(float)))
    rows = []
    for r in runs:
        if r.is_apo or r.protein_rmsf is None or r.residue_frequency is None or r.residue_frequency.empty:
            continue
        bound_map = dict(zip(r.protein_rmsf["ResidueKey"].astype(str), r.protein_rmsf["RMSF_A"].astype(float)))
        for _, contact in r.residue_frequency.iterrows():
            key = str(contact["ResidueKey"])
            if key not in apo_map or key not in bound_map:
                continue
            apo_val, bound_val = float(apo_map[key]), float(bound_map[key])
            rows.append({
                "Run": r.name,
                "Ligand": r.ligand_display,
                "Residue": contact["Residue"],
                "ResidueKey": key,
                "Chain": contact.get("Chain", ""),
                "ResNum": contact.get("ResNum", np.nan),
                "ContactFraction": contact.get("contact_fraction", np.nan),
                "ContactCount": contact.get("count", np.nan),
                "Apo_RMSF_A": apo_val,
                "Bound_RMSF_A": bound_val,
                "Delta_RMSF_A": bound_val - apo_val,
                "Direction": "stabilized / rigidified" if bound_val < apo_val else "loosened / more flexible",
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["Abs_Delta_RMSF_A"] = pd.to_numeric(out["Delta_RMSF_A"], errors="coerce").abs()
    out["__order"] = out.apply(parse_residue_sort_tuple_from_row, axis=1)
    out = out.sort_values(["Run", "__order"]).drop(columns=["__order"]).reset_index(drop=True)
    return out


def build_binding_site_summary(binding_delta: pd.DataFrame) -> pd.DataFrame:
    if binding_delta is None or binding_delta.empty:
        return pd.DataFrame()
    work = binding_delta.dropna(subset=["Delta_RMSF_A"]).copy()
    if work.empty:
        return pd.DataFrame()
    rows = []
    for (run, lig), grp in work.groupby(["Run", "Ligand"]):
        mean_delta = float(grp["Delta_RMSF_A"].mean())
        weights = pd.to_numeric(grp.get("ContactFraction"), errors="coerce") if "ContactFraction" in grp else None
        weighted = None
        if weights is not None and weights.notna().any() and weights.fillna(0).sum() > 0:
            weighted = float(np.average(pd.to_numeric(grp["Delta_RMSF_A"], errors="coerce").fillna(0), weights=weights.fillna(0)))
        rows.append({"Run": run, "Ligand": lig, "MatchedContactResidues": int(len(grp)), "Mean_BindingSite_Delta_RMSF_A": mean_delta, "Weighted_BindingSite_Delta_RMSF_A": weighted, "MeanAbs_BindingSite_Delta_RMSF_A": float(grp["Delta_RMSF_A"].abs().mean()), "RigidifiedResidues_N": int((grp["Delta_RMSF_A"] < 0).sum()), "MoreFlexibleResidues_N": int((grp["Delta_RMSF_A"] > 0).sum()), "BindingSiteCall": "binding site rigidified/stabilized" if mean_delta < -0.05 else "binding site loosened/more flexible" if mean_delta > 0.05 else "binding site similar to apo"})
    return pd.DataFrame(rows).sort_values("Run").reset_index(drop=True)


def build_contact_matrix(runs: List[RunData], topn: int = 40, min_contact_fraction: float = 0.0) -> pd.DataFrame:
    """
    Build the ligand-contact residue matrix used by the contact heatmap.

    Important dashboard behavior:
      APO/reference runs are intentionally excluded because they do not have a
      ligand-contact profile. Keeping APO in this heatmap makes a zero row that
      visually looks like a biological result, so ligand-only rows are safer and
      easier to interpret.
    """
    ligand_runs = [r for r in runs if not r.is_apo]
    score: Dict[str, float] = {}
    for r in ligand_runs:
        if r.residue_frequency is None or r.residue_frequency.empty:
            continue
        for _, row in r.residue_frequency.iterrows():
            frac = pd.to_numeric(pd.Series([row.get("contact_fraction")]), errors="coerce").iloc[0]
            frac = 0.0 if pd.isna(frac) else float(frac)
            if frac < min_contact_fraction:
                continue
            res = str(row["Residue"])
            score[res] = score.get(res, 0.0) + frac
    residues = sorted(score, key=score.get, reverse=True)
    if topn > 0:
        residues = residues[:topn]
    residues = sorted(residues, key=residue_sort_key)
    mat = pd.DataFrame(index=[r.name for r in ligand_runs], columns=residues, dtype=float).fillna(0.0)
    for r in ligand_runs:
        if r.residue_frequency is None:
            continue
        freq = dict(zip(r.residue_frequency["Residue"], pd.to_numeric(r.residue_frequency["contact_fraction"], errors="coerce").fillna(0.0)))
        for res in residues:
            mat.loc[r.name, res] = float(freq.get(res, 0.0))
    return mat

def build_interaction_type_matrix(runs: List[RunData]) -> pd.DataFrame:
    types = sorted({str(t) for r in runs if r.type_fraction is not None for t in r.type_fraction["InteractionType"].tolist() if not is_excluded_interaction_type(t)})
    mat = pd.DataFrame(index=[r.name for r in runs], columns=types, dtype=float).fillna(0.0)
    for r in runs:
        if r.type_fraction is None:
            continue
        freq = {str(k): float(v) for k, v in zip(r.type_fraction["InteractionType"], pd.to_numeric(r.type_fraction["fraction"], errors="coerce").fillna(0.0)) if not is_excluded_interaction_type(k)}
        for t in types:
            mat.loc[r.name, t] = float(freq.get(t, 0.0))
    return mat


def jaccard(a: set, b: set) -> float:
    return float(len(a & b) / len(a | b)) if (a | b) else 0.0


def build_pairwise_ligand_table(runs: List[RunData]) -> pd.DataFrame:
    lig_runs = [r for r in runs if not r.is_apo]
    rows = []
    for i, a in enumerate(lig_runs):
        for b in lig_runs[i + 1:]:
            a_res = set(a.residue_frequency["ResidueKey"].astype(str)) if a.residue_frequency is not None else set()
            b_res = set(b.residue_frequency["ResidueKey"].astype(str)) if b.residue_frequency is not None else set()
            a_types = {str(k): float(v) for k, v in zip(a.type_fraction["InteractionType"], a.type_fraction["fraction"]) if not is_excluded_interaction_type(k)} if a.type_fraction is not None else {}
            b_types = {str(k): float(v) for k, v in zip(b.type_fraction["InteractionType"], b.type_fraction["fraction"]) if not is_excluded_interaction_type(k)} if b.type_fraction is not None else {}
            all_types = sorted(set(a_types) | set(b_types))
            avec = np.array([float(a_types.get(t, 0.0) or 0.0) for t in all_types])
            bvec = np.array([float(b_types.get(t, 0.0) or 0.0) for t in all_types])
            denom = np.linalg.norm(avec) * np.linalg.norm(bvec)
            rows.append({"Ligand_A_Run": a.name, "Ligand_A": a.ligand_display, "Ligand_B_Run": b.name, "Ligand_B": b.ligand_display, "ContactResidue_Jaccard": jaccard(a_res, b_res), "SharedContactResidues_N": len(a_res & b_res), "UniqueContactResidues_A_N": len(a_res - b_res), "UniqueContactResidues_B_N": len(b_res - a_res), "InteractionType_CosineSimilarity": float(np.dot(avec, bvec) / denom) if denom else 0.0, "Delta_ProteinRMSDMean_A_minus_B": safe_delta(a.metrics.get("protein_rmsd_mean_A"), b.metrics.get("protein_rmsd_mean_A")), "Delta_RgMean_A_minus_B": safe_delta(a.metrics.get("rg_mean_A"), b.metrics.get("rg_mean_A")), "Delta_LigandRMSDMean_A_minus_B": safe_delta(a.metrics.get("ligand_rmsd_mean_A"), b.metrics.get("ligand_rmsd_mean_A"))})
    return pd.DataFrame(rows)


def build_stability_flags(runs: List[RunData], apo_delta: pd.DataFrame, binding_summary: pd.DataFrame) -> pd.DataFrame:
    delta_lookup = {row["Run"]: row for _, row in apo_delta.iterrows()} if apo_delta is not None and not apo_delta.empty else {}
    bind_lookup = {row["Run"]: row for _, row in binding_summary.iterrows()} if binding_summary is not None and not binding_summary.empty else {}
    rows = []
    for r in runs:
        if r.is_apo:
            continue
        possible, reasons = False, []
        ligand_mean, ligand_final = r.metrics.get("ligand_rmsd_mean_A"), r.metrics.get("ligand_rmsd_final_A")
        load_ratio = r.metrics.get("interaction_load_ratio_last_over_first")
        try:
            if ligand_final is not None and float(ligand_final) >= 8.0:
                possible = True; reasons.append("high final ligand RMSD")
        except Exception: pass
        try:
            if ligand_mean is not None and float(ligand_mean) >= 5.0:
                possible = True; reasons.append("high mean ligand RMSD")
        except Exception: pass
        try:
            if load_ratio is not None and float(load_ratio) <= 0.50:
                possible = True; reasons.append("interaction load collapsed")
        except Exception: pass
        protein_call = delta_lookup.get(r.name, {}).get("protein_state_call", "no apo delta")
        binding_call = bind_lookup.get(r.name, {}).get("BindingSiteCall", "no binding-site RMSF summary")
        if possible:
            overall = "possible ligand dissociation / unstable binding"
        elif "rigidified" in str(binding_call):
            overall = "stable binding-site stabilizer"
        elif "loosened" in str(binding_call):
            overall = "binding-site flexibility increaser"
        elif "more compact" in str(protein_call):
            overall = "protein-compacting binder"
        else:
            overall = "no obvious instability flag"
        rows.append({"Run": r.name, "Ligand": r.ligand_display, "OverallCall": overall, "ProteinStateCall": protein_call, "BindingSiteCall": binding_call, "PossibleDissociation": possible, "FlagReasons": "; ".join(reasons), "LigandRMSDMean_A": ligand_mean, "LigandRMSDFinal_A": ligand_final, "InteractionLoadRatio_LastOverFirst": load_ratio})
    return pd.DataFrame(rows)


def apply_tron_layout(fig: go.Figure, title: str = "", height: int = 520) -> go.Figure:
    fig.update_layout(title=dict(text=title, font=dict(color=TRON_TEXT, size=19)), paper_bgcolor=TRON_BG, plot_bgcolor=TRON_BG, font=dict(color=TRON_TEXT, size=13), margin=dict(l=70, r=36, t=76, b=64), height=height, legend=dict(bgcolor="rgba(0,0,0,0)", borderwidth=0, font=dict(color=TRON_MUTED)), xaxis=dict(showgrid=True, gridcolor=TRON_GRID, zeroline=False, linecolor=TRON_GRID, tickcolor=TRON_GRID, title_font=dict(color=TRON_TEXT), tickfont=dict(color=TRON_MUTED)), yaxis=dict(showgrid=True, gridcolor=TRON_GRID, zeroline=False, linecolor=TRON_GRID, tickcolor=TRON_GRID, title_font=dict(color=TRON_TEXT), tickfont=dict(color=TRON_MUTED)))
    return fig


def fig_to_div(fig: Optional[go.Figure]) -> str:
    return "" if fig is None else pio.to_html(fig, include_plotlyjs=False, full_html=False, config=PLOTLY_CONFIG)


def fig_overlay_time_series(runs: List[RunData], attr: str, ycol: str, title: str, ytitle: str, apo_name: Optional[str] = None) -> Optional[go.Figure]:
    fig = go.Figure(); added = 0
    palette = build_run_palette(runs, apo_name=apo_name)
    for r in runs:
        df = getattr(r, attr)
        if df is None or df.empty or ycol not in df.columns:
            continue
        color = palette.get(r.name, TRON_CYAN)
        width = 4.2 if apo_name and r.name == apo_name else 2.5
        dash = "dashdot" if apo_name and r.name == apo_name else ("dash" if r.is_apo else "solid")
        fig.add_trace(go.Scatter(x=df["time_ns"], y=df[ycol], mode="lines", name=r.name, line=dict(color=color, width=width, dash=dash), hovertemplate=f"Run: {r.name}<br>Time/frame: %{{x:.3f}}<br>{ytitle}: %{{y:.3f}}<extra></extra>"))
        added += 1
    if not added:
        return None
    apply_tron_layout(fig, title); fig.update_xaxes(title="Time (ns) / frame index"); fig.update_yaxes(title=ytitle)
    return fig


def fig_protein_rmsf_overlay(runs: List[RunData], apo_name: Optional[str]) -> Optional[go.Figure]:
    fig = go.Figure(); added = 0
    palette = build_run_palette(runs, apo_name=apo_name)
    ymax = 0.0
    for r in runs:
        if r.protein_rmsf is not None and not r.protein_rmsf.empty:
            ymax = max(ymax, float(pd.to_numeric(r.protein_rmsf["RMSF_A"], errors="coerce").max()))
    ymax = ymax * 1.05 if ymax > 0 else 1.0

    for r in runs:
        df = r.protein_rmsf
        if df is None or df.empty:
            continue
        work = df.copy(); work["__order"] = work["Residue"].map(residue_sort_key); work = work.sort_values("__order")
        x = pd.to_numeric(work["ResNum"], errors="coerce")
        if x.isna().all():
            x = pd.Series(np.arange(len(work)), index=work.index)
        color = palette.get(r.name, TRON_CYAN)
        width = 4.2 if apo_name and r.name == apo_name else 2.1
        dash = "dashdot" if apo_name and r.name == apo_name else ("dash" if r.is_apo else "solid")
        fig.add_trace(go.Scatter(x=x, y=work["RMSF_A"], mode="lines", name=r.name, line=dict(color=color, width=width, dash=dash), customdata=np.stack([work["Residue"].astype(str)], axis=1), hovertemplate="Run: " + r.name + "<br>Residue: %{customdata[0]}<br>RMSF: %{y:.3f} Å<extra></extra>"))
        added += 1

    # Add per-ligand interacting-residue vertical guides, hidden by default but toggleable from the legend.
    for r in runs:
        if r.is_apo or r.residue_frequency is None or r.residue_frequency.empty:
            continue
        tmp = r.residue_frequency.copy()
        xs = pd.to_numeric(tmp.get("ResNum"), errors="coerce")
        if xs.isna().all():
            parsed = tmp["Residue"].map(lambda v: parse_residue_label(v).get("resnum"))
            xs = pd.to_numeric(parsed, errors="coerce")
        xs = sorted({float(v) for v in xs.dropna().tolist()})
        if not xs:
            continue
        fig.add_trace(make_vertical_line_trace(xs, 0.0, ymax, palette.get(r.name, TRON_CYAN), f"{r.name} interacting residues"))

    if not added: return None
    apply_tron_layout(fig, "Protein RMSF by residue")
    fig.update_xaxes(title="Residue number / residue order")
    fig.update_yaxes(title="RMSF (Å)", range=[0, ymax])
    fig.update_layout(legend=dict(bgcolor="rgba(0,0,0,0)", borderwidth=0, font=dict(color=TRON_MUTED), groupclick="toggleitem"))
    return fig


def fig_metric_delta_bars(apo_delta: pd.DataFrame) -> Optional[go.Figure]:
    if apo_delta is None or apo_delta.empty: return None
    metrics = [("delta_protein_rmsd_last_half_mean_A", "Δ protein RMSD"), ("delta_rg_last_half_mean_A", "Δ Rg"), ("delta_protein_rmsf_mean_A", "Δ protein RMSF")]
    fig = go.Figure(); added = 0
    for i, (col, label) in enumerate(metrics):
        if col not in apo_delta.columns: continue
        vals = pd.to_numeric(apo_delta[col], errors="coerce")
        if vals.dropna().empty: continue
        fig.add_trace(go.Bar(x=apo_delta["Run"], y=vals, name=label, marker_color=COLOR_CYCLE[i % len(COLOR_CYCLE)], hovertemplate=f"{label}<br>Run: %{{x}}<br>Bound - Apo: %{{y:.4f}}<extra></extra>")); added += 1
    if not added: return None
    apply_tron_layout(fig, "Apo-relative global protein-state deltas"); fig.update_layout(barmode="group"); fig.update_xaxes(title="Ligand-bound run"); fig.update_yaxes(title="Bound - Apo"); fig.add_hline(y=0, line_color=TRON_MUTED, line_width=1)
    return fig


def fig_binding_site_delta(binding_delta: pd.DataFrame, topn: int = 60) -> Optional[go.Figure]:
    if binding_delta is None or binding_delta.empty: return None
    work = binding_delta.dropna(subset=["Delta_RMSF_A"]).copy()
    if work.empty: return None
    if len(work) > topn:
        work = work.head(topn).copy()
    work["__order"] = work.apply(parse_residue_sort_tuple_from_row, axis=1)
    work = work.sort_values(["Run", "__order"]).copy()
    colors = [TRON_RED if v > 0 else TRON_GREEN for v in work["Delta_RMSF_A"]]
    fig = go.Figure(go.Bar(x=work["Delta_RMSF_A"], y=work["Run"].astype(str) + " | " + work["Residue"].astype(str), orientation="h", marker_color=colors, hovertemplate="%{y}<br>Δ RMSF: %{x:.3f} Å<extra></extra>"))
    apply_tron_layout(fig, "Ligand-contact residues: RMSF change vs apo", height=max(560, 24 * len(work) + 170)); fig.update_xaxes(title="Δ RMSF at ligand-contact residues (bound - apo, Å)"); fig.update_yaxes(title=""); fig.add_vline(x=0, line_color=TRON_MUTED, line_width=1)
    return fig


def fig_contact_heatmap(mat: pd.DataFrame) -> Optional[go.Figure]:
    if mat is None or mat.empty: return None
    fig = go.Figure(go.Heatmap(z=mat.values, x=list(mat.columns), y=list(mat.index), colorscale="Turbo", colorbar=dict(title="Contact fraction", tickfont=dict(color=TRON_MUTED), title_font=dict(color=TRON_TEXT)), hovertemplate="Run: %{y}<br>Residue: %{x}<br>Contact fraction: %{z:.4f}<extra></extra>"))
    apply_tron_layout(fig, "Ligand-contact residue summary", height=max(460, 42 * len(mat.index) + 190))
    return fig


def fig_type_fraction_stack(type_mat: pd.DataFrame) -> Optional[go.Figure]:
    if type_mat is None or type_mat.empty: return None
    clean_cols = [c for c in type_mat.columns if not is_excluded_interaction_type(c)]
    type_mat = type_mat[clean_cols]
    if type_mat.empty: return None
    fig = go.Figure()
    for i, typ in enumerate(type_mat.columns):
        fig.add_trace(go.Bar(x=type_mat.index, y=type_mat[typ], name=typ, marker_color=TYPE_COLORS.get(typ, COLOR_CYCLE[i % len(COLOR_CYCLE)]), hovertemplate=f"Type: {typ}<br>Run: %{{x}}<br>Fraction: %{{y:.4f}}<extra></extra>"))
    apply_tron_layout(fig, "Interaction type composition"); fig.update_layout(barmode="stack"); fig.update_xaxes(title="Run"); fig.update_yaxes(title="Fraction")
    return fig


def fig_ligand_summary_scatter(summary: pd.DataFrame) -> Optional[go.Figure]:
    if summary is None or summary.empty or "is_apo" not in summary.columns: return None
    work = summary[summary["is_apo"] == False].copy()
    if work.empty: return None
    x = pd.to_numeric(work.get("ligand_rmsd_last_half_mean_A"), errors="coerce")
    if x.dropna().empty: x = pd.to_numeric(work.get("ligand_rmsd_mean_A"), errors="coerce")
    y = pd.to_numeric(work.get("protein_rmsd_last_half_mean_A"), errors="coerce")
    if y.dropna().empty: y = pd.to_numeric(work.get("protein_rmsd_mean_A"), errors="coerce")
    if x.dropna().empty or y.dropna().empty: return None
    size = pd.to_numeric(work.get("unique_contact_residues"), errors="coerce").fillna(5) * 3 + 14
    color = pd.to_numeric(work.get("rg_last_half_mean_A"), errors="coerce")
    fig = go.Figure(go.Scatter(x=x, y=y, mode="markers+text", text=work["Run"], textposition="top center", marker=dict(size=size, color=color, colorscale="Turbo", showscale=True, colorbar=dict(title="Rg", tickfont=dict(color=TRON_MUTED), title_font=dict(color=TRON_TEXT)), line=dict(color=TRON_WHITE, width=1)), hovertemplate="Run: %{text}<br>Ligand RMSD: %{x:.3f}<br>Protein RMSD: %{y:.3f}<extra></extra>"))
    apply_tron_layout(fig, "Ligand behavior vs protein-state behavior"); fig.update_xaxes(title="Ligand RMSD, last-half mean if available (Å)"); fig.update_yaxes(title="Protein RMSD, last-half mean if available (Å)")
    return fig


def fig_similarity_heatmap(mat: pd.DataFrame, title: str) -> Optional[go.Figure]:
    sim = cosine_similarity_matrix(mat)
    if sim is None or sim.empty: return None
    fig = go.Figure(go.Heatmap(z=sim.values, x=list(sim.columns), y=list(sim.index), colorscale="Viridis", colorbar=dict(title="Cosine", tickfont=dict(color=TRON_MUTED), title_font=dict(color=TRON_TEXT)), hovertemplate="Run A: %{y}<br>Run B: %{x}<br>Similarity: %{z:.4f}<extra></extra>"))
    apply_tron_layout(fig, title, height=max(440, 42 * len(sim.index) + 190))
    return fig


def fig_cross_overlay_time_series(runs: List[RunData], long_df: pd.DataFrame, title: str, ytitle: str) -> Optional[go.Figure]:
    if long_df is None or long_df.empty:
        return None
    fig = go.Figure()
    palette = build_cross_species_palette(runs)
    added = 0
    for run in runs:
        sub = long_df[long_df["run_name"] == run.name].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("time_ns")
        fig.add_trace(go.Scatter(
            x=pd.to_numeric(sub["time_ns"], errors="coerce"),
            y=pd.to_numeric(sub["value"], errors="coerce"),
            mode="lines",
            name=cross_label(run),
            line=dict(color=palette.get(run.name, TRON_CYAN), width=3.0 if run.is_apo else 2.4, dash="dash" if run.is_apo else "solid"),
            opacity=0.95,
            customdata=np.column_stack([
                sub["species"].astype(str),
                sub["state"].astype(str),
                sub["run_name"].astype(str),
                sub["ligand_display"].astype(str),
            ]),
            hovertemplate="Label: %{fullData.name}<br>Run: %{customdata[2]}<br>Species: %{customdata[0]}<br>State: %{customdata[1]}<br>Ligand: %{customdata[3]}<br>Time/frame: %{x:.3f}<br>" + ytitle + ": %{y:.3f}<extra></extra>",
        ))
        added += 1
    if not added:
        return None
    apply_tron_layout(fig, title)
    fig.update_xaxes(title="Time (ns) / frame index")
    fig.update_yaxes(title=ytitle)
    fig.update_layout(legend=cross_species_legend_layout())
    return fig


def fig_cross_rmsf_overlay(runs: List[RunData], long_df: pd.DataFrame, title: str, x_mode: str = "residue", add_contact_guides: bool = False) -> Optional[go.Figure]:
    if long_df is None or long_df.empty:
        return None
    fig = go.Figure()
    palette = build_cross_species_palette(runs)
    added = 0
    for run in runs:
        sub = long_df[long_df["run_name"] == run.name].copy()
        if sub.empty:
            continue
        if x_mode == "normalized":
            sub["__x"] = pd.to_numeric(sub["normalized_position"], errors="coerce")
        else:
            resnum = pd.to_numeric(sub["residue_num"], errors="coerce")
            residx = pd.to_numeric(sub["residue_index"], errors="coerce")
            sub["__x"] = resnum.where(resnum.notna(), residx)
        sub["__y"] = pd.to_numeric(sub["rmsf_A"], errors="coerce")
        sub = sub.dropna(subset=["__x", "__y"]).sort_values("__x")
        if sub.empty:
            continue
        fig.add_trace(go.Scatter(
            x=sub["__x"],
            y=sub["__y"],
            mode="lines",
            name=cross_label(run),
            line=dict(color=palette.get(run.name, TRON_CYAN), width=3.0 if run.is_apo else 2.4, dash="dash" if run.is_apo else "solid"),
            opacity=0.95,
            customdata=np.column_stack([
                sub["run_name"].astype(str),
                sub["species"].astype(str),
                sub["state"].astype(str),
                sub["residue_raw"].astype(str),
            ]),
            hovertemplate="Label: %{fullData.name}<br>Run: %{customdata[0]}<br>Species: %{customdata[1]}<br>State: %{customdata[2]}<br>Residue: %{customdata[3]}<br>RMSF: %{y:.3f} Å<extra></extra>",
        ))
        added += 1
    if not added:
        return None
    ymax = 0.0
    for trace in fig.data:
        yvals = pd.to_numeric(pd.Series(trace.y), errors="coerce")
        if yvals.notna().any():
            ymax = max(ymax, float(yvals.max()))
    if x_mode == "residue" and add_contact_guides:
        for run in runs:
            sub = long_df[long_df["run_name"] == run.name].copy()
            records = build_contact_hover_records_for_run(run, sub, x_mode=x_mode)
            if records is None or records.empty:
                continue
            fig.add_trace(make_contact_residue_hover_trace(records, palette.get(run.name, TRON_CYAN), f"{cross_label(run)} interacting residues", ymax if ymax > 0 else 1.0))
    apply_tron_layout(fig, title)
    fig.update_xaxes(title="Normalized residue position (0-1)" if x_mode == "normalized" else "Residue number when available, otherwise residue index")
    fig.update_yaxes(title="RMSF (Å)")
    fig.update_layout(legend=cross_species_legend_layout())
    return fig


def fig_cross_metric_bars(summary: pd.DataFrame, runs: List[RunData]) -> Optional[go.Figure]:
    if summary is None or summary.empty:
        return None
    metrics = [
        ("protein_rmsd_last_half_mean_A", "Protein RMSD, last-half mean (Å)"),
        ("protein_rmsf_mean_A", "Protein RMSF mean (Å)"),
        ("rg_last_half_mean_A", "Radius of gyration, last-half mean (Å)"),
        ("ligand_rmsd_last_half_mean_A", "Ligand RMSD, last-half mean (Å)"),
        ("contact_residue_count", "Contact residue count"),
    ]
    lookup = {r.name: r for r in runs}
    palette = build_cross_species_palette(runs)
    work = summary.copy()
    work["run_name"] = work["run_name"].astype(str)
    work["label"] = work["run_name"].map(lambda name: cross_label(lookup[name]) if name in lookup else str(name))
    fig = go.Figure()
    traces_meta: List[Tuple[str, str]] = []
    for col, label in metrics:
        if col not in work.columns:
            continue
        sub = work.dropna(subset=[col]).copy()
        if sub.empty:
            continue
        sub = sub.sort_values(["species", "is_apo", "label"])
        fig.add_trace(go.Bar(
            x=sub["label"],
            y=pd.to_numeric(sub[col], errors="coerce"),
            name=label,
            visible=(len(traces_meta) == 0),
            marker_color=[palette.get(name, TRON_CYAN) for name in sub["run_name"]],
            opacity=0.95,
            customdata=np.column_stack([
                sub["run_name"].astype(str),
                sub["species"].astype(str),
                sub["state"].astype(str),
                sub["ligand_display"].astype(str),
            ]),
            hovertemplate="Label: %{x}<br>Run: %{customdata[0]}<br>Species: %{customdata[1]}<br>State: %{customdata[2]}<br>Ligand: %{customdata[3]}<br>Value: %{y:.3f}<extra></extra>",
        ))
        traces_meta.append((col, label))
    if not traces_meta:
        return None
    buttons = []
    for i, (col, label) in enumerate(traces_meta):
        visible = [j == i for j in range(len(traces_meta))]
        buttons.append(dict(label=label, method="update", args=[{"visible": visible}, {"title": f"Cross-species metric comparison: {label}", "yaxis": {"title": label}}]))
    apply_tron_layout(fig, f"Cross-species metric comparison: {traces_meta[0][1]}", height=max(520, 28 * len(work) + 220))
    fig.update_xaxes(title="Run")
    fig.update_yaxes(title=traces_meta[0][1])
    fig.update_layout(
        showlegend=False,
        updatemenus=[dict(type="dropdown", direction="down", x=1.0, y=1.18, xanchor="right", yanchor="top", buttons=buttons, bgcolor=TRON_PANEL, bordercolor=TRON_GRID, font=dict(color=TRON_TEXT))],
    )
    return fig


def fig_cross_apo_delta_bars(apo_delta_all: pd.DataFrame, runs: List[RunData]) -> Optional[go.Figure]:
    if apo_delta_all is None or apo_delta_all.empty:
        return None
    metrics = [
        ("delta_protein_rmsd_last_half_mean_A", "Δ protein RMSD vs same-species apo (Å)"),
        ("delta_rg_last_half_mean_A", "Δ Rg vs same-species apo (Å)"),
        ("delta_protein_rmsf_mean_A", "Δ protein RMSF vs same-species apo (Å)"),
    ]
    lookup = {r.name: r for r in runs}
    palette = build_cross_species_palette(runs)
    work = apo_delta_all.copy()
    work["Run"] = work["Run"].astype(str)
    work["label"] = work["Run"].map(lambda name: cross_label(lookup[name]) if name in lookup else str(name))
    fig = go.Figure()
    traces_meta: List[Tuple[str, str]] = []
    for col, label in metrics:
        if col not in work.columns:
            continue
        sub = work.dropna(subset=[col]).copy()
        if sub.empty:
            continue
        sub = sub.sort_values(["species", "label"])
        fig.add_trace(go.Bar(
            x=sub["label"],
            y=pd.to_numeric(sub[col], errors="coerce"),
            name=label,
            visible=(len(traces_meta) == 0),
            marker_color=[palette.get(name, TRON_CYAN) for name in sub["Run"]],
            opacity=0.95,
            customdata=np.column_stack([
                sub["Run"].astype(str),
                sub["species"].astype(str),
                sub["Ligand"].astype(str),
                sub["ApoRun"].astype(str),
            ]),
            hovertemplate="Label: %{x}<br>Run: %{customdata[0]}<br>Species: %{customdata[1]}<br>Ligand: %{customdata[2]}<br>Apo baseline: %{customdata[3]}<br>Delta: %{y:.3f}<extra></extra>",
        ))
        traces_meta.append((col, label))
    if not traces_meta:
        return None
    buttons = []
    for i, (col, label) in enumerate(traces_meta):
        visible = [j == i for j in range(len(traces_meta))]
        buttons.append(dict(label=label, method="update", args=[{"visible": visible}, {"title": f"Cross-species apo-relative comparison: {label}", "yaxis": {"title": label}}]))
    apply_tron_layout(fig, f"Cross-species apo-relative comparison: {traces_meta[0][1]}", height=max(520, 28 * len(work) + 220))
    fig.update_xaxes(title="Ligand-bound run")
    fig.update_yaxes(title=traces_meta[0][1])
    fig.add_hline(y=0, line_color=TRON_MUTED, line_width=1)
    fig.update_layout(
        showlegend=False,
        updatemenus=[dict(type="dropdown", direction="down", x=1.0, y=1.18, xanchor="right", yanchor="top", buttons=buttons, bgcolor=TRON_PANEL, bordercolor=TRON_GRID, font=dict(color=TRON_TEXT))],
    )
    return fig


def fig_cross_contact_count_bar(summary: pd.DataFrame, runs: List[RunData]) -> Optional[go.Figure]:
    if summary is None or summary.empty or "contact_residue_count" not in summary.columns:
        return None
    lookup = {r.name: r for r in runs}
    palette = build_cross_species_palette(runs)
    work = summary[summary["is_apo"] == False].dropna(subset=["contact_residue_count"]).copy()
    if work.empty:
        return None
    work = work.sort_values(["species", "ligand_display", "run_name"])
    work["label"] = work["run_name"].map(lambda name: cross_label(lookup[str(name)]) if str(name) in lookup else str(name))
    fig = go.Figure(go.Bar(
        x=work["label"],
        y=pd.to_numeric(work["contact_residue_count"], errors="coerce"),
        marker_color=[palette.get(str(name), TRON_CYAN) for name in work["run_name"]],
        opacity=0.95,
        customdata=np.column_stack([
            work["run_name"].astype(str),
            work["species"].astype(str),
            work["ligand_display"].astype(str),
        ]),
        hovertemplate="Label: %{x}<br>Run: %{customdata[0]}<br>Species: %{customdata[1]}<br>Ligand: %{customdata[2]}<br>Contact residues: %{y:.0f}<extra></extra>",
    ))
    apply_tron_layout(fig, "Cross-species ligand-contact coverage", height=max(500, 28 * len(work) + 220))
    fig.update_xaxes(title="Ligand-bound run")
    fig.update_yaxes(title="Contact residue count")
    fig.update_layout(showlegend=False)
    return fig


def fig_cross_stability_scatter(summary: pd.DataFrame, runs: List[RunData]) -> Optional[go.Figure]:
    if summary is None or summary.empty:
        return None
    lookup = {r.name: r for r in runs}
    palette = build_cross_species_palette(runs)
    work = summary[summary["is_apo"] == False].copy()
    if work.empty:
        return None
    x = pd.to_numeric(work.get("ligand_rmsd_last_half_mean_A"), errors="coerce")
    if x.dropna().empty:
        x = pd.to_numeric(work.get("ligand_rmsd_mean_A"), errors="coerce")
    y = pd.to_numeric(work.get("protein_rmsd_last_half_mean_A"), errors="coerce")
    if y.dropna().empty:
        y = pd.to_numeric(work.get("protein_rmsd_mean_A"), errors="coerce")
    if x.dropna().empty or y.dropna().empty:
        return None
    size = pd.to_numeric(work.get("contact_residue_count"), errors="coerce").fillna(0)
    size = size.clip(lower=0) * 2.5 + 14
    work["label"] = work["run_name"].map(lambda name: cross_label(lookup[str(name)]) if str(name) in lookup else str(name))
    fig = go.Figure(go.Scatter(
        x=x,
        y=y,
        mode="markers+text",
        text=work["label"],
        textposition="top center",
        marker=dict(size=size, color=[palette.get(str(name), TRON_CYAN) for name in work["run_name"]], opacity=0.95, line=dict(color=TRON_WHITE, width=1)),
        customdata=np.column_stack([
            work["species"].astype(str),
            work["run_name"].astype(str),
            work["ligand_display"].astype(str),
            pd.to_numeric(work.get("rg_last_half_mean_A"), errors="coerce").fillna(np.nan),
            pd.to_numeric(work.get("protein_rmsf_mean_A"), errors="coerce").fillna(np.nan),
            pd.to_numeric(work.get("contact_residue_count"), errors="coerce").fillna(np.nan),
        ]),
        hovertemplate="Label: %{text}<br>Species: %{customdata[0]}<br>Run: %{customdata[1]}<br>Ligand: %{customdata[2]}<br>Ligand RMSD: %{x:.3f} Å<br>Protein RMSD: %{y:.3f} Å<br>Rg: %{customdata[3]:.3f} Å<br>Protein RMSF: %{customdata[4]:.3f} Å<br>Contact count: %{customdata[5]:.0f}<extra></extra>",
    ))
    apply_tron_layout(fig, "Cross-species stability triage view")
    fig.update_xaxes(title="Ligand RMSD, last-half mean if available (Å)")
    fig.update_yaxes(title="Protein RMSD, last-half mean if available (Å)")
    fig.update_layout(showlegend=False)
    return fig


def fig_cross_interaction_type_stack(type_mat: pd.DataFrame) -> Optional[go.Figure]:
    if type_mat is None or type_mat.empty:
        return None
    clean_cols = [c for c in type_mat.columns if not is_excluded_interaction_type(c)]
    work = type_mat[clean_cols].copy()
    if work.empty:
        return None
    fig = go.Figure()
    for i, typ in enumerate(work.columns):
        fig.add_trace(go.Bar(
            x=work.index.tolist(),
            y=pd.to_numeric(work[typ], errors="coerce"),
            name=str(typ),
            marker_color=CROSS_SPECIES_TRACE_COLORS[i % len(CROSS_SPECIES_TRACE_COLORS)],
            opacity=0.95,
            hovertemplate="Type: " + str(typ) + "<br>Run: %{x}<br>Fraction: %{y:.4f}<extra></extra>",
        ))
    apply_tron_layout(fig, "Cross-species interaction-type composition", height=max(500, 26 * len(work.index) + 220))
    fig.update_xaxes(title="Ligand-bound run")
    fig.update_yaxes(title="Fraction")
    fig.update_layout(barmode="stack", legend=cross_species_legend_layout())
    return fig


def make_cross_species_overlay_figures(runs: List[RunData], comps: Dict[str, Any]) -> Dict[str, Optional[go.Figure]]:
    return {
        "protein_rmsd": fig_cross_overlay_time_series(runs, comps.get("rmsd_long"), "Protein RMSD overlay", "Protein RMSD (Å)"),
        "rg": fig_cross_overlay_time_series(runs, comps.get("rg_long"), "Radius of gyration overlay", "Rg (Å)"),
        "protein_rmsf_residue": fig_cross_rmsf_overlay(runs, comps.get("rmsf_long"), "Protein RMSF overlay by residue index / number", x_mode="residue", add_contact_guides=True),
        "cross_apo_delta_bars": fig_cross_apo_delta_bars(comps.get("apo_delta_all"), runs),
        "ligand_rmsd": fig_cross_overlay_time_series(runs, comps.get("ligand_rmsd_long"), "Ligand RMSD overlay", "Ligand RMSD (Å)"),
        "ligand_rmsf": fig_cross_rmsf_overlay(runs, comps.get("ligand_rmsf_long"), "Ligand RMSF overlay", x_mode="residue", add_contact_guides=False),
        "cross_metric_bars": fig_cross_metric_bars(comps.get("summary"), runs),
        "cross_contact_count_bar": fig_cross_contact_count_bar(comps.get("summary"), runs),
        "cross_stability_scatter": fig_cross_stability_scatter(comps.get("summary"), runs),
        "cross_interaction_type_stack": fig_cross_interaction_type_stack(comps.get("cross_interaction_type_matrix")),
    }


GLOBAL_CSS = r'''
:root{--bg:#05060a;--panel:#0b1020;--grid:#1c2a4a;--text:#cfe8ff;--muted:#7aa6d9;--cyan:#35d0ff;--green:#18ff8b;--yellow:#ffe66d;--magenta:#ff2bd6;--purple:#8a5cff;--red:#ff3b5c;--orange:#ff8a3d;--white:#f5fbff}*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;padding:22px;background:radial-gradient(circle at top right,rgba(53,208,255,.13),transparent 30%),radial-gradient(circle at top left,rgba(255,43,214,.09),transparent 28%),linear-gradient(180deg,#05060a,#070a12 48%,#05060a);color:var(--text);font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}a{color:var(--cyan);text-decoration:none}a:hover{text-decoration:underline}.page{max-width:1580px;margin:0 auto}.card{background:linear-gradient(180deg,rgba(11,16,32,.98),rgba(7,11,22,.98));border:1px solid rgba(28,42,74,.95);border-radius:22px;padding:18px;box-shadow:0 0 0 1px rgba(53,208,255,.025),0 16px 42px rgba(0,0,0,.24);margin-bottom:16px}.hero{position:relative;overflow:hidden;border:1px solid rgba(53,208,255,.28);box-shadow:0 0 46px rgba(53,208,255,.10),inset 0 0 58px rgba(53,208,255,.04);padding:24px}.hero:before{content:"";position:absolute;inset:-45%;background:conic-gradient(from 180deg,transparent,rgba(53,208,255,.12),transparent,rgba(255,43,214,.10),transparent);animation:slowSpin 28s linear infinite;opacity:.7;pointer-events:none}.hero>*{position:relative;z-index:1}@keyframes slowSpin{to{transform:rotate(360deg)}}h1,h2,h3,h4{margin:0 0 10px 0}h1{font-size:clamp(30px,4.4vw,60px);line-height:.98;letter-spacing:-.055em}h2{font-size:22px;letter-spacing:-.02em}.sub,.muted{color:var(--muted)}.eyebrow{color:var(--green);font-weight:800;letter-spacing:.16em;text-transform:uppercase;font-size:12px;margin-bottom:8px}.heroGrid{display:grid;grid-template-columns:minmax(0,1fr) 270px;gap:20px;align-items:center}.brandOrb{display:flex;flex-direction:column;align-items:center;justify-content:center;gap:10px;min-height:220px;border-radius:28px;background:radial-gradient(circle at center,rgba(53,208,255,.16),rgba(138,92,255,.07),rgba(0,0,0,.10) 68%);border:1px solid rgba(53,208,255,.24);box-shadow:inset 0 0 55px rgba(53,208,255,.08),0 0 42px rgba(53,208,255,.06)}.brandOrb img{width:min(185px,82%);height:auto;filter:drop-shadow(0 0 18px rgba(53,208,255,.55)) drop-shadow(0 0 34px rgba(255,43,214,.28));animation:logoPulse 4.8s ease-in-out infinite}.brandOrb .orbText{font-size:12px;color:var(--muted);letter-spacing:.10em;text-transform:uppercase;font-weight:800}@keyframes logoPulse{0%,100%{transform:translateY(0) scale(1);filter:drop-shadow(0 0 16px rgba(53,208,255,.50)) drop-shadow(0 0 30px rgba(255,43,214,.22))}50%{transform:translateY(-3px) scale(1.025);filter:drop-shadow(0 0 26px rgba(53,208,255,.72)) drop-shadow(0 0 48px rgba(255,43,214,.34))}}.heroMeta{display:flex;gap:10px;flex-wrap:wrap;margin:14px 0}.metaPill{display:inline-flex;gap:8px;align-items:center;padding:8px 11px;border-radius:999px;background:rgba(53,208,255,.08);border:1px solid rgba(53,208,255,.22);font-size:12px;color:var(--text)}.grid2{display:grid;grid-template-columns:repeat(2,minmax(260px,1fr));gap:16px}.grid3{display:grid;grid-template-columns:repeat(3,minmax(240px,1fr));gap:16px}.grid4{display:grid;grid-template-columns:repeat(4,minmax(190px,1fr));gap:16px}.metricCard{background:linear-gradient(180deg,rgba(16,23,45,.92),rgba(10,15,28,.96));border:1px solid rgba(53,208,255,.18);border-radius:18px;padding:16px;position:relative;overflow:hidden}.metricCard:after{content:"";position:absolute;inset:auto -20% -50% -20%;height:70%;background:radial-gradient(circle at center,rgba(53,208,255,.09),transparent 70%);pointer-events:none}.metricLabel{color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.10em}.metricValue{font-size:26px;font-weight:850;margin-top:6px}.metricHint{color:var(--muted);font-size:12px;margin-top:6px}.brandLinks,.sectionNav,.toggleBar{display:flex;gap:10px;flex-wrap:wrap;margin-top:14px}.buttonLink,.pill{display:inline-flex;align-items:center;gap:8px;background:rgba(53,208,255,.08);border:1px solid rgba(53,208,255,.24);color:var(--text);border-radius:999px;padding:9px 14px;font-size:13px}.buttonLink.primary{background:linear-gradient(90deg,rgba(53,208,255,.18),rgba(138,92,255,.16))}.buttonLink.green{border-color:rgba(24,255,139,.28);background:rgba(24,255,139,.08)}.buttonLink.magenta{border-color:rgba(255,43,214,.28);background:rgba(255,43,214,.08)}.buttonLink.yellow{border-color:rgba(255,230,109,.30);background:rgba(255,230,109,.08)}.toggleBtn{appearance:none;border:1px solid rgba(53,208,255,.22);background:rgba(53,208,255,.08);color:var(--text);padding:9px 12px;border-radius:999px;cursor:pointer}.toggleBtn:hover{background:rgba(53,208,255,.16)}.readerGuide{border:1px solid rgba(24,255,139,.20);background:linear-gradient(180deg,rgba(24,255,139,.06),rgba(53,208,255,.035));}.guideList{display:grid;grid-template-columns:repeat(3,minmax(220px,1fr));gap:12px;margin-top:12px}.guideItem{border:1px solid rgba(53,208,255,.14);border-radius:16px;padding:13px;background:rgba(0,0,0,.18)}.guideItem b{color:var(--white)}.graphCard{padding:0;overflow:hidden}.graphHeader{display:flex;align-items:flex-start;justify-content:space-between;gap:16px;padding:18px 18px 0 18px}.graphKicker{color:var(--green);font-weight:850;letter-spacing:.15em;text-transform:uppercase;font-size:11px;margin-bottom:5px}.graphTitleLine{display:flex;align-items:center;gap:9px;flex-wrap:wrap}.graphHeader p{margin:0 0 6px 0;max-width:980px}.infoTip{position:relative;display:inline-flex;align-items:center;justify-content:center;width:20px;height:20px;border-radius:50%;border:1px solid rgba(53,208,255,.42);background:rgba(53,208,255,.10);color:var(--cyan);font-size:12px;font-weight:900;cursor:help;line-height:1;vertical-align:middle}.tipPanel{position:absolute;left:50%;top:28px;transform:translateX(-50%) translateY(-4px);width:min(430px,calc(100vw - 46px));padding:13px 14px;border-radius:16px;border:1px solid rgba(53,208,255,.28);background:rgba(5,6,10,.98);color:var(--text);box-shadow:0 18px 55px rgba(0,0,0,.55),0 0 32px rgba(53,208,255,.12);font-size:12px;line-height:1.42;text-align:left;opacity:0;visibility:hidden;transition:.16s ease;z-index:50}.infoTip:hover .tipPanel,.infoTip:focus .tipPanel,.infoTip:focus-within .tipPanel{opacity:1;visibility:visible;transform:translateX(-50%) translateY(0)}.tipPanel strong{color:var(--green)}.tipPanel .tipBlock{margin-top:7px}.tableWrap{overflow-x:auto}.search{width:100%;padding:10px 12px;background:rgba(0,0,0,.24);color:var(--text);border:1px solid rgba(28,42,74,.95);border-radius:12px;margin:8px 0 12px 0}table{width:100%;border-collapse:collapse;font-size:13px}th,td{border-bottom:1px solid rgba(28,42,74,.78);padding:8px 10px;text-align:left;vertical-align:top}th{cursor:pointer;color:var(--cyan);position:sticky;top:0;background:#0b1020;z-index:2}.small{font-size:12px}.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,"Liberation Mono",monospace}.monoBlock{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,"Liberation Mono",monospace;white-space:pre-wrap;overflow-wrap:anywhere;background:rgba(0,0,0,.20);padding:12px;border-radius:14px;border:1px solid rgba(28,42,74,.95)}.callout{border-left:3px solid var(--green);padding:12px 14px;background:rgba(24,255,139,.06);border-radius:14px}.warnout{border-left:3px solid var(--yellow);padding:12px 14px;background:rgba(255,230,109,.06);border-radius:14px}hr.sep{border:none;border-top:1px solid rgba(28,42,74,.9);margin:16px 0}details summary{cursor:pointer}.ligandBox{display:flex;flex-direction:column;gap:10px;min-height:100%;}.ligand2dLabel,.ligand3dHead{display:flex;align-items:center;justify-content:space-between;font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);font-weight:800;}.ligandSvgWrap{min-height:330px;display:flex;align-items:center;justify-content:center;padding:8px 4px;border-radius:16px;background:radial-gradient(circle at center,rgba(53,208,255,.045),transparent 62%);border:1px solid rgba(53,208,255,.08);}.pymacs-ligand-svg{width:100%;max-width:660px;height:auto;display:block;overflow:visible;}.lig3dViewer{position:relative;display:block;width:100%;height:320px;border-radius:16px;overflow:hidden;border:1px solid rgba(53,208,255,.18);background:radial-gradient(circle at center,rgba(53,208,255,.05),rgba(0,0,0,.18));box-shadow:inset 0 0 28px rgba(53,208,255,.035);}.lig3dViewer canvas{position:absolute!important;inset:0!important;width:100%!important;height:100%!important;display:block!important;}.ligandMeta{line-height:1.25;}.footerNote{color:var(--muted);font-size:12px;margin-top:14px}@media(max-width:1200px){.grid4{grid-template-columns:repeat(2,minmax(190px,1fr))}.grid3,.guideList{grid-template-columns:repeat(2,minmax(240px,1fr))}.heroGrid{grid-template-columns:1fr}.brandOrb{min-height:180px}}@media(max-width:820px){body{padding:12px}.grid2,.grid3,.grid4,.guideList{grid-template-columns:1fr}.hero{padding:18px}.tipPanel{left:auto;right:-10px;transform:translateY(-4px)}.infoTip:hover .tipPanel,.infoTip:focus .tipPanel,.infoTip:focus-within .tipPanel{transform:translateY(0)}}
/* PyMACS polished dashboard overrides */
.readerGuide,.readerGuide{border:1px solid rgba(24,255,139,.20);background:linear-gradient(180deg,rgba(24,255,139,.06),rgba(53,208,255,.035));}.ligandGrid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px;margin-top:16px}.ligandBox{position:relative;display:flex;flex-direction:column;gap:12px;min-height:100%;padding:16px;border-radius:22px;background:linear-gradient(180deg,rgba(13,21,42,.98),rgba(7,10,20,.98));border:1px solid rgba(53,208,255,.18);box-shadow:inset 0 0 38px rgba(53,208,255,.025),0 18px 42px rgba(0,0,0,.20);overflow:hidden}.ligandBox:before{content:"";position:absolute;inset:0 0 auto 0;height:2px;background:linear-gradient(90deg,var(--cyan),var(--green),var(--magenta));opacity:.72}.ligandHeader{display:flex;align-items:flex-start;justify-content:space-between;gap:12px}.ligandName{font-size:22px;letter-spacing:-.025em;margin-bottom:4px}.ligandRun{letter-spacing:.08em;text-transform:uppercase}.ligandSvgWrap{min-height:285px;background:radial-gradient(circle at center,rgba(53,208,255,.055),transparent 64%),rgba(0,0,0,.10);border:1px solid rgba(53,208,255,.12);box-shadow:inset 0 0 30px rgba(53,208,255,.025)}.pymacs-ligand-svg{max-width:100%;height:auto}.ligandMetaPills{display:flex;flex-wrap:wrap;gap:8px}.ligandPill{display:inline-flex;gap:6px;align-items:center;padding:7px 10px;border-radius:999px;background:rgba(53,208,255,.08);border:1px solid rgba(53,208,255,.16);font-size:12px;color:var(--text)}.ligandPill b{color:var(--green);letter-spacing:.08em;text-transform:uppercase;font-size:10px}.lig3dViewer{height:265px}.ligandDetails{border:1px solid rgba(28,42,74,.95);border-radius:14px;padding:9px 11px;background:rgba(0,0,0,.14)}.ligandDetails summary{font-size:12px;color:var(--muted);font-weight:800;text-transform:uppercase;letter-spacing:.08em}.smilesText{margin-top:8px;overflow-wrap:anywhere}.footerNote{line-height:1.45}.hero .sub{max-width:1050px;font-size:16px;line-height:1.5}.card>p.muted{line-height:1.45;max-width:1120px}
@media(max-width:1320px){.ligandGrid{grid-template-columns:repeat(2,minmax(0,1fr))}}@media(max-width:820px){.ligandGrid{grid-template-columns:1fr}.lig3dViewer{height:255px}.ligandSvgWrap{min-height:245px}}

'''

GLOBAL_JS = r'''
function hookSearch(tableId){const search=document.getElementById(tableId+'_search');const table=document.getElementById(tableId);if(!search||!table)return;search.addEventListener('input',function(){const q=(this.value||'').toLowerCase().trim();table.querySelectorAll('tbody tr').forEach(tr=>{tr.style.display=tr.innerText.toLowerCase().includes(q)?'':'none';});});}function hookSort(tableId){const table=document.getElementById(tableId);if(!table||!table.tHead)return;const ths=Array.from(table.tHead.rows[0].cells);ths.forEach((th,colIndex)=>{th.addEventListener('click',()=>{const tbody=table.tBodies[0];const rows=Array.from(tbody.querySelectorAll('tr'));const asc=!(th.dataset.asc==='1');ths.forEach(x=>x.dataset.asc='');th.dataset.asc=asc?'1':'0';rows.sort((a,b)=>{const ta=a.cells[colIndex].innerText.trim();const tb=b.cells[colIndex].innerText.trim();const na=parseFloat(ta);const nb=parseFloat(tb);const aNum=!Number.isNaN(na)&&ta.match(/^-?\d+(\.\d+)?/);const bNum=!Number.isNaN(nb)&&tb.match(/^-?\d+(\.\d+)?/);if(aNum&&bNum)return asc?na-nb:nb-na;return asc?ta.localeCompare(tb):tb.localeCompare(ta);});rows.forEach(r=>tbody.appendChild(r));});});}function wireTables(ids){ids.forEach(id=>{hookSearch(id);hookSort(id);});}function togglePanel(kind){document.querySelectorAll('[data-panel="'+kind+'"]').forEach(n=>{n.hidden=!n.hidden;});}


function initLigand3DViewers(){
  if (!window.$3Dmol) {
    document.querySelectorAll('.js-lig3d').forEach(el => {
      el.innerHTML = '<div style="padding:12px;color:#7aa6d9;font-size:12px;">3Dmol.js did not load. Check internet access or CDN blocking.</div>';
    });
    return;
  }
  document.querySelectorAll('.js-lig3d').forEach((el) => {
    if (el.dataset.initialized === '1') return;
    el.dataset.initialized = '1';
    const fmt = (el.dataset.format || 'sdf').toLowerCase();
    const b64 = el.dataset.structureB64 || '';
    if (!b64) {
      el.innerHTML = '<div style="padding:12px;color:#7aa6d9;font-size:12px;">No 3D structure payload available.</div>';
      return;
    }
    let txt = '';
    try { txt = decodeURIComponent(escape(window.atob(b64))); }
    catch (e) { try { txt = window.atob(b64); } catch (_) { txt = ''; } }
    if (!txt) {
      el.innerHTML = '<div style="padding:12px;color:#ff8a3d;font-size:12px;">Could not decode 3D structure payload.</div>';
      return;
    }
    try {
      // 3Dmol uses an absolutely positioned canvas. Force the container to be
      // the positioning context so the canvas lands inside the card, not at
      // the top-left of the page.
      el.style.position = 'relative';
      el.style.overflow = 'hidden';
      el.style.minHeight = el.style.minHeight || '320px';
      const viewer = $3Dmol.createViewer(el, {backgroundColor:'#05060a', antialias:true});
      viewer.addModel(txt, fmt);
      viewer.setStyle({}, {stick:{radius:0.16,color:'#ffe66d'}, sphere:{scale:0.22,color:'#ffe66d'}});
      viewer.setStyle({elem:'H'},  {stick:{radius:0.08,color:'#9fb7d5'}, sphere:{scale:0.13,color:'#9fb7d5'}});
      viewer.setStyle({elem:'N'},  {stick:{radius:0.16,color:'#35d0ff'}, sphere:{scale:0.26,color:'#35d0ff'}});
      viewer.setStyle({elem:'O'},  {stick:{radius:0.16,color:'#ff2b2b'}, sphere:{scale:0.26,color:'#ff2b2b'}});
      viewer.setStyle({elem:'F'},  {stick:{radius:0.16,color:'#7cff9b'}, sphere:{scale:0.26,color:'#7cff9b'}});
      viewer.setStyle({elem:'Cl'}, {stick:{radius:0.17,color:'#00ff66'}, sphere:{scale:0.34,color:'#00ff66'}});
      viewer.setStyle({elem:'BR'}, {stick:{radius:0.17,color:'#ff8a3d'}, sphere:{scale:0.36,color:'#ff8a3d'}});
      viewer.setStyle({elem:'Br'}, {stick:{radius:0.17,color:'#ff8a3d'}, sphere:{scale:0.36,color:'#ff8a3d'}});
      viewer.setStyle({elem:'I'},  {stick:{radius:0.17,color:'#d36bff'}, sphere:{scale:0.38,color:'#d36bff'}});
      viewer.setStyle({elem:'P'},  {stick:{radius:0.17,color:'#ffb347'}, sphere:{scale:0.32,color:'#ffb347'}});
      viewer.setStyle({elem:'S'},  {stick:{radius:0.17,color:'#ffd84d'}, sphere:{scale:0.32,color:'#ffd84d'}});
      viewer.zoomTo(); viewer.rotate(18,'y'); viewer.rotate(-8,'x'); viewer.render();
      setTimeout(() => { try { viewer.resize(); viewer.zoomTo(); viewer.render(); } catch(e) {} }, 250);
    } catch (e) {
      console.warn('3D viewer failed:', e);
      el.innerHTML = '<div style="padding:12px;color:#ff8a3d;font-size:12px;">3D viewer failed to parse this ligand. Try SDF/MOL2 input.</div>';
    }
  });
}
document.addEventListener('DOMContentLoaded', initLigand3DViewers);
window.addEventListener('load', initLigand3DViewers);

'''

GLOBAL_CSS += r'''
.miniNavWrap{position:sticky;top:10px;z-index:60;margin-bottom:16px}.miniNav{display:flex;gap:8px;flex-wrap:wrap;padding:10px 12px;border-radius:16px;background:rgba(5,6,10,.88);backdrop-filter:blur(10px);border:1px solid rgba(53,208,255,.18);box-shadow:0 12px 30px rgba(0,0,0,.24)}.miniNav a{display:inline-flex;align-items:center;padding:8px 10px;border-radius:999px;background:rgba(53,208,255,.08);border:1px solid rgba(53,208,255,.16);color:var(--text);font-size:12px}.miniNav a.active{background:rgba(24,255,139,.12);border-color:rgba(24,255,139,.34);color:var(--white)}.graphCard[id],.card[id]{scroll-margin-top:90px}.backTop{position:fixed;right:16px;bottom:16px;z-index:65;display:inline-flex;align-items:center;justify-content:center;padding:11px 13px;border-radius:999px;background:rgba(53,208,255,.14);border:1px solid rgba(53,208,255,.28);color:var(--text);box-shadow:0 10px 28px rgba(0,0,0,.24)}.notePill{margin:0 0 12px 0;padding:10px 12px;border-radius:14px;background:rgba(53,208,255,.06);border:1px solid rgba(53,208,255,.14);font-size:12px;color:var(--muted)}
@media(max-width:820px){.miniNav{padding:8px}.miniNav a{font-size:11px;padding:7px 9px}.backTop{right:10px;bottom:10px}}
'''

GLOBAL_JS += r'''
function initActiveMiniNav(){const links=Array.from(document.querySelectorAll('.miniNav a[href^="#"]'));if(!links.length||!('IntersectionObserver'in window))return;const byId=new Map();links.forEach(link=>{const id=link.getAttribute('href').slice(1);const el=document.getElementById(id);if(el)byId.set(id,{link,el});});const observer=new IntersectionObserver((entries)=>{entries.forEach(entry=>{const id=entry.target.id;const meta=byId.get(id);if(!meta)return;if(entry.isIntersecting){links.forEach(l=>l.classList.remove('active'));meta.link.classList.add('active');}});},{rootMargin:'-25% 0px -60% 0px',threshold:[0,0.2,0.6]});byId.forEach(meta=>observer.observe(meta.el));}
document.addEventListener('DOMContentLoaded', initActiveMiniNav);
window.addEventListener('load', initActiveMiniNav);
'''

def load_time_series_csv(path: Optional[Path], value_candidates: Sequence[str], out_col: str) -> Optional[pd.DataFrame]:
    df = safe_read_csv(path)
    if df is None or df.empty:
        return None
    t = detect_time_ns(df)
    val_col = find_first_col(df, value_candidates)
    if val_col is None:
        excluded = {"time_ns", "time_ps", "time", "frame", "frames", "t_ns", "t_ps"}
        numeric_cols = []
        for c in df.columns:
            if any(k in c.lower() for k in excluded):
                continue
            if pd.to_numeric(df[c], errors="coerce").notna().sum() > 0:
                numeric_cols.append(c)
        val_col = numeric_cols[0] if numeric_cols else None
    if t is None or val_col is None:
        return None
    out = pd.DataFrame({"time_ns": pd.to_numeric(t, errors="coerce"), out_col: pd.to_numeric(df[val_col], errors="coerce")}).dropna()
    return out if not out.empty else None



def load_combined_protein_rmsf(paths: List[Path]) -> Optional[pd.DataFrame]:
    frames = []
    for p in paths:
        sub = load_rmsf_csv(p, source_label=p.stem)
        if sub is not None and not sub.empty:
            frames.append(sub)
    if not frames:
        return None
    combined = pd.concat(frames, ignore_index=True)
    combined = combined.groupby("ResidueKey", as_index=False).agg({"Residue": "first", "Chain": "first", "ResNum": "first", "RMSF_A": "mean", "Source": lambda x: ",".join(sorted(set(map(str, x))))})
    return combined.sort_values("Residue", key=lambda s: s.map(residue_sort_key)).reset_index(drop=True)



def load_rmsf_csv(path: Path, source_label: Optional[str] = None) -> Optional[pd.DataFrame]:
    df = safe_read_csv(path)
    if df is None or df.empty:
        return None
    val_col = find_first_col(df, ["RMSF_A", "RMSF_Å", "RMSF", "rmsf", "CA_RMSF", "Backbone_RMSF", "Protein_RMSF", "Fluctuation", "fluctuation"])
    if val_col is None:
        numeric_cols = [c for c in df.columns if pd.to_numeric(df[c], errors="coerce").notna().sum() > 0]
        if not numeric_cols:
            return None
        val_col = numeric_cols[-1]
    residue_col = find_first_col(df, ["Residue", "residue", "Residue_Label", "residue_label", "Protein_Residue", "Resid", "resid", "ResID", "resSeq", "ResidueID", "Residue_Number", "ResNum", "resnum"])
    chain_col = find_first_col(df, ["Chain", "chain", "Protein_Chain", "protein_chain", "chainID", "segid"])
    resname_col = find_first_col(df, ["ResName", "resname", "Protein_Resname", "protein_resname", "AA", "ResidueName"])
    resnum_col = find_first_col(df, ["ResNum", "resnum", "Residue_Number", "residue_number", "Protein_Resnum", "ResidueID", "Resid", "resid", "resSeq"])
    rows = []
    for _, row in df.iterrows():
        key, label, chain, resnum = residue_key_from_values(row.get(residue_col) if residue_col else None, row.get(chain_col) if chain_col else None, row.get(resname_col) if resname_col else None, row.get(resnum_col) if resnum_col else None)
        rmsf = pd.to_numeric(pd.Series([row.get(val_col)]), errors="coerce").iloc[0]
        if pd.isna(rmsf) or not key:
            continue
        rows.append({"Residue": label, "ResidueKey": key, "Chain": chain, "ResNum": resnum, "RMSF_A": float(rmsf), "Source": source_label or path.name})
    out = pd.DataFrame(rows)
    if out.empty:
        return None
    out = out.groupby("ResidueKey", as_index=False).agg({"Residue": "first", "Chain": "first", "ResNum": "first", "RMSF_A": "mean", "Source": lambda x: ",".join(sorted(set(map(str, x))))})
    return out.sort_values("Residue", key=lambda s: s.map(residue_sort_key)).reset_index(drop=True)


def load_residue_summary(path: Optional[Path]) -> Optional[pd.DataFrame]:
    df = safe_read_csv(path)
    if df is None or df.empty:
        return None
    residue_col = find_first_col(df, ["Residue", "residue", "Residue_Label", "residue_label", "Protein_Residue", "Resid", "resid", "ResID", "ResidueID"])
    chain_col = find_first_col(df, ["Chain", "chain", "Protein_Chain", "protein_chain", "chainID", "segid"])
    resname_col = find_first_col(df, ["ResName", "resname", "Protein_Resname", "protein_resname", "AA", "ResidueName"])
    resnum_col = find_first_col(df, ["ResNum", "resnum", "Residue_Number", "residue_number", "Protein_Resnum", "ResidueID", "Resid", "resid", "resSeq"])
    if residue_col is None and resnum_col is None:
        return None
    frac_col = find_first_col(df, ["contact_fraction", "ContactFraction", "fraction", "Fraction", "Contact_Fraction", "Frequency", "frequency", "ContactFrequency", "contact_frequency", "Percent_Frames", "percent_frames", "Percent", "percent", "%Frames", "% Frames"])
    count_col = find_first_col(df, ["count", "Count", "ContactCount", "contact_count", "frames_contacted", "FramesContacted", "Frames", "frames", "n_contacts", "N_contacts", "Contacts", "contacts"])
    rows = []
    for _, row in df.iterrows():
        key, label, chain, resnum = residue_key_from_values(row.get(residue_col) if residue_col else None, row.get(chain_col) if chain_col else None, row.get(resname_col) if resname_col else None, row.get(resnum_col) if resnum_col else None)
        if not key:
            continue
        frac = None
        if frac_col:
            val = pd.to_numeric(pd.Series([row.get(frac_col)]), errors="coerce").iloc[0]
            if pd.notna(val):
                frac = float(val)
                if frac > 1.5:
                    frac /= 100.0
        count = None
        if count_col:
            val = pd.to_numeric(pd.Series([row.get(count_col)]), errors="coerce").iloc[0]
            if pd.notna(val):
                count = float(val)
        rows.append({"Residue": label, "ResidueKey": key, "Chain": chain, "ResNum": resnum, "contact_fraction": frac, "count": count})
    out = pd.DataFrame(rows)
    if out.empty:
        return None
    if out["contact_fraction"].isna().all():
        if out["count"].notna().any():
            total = out["count"].fillna(0).sum()
            out["contact_fraction"] = out["count"].fillna(0) / total if total else 0.0
        else:
            out["contact_fraction"] = np.nan
    out = out.groupby("ResidueKey", as_index=False).agg({"Residue": "first", "Chain": "first", "ResNum": "first", "contact_fraction": "max", "count": "sum"})
    return out.sort_values(["contact_fraction", "count", "Residue"], ascending=[False, False, True]).reset_index(drop=True)


def frame_id_series_from_df(df: pd.DataFrame) -> Tuple[Optional[pd.Series], str]:
    if df is None or df.empty:
        return None, ""
    for candidates in [
        ["Time_ns", "time_ns", "t_ns", "Time (ns)", "Time(ns)", "ns"],
        ["Time_ps", "time_ps", "t_ps", "Time (ps)", "Time(ps)", "ps"],
        ["Frame", "frame", "Frame_ID", "frame_id", "Frames", "frames"],
    ]:
        col = find_first_col(df, candidates)
        if col:
            return pd.to_numeric(df[col], errors="coerce"), col
    return None, ""


def canonical_residue_label(label: Any, chain: Any = None, resnum: Any = None) -> Tuple[str, str, str, Optional[int]]:
    key, clean_label, clean_chain, clean_resnum = residue_key_from_values(label, chain, None, resnum)
    if not clean_label:
        clean_label = "" if label is None else str(label)
    return clean_label, key, clean_chain, clean_resnum


def standardize_hover_summary_from_residue_frequency(run: RunData, source_label: str = "residue summary") -> Optional[pd.DataFrame]:
    if run.residue_frequency is None or run.residue_frequency.empty:
        return None
    rows = []
    for _, row in run.residue_frequency.iterrows():
        residue_label, residue_key, chain, resnum = canonical_residue_label(row.get("Residue"), row.get("Chain"), row.get("ResNum"))
        frac = pd.to_numeric(pd.Series([row.get("contact_fraction")]), errors="coerce").iloc[0]
        count = pd.to_numeric(pd.Series([row.get("count")]), errors="coerce").iloc[0]
        rows.append({
            "residue_label": residue_label,
            "ResidueKey": residue_key,
            "Chain": chain,
            "ResNum": resnum,
            "contact_fraction": float(frac) if pd.notna(frac) else np.nan,
            "contact_percent": float(frac) * 100.0 if pd.notna(frac) else np.nan,
            "contact_count": float(count) if pd.notna(count) else np.nan,
            "contact_frames": np.nan,
            "total_frames": np.nan,
            "mean_distance_A": np.nan,
            "min_distance_A": np.nan,
            "top_ligand_atom": None,
            "top_protein_atom": None,
            "source": source_label,
            "source_method": "summary_table",
        })
    out = pd.DataFrame(rows)
    return out.drop_duplicates(subset=["ResidueKey"]).reset_index(drop=True) if not out.empty else None


def load_residue_hover_summary_from_framewise(path: Optional[Path], max_mb: float = 250.0) -> Optional[pd.DataFrame]:
    size = file_size_mb(path)
    if path is None or size is None or size > max_mb:
        return None
    names = read_header(path)
    if not names:
        return None
    wanted = [
        "Time_ns", "time_ns", "Time_ps", "time_ps", "Frame", "frame",
        "Residue", "residue", "ResidNum", "Residue_Number", "ResNum", "resnum",
        "Chain", "chain", "LigAtomName", "ligatomname", "ProtAtomName", "protatomname",
        "Distance", "distance", "Type", "type",
    ]
    usecols = [c for c in names if c in wanted]
    df = safe_read_csv(path, usecols=usecols if usecols else None)
    if df is None or df.empty:
        return None
    frame_ids, frame_col = frame_id_series_from_df(df)
    if frame_ids is not None:
        df["__frame_id"] = frame_ids
        total_frames = int(df["__frame_id"].dropna().nunique())
        source_method = "unique_time"
    else:
        df["__frame_id"] = np.arange(len(df))
        total_frames = int(len(df))
        source_method = "count_estimate"
    residue_col = find_first_col(df, ["Residue", "residue"])
    chain_col = find_first_col(df, ["Chain", "chain"])
    resnum_col = find_first_col(df, ["ResidNum", "Residue_Number", "ResNum", "resnum"])
    dist_col = find_first_col(df, ["Distance", "distance"])
    lig_atom_col = find_first_col(df, ["LigAtomName", "ligatomname"])
    prot_atom_col = find_first_col(df, ["ProtAtomName", "protatomname"])
    if residue_col is None and resnum_col is None:
        return None

    def top_value(series: pd.Series) -> Optional[str]:
        ss = series.dropna().astype(str).str.strip()
        ss = ss[ss != ""]
        return None if ss.empty else str(ss.value_counts().idxmax())

    rows = []
    for _, grp in df.groupby(df.apply(lambda row: residue_key_from_values(row.get(residue_col) if residue_col else None, row.get(chain_col) if chain_col else None, None, row.get(resnum_col) if resnum_col else None)[0], axis=1)):
        first = grp.iloc[0]
        residue_label, residue_key, chain, resnum = canonical_residue_label(first.get(residue_col) if residue_col else None, first.get(chain_col) if chain_col else None, first.get(resnum_col) if resnum_col else None)
        contact_frames = int(grp["__frame_id"].dropna().nunique()) if "__frame_id" in grp else int(len(grp))
        contact_fraction = (float(contact_frames) / float(total_frames)) if total_frames else np.nan
        dist_vals = pd.to_numeric(grp[dist_col], errors="coerce") if dist_col else pd.Series(dtype=float)
        rows.append({
            "residue_label": residue_label,
            "ResidueKey": residue_key,
            "Chain": chain,
            "ResNum": resnum,
            "contact_frames": contact_frames,
            "total_frames": total_frames,
            "contact_count": int(len(grp)),
            "contact_fraction": contact_fraction,
            "contact_percent": contact_fraction * 100.0 if pd.notna(contact_fraction) else np.nan,
            "mean_distance_A": float(dist_vals.mean()) if not dist_vals.dropna().empty else np.nan,
            "min_distance_A": float(dist_vals.min()) if not dist_vals.dropna().empty else np.nan,
            "top_ligand_atom": top_value(grp[lig_atom_col]) if lig_atom_col else None,
            "top_protein_atom": top_value(grp[prot_atom_col]) if prot_atom_col else None,
            "source": path.name,
            "source_method": source_method if frame_col else "count_estimate",
        })
    out = pd.DataFrame(rows)
    return out.sort_values(["contact_fraction", "contact_count", "residue_label"], ascending=[False, False, True]).reset_index(drop=True) if not out.empty else None


def interaction_fraction_column_name(interaction_type: str) -> Optional[str]:
    mapping = {
        "H-bond": "hbond_fraction",
        "Hydrophobic": "hydrophobic_fraction",
        "Ionic": "ionic_fraction",
        "Polar": "polar_fraction",
        "Water bridge": "water_bridge_fraction",
        "Van der Waals": "vdw_fraction",
    }
    return mapping.get(interaction_type)


def load_residue_interaction_summary(path: Optional[Path], framewise_path: Optional[Path] = None, max_mb: float = 250.0) -> Optional[pd.DataFrame]:
    frac_cols = ["hbond_fraction", "hydrophobic_fraction", "ionic_fraction", "polar_fraction", "water_bridge_fraction", "vdw_fraction"]
    size = file_size_mb(framewise_path)
    if framewise_path is not None and size is not None and size <= max_mb:
        names = read_header(framewise_path)
        usecols = [c for c in names if c in {"Time_ns", "time_ns", "Time_ps", "time_ps", "Frame", "frame", "Residue", "residue", "ResNum", "resnum", "ResidNum", "Chain", "chain", "Type", "type"}]
        df = safe_read_csv(framewise_path, usecols=usecols if usecols else None)
        if df is not None and not df.empty:
            frame_ids, frame_col = frame_id_series_from_df(df)
            if frame_ids is not None:
                df["__frame_id"] = frame_ids
                total_frames = int(df["__frame_id"].dropna().nunique())
            else:
                df["__frame_id"] = np.arange(len(df))
                total_frames = int(len(df))
            residue_col = find_first_col(df, ["Residue", "residue"])
            type_col = find_first_col(df, ["Type", "type", "InteractionType", "interaction_type"])
            chain_col = find_first_col(df, ["Chain", "chain"])
            resnum_col = find_first_col(df, ["ResNum", "resnum", "ResidNum"])
            if residue_col and type_col:
                rows = []
                for res_key, grp in df.groupby(df.apply(lambda row: residue_key_from_values(row.get(residue_col), row.get(chain_col) if chain_col else None, None, row.get(resnum_col) if resnum_col else None)[0], axis=1)):
                    first = grp.iloc[0]
                    residue_label, residue_key, chain, resnum = canonical_residue_label(first.get(residue_col), first.get(chain_col) if chain_col else None, first.get(resnum_col) if resnum_col else None)
                    item: Dict[str, Any] = {
                        "Residue": residue_label,
                        "ResidueKey": residue_key,
                        "Chain": chain,
                        "ResNum": resnum,
                        "dominant_interaction_type": None,
                        "interaction_type_list": "",
                        "source": framewise_path.name,
                        "source_method": "framewise_unique_time" if frame_col else "count_estimate",
                    }
                    for c in frac_cols:
                        item[c] = np.nan
                    seen_types = []
                    best_type = None
                    best_frames = -1
                    for typ_raw, typ_grp in grp.groupby(grp[type_col].astype(str).map(normalize_interaction_type)):
                        if is_excluded_interaction_type(typ_raw):
                            continue
                        type_frames = int(typ_grp["__frame_id"].dropna().nunique())
                        frac = (float(type_frames) / float(total_frames)) if total_frames else np.nan
                        out_col = interaction_fraction_column_name(str(typ_raw))
                        if out_col:
                            item[out_col] = frac
                        seen_types.append(str(typ_raw))
                        if type_frames > best_frames:
                            best_frames = type_frames
                            best_type = str(typ_raw)
                    item["dominant_interaction_type"] = best_type
                    item["interaction_type_list"] = ", ".join(sorted(set(seen_types)))
                    rows.append(item)
                out = pd.DataFrame(rows)
                return out.sort_values(["dominant_interaction_type", "Residue"], na_position="last").reset_index(drop=True) if not out.empty else None

    df = safe_read_csv(path)
    if df is None or df.empty:
        return None
    residue_col = find_first_col(df, ["Residue", "residue", "Residue_Label", "residue_label", "Protein_Residue", "Resid", "resid", "ResID", "ResidueID"])
    chain_col = find_first_col(df, ["Chain", "chain"])
    resnum_col = find_first_col(df, ["ResNum", "resnum", "ResidNum", "Residue_Number", "residue_number"])
    if residue_col is None:
        return None
    total_col = find_first_col(df, ["Total", "total", "Count", "count", "Rows", "rows"])
    rows = []
    reserved = {c.lower() for c in [residue_col, chain_col or "", resnum_col or "", total_col or ""]}
    for _, row in df.iterrows():
        residue_label, residue_key, chain, resnum = canonical_residue_label(row.get(residue_col), row.get(chain_col) if chain_col else None, row.get(resnum_col) if resnum_col else None)
        item: Dict[str, Any] = {
            "Residue": residue_label,
            "ResidueKey": residue_key,
            "Chain": chain,
            "ResNum": resnum,
            "dominant_interaction_type": None,
            "interaction_type_list": "",
            "source": path.name if path else "interaction summary",
            "source_method": "summary_table",
        }
        for c in frac_cols:
            item[c] = np.nan
        total = pd.to_numeric(pd.Series([row.get(total_col)]), errors="coerce").iloc[0] if total_col else np.nan
        seen_types = []
        best_val = -np.inf
        best_type = None
        for col in df.columns:
            if col.lower() in reserved:
                continue
            val = pd.to_numeric(pd.Series([row.get(col)]), errors="coerce").iloc[0]
            if pd.isna(val):
                continue
            typ = normalize_interaction_type(col)
            if is_excluded_interaction_type(typ):
                continue
            out_col = interaction_fraction_column_name(typ)
            frac = float(val / total) if pd.notna(total) and float(total) != 0 else float(val)
            if out_col:
                item[out_col] = frac
            seen_types.append(typ)
            if float(val) > best_val:
                best_val = float(val)
                best_type = typ
        item["dominant_interaction_type"] = best_type
        item["interaction_type_list"] = ", ".join(sorted(set(seen_types)))
        rows.append(item)
    out = pd.DataFrame(rows)
    return out.drop_duplicates(subset=["ResidueKey"]).reset_index(drop=True) if not out.empty else None


def load_type_summary(path: Optional[Path]) -> Optional[pd.DataFrame]:
    df = safe_read_csv(path)
    if df is None or df.empty:
        return None
    type_col = find_first_col(df, ["InteractionType", "interaction_type", "Type", "type", "Interaction", "interaction", "ContactType", "contact_type"])
    count_col = find_first_col(df, ["count", "Count", "N", "n", "Rows", "rows", "Contacts", "contacts", "n_contacts", "N_contacts"])
    frac_col = find_first_col(df, ["fraction", "Fraction", "Frac", "frac", "Percent", "percent", "Percentage", "percentage", "Normalized", "normalized"])
    if type_col:
        rows = []
        for _, row in df.iterrows():
            typ = normalize_interaction_type(row.get(type_col))
            count = None
            frac = None
            if count_col:
                val = pd.to_numeric(pd.Series([row.get(count_col)]), errors="coerce").iloc[0]
                if pd.notna(val):
                    count = float(val)
            if frac_col:
                val = pd.to_numeric(pd.Series([row.get(frac_col)]), errors="coerce").iloc[0]
                if pd.notna(val):
                    frac = float(val)
                    if frac > 1.5:
                        frac /= 100.0
            rows.append({"InteractionType": typ, "count": count, "fraction": frac})
        out = pd.DataFrame(rows)
        if out.empty:
            return None
        if out["fraction"].isna().all() and out["count"].notna().any():
            total = out["count"].fillna(0).sum()
            out["fraction"] = out["count"].fillna(0) / total if total else 0.0
        return out.groupby("InteractionType", as_index=False).agg({"count": "sum", "fraction": "sum"}).sort_values("fraction", ascending=False).reset_index(drop=True)
    numeric_cols = [c for c in df.columns if pd.to_numeric(df[c], errors="coerce").notna().sum() > 0]
    rows = []
    for c in numeric_cols:
        low = c.lower()
        if any(k in low for k in ["frame", "time", "residue", "resnum", "resid"]):
            continue
        val = pd.to_numeric(df[c], errors="coerce").dropna()
        if val.empty:
            continue
        rows.append({"InteractionType": normalize_interaction_type(c), "count": float(val.sum()), "fraction": float(val.mean())})
    out = pd.DataFrame(rows)
    if out.empty:
        return None
    total = out["count"].sum()
    if total:
        out["fraction"] = out["count"] / total
    return out.sort_values("fraction", ascending=False).reset_index(drop=True)




GRAPH_HELP: Dict[str, Dict[str, str]] = {
    "protein_rmsd": {"title": "Protein RMSD overlay", "caption": "Global protein structural drift over time for apo/reference and ligand-bound runs.", "means": "RMSD measures how far the protein conformation has moved from the starting/reference structure.", "expected": "A stable simulation usually rises early, then plateaus. A ligand-bound trace lower than apo can suggest global stabilization; higher or continuously rising traces may indicate larger rearrangement or incomplete equilibration.", "read": "Compare plateau height, late-time drift, and whether ligand-bound runs separate from apo. Use this as a global signal, not as direct binding proof."},
    "rg": {"title": "Radius of gyration overlay", "caption": "Protein compactness over time. Lower Rg generally means a more compact protein state.", "means": "Rg summarizes how spread out the protein mass is around its center of mass.", "expected": "Flat Rg suggests stable compactness. A ligand that lowers Rg versus apo may compact the protein; a higher Rg may indicate expansion/opening.", "read": "Interpret with RMSD/RMSF: Rg changes tell you compactness, not necessarily whether the pose is good."},
    "protein_rmsf": {"title": "Protein RMSF overlay", "caption": "Per-residue flexibility profile, with optional interacting-residue guide traces in the legend.", "means": "RMSF measures how much each residue fluctuates across the trajectory.", "expected": "Ligand-stabilized local regions often show reduced RMSF near contact residues. Flexible loops/termini naturally show peaks.", "read": "Focus on contact-adjacent regions and ligand-specific peak changes. Toggle guide traces to see where contacts occur."},
    "protein_rmsf_residue": {"title": "Protein RMSF overlay by residue index / number", "caption": "Cross-species RMSF is shown only in native residue-number/index space. This avoids implying residue equivalence where no ortholog/sequence alignment has been supplied.", "means": "Each trace shows residue-level flexibility for one run using native numbering when available, otherwise residue order.", "expected": "Use this to inspect broad within-run flexibility patterns and compare traces cautiously across species without claiming residue equivalence.", "read": "Toggle interacting-residue guide traces from the legend, then hover guide markers to inspect residue identity and contact fraction without assuming ortholog equivalence."},
    "apo_delta": {"title": "Apo-relative global protein deltas", "caption": "Ligand-bound global metrics minus the apo/reference metric.", "means": "Positive bars mean the ligand-bound run is higher than apo; negative bars mean lower than apo.", "expected": "Negative RMSD/RMSF can imply stabilization; negative Rg suggests compaction. Large positive values may reflect opening, flexibility gain, or drift.", "read": "This is a compact view of how each ligand shifts the global protein state."},
    "binding_delta": {"title": "Ligand-contact residue RMSF deltas vs apo", "caption": "Local flexibility change at residues contacted by ligand, relative to apo.", "means": "Each point/bar is bound-state RMSF minus apo RMSF for a ligand-contacting residue.", "expected": "Negative deltas indicate contact residues became more rigid/stabilized; positive deltas indicate they became more flexible.", "read": "Prioritize high-contact residues with consistent negative deltas when arguing local stabilization."},
    "ligand_rmsd": {"title": "Ligand RMSD overlay", "caption": "Ligand pose movement over time for bound simulations.", "means": "Ligand RMSD tracks how much the ligand changes position/conformation relative to its starting pose.", "expected": "A stable binder often has low-to-moderate RMSD and a plateau. Large final RMSD or monotonic increase can flag pose drift or possible dissociation.", "read": "Use together with contact persistence. A ligand can move internally while maintaining key contacts."},
    "ligand_rmsf": {"title": "Ligand RMSF overlay", "caption": "Ligand-atom or ligand-fragment flexibility overlay when ligand RMSF tables are available.", "means": "This readout summarizes how different ligand positions fluctuate across each trajectory.", "expected": "Stable ligands often show lower overall fluctuations, though flexible substituents can remain mobile even in a retained pose.", "read": "Interpret this together with ligand RMSD and contact persistence rather than as a standalone stability score."},
    "bound_complex_rmsd": {"title": "Bound-complex RMSD overlay", "caption": "Complex-level drift for ligand-bound systems when available.", "means": "This captures movement of the protein-ligand complex as a combined object.", "expected": "Plateau behavior supports equilibration. Divergent traces suggest ligand-specific complex rearrangements.", "read": "Helpful as a bridge between global protein RMSD and ligand RMSD."},
    "contact_heatmap": {"title": "Ligand-contact residue summary heatmap", "caption": "Ligand-only contact fraction matrix. APO is intentionally excluded because it has no ligand-contact profile.", "means": "Color intensity is the fraction of frames in which a residue contacted the ligand.", "expected": "Persistent contacts appear warm/high. Ligands with different binding modes show different residue fingerprints.", "read": "Rows are ligand-bound runs; columns are top contact residues. Use this as a pocket fingerprint, not a direct affinity score."},
    "type_stack": {"title": "Interaction type composition", "caption": "Fractional makeup of contact/interaction categories by run.", "means": "Shows whether interactions are mostly hydrophobic, polar, hydrogen bonding, ionic, water-mediated, or van der Waals-like.", "expected": "Better-behaved poses often preserve a coherent interaction signature instead of collapsing to nonspecific weak contacts.", "read": "Compare composition across ligands to explain why two ligands with similar RMSD may behave differently."},
    "ligand_scatter": {"title": "Ligand behavior vs protein-state behavior", "caption": "Compact ligand-vs-protein behavior map using ligand RMSD, protein RMSD, contacts, and Rg.", "means": "Each point is a ligand-bound run. Position shows ligand movement and protein movement; marker size reflects contact-residue count; color reflects Rg.", "expected": "Promising stable systems often cluster toward lower ligand RMSD and lower/plateaued protein RMSD, with meaningful contact persistence.", "read": "Use as a triage plot to select systems worth deeper structural inspection."},
    "contact_similarity": {"title": "Contact-profile similarity", "caption": "Ligand-vs-ligand cosine similarity based on residue contact fingerprints.", "means": "High similarity means two ligands contacted a similar residue set with similar frequencies.", "expected": "Ligands sharing a binding mode should cluster high; divergent binding modes should be lower.", "read": "Useful for grouping ligands by pocket behavior before selecting representative cases."},
    "type_similarity": {"title": "Interaction-type similarity", "caption": "Ligand-vs-ligand similarity based on interaction type composition.", "means": "High similarity means two ligands rely on similar classes of interactions.", "expected": "Similar scaffolds or poses may share interaction signatures; outliers may reveal different binding chemistry.", "read": "Pair with contact-profile similarity to separate same-pocket/different-chemistry cases from different-pocket cases."},
    "cross_metric_bars": {"title": "Cross-species metric comparison", "caption": "Compares one run-level metric at a time across all detected systems.", "means": "The dropdown switches between run-level metrics without mixing them on one y-axis.", "expected": "This helps compare broad system-level behavior across species without implying residue-level equivalence.", "read": "Use the selector to compare one metric at a time and watch for species- or ligand-specific outliers."},
    "cross_apo_delta_bars": {"title": "Cross-species apo-relative comparison", "caption": "Shows ligand-bound minus same-species APO baseline.", "means": "Each bar compares a ligand-bound run only against its own species-matched apo/reference baseline.", "expected": "Negative protein RMSD/RMSF shifts can suggest stabilization; negative Rg shifts can suggest compaction within that species.", "read": "Interpret these as within-species baseline shifts, not as direct residue-level cross-species equivalence."},
    "cross_contact_count_bar": {"title": "Cross-species ligand-contact coverage", "caption": "Compares ligand-contact coverage, not residue identity.", "means": "Bar height is the number of contact residues summarized for each ligand-bound run.", "expected": "Higher counts can reflect broader pocket engagement, but not necessarily stronger or better binding.", "read": "Use this as a coverage summary rather than as an affinity proxy."},
    "cross_stability_scatter": {"title": "Cross-species stability triage view", "caption": "Triage view of ligand pose movement versus protein drift; not a binding-affinity estimate.", "means": "Each point summarizes one ligand-bound run using ligand RMSD, protein RMSD, and contact coverage.", "expected": "Lower ligand and protein RMSD values can flag stable-looking systems worth deeper inspection, especially with nontrivial contact coverage.", "read": "Use this to prioritize systems for follow-up, not to rank binding affinity."},
    "cross_interaction_type_stack": {"title": "Cross-species interaction-type composition", "caption": "Compares interaction-category composition across ligand-bound runs.", "means": "Each stacked bar shows how contact categories are distributed within one ligand-bound simulation.", "expected": "Runs with similar chemotypes or binding modes may share broad interaction compositions even across species.", "read": "Interpret this as category composition only; it does not imply residue-level conservation."},
    "contact_hover_guides": {"title": "Interacting-residue hover guides", "caption": "Hoverable RMSF guide traces for ligand-contact residues.", "means": "Each dotted guide trace marks residues that contact one ligand in that run, with hover metadata describing contact fraction, residue identity, and dominant interaction type when available.", "expected": "High-contact residues often align with local RMSF features, but cross-species interpretation should remain conservative without a residue mapping.", "read": "Toggle the guide from the legend, then hover the diamond markers to inspect residue-level contact context."},
}


def info_tip_html(key: str) -> str:
    meta = GRAPH_HELP.get(key, {})
    if not meta:
        return ""
    parts = []
    for label, field in [("What it means", "means"), ("Expected trend", "expected"), ("How to read it", "read")]:
        val = meta.get(field, "")
        if val:
            parts.append(f'<div class="tipBlock"><strong>{html_escape(label)}:</strong> {html_escape(val)}</div>')
    return f'<span class="infoTip" tabindex="0" aria-label="Graph explainer">i<span class="tipPanel">{"".join(parts)}</span></span>'


def graph_card(fig: Optional[go.Figure], key: str, panel: str) -> str:
    if fig is None:
        return ""
    meta = GRAPH_HELP.get(key, {"title": key.replace("_", " ").title(), "caption": ""})
    title = meta.get("title", key.replace("_", " ").title())
    caption = meta.get("caption", "")
    return (f'<div id="fig-{html_escape(key)}" class="card graphCard" data-panel="{html_escape(panel)}">'
            f'<div class="graphHeader"><div><div class="graphKicker">Interactive figure</div>'
            f'<div class="graphTitleLine"><h2>{html_escape(title)}</h2>{info_tip_html(key)}</div>'
            f'<p class="small muted">{html_escape(caption)}</p></div></div>'
            f'{fig_to_div(fig)}</div>')


def dashboard_hero_html(title: str, subtitle: str, logo_rel: str, github_url: str, preprint_url: str, contact_email: str, extra_meta: Optional[List[Tuple[str, str]]] = None, home_link: bool = False, data_dir_rel: Optional[str] = None) -> str:
    meta = extra_meta or []
    meta_html = "".join(f'<span class="metaPill"><b>{html_escape(k)}:</b> {html_escape(v)}</span>' for k, v in meta if v is not None)
    nav = '<div class="sectionNav">'
    if home_link:
        nav += '<a class="pill" href="index.html">Home</a>'
    nav += '<a class="pill" href="#overview">Overview</a><a class="pill" href="#guide">How to read</a><a class="pill" href="#tables">Tables</a>'
    if data_dir_rel:
        nav += f'<a class="pill" href="{html_escape(data_dir_rel)}/summary_metrics.csv">CSV exports</a>'
    nav += '</div>'
    return ('<div class="card hero"><div class="heroGrid"><div>'
            '<div class="eyebrow">PyMACS dashboard</div>'
            f'<h1>{html_escape(title)}</h1><p class="sub">{html_escape(subtitle)}</p>'
            f'<div class="heroMeta">{meta_html}</div>{pymacs_brand_links(github_url, preprint_url, contact_email)}{nav}'
            '</div>'
            f'<div class="brandOrb"><img src="{html_escape(logo_rel)}" alt="PyMACS logo"/><div class="orbText">PyMACS</div></div>'
            '</div></div>')


def reading_guide_html(use_framewise: bool) -> str:
    note = "Framewise parsing was enabled to recover missing contact summaries." if use_framewise else "Fast mode uses existing summary files and skips large framewise contact tables unless --use-framewise is requested."
    return f'''
    <div id="guide" class="card readerGuide">
      <div class="eyebrow">Reading guide</div>
      <h2>How to read this dashboard</h2>
      <p class="muted">Use the panels as a linked set of molecular-dynamics readouts: global protein motion, local residue flexibility, ligand pose stability, and contact persistence are most informative when interpreted together.</p>
      <div class="guideList">
        <div class="guideItem"><b>Global protein state</b><p class="small muted">RMSD tracks overall conformational drift, Rg tracks compactness, and RMSF identifies flexible or stabilized residue regions.</p></div>
        <div class="guideItem"><b>Ligand behavior</b><p class="small muted">Ligand RMSD is strongest when interpreted with contact persistence. A stable pose should retain meaningful residue contacts over time.</p></div>
        <div class="guideItem"><b>Contact fingerprints</b><p class="small muted">The contact heatmap includes ligand-bound systems only and summarizes which pocket residues are repeatedly contacted.</p></div>
      </div>
      <p class="footerNote">{html_escape(note)}</p>
    </div>
    '''


def cross_species_nav_html(data_dir_rel: str) -> str:
    return f'''
    <div class="miniNavWrap">
      <nav class="miniNav" aria-label="Cross-species navigation">
        <a href="#overview">Overview</a>
        <a href="#fig-protein_rmsd">Protein RMSD</a>
        <a href="#fig-rg">Rg</a>
        <a href="#fig-protein_rmsf_residue">Protein RMSF</a>
        <a href="#fig-cross_apo_delta_bars">Apo deltas</a>
        <a href="#fig-ligand_rmsd">Ligand RMSD</a>
        <a href="#fig-cross_stability_scatter">Stability scatter</a>
        <a href="#fig-cross_metric_bars">Summary figures</a>
        <a href="#tables">Tables</a>
        <a href="{html_escape(data_dir_rel)}/summary_metrics.csv">Exports</a>
      </nav>
    </div>
    '''

def dataframe_to_html_table(df: pd.DataFrame, table_id: str, max_rows: Optional[int] = None) -> str:
    if df is None or df.empty:
        return '<div class="muted">No data available.</div>'
    work = df.copy()
    more = ""
    if max_rows is not None and len(work) > max_rows:
        work = work.head(max_rows)
        more = f'<div class="small muted">Showing first {max_rows} of {len(df)} rows. Full table exported as CSV.</div>'
    cols = list(work.columns)
    head = "".join(f"<th>{html_escape(c)}</th>" for c in cols)
    rows = []
    for _, row in work.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            cells.append(f"<td>{nice(v,4)}</td>" if isinstance(v, (float, np.floating)) else f"<td>{html_escape(v)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f'<div class="tableWrap"><input id="{table_id}_search" class="search" placeholder="Filter table..." />{more}<table id="{table_id}"><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'


def page_shell(title: str, body: str, bundle_root: Path, current_file: Path, table_ids: List[str]) -> str:
    plotly_rel = rel_href(current_file, bundle_root / "assets" / "plotly.min.js")
    favicon_rel = rel_href(current_file, bundle_root / PYMACS_LOGO_OUTPUT_REL)
    ids_js = json.dumps(table_ids)
    return f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1"/><title>{html_escape(title)}</title><link rel="icon" type="image/png" href="{html_escape(favicon_rel)}"/><link rel="apple-touch-icon" href="{html_escape(favicon_rel)}"/><meta name="theme-color" content="#05060a"/><meta name="description" content="PyMACS comparative molecular dynamics dashboard export."/><script src="{html_escape(plotly_rel)}"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script><style>{GLOBAL_CSS}</style></head><body><div id="top" class="page">{body}</div><script>{GLOBAL_JS}</script><script>wireTables({ids_js});</script></body></html>'

def pymacs_brand_links(github_url: str, preprint_url: str, contact_email: str) -> str:
    subject = "PyMACS comparative MD dashboard"
    body = "Hello PyMACS authors,%0D%0A%0D%0AI have a question about the PyMACS comparative MD dashboard.%0D%0A"
    mailto = f"mailto:{contact_email}?subject={subject.replace(' ', '%20')}&body={body}"
    return f'<div class="brandLinks"><a class="buttonLink primary" href="{html_escape(github_url)}" target="_blank" rel="noopener">GitHub: PyMACS</a><a class="buttonLink green" href="{html_escape(preprint_url)}" target="_blank" rel="noopener">Preprint / Paper</a><a class="buttonLink magenta" href="{html_escape(mailto)}">Contact Authors</a></div>'



def ligand_card_html(run: RunData) -> str:
    lv = run.ligand_visual
    svg = lv.svg if lv.svg else '<div class="muted">No 2D ligand depiction available.</div>'

    pills = []
    if lv.formula:
        pills.append(f"<span class='ligandPill'><b>Formula</b> {html_escape(lv.formula)}</span>")
    if lv.mw is not None:
        pills.append(f"<span class='ligandPill'><b>MW</b> {nice(lv.mw, 2)}</span>")
    if not pills:
        pills.append("<span class='ligandPill muted'>Structure metadata unavailable</span>")

    smiles_html = ""
    if lv.smiles:
        smiles_html = (
            "<details class='ligandDetails'><summary>Chemical identifier</summary>"
            f"<div class='mono small smilesText'>{html_escape(lv.smiles)}</div></details>"
        )

    viewer_id = 'viewer3d_' + slugify(run.name)
    if lv.structure_text_b64 and lv.structure_format:
        viewer = (
            f'<div class="lig3dViewer js-lig3d" id="{viewer_id}" '
            f'data-format="{html_escape(lv.structure_format)}" '
            f'data-structure-b64="{html_escape(lv.structure_text_b64)}"></div>'
        )
    else:
        viewer = '<div class="lig3dViewer muted">No 3D ligand payload available.</div>'

    return f"""
    <article class="ligandBox">
      <div class="ligandHeader">
        <div>
          <h3 class="ligandName">{html_escape(run.ligand_display or run.state)}</h3>
          <div class="ligandRun muted small">{html_escape(run.name)}</div>
        </div>
      </div>
      <div class="ligandSvgWrap">{svg}</div>
      <div class="ligandMetaPills">{''.join(pills)}</div>
      <div class="ligand3dHead"><span>3D view</span><span class="small muted">drag · zoom</span></div>
      {viewer}
      {smiles_html}
    </article>
    """

def metric_cards(summary: pd.DataFrame, apo: Optional[RunData], contact_mat: pd.DataFrame, binding_summary: pd.DataFrame, pairwise: pd.DataFrame) -> str:
    n_runs = len(summary) if summary is not None else 0
    n_lig = int((summary["is_apo"] == False).sum()) if summary is not None and "is_apo" in summary else 0
    mean_pr = numeric_mean(summary["protein_rmsd_mean_A"]) if summary is not None and "protein_rmsd_mean_A" in summary else None
    mean_rg = numeric_mean(summary["rg_mean_A"]) if summary is not None and "rg_mean_A" in summary else None
    mean_lr = numeric_mean(summary["ligand_rmsd_mean_A"]) if summary is not None and "ligand_rmsd_mean_A" in summary else None
    contact_count = len(contact_mat.columns) if contact_mat is not None and not contact_mat.empty else 0
    return f'<div class="grid4"><div class="metricCard"><div class="metricLabel">Runs compared</div><div class="metricValue">{n_runs}</div><div class="metricHint">{n_lig} ligand-bound + apo/reference</div></div><div class="metricCard"><div class="metricLabel">Apo reference</div><div class="metricValue">{html_escape(apo.name if apo else "None")}</div><div class="metricHint">Baseline protein state</div></div><div class="metricCard"><div class="metricLabel">Mean protein RMSD</div><div class="metricValue">{nice(mean_pr,3)}</div><div class="metricHint">Å, across available runs</div></div><div class="metricCard"><div class="metricLabel">Mean Rg</div><div class="metricValue">{nice(mean_rg,3)}</div><div class="metricHint">Å, compactness signal</div></div></div><div class="grid4" style="margin-top:16px;"><div class="metricCard"><div class="metricLabel">Mean ligand RMSD</div><div class="metricValue">{nice(mean_lr,3)}</div><div class="metricHint">Å, ligand pose stability</div></div><div class="metricCard"><div class="metricLabel">Contact residues tracked</div><div class="metricValue">{contact_count}</div><div class="metricHint">From summary files unless --use-framewise</div></div><div class="metricCard"><div class="metricLabel">Binding-site summaries</div><div class="metricValue">{len(binding_summary) if binding_summary is not None else 0}</div><div class="metricHint">Contact residue RMSF vs apo</div></div><div class="metricCard"><div class="metricLabel">Pairwise ligand comparisons</div><div class="metricValue">{len(pairwise) if pairwise is not None else 0}</div><div class="metricHint">Ligand-vs-ligand table rows</div></div></div>'


def cross_metric_cards(summary: pd.DataFrame, contact_mat: pd.DataFrame, stability: pd.DataFrame, pairwise: pd.DataFrame) -> str:
    n_runs = len(summary) if summary is not None else 0
    n_species = int(summary["species"].nunique()) if summary is not None and "species" in summary.columns else 0
    n_lig = int((summary["is_apo"] == False).sum()) if summary is not None and "is_apo" in summary.columns else 0
    n_apo = int((summary["is_apo"] == True).sum()) if summary is not None and "is_apo" in summary.columns else 0
    mean_pr = numeric_mean(summary["protein_rmsd_mean_A"]) if summary is not None and "protein_rmsd_mean_A" in summary.columns else None
    mean_rg = numeric_mean(summary["rg_mean_A"]) if summary is not None and "rg_mean_A" in summary.columns else None
    mean_lr = numeric_mean(summary["ligand_rmsd_mean_A"]) if summary is not None and "ligand_rmsd_mean_A" in summary.columns else None
    contact_count = len(contact_mat.columns) if contact_mat is not None and not contact_mat.empty else 0
    flagged = int(stability["PossibleDissociation"].fillna(False).sum()) if stability is not None and not stability.empty and "PossibleDissociation" in stability.columns else 0
    return f'<div class="grid4"><div class="metricCard"><div class="metricLabel">Species aggregated</div><div class="metricValue">{n_species}</div><div class="metricHint">Combined in one overlay view</div></div><div class="metricCard"><div class="metricLabel">Runs compared</div><div class="metricValue">{n_runs}</div><div class="metricHint">{n_lig} ligand-bound + {n_apo} apo/reference</div></div><div class="metricCard"><div class="metricLabel">Mean protein RMSD</div><div class="metricValue">{nice(mean_pr,3)}</div><div class="metricHint">Å, across available runs</div></div><div class="metricCard"><div class="metricLabel">Mean Rg</div><div class="metricValue">{nice(mean_rg,3)}</div><div class="metricHint">Å, compactness signal</div></div></div><div class="grid4" style="margin-top:16px;"><div class="metricCard"><div class="metricLabel">Mean ligand RMSD</div><div class="metricValue">{nice(mean_lr,3)}</div><div class="metricHint">Å, ligand pose stability</div></div><div class="metricCard"><div class="metricLabel">Contact residues tracked</div><div class="metricValue">{contact_count}</div><div class="metricHint">Union of top ligand-contact residues</div></div><div class="metricCard"><div class="metricLabel">Instability flags</div><div class="metricValue">{flagged}</div><div class="metricHint">Runs flagged as possible dissociation</div></div><div class="metricCard"><div class="metricLabel">Ligand comparisons</div><div class="metricValue">{len(pairwise) if pairwise is not None else 0}</div><div class="metricHint">Cross-run pairwise table rows</div></div></div>'


def build_family_components(runs: List[RunData], apo: Optional[RunData], top_residues: int, min_contact_fraction: float) -> Dict[str, Any]:
    summary = build_summary_table(runs)
    availability = build_availability_table(runs)
    apo_delta = build_apo_delta_summary(runs, apo)
    rmsf_delta = build_rmsf_delta_long(runs, apo)
    binding_delta = build_binding_site_rmsf_delta(runs, apo)
    binding_summary = build_binding_site_summary(binding_delta)
    contact_mat = build_contact_matrix(runs, topn=top_residues, min_contact_fraction=min_contact_fraction)
    type_mat = build_interaction_type_matrix(runs)
    ligand_names = [r.name for r in runs if not r.is_apo]
    contact_mat_lig = contact_mat.loc[[n for n in ligand_names if n in contact_mat.index]] if contact_mat is not None and not contact_mat.empty else pd.DataFrame()
    type_mat_lig = type_mat.loc[[n for n in ligand_names if n in type_mat.index]] if type_mat is not None and not type_mat.empty else pd.DataFrame()
    pairwise = build_pairwise_ligand_table(runs)
    stability = build_stability_flags(runs, apo_delta, binding_summary)
    figs = {
        "protein_rmsd": fig_overlay_time_series(runs, "protein_rmsd", "protein_rmsd_A", "Protein RMSD overlay", "Protein RMSD (Å)", apo.name if apo else None),
        "rg": fig_overlay_time_series(runs, "radius_gyration", "rg_A", "Radius of gyration overlay", "Rg (Å)", apo.name if apo else None),
        "protein_rmsf": fig_protein_rmsf_overlay(runs, apo.name if apo else None),
        "apo_delta": fig_metric_delta_bars(apo_delta),
        "binding_delta": fig_binding_site_delta(binding_delta),
        "ligand_rmsd": fig_overlay_time_series(runs, "ligand_rmsd", "ligand_rmsd_A", "Ligand RMSD overlay", "Ligand RMSD (Å)", None),
        "bound_complex_rmsd": fig_overlay_time_series(runs, "bound_complex_rmsd", "bound_complex_rmsd_A", "Bound-complex RMSD overlay", "Bound-complex RMSD (Å)", None),
        "contact_heatmap": fig_contact_heatmap(contact_mat),
        "type_stack": fig_type_fraction_stack(type_mat),
        "ligand_scatter": fig_ligand_summary_scatter(summary),
        "contact_similarity": fig_similarity_heatmap(contact_mat_lig, "Contact-profile similarity (ligands only)"),
        "type_similarity": fig_similarity_heatmap(type_mat_lig, "Interaction-type similarity (ligands only)"),
    }
    return {"summary": summary, "availability": availability, "apo_delta": apo_delta, "rmsf_delta": rmsf_delta, "binding_delta": binding_delta, "binding_summary": binding_summary, "contact_mat": contact_mat, "type_mat": type_mat, "pairwise": pairwise, "stability": stability, "figs": figs}


def save_df(df: pd.DataFrame, path: Path, index: bool = False) -> None:
    if df is None or df.empty:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=index)


def save_family_exports(data_dir: Path, comps: Dict[str, Any]) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    save_df(comps["summary"], data_dir / "summary_metrics.csv")
    save_df(comps["availability"], data_dir / "data_availability.csv")
    save_df(comps["apo_delta"], data_dir / "apo_relative_global_deltas.csv")
    save_df(comps["rmsf_delta"], data_dir / "all_residue_rmsf_delta_vs_apo.csv")
    save_df(comps["binding_delta"], data_dir / "binding_site_contact_residue_rmsf_delta_vs_apo.csv")
    save_df(comps["binding_summary"], data_dir / "binding_site_summary.csv")
    save_df(comps["pairwise"], data_dir / "ligand_vs_ligand_pairwise_comparison.csv")
    save_df(comps["stability"], data_dir / "ligand_stability_flags.csv")
    save_df(comps["contact_mat"], data_dir / "contact_residue_matrix.csv", index=True)
    save_df(comps["type_mat"], data_dir / "interaction_type_matrix.csv", index=True)


def save_cross_species_exports(data_dir: Path, comps: Dict[str, Any]) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    save_df(comps.get("summary"), data_dir / "summary_metrics.csv")
    save_df(comps.get("availability"), data_dir / "data_availability.csv")
    save_df(comps.get("rmsd_long"), data_dir / "rmsd_long.csv")
    save_df(comps.get("rmsf_long"), data_dir / "rmsf_long.csv")
    save_df(comps.get("rg_long"), data_dir / "rg_long.csv")
    save_df(comps.get("ligand_rmsd_long"), data_dir / "ligand_rmsd_long.csv")
    save_df(comps.get("cross_index"), data_dir / "cross_species_index.csv")
    save_df(comps.get("apo_delta_all"), data_dir / "apo_delta_all.csv")
    save_df(comps.get("binding_summary_all"), data_dir / "binding_site_summary_all.csv")
    save_df(comps.get("cross_metric_summary"), data_dir / "cross_metric_summary.csv")
    save_df(comps.get("cross_interaction_type_matrix"), data_dir / "cross_interaction_type_matrix.csv", index=True)
    save_df(comps.get("contact_hover_long"), data_dir / "contact_hover_long.csv")


def build_family_page(out_html: Path, bundle_root: Path, family: str, runs: List[RunData], apo: Optional[RunData], comps: Dict[str, Any], github_url: str, preprint_url: str, contact_email: str, use_framewise: bool) -> None:
    figs = comps["figs"]
    fig_groups = [
        ("protein", "Protein-state figures", ["protein_rmsd", "rg", "protein_rmsf", "apo_delta", "binding_delta"]),
        ("ligand", "Ligand-behavior figures", ["ligand_rmsd", "bound_complex_rmsd", "contact_heatmap", "type_stack"]),
        ("compare", "Ligand comparison figures", ["ligand_scatter", "contact_similarity", "type_similarity"]),
    ]
    fig_html = ""
    for panel, heading, keys in fig_groups:
        valid_keys = [k for k in keys if figs.get(k) is not None]
        if not valid_keys:
            continue
        fig_html += f'<div class="card" data-panel="{html_escape(panel)}"><div class="eyebrow">{html_escape(panel)} panels</div><h2>{html_escape(heading)}</h2><p class="small muted">Hover over the circular <b>i</b> icon on each graph for interpretation notes and expected trends.</p></div>'
        for key in valid_keys:
            fig_html += graph_card(figs.get(key), key, panel)

    ligand_cards = "".join(ligand_card_html(r) for r in runs if not r.is_apo)
    if not ligand_cards:
        ligand_cards = '<div class="muted">No ligand-bound runs detected or ligand visuals were disabled.</div>'
    ligand_cards_html = f'<div class="card" data-panel="ligand"><div class="eyebrow">Ligand identity</div><h2>Ligand structures</h2><p class="muted">Visual identity and embedded 3D geometry checks for each ligand-bound simulation.</p><div class="ligandGrid">{ligand_cards}</div></div>'

    table_ids = ["summaryTable", "availabilityTable", "apoDeltaTable", "bindingSummaryTable", "bindingDeltaTable", "pairwiseTable", "stabilityTable"]
    run_details = []
    for r in runs:
        metrics_df = pd.DataFrame([{"Metric": k, "Value": v} for k, v in r.metrics.items()])
        tid = "metrics_" + slugify(r.name); table_ids.append(tid)
        run_details.append(f'<details class="card"><summary><b>{html_escape(r.name)}</b> — state={html_escape(r.state)} — ligand={html_escape(r.ligand_display)}</summary><div class="grid2" style="margin-top:14px;"><div><h3>Metrics</h3>{dataframe_to_html_table(metrics_df, tid)}</div><div><h3>Detected files</h3><pre class="monoBlock">{html_escape(json.dumps(r.files, indent=2))}</pre></div></div></details>')

    data_dir_rel = rel_href(out_html, bundle_root / "data" / family)
    logo_rel = rel_href(out_html, bundle_root / PYMACS_LOGO_OUTPUT_REL)
    n_lig = len([r for r in runs if not r.is_apo])
    hero = dashboard_hero_html(
        f"{family} apo-vs-ligand comparison",
        "PyMACS molecular-dynamics dashboard summarizing global protein behavior, local flexibility, ligand pose stability, contact persistence, and ligand-to-ligand differences.",
        logo_rel,
        github_url,
        preprint_url,
        contact_email,
        extra_meta=[("Runs", str(len(runs))), ("Ligand-bound", str(n_lig)), ("Apo/reference", apo.name if apo else "None"), ("Mode", "Framewise" if use_framewise else "Fast")],
        home_link=True,
        data_dir_rel=data_dir_rel,
    )

    body = f'''
    {hero}
    <div class="toggleBar"><button class="toggleBtn" onclick="togglePanel('protein')">Toggle protein panels</button><button class="toggleBtn" onclick="togglePanel('ligand')">Toggle ligand panels</button><button class="toggleBtn" onclick="togglePanel('compare')">Toggle comparison panels</button><button class="toggleBtn" onclick="togglePanel('tables')">Toggle tables</button></div>
    <div id="overview" class="card"><div class="eyebrow">Overview</div><h2>At-a-glance summary</h2>{metric_cards(comps["summary"], apo, comps["contact_mat"], comps["binding_summary"], comps["pairwise"])}</div>
    {reading_guide_html(use_framewise)}
    <div class="grid2"><div class="card callout"><h2>Analysis scope</h2><p>This dashboard compares the apo/reference protein state with ligand-bound simulations using global dynamics, local residue flexibility, and ligand-contact behavior.</p><p class="small muted">In the protein RMSF plot, each ligand includes a hidden-by-default interacting-residue guide trace. Toggle those traces from the legend to connect flexibility changes with contacting residues.</p></div><div class="card warnout"><h2>Data exports</h2><p>Interactive figures are shown in the HTML. Core numeric tables are exported under <span class="mono">{html_escape(data_dir_rel)}</span>, including summary metrics, apo-relative deltas, contact matrices, and ligand comparison tables.</p></div></div>
    {ligand_cards_html}
    {fig_html}
    <div id="tables" class="card" data-panel="tables"><div class="eyebrow">Data tables</div><h2>Summary metrics</h2>{dataframe_to_html_table(comps["summary"], "summaryTable")}<hr class="sep"/><h2>Data availability</h2>{dataframe_to_html_table(comps["availability"], "availabilityTable")}</div>
    <div class="card" data-panel="tables"><h2>Apo-relative global protein deltas</h2>{dataframe_to_html_table(comps["apo_delta"], "apoDeltaTable")}</div>
    <div class="card" data-panel="tables"><h2>Binding-site RMSF summary</h2>{dataframe_to_html_table(comps["binding_summary"], "bindingSummaryTable")}</div>
    <div class="card" data-panel="tables"><h2>Ligand-contact residue RMSF deltas vs apo</h2>{dataframe_to_html_table(comps["binding_delta"], "bindingDeltaTable", max_rows=250)}</div>
    <div class="card" data-panel="tables"><h2>Ligand-vs-ligand comparison</h2>{dataframe_to_html_table(comps["pairwise"], "pairwiseTable")}</div>
    <div class="card" data-panel="tables"><h2>Ligand stability flags</h2>{dataframe_to_html_table(comps["stability"], "stabilityTable")}</div>
    <div class="card" data-panel="tables"><h2>Run details</h2>{''.join(run_details)}</div>'''
    out_html.write_text(page_shell(f"PyMACS Dashboard — {family}", body, bundle_root, out_html, table_ids), encoding="utf-8")


def cross_species_rmsf_warning_html(rmsf_long: pd.DataFrame) -> str:
    if rmsf_long is None or rmsf_long.empty:
        return ""
    grouped = rmsf_long.groupby("run_name")["rmsf_A"].agg(["max", "mean"]).reset_index()
    if grouped.empty:
        return ""
    if bool(((pd.to_numeric(grouped["max"], errors="coerce") > 20.0) | (pd.to_numeric(grouped["mean"], errors="coerce") > 8.0)).any()):
        return '<div class="card warnout" data-panel="protein"><p>One or more RMSF traces are high relative to the rest of the dataset. This may reflect real flexibility, incomplete equilibration, alignment/reference issues, or preprocessing differences. The plot is shown on the full scale for transparency.</p></div>'
    return ""


def build_cross_species_page(out_html: Path, bundle_root: Path, runs: List[RunData], comps: Dict[str, Any], github_url: str, preprint_url: str, contact_email: str, use_framewise: bool) -> None:
    figs = comps["figs"]
    fig_groups = [
        ("protein", "Protein-state figures", ["protein_rmsd", "rg", "protein_rmsf_residue", "cross_apo_delta_bars"]),
        ("ligand", "Ligand-behavior figures", ["ligand_rmsd", "ligand_rmsf", "cross_stability_scatter"]),
        ("summary", "Cross-species summary figures", ["cross_metric_bars", "cross_contact_count_bar", "cross_interaction_type_stack"]),
    ]
    fig_html = ""
    for panel, heading, keys in fig_groups:
        valid_keys = [k for k in keys if figs.get(k) is not None]
        if not valid_keys:
            continue
        panel_note = "Cross-species RMSF is shown only in native residue-number/index space. This avoids implying residue equivalence where no ortholog/sequence alignment has been supplied." if panel == "protein" else "Hover over the circular <b>i</b> icon on each graph for interpretation notes and expected trends."
        fig_html += f'<div class="card" data-panel="{html_escape(panel)}"><div class="eyebrow">{html_escape(panel)} panels</div><h2>{html_escape(heading)}</h2><p class="small muted">{panel_note}</p></div>'
        for key in valid_keys:
            if key == "protein_rmsf_residue":
                fig_html += cross_species_rmsf_warning_html(comps.get("rmsf_long"))
                fig_html += '<div class="card" data-panel="protein"><p class="notePill">Toggle each ligand\'s interacting-residue guide from the legend. Hover guide markers to see residue identity and contact fraction.</p></div>'
            fig_html += graph_card(figs.get(key), key, panel)

    species = sorted({r.family for r in runs})
    scope_note = "Single-species aggregate" if len(species) == 1 else "Cross-species aggregate"
    data_dir_rel = rel_href(out_html, bundle_root / "data" / "cross_species_ALL")
    logo_rel = rel_href(out_html, bundle_root / PYMACS_LOGO_OUTPUT_REL)
    hero = dashboard_hero_html(
        "ALL species comparison",
        "PyMACS molecular-dynamics dashboard overlaying all detected runs across species/families using shared run-level summaries and already-loaded analysis tables.",
        logo_rel,
        github_url,
        preprint_url,
        contact_email,
        extra_meta=[("Scope", scope_note), ("Species", str(len(species))), ("Runs", str(len(runs))), ("Mode", "Framewise" if use_framewise else "Fast")],
        home_link=True,
        data_dir_rel=data_dir_rel,
    )

    table_ids = ["crossSummaryTable", "crossAvailabilityTable", "crossIndexTable", "crossStabilityTable", "crossPairwiseTable"]
    body = f'''
    {hero}
    {cross_species_nav_html(data_dir_rel)}
    <div class="toggleBar"><button class="toggleBtn" onclick="togglePanel('protein')">Toggle protein panels</button><button class="toggleBtn" onclick="togglePanel('ligand')">Toggle ligand panels</button><button class="toggleBtn" onclick="togglePanel('tables')">Toggle tables</button></div>
    <div id="overview" class="card"><div class="eyebrow">Overview</div><h2>At-a-glance summary</h2>{cross_metric_cards(comps["summary"], comps["contact_mat"], comps["stability"], comps["pairwise"])}</div>
    {reading_guide_html(use_framewise)}
    <div class="grid2"><div class="card callout"><h2>Analysis scope</h2><p>This aggregate view overlays all detected runs across species/families without reloading trajectories or forcing framewise contact parsing.</p><p class="small muted">Trace labels use the explicit <span class="mono">SPECIES | LIGAND</span> format, for example <span class="mono">HUMAN | APO</span> and <span class="mono">RAT | UMF504</span>. APO traces are dashed and ligand-bound traces are solid.</p></div><div class="card warnout"><h2>Cross-species interpretation</h2><p>Residue-level cross-species interpretation is conservative unless an explicit alignment or residue-mapping file is supplied. The RMSF plot stays in native residue-number/index space, and interacting-residue guides can be toggled from the legend.</p></div></div>
    {fig_html}
    <div id="tables" class="card" data-panel="tables"><div class="eyebrow">Data tables</div><h2 id="table-summary">Summary metrics</h2>{dataframe_to_html_table(comps["summary"], "crossSummaryTable")}<hr class="sep"/><h2 id="table-availability">Data availability</h2>{dataframe_to_html_table(comps["availability"], "crossAvailabilityTable")}</div>
    <div id="table-index" class="card" data-panel="tables"><h2>Cross-species index</h2>{dataframe_to_html_table(comps["cross_index"], "crossIndexTable")}</div>
    <div class="card" data-panel="tables"><h2>Ligand / protein stability flags</h2>{dataframe_to_html_table(comps["stability"], "crossStabilityTable")}</div>
    <div class="card" data-panel="tables"><h2>Cross-species ligand comparison</h2>{dataframe_to_html_table(comps["pairwise"], "crossPairwiseTable", max_rows=300)}</div>
    <a class="backTop" href="#top">Back to top</a>
    '''
    out_html.write_text(page_shell("PyMACS Dashboard — ALL species comparison", body, bundle_root, out_html, table_ids), encoding="utf-8")


def build_index_page(out_html: Path, bundle_root: Path, payloads: Dict[str, Tuple[List[RunData], Optional[RunData], Dict[str, Any]]], github_url: str, preprint_url: str, contact_email: str, cross_species_payload: Optional[Tuple[List[RunData], Dict[str, Any]]] = None) -> None:
    rows, cards = [], []
    if cross_species_payload is not None:
        cross_runs, cross_comps = cross_species_payload
        cross_summary = cross_comps["summary"]
        n_lig = int((cross_summary["is_apo"] == False).sum()) if "is_apo" in cross_summary.columns else 0
        species_count = int(cross_summary["species"].nunique()) if "species" in cross_summary.columns else len({r.family for r in cross_runs})
        rows.append({"Family": "ALL_SPECIES", "Runs": len(cross_runs), "LigandRuns": n_lig, "ApoReference": "Per-species apo where available", "MeanProteinRMSD_A": numeric_mean(cross_summary["protein_rmsd_mean_A"]) if "protein_rmsd_mean_A" in cross_summary.columns else None, "MeanRg_A": numeric_mean(cross_summary["rg_mean_A"]) if "rg_mean_A" in cross_summary.columns else None, "Dashboard": "cross_species_ALL.html"})
        cards.append(f'<div class="metricCard" style="border-color:rgba(255,230,109,.32);box-shadow:0 0 0 1px rgba(255,230,109,.06),0 16px 42px rgba(0,0,0,.24);"><div class="metricLabel">ALL species comparison</div><div class="metricValue">Aggregate view</div><div class="metricHint">{species_count} species • {len(cross_runs)} runs • {n_lig} ligand-bound</div><div style="margin-top:12px;"><a class="buttonLink yellow" href="cross_species_ALL.html">Open ALL dashboard</a></div></div>')
    for family, (runs, apo, comps) in payloads.items():
        summary = comps["summary"]
        n_lig = int((summary["is_apo"] == False).sum()) if "is_apo" in summary else 0
        rows.append({"Family": family, "Runs": len(runs), "LigandRuns": n_lig, "ApoReference": apo.name if apo else "", "MeanProteinRMSD_A": numeric_mean(summary["protein_rmsd_mean_A"]) if "protein_rmsd_mean_A" in summary else None, "MeanRg_A": numeric_mean(summary["rg_mean_A"]) if "rg_mean_A" in summary else None, "Dashboard": f"family_{slugify(family)}.html"})
        cards.append(f'<div class="metricCard"><div class="metricLabel">Family / species</div><div class="metricValue">{html_escape(family)}</div><div class="metricHint">{len(runs)} runs • {n_lig} ligand-bound • apo={html_escape(apo.name if apo else "none")}</div><div style="margin-top:12px;"><a class="buttonLink primary" href="family_{slugify(family)}.html">Open dashboard</a></div></div>')
    table = pd.DataFrame(rows)
    logo_rel = rel_href(out_html, bundle_root / PYMACS_LOGO_OUTPUT_REL)
    hero = dashboard_hero_html(
        "PyMACS comparative MD dashboard",
        "Interactive molecular-dynamics summaries generated from existing PyMACS analysis outputs. Select a family or species below to inspect apo-vs-ligand dynamics, ligand behavior, contact fingerprints, and exported data tables.",
        logo_rel,
        github_url,
        preprint_url,
        contact_email,
        extra_meta=[("Families", str(len(payloads))), ("Output", "HTML + CSV"), ("Assets", "bundled locally")],
        home_link=False,
        data_dir_rel="data",
    )
    body = f'''
    {hero}
    <div id="overview" class="card"><div class="eyebrow">Available dashboards</div><h2>Families / species</h2><div class="grid3">{"".join(cards)}</div></div>
    <div id="guide" class="card readerGuide"><div class="eyebrow">Reading guide</div><h2>Dashboard contents</h2><div class="guideList"><div class="guideItem"><b>Protein dynamics</b><p class="small muted">RMSD, radius of gyration, and RMSF summarize global motion, compactness, and residue-level flexibility.</p></div><div class="guideItem"><b>Ligand behavior</b><p class="small muted">Ligand RMSD, 2D/3D ligand cards, contact heatmaps, and interaction-type panels summarize pose stability and pocket engagement.</p></div><div class="guideItem"><b>Exported tables</b><p class="small muted">CSV tables provide the numeric values behind the figures, including apo-relative deltas, contact matrices, and ligand comparisons.</p></div></div></div>
    <div class="card"><h2>Bundle summary</h2>{dataframe_to_html_table(table, "familySummaryTable")}</div>
    '''
    out_html.write_text(page_shell("PyMACS Comparative Dashboard", body, bundle_root, out_html, ["familySummaryTable"]), encoding="utf-8")

def ensure_assets(bundle_root: Path, source_root: Optional[Path] = None, logo_url: str = PYMACS_LOGO_URL, logo_file: Optional[str] = None, skip_logo_download: bool = False) -> None:
    assets = bundle_root / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    plotly_js = assets / "plotly.min.js"
    if not plotly_js.exists():
        plotly_js.write_text(get_plotlyjs(), encoding="utf-8")

    images_dir = bundle_root / "images"
    images_dir.mkdir(parents=True, exist_ok=True)
    favicon_path = bundle_root / PYMACS_LOGO_OUTPUT_REL
    canonical_logo_path = images_dir / "Pymacs_Logo.png"

    candidates: List[Path] = []
    if logo_file:
        candidates.append(Path(logo_file).expanduser())
    if source_root:
        candidates.extend([
            source_root / "images" / "Pymacs_Logo.png",
            source_root / "images" / "favicon.png",
            source_root / "Pymacs_Logo.png",
        ])

    copied = False
    for candidate in candidates:
        try:
            if candidate.exists() and candidate.is_file():
                shutil.copyfile(candidate, favicon_path)
                shutil.copyfile(candidate, canonical_logo_path)
                copied = True
                info(f"Bundled PyMACS logo from local file: {candidate}")
                break
        except Exception as e:
            warn(f"Could not copy logo from {candidate}: {e}")

    if not copied and not favicon_path.exists() and not skip_logo_download and logo_url:
        try:
            req = Request(logo_url, headers={"User-Agent": "PyMACS-dashboard-builder/1.0"})
            with urlopen(req, timeout=25) as resp:
                data = resp.read()
            if data:
                favicon_path.write_bytes(data)
                canonical_logo_path.write_bytes(data)
                copied = True
                info(f"Downloaded PyMACS logo to: {favicon_path}")
        except Exception as e:
            warn(f"Could not download PyMACS logo from {logo_url}: {e}")

    if favicon_path.exists() and not canonical_logo_path.exists():
        try:
            shutil.copyfile(favicon_path, canonical_logo_path)
        except Exception:
            pass

    if not favicon_path.exists():
        fallback_png_b64 = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAFgwJ/l4tQzgAAAABJRU5ErkJggg=="
        try:
            favicon_path.write_bytes(base64.b64decode(fallback_png_b64))
            canonical_logo_path.write_bytes(base64.b64decode(fallback_png_b64))
            warn("Using fallback tiny favicon because the PyMACS logo was not available locally and could not be downloaded.")
        except Exception as e:
            warn(f"Could not create fallback favicon: {e}")

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fast PyMACS comparative dashboard builder.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--root", default=".", help="Root directory containing run folders.")
    p.add_argument("--outdir", default=None, help="Output dashboard directory. Default uses timestamp.")
    p.add_argument("--family", default=None, help="Only build one family/species, e.g. HUMAN, RAT, MOUSE.")
    p.add_argument("--cross-species", "--all-species", dest="cross_species", action="store_true", help="Build an additional aggregate dashboard across all detected families/species.")
    p.add_argument("--apo", default=None, help="Explicit apo reference run name, e.g. HUMAN_APO.")
    p.add_argument("--top-residues", type=int, default=40, help="Top contact residues to include in heatmaps.")
    p.add_argument("--min-contact-fraction", type=float, default=0.0, help="Minimum contact fraction for contact heatmap.")
    p.add_argument("--use-framewise", action="store_true", help="Optional: parse framewise contact files if summary files are missing.")
    p.add_argument("--contact-hover-framewise", action="store_true", help="Use framewise contact files to calculate exact residue hover fractions for RMSF guide traces.")
    p.add_argument("--max-hover-framewise-mb", type=float, default=250.0, help="Maximum framewise contact file size to parse for RMSF hover metadata.")
    p.add_argument("--no-ligand-visuals", action="store_true", help="Skip RDKit/OpenBabel ligand 2D rendering from .cgenff.mol2/.mol2 files.")
    p.add_argument("--max-framewise-mb", type=float, default=500.0, help="Safety limit for optional framewise parsing.")
    p.add_argument("--open", action="store_true", help="Open index.html in browser after building.")
    p.add_argument("--github-url", default="https://github.com/schurerlab/Pymacs", help="PyMACS GitHub URL.")
    p.add_argument("--preprint-url", default="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6584500", help="PyMACS preprint/paper URL.")
    p.add_argument("--contact-email", default="jmschulz@med.miami.edu", help="Contact email for dashboard mailto link.")
    p.add_argument("--logo-url", default=PYMACS_LOGO_URL, help="URL for PyMACS logo to bundle as images/favicon.png when no local logo is found.")
    p.add_argument("--logo-file", default=None, help="Optional local logo PNG to copy into each dashboard export as images/favicon.png.")
    p.add_argument("--skip-logo-download", action="store_true", help="Do not attempt to download the PyMACS logo if a local file is not found.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    if not root.exists(): die(f"Root does not exist: {root}")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir).resolve() if args.outdir else (root / f"PyMACS_FAST_Dashboard_{timestamp}")
    outdir.mkdir(parents=True, exist_ok=True)
    ensure_assets(outdir, source_root=root, logo_url=args.logo_url, logo_file=args.logo_file, skip_logo_download=args.skip_logo_download)
    phase("Scanning run folders")
    run_dirs = find_run_dirs(root)
    if args.family:
        ff = args.family.upper()
        run_dirs = [rd for rd in run_dirs if rd.name.upper().startswith(ff + "_")]
    if not run_dirs: die(f"No matching run folders with Analysis_Results/ found under {root}")
    info(f"Detected {len(run_dirs)} run folders")
    phase("Loading runs")
    runs: List[RunData] = []
    for rd in run_dirs:
        try:
            runs.append(load_one_run(rd, use_framewise=args.use_framewise, max_framewise_mb=args.max_framewise_mb, make_ligand_visual=not args.no_ligand_visuals, contact_hover_framewise=args.contact_hover_framewise, max_hover_framewise_mb=args.max_hover_framewise_mb))
        except Exception as e:
            warn(f"Failed loading {rd.name}: {e}")
    if not runs: die("No runs loaded successfully.")
    families = sorted({r.family for r in runs})
    payloads: Dict[str, Tuple[List[RunData], Optional[RunData], Dict[str, Any]]] = {}
    for family in families:
        phase(f"Building family dashboard: {family}")
        fam_runs = sorted([r for r in runs if r.family == family], key=lambda r: (not r.is_apo, r.state, r.name))
        apo = choose_apo(fam_runs, args.apo)
        info(f"Apo reference: {apo.name}" if apo else f"No apo reference detected for {family}; apo-relative deltas skipped.")
        comps = build_family_components(fam_runs, apo=apo, top_residues=args.top_residues, min_contact_fraction=args.min_contact_fraction)
        data_dir = outdir / "data" / family; save_family_exports(data_dir, comps)
        save_df(build_contact_hover_export(fam_runs), data_dir / "contact_hover_summary.csv")
        family_html = outdir / f"family_{slugify(family)}.html"
        build_family_page(family_html, outdir, family, fam_runs, apo, comps, github_url=args.github_url, preprint_url=args.preprint_url, contact_email=args.contact_email, use_framewise=args.use_framewise)
        info(f"Wrote {family_html}")
        payloads[family] = (fam_runs, apo, comps)
    cross_species_payload: Optional[Tuple[List[RunData], Dict[str, Any]]] = None
    if args.cross_species:
        phase("Building cross-species dashboard: ALL species")
        cross_runs = sorted(runs, key=lambda r: (r.family, not r.is_apo, r.state, r.name))
        cross_comps = build_cross_species_components(cross_runs, top_residues=args.top_residues, min_contact_fraction=args.min_contact_fraction)
        cross_data_dir = outdir / "data" / "cross_species_ALL"
        save_cross_species_exports(cross_data_dir, cross_comps)
        cross_html = outdir / "cross_species_ALL.html"
        build_cross_species_page(cross_html, outdir, cross_runs, cross_comps, github_url=args.github_url, preprint_url=args.preprint_url, contact_email=args.contact_email, use_framewise=args.use_framewise)
        info(f"Wrote {cross_html}")
        cross_species_payload = (cross_runs, cross_comps)
    index_html = outdir / "index.html"
    build_index_page(index_html, outdir, payloads, github_url=args.github_url, preprint_url=args.preprint_url, contact_email=args.contact_email, cross_species_payload=cross_species_payload)
    phase("Done")
    print(f"✅ Dashboard index: {index_html}", flush=True)
    print(f"✅ CSV exports:      {outdir / 'data'}", flush=True)
    if args.open: webbrowser.open(index_html.as_uri())


if __name__ == "__main__":
    main()
