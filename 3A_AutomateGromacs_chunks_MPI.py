#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=======================================================================
💠 Advanced MD Trajectory Analyzer — RMSD/RMSF + Pocket Interaction Maps
=======================================================================

Author:  Joseph-Michael Schulz (PhD Candidate, UM Biochemistry & Molecular Biology)
Version: 2.2 — October 2025
Dependencies: MDAnalysis, MDTraj, CuPy (optional), NumPy, Pandas, Matplotlib, Seaborn, NetworkX, tqdm

📘 Overview
-----------
This script performs *comprehensive structural trajectory analysis* of molecular
dynamics simulations, integrating **MDTraj**, **MDAnalysis**, and **CuPy** for GPU acceleration.

It automatically detects the ligand, extracts the binding pocket, and computes
high-resolution RMSD/RMSF profiles, residue contact maps, and protein–ligand
interaction networks — all in a reproducible, step-wise workflow with full multithreading
support and detailed CLI logging.

💡 Key Features
---------------
1️⃣ **Automatic Ligand Detection**
    • Detects the largest non-protein, non-water residue in the system (or use --ligand)
    • Optional fallback via atom indices (--lig-atom-start / --lig-atom-end)

2️⃣ **OpenMP / Thread Optimization**
    • Interactive detection of CPU threads with auto-scaling to 66 % of total cores
    • Environment variable propagation for MKL / OpenBLAS / NumExpr
    • Compatible with headless or bash-pipeline execution

3️⃣ **GPU Acceleration (CuPy-native)**
    • Vectorized ligand–protein distance and contact calculations on GPU
    • Automatic fallback to NumPy on CPU when GPU unavailable

4️⃣ **Dynamic Binding Pocket Refinement**
    • Chain-aware residue selection within a configurable ligand pocket cutoff
    • Writes reduced pocket-only trajectory and topology for downstream analysis

5️⃣ **Comprehensive RMSD/RMSF Analysis**
    • Global backbone RMSD
    • Ligand RMSD and per-atom RMSF
    • Per-protein RMSD/RMSF from atomIndex.txt
    • Bound-complex and full-complex RMSD time series

6️⃣ **Residue Contact & Interaction Analysis**
    • Per-frame ligand–protein contact mapping (with tqdm progress bars)
    • Heatmap of contact persistence over simulation time
    • Contact frequency histogram (residue-wise)
    • Classification of H-bond, hydrophobic, and ionic interactions
    • Automatically colored network graph (percent of frames per residue)

7️⃣ **Robust CLI & Logging**
    • Modular `phase()` / `tick()` / `done()` functions for live status feedback
    • Auto-timestamps, elapsed phase timing, and parallel progress visualization
    • Clean failure handling for missing files, ligand misnaming, or topology issues

📊 Output Summary
----------------
All results are written to the output directory (`--outdir`, default: `Analysis_Results/`):

│── RMSD / RMSF Data
│   ├── rmsd_over_time.csv / .png
│   ├── {LIG}_rmsf.csv / _rmsf_plot.png
│   ├── RMSD_{CHAIN}.csv / .png  (per chain)
│   ├── RMSF_{CHAIN}.csv / .png
│   ├── {LIG}_BoundComplex_RMSD.csv / .png
│   ├── GlobalComplex_RMSD.csv / .png
│
│── Pocket Contact & Interaction Analysis
│   ├── {LIG}_contact_map_residue_time.png
│   ├── residue_contact_frequency.png
│   ├── interaction_stacked_normalized.png
│   └── {LIG}_interaction_network_percent_in_nodes_pastel.png
│
│── Pocket & Subset Trajectories
│   ├── Final_Trajectory.xtc / .pdb (Protein|Ligand subset)
│   ├── binding_pocket_only.xtc / .pdb
│
└── Log Output
    └── Detailed phase-by-phase progress with timestamps and thread summary

⚙️ Example Usage
----------------
Interactive (auto-detect threads):
    $ python analyze_md.py --topo md_0_1.gro --traj md_0_1.xtc -l GDP

Headless (predefined threading):
    $ export OMP_NUM_THREADS=16
    $ python analyze_md.py --threads 16 --topo md_0_1.gro --traj md_0_1.xtc

🧩 Workflow Summary
-------------------
STEP A  — Recenter & Subset (Protein|Ligand)
STEP B  — Dynamic Pocket Refinement (configurable pocket cutoff)
STEP C  — Identify Binding Chains (via atomIndex.txt)
STEP D1 — Global Protein RMSD
STEP D2 — Ligand RMSD / RMSF
STEP D3 — Per-Protein RMSD/RMSF
STEP D4 — Bound Complex RMSD
STEP D5 — Full Complex RMSD
STEP E1 — Pocket Contact Analysis (MDAnalysis)
STEP E2 — Interaction Classification (H-bond / Hydrophobic / Ionic)
STEP E3 — Network Graph Generation

✅ Final Deliverable
-------------------
A complete post-processing suite producing quantitative and visual analyses
of ligand binding, residue dynamics, and interaction persistence — ideal for
supplementary figures, QC validation, and comparative complex evaluation.

=======================================================================
"""

import os
import sys
import shutil
import argparse
import subprocess
import shlex
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import gc, MDAnalysis as mda
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import networkx as nx
import multiprocessing
from tqdm import tqdm
import time
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pymacs_component_utils import (
    build_polymer_group_name,
    get_all_retained_components,
    index_has_group,
    load_system_registry,
    load_component_registry,
    merged_component_group_name,
    parse_component_list,
    polymer_type_label,
)
from pymacs_analysis_cache import (
    compute_signature,
    fingerprint_files,
    record_stage_completion,
    validate_stage_cache,
)

GMX_BIN = None


def resolve_gmx_binary(requested="auto"):
    explicit_requested = requested not in (None, "", "auto")
    env_requested = os.environ.get("PYMACS_GMX_BIN")

    if explicit_requested:
        candidates = [requested]
    elif env_requested:
        candidates = [env_requested]
    else:
        candidates = ["gmx_mpi", "gmx-mpi", "gmx"]

    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved

    if explicit_requested:
        raise FileNotFoundError(f"Requested GROMACS executable was not found: {requested}")
    if env_requested:
        raise FileNotFoundError(f"PYMACS_GMX_BIN is set but not callable: {env_requested}")
    raise FileNotFoundError(
        "No GROMACS executable found. Tried: gmx_mpi, gmx-mpi, gmx. "
        "Use --gmx-bin or set PYMACS_GMX_BIN."
    )


def gmx_cmd(subcommand, args=""):
    base = f"{shlex.quote(GMX_BIN)} {subcommand}"
    return f"{base} {args}".strip()


def gmx_argv(subcommand, args=None):
    return [GMX_BIN, subcommand] + list(args or [])


def print_gmx_preflight():
    print(f"🧬 Using GROMACS executable: {GMX_BIN}")
    try:
        result = subprocess.run([GMX_BIN, "--version"], text=True, capture_output=True)
    except Exception as exc:
        raise RuntimeError(f"Failed to execute '{GMX_BIN} --version': {exc}") from exc

    if result.returncode != 0:
        stderr = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(f"GROMACS executable is not callable: {GMX_BIN}\n{stderr}")

    version_line = next((line.strip() for line in result.stdout.splitlines() if line.strip()), None)
    if version_line:
        print(f"🧪 GROMACS version: {version_line}")




def phase(msg):
    print("\n" + "="*70)
    print(f"🧩 {msg.upper()} — {datetime.now().strftime('%H:%M:%S')}")
    print("="*70)
    time.sleep(0.3)

def tick(msg):
    print(f"   ⏳ {msg}...", end="", flush=True)

def done(msg="Done"):
    print(f" ✅ {msg}")

# ------------------------- CLI -------------------------
parser = argparse.ArgumentParser(
    description="Generate RMSD/RMSF and pocket interaction plots (robust ligand detection)."
)

parser.add_argument(
    "-l", "--ligand", type=str, default=None,
    help="Ligand resname (e.g., PTC, GDP). If omitted, auto-detects largest non-protein residue."
)
parser.add_argument(
    "--cofactors", "--context-ligands", dest="cofactors", type=str, default=None,
    help="Comma-separated retained cofactors/context ligands kept in the trajectory but excluded from ligand-centered analysis by default."
)

parser.add_argument(
    "--lig-atom-start", type=int, default=None,
    help="Fallback: starting atom index (global, 0-based) for ligand in MDTraj topology."
)

parser.add_argument(
    "--lig-atom-end", type=int, default=None,
    help="Fallback: ending atom index (inclusive) for ligand in MDTraj topology."
)

# 👇 Combined alias: --topo OR --gro (both map to args.topo)
parser.add_argument(
    "--topo", "--gro", default="Final_Trajectory.pdb",
    help="Topology for MDTraj (e.g., .gro or .pdb). (--gro is an alias for convenience.)"
)

parser.add_argument(
    "--traj", default="Final_Trajectory.xtc",
    help="Trajectory for MDTraj (e.g., .xtc)."
)

parser.add_argument(
    "--pocket_pdb", default="binding_pocket_only.pdb",
    help="Pocket-only PDB for MDAnalysis."
)

parser.add_argument(
    "--pocket_xtc", default="binding_pocket_only.xtc",
    help="Pocket-only XTC for MDAnalysis."
)

parser.add_argument(
    "--outdir", default="Analysis_Results",
    help="Output directory for results and plots."
)

parser.add_argument(
    "--contact_cutoff", type=float, default=4.0,
    help="Contact cutoff distance in Å (default 4.0)."
)

parser.add_argument(
    "--pocket-cutoff", type=float, default=5.0,
    help="Pocket extraction cutoff in Å used when building binding_pocket_only.* around the ligand."
)

parser.add_argument(
    "--hbond-distance-cutoff", type=float, default=3.5,
    help="Distance cutoff in Å for hydrogen-bond-like interaction classification."
)

parser.add_argument(
    "--hbond-angle-cutoff", type=float, default=135.0,
    help="Reserved donor-hydrogen-acceptor angle cutoff in degrees for documenting hydrogen-bond geometry expectations."
)

parser.add_argument(
    "--min_contact_frac", type=float, default=0.10,
    help="Minimum fraction of frames a residue must contact ligand to appear in plots (default 0.10 = 10%)."
)

parser.add_argument(
    "--threads", type=int, default=None,
    help="Number of CPU threads to use (default: auto-detect 66%% of available cores). "
         "If empty in interactive mode, it will auto-fill."
)

parser.add_argument(
    "--chunk-size",
    type=int,
    default=1000,
    help="Number of trajectory frames to process at once for memory-safe analyses."
)

parser.add_argument(
    "--compound-name",
    type=str,
    default=None,
    help="Readable compound name to use in output (e.g., Paprotrain or CPD32_Inhibitor)"
)

parser.add_argument(
    "--numbering-template",
    type=str,
    default=None,
    help="Optional PDB/CIF/mmCIF file to use as residue-numbering template."
)

parser.add_argument(
    "--mode",
    choices=["ligand", "protein", "protac", "peptide", "biological"],
    default=None,
    help="Analysis mode: ligand, protein, peptide, protac, or biological. If omitted, prompts interactively like Step 2."
)

parser.add_argument(
    "--headless",
    action="store_true",
    help="Non-interactive mode; do not prompt for simulation type or optional ligand-only rendering settings."
)
parser.add_argument(
    "--gmx-bin",
    default="auto",
    help="GROMACS executable to use: auto, gmx, gmx_mpi, gmx-mpi, or an explicit path.",
)
parser.add_argument(
    "--skip-interface-rin",
    action="store_true",
    help="Do not dispatch to 3C_Interface_RIN.py even if a ligandless multi-chain interface is detected.",
)
parser.add_argument(
    "--interface-contact-cutoff",
    type=float,
    default=None,
    help="Distance cutoff in Å passed to 3C_Interface_RIN.py (default: --contact_cutoff).",
)
parser.add_argument(
    "--interface-min-contact-frac",
    type=float,
    default=None,
    help="Minimum residue-residue contact fraction passed to 3C_Interface_RIN.py (default: --min_contact_frac).",
)
parser.add_argument(
    "--interface-frame-step",
    type=int,
    default=1,
    help="Analyze every Nth frame in 3C_Interface_RIN.py (default 1).",
)
parser.add_argument(
    "--interface-max-edges",
    type=int,
    default=100,
    help="Maximum number of residue-residue edges shown in the final interface network (default 100).",
)
parser.add_argument(
    "--interface-chain-pairs",
    type=str,
    default=None,
    help="Optional comma-separated chain pairs for 3C_Interface_RIN.py (example: A:B,A:C).",
)

# ------------------ NETWORX-compatible flags ------------------
parser.add_argument(
    "--minfrac", type=float, default=None,
    help="Minimum fraction of frames a residue must contact ligand (NETWORX threshold, default 0.10 = 10%)."
)

parser.add_argument(
    "--net-out", type=str, default=None,
    help="Optional output base name for NETWORX figure (e.g., A1D_network)."
)

parser.add_argument(
    "--ellipse-rx", type=float, default=None,
    help="Override ellipse horizontal radius (NETWORX)."
)

parser.add_argument(
    "--ellipse-ry", type=float, default=None,
    help="Override ellipse vertical radius (NETWORX)."
)

parser.add_argument(
    "--no-edges", action="store_true",
    help="If set, NETWORX expanded panel will hide all edges (node-only mode)."
)


args = parser.parse_args()
if args.interface_contact_cutoff is None:
    args.interface_contact_cutoff = args.contact_cutoff
if args.interface_min_contact_frac is None:
    args.interface_min_contact_frac = args.min_contact_frac
try:
    GMX_BIN = resolve_gmx_binary(args.gmx_bin)
    print_gmx_preflight()
except (FileNotFoundError, RuntimeError) as exc:
    raise SystemExit(f"❌ {exc}")






import os
import json

def load_peptide_chainmap(run_dir="."):
    """
    Auto-detect peptide metadata from sidecar JSON.
    Returns dict or None.
    """
    candidates = [
        os.path.join(run_dir, ".chainmap.peptide.json"),
        os.path.join(run_dir, "Analysis_Results", ".chainmap.peptide.json"),
    ]

    for path in candidates:
        if os.path.isfile(path):
            with open(path, "r") as f:
                data = json.load(f)

            chain_id = str(data.get("chain_id", "")).strip()
            name = str(data.get("name", "Peptide")).strip()
            group_name = str(data.get("group_name", name.replace(" ", "_"))).strip()
            atom_range = str(data.get("atom_range", "")).strip()

            if not chain_id:
                raise ValueError(f"Peptide chainmap exists but has no chain_id: {path}")

            return {
                "path": path,
                "chain_id": chain_id,
                "name": name,
                "group_name": group_name,
                "atom_range": atom_range,
            }


    return None

def parse_atomindex_file(path="atomIndex.txt"):
    """
    Parse atomIndex.txt robustly.

    Accepts lines like:
      A 1-1975
      chainA 1-1975
      a 1-1975
      a 1-1975 # comment

    Returns dict[name] = (start0, end0)
    """
    chain_map = {}
    if not os.path.exists(path):
        return chain_map

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            line = line.split("#", 1)[0].strip()
            if not line:
                continue

            m = re.match(r"^(\S+)\s+(\d+)\s*-\s*(\d+)$", line)
            if not m:
                toks = line.split()
                rng_tok = None
                name_tok = None
                for tok in toks:
                    if re.fullmatch(r"\d+\s*-\s*\d+", tok):
                        rng_tok = tok.replace(" ", "")
                        break
                if rng_tok:
                    for tok in toks:
                        if tok != rng_tok:
                            name_tok = tok
                            break
                if not (name_tok and rng_tok):
                    print(f"⚠️ Skipping unparsable atomIndex line: {raw.rstrip()}")
                    continue
                name = name_tok
                start_s, end_s = rng_tok.split("-", 1)
            else:
                name = m.group(1)
                start_s, end_s = m.group(2), m.group(3)

            start = int(start_s) - 1
            end = int(end_s) - 1
            chain_map[name] = (start, end)

    return chain_map


def parse_atomindex_entries(path="atomIndex.txt"):
    """
    Parse atomIndex.txt and preserve inline provenance comments like:
      Label 1-100 # current_chain=A chain_type=rna ...
    """
    entries = []
    if not os.path.exists(path):
        return entries

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            comment = ""
            if "#" in line:
                line, comment = line.split("#", 1)
                line = line.strip()
                comment = comment.strip()

            if not line:
                continue

            m = re.match(r"^(\S+)\s+(\d+)\s*-\s*(\d+)$", line)
            if not m:
                print(f"⚠️ Skipping unparsable atomIndex line: {raw.rstrip()}")
                continue

            label = m.group(1)
            start = int(m.group(2)) - 1
            end = int(m.group(3)) - 1

            meta = {}
            for key, value in re.findall(r"(\w+)=([^\s#]+)", comment):
                meta[key] = value

            entries.append(
                {
                    "label": label,
                    "start": start,
                    "end": end,
                    "current_chain": meta.get("current_chain"),
                    "chain_type": (meta.get("chain_type") or "mixed_or_unknown").strip().lower(),
                    "current_residues": meta.get("current_residues"),
                    "source_chain": meta.get("source_chain"),
                    "source_residues": meta.get("source_residues"),
                }
            )
    return entries


def normalize_chain_type(chain_type):
    return (chain_type or "mixed_or_unknown").strip().lower()


def is_polymer_chain_type(chain_type):
    normalized = normalize_chain_type(chain_type)
    if not normalized or normalized == "mixed_or_unknown":
        return True
    return normalized not in {
        "ligand",
        "small_molecule",
        "cofactor",
        "ion",
        "water",
        "solvent",
        "lipid",
    }


def chain_display_label(entry, index):
    return (
        entry.get("current_chain")
        or entry.get("source_chain")
        or entry.get("display_label")
        or entry.get("label")
        or f"Chain_{index}"
    )


def collect_polymer_chain_entries(chain_map, chain_entry_meta, system_registry=None):
    polymer_entries = []
    registry_chain_types = (system_registry or {}).get("chain_types", {})

    for idx, (label, (start, end)) in enumerate(chain_map.items(), start=1):
        meta = chain_entry_meta.get(label, {})
        current_chain = meta.get("current_chain")
        source_chain = meta.get("source_chain")
        chain_type = normalize_chain_type(meta.get("chain_type"))
        if (not chain_type or chain_type == "mixed_or_unknown") and current_chain:
            chain_type = normalize_chain_type(registry_chain_types.get(current_chain))
        if not chain_type:
            chain_type = "mixed_or_unknown"
        if not is_polymer_chain_type(chain_type):
            continue
        entry = {
            "label": label,
            "start": start,
            "end": end,
            "current_chain": current_chain,
            "source_chain": source_chain,
            "chain_type": chain_type,
        }
        if meta.get("display_label"):
            entry["display_label"] = meta.get("display_label")
        entry["display_label"] = chain_display_label(entry, idx)
        polymer_entries.append(entry)

    return polymer_entries


def anchor_atom_candidates_for_chain(chain_type):
    normalized = (chain_type or "").strip().lower()
    if normalized == "protein":
        return ["CA"]
    if normalized in {"rna", "dna"}:
        return ["P", "C4'", "C4*", "C3'", "C3*", "O5'", "O5*"]
    return ["CA", "P", "C4'", "C4*", "C3'", "C3*", "O5'", "O5*"]


def select_chain_anchor_atoms(traj, start, end, chain_type="mixed_or_unknown"):
    """
    Pick one representative anchor atom per residue:
    - protein: CA
    - RNA/DNA: prefer P, then sugar anchors
    Returns (anchor_indices, residue_numbers, anchor_name)
    """
    residue_map = []
    for residue in traj.topology.residues:
        atom_indices = [atom.index for atom in residue.atoms if start <= atom.index <= end]
        if atom_indices:
            residue_map.append((residue, atom_indices))

    candidates = anchor_atom_candidates_for_chain(chain_type)
    selected_indices = []
    residue_numbers = []
    used_anchor_name = None

    for residue, atom_indices in residue_map:
        residue_atoms = {traj.topology.atom(idx).name: idx for idx in atom_indices}
        chosen = None
        chosen_name = None
        for atom_name in candidates:
            if atom_name in residue_atoms:
                chosen = residue_atoms[atom_name]
                chosen_name = atom_name
                break
        if chosen is None:
            continue
        if used_anchor_name is None:
            used_anchor_name = chosen_name
        selected_indices.append(chosen)
        residue_numbers.append(residue.resSeq)

    return np.array(selected_indices, dtype=int), residue_numbers, used_anchor_name or "mixed"

def peptide_atom_bounds(peptide_info):
    if not peptide_info:
        return None
    atom_range = str(peptide_info.get("atom_range", "")).strip()
    if "-" not in atom_range:
        return None
    try:
        s, e = atom_range.split("-", 1)
        return int(s) - 1, int(e) - 1
    except Exception:
        return None

def mda_partner_selection(peptide_info=None, ligand_code=None):
    """
    Return MDAnalysis selection string for ligand/peptide partner.
    """
    if peptide_info:
        bounds = peptide_atom_bounds(peptide_info)
        if bounds is not None:
            s, e = bounds
            return f"index {s}:{e}"
        chain_id = str(peptide_info.get("chain_id", "")).strip()
        if chain_id:
            return f"segid {chain_id} or chainID {chain_id}"
        group_name = str(peptide_info.get("group_name", "")).strip()
        if group_name:
            return f"resname {group_name}"
        raise ValueError("Peptide mode selected but no usable peptide selector could be built.")
    if ligand_code:
        return f"resname {ligand_code}"
    return None

def mdtraj_partner_indices(traj, peptide_info=None, ligand_code=None):
    """
    Return MDTraj atom indices for ligand/peptide partner.
    """
    if peptide_info:
        bounds = peptide_atom_bounds(peptide_info)
        if bounds is not None:
            s, e = bounds
            s = max(0, s)
            e = min(e, traj.n_atoms - 1)
            if s <= e:
                return np.arange(s, e + 1, dtype=int)
        chain_id = str(peptide_info.get("chain_id", "")).strip()
        if chain_id:
            atoms = []
            for atom in traj.topology.atoms:
                cid = getattr(atom.residue.chain, "chain_id", None)
                if cid is None:
                    cid = getattr(atom.residue.chain, "id", None)
                if cid is not None and str(cid).strip() == chain_id:
                    atoms.append(atom.index)
            if atoms:
                return np.array(sorted(atoms), dtype=int)
        group_name = str(peptide_info.get("group_name", "")).strip()
        if group_name:
            atoms = traj.topology.select("resname == " + repr(group_name))
            if len(atoms) > 0:
                return atoms
        name = str(peptide_info.get("name", "")).strip()
        if name:
            # Last-ditch fallback: search by residue name tokens
            token = name.replace(" ", "_")
            atoms = traj.topology.select("resname == " + repr(token))
            if len(atoms) > 0:
                return atoms
        raise RuntimeError("❌ Peptide partner atoms could not be resolved from .chainmap.peptide.json.")
    if ligand_code:
        return traj.topology.select("resname == " + repr(ligand_code.upper()))
    return np.array([], dtype=int)

def label_from_mode(ligand_code=None, peptide_info=None):
    if peptide_info:
        return peptide_info.get("group_name") or peptide_info.get("name") or "Peptide"
    return ligand_code

def select_simulation_type(args):
    sim_map = {
        "1": "Ligand:Protein",
        "2": "Protein:Protein/Peptide",
        "3": "PROTAC",
        "4": "Just Protein",
        "5": "Biological System",
    }
    mode_map = {
        "ligand": "Ligand:Protein",
        "peptide": "Protein:Protein/Peptide",
        "protac": "PROTAC",
        "protein": "Just Protein",
        "biological": "Biological System",
    }

    if args.mode:
        return mode_map.get(args.mode.lower(), args.mode)

    if args.headless:
        return "Just Protein"

    print("\nPlease select the simulation type:")
    print("  1. Ligand:Protein")
    print("  2. Protein:Protein/Peptide")
    print("  3. PROTAC")
    print("  4. Just Protein")
    print("  5. Biological System (Protein:RNA/DNA or RNA/DNA)")

    sim_choice = input("Enter choice (1-5): ").strip()
    if sim_choice not in sim_map:
        print("❌ Invalid selection. Defaulting to Just Protein.")
        return "Just Protein"
    return sim_map[sim_choice]


def dispatch_to_protac_analysis():
    protac_script_mpi = os.path.join(os.getcwd(), "3_PROTAC_Analysis_MPI.py")
    protac_script_legacy = os.path.join(os.getcwd(), "3_PROTAC_Analysis.py")

    if os.path.exists(protac_script_mpi):
        protac_script = protac_script_mpi
    elif os.path.exists(protac_script_legacy):
        protac_script = protac_script_legacy
        print("⚠️ 3_PROTAC_Analysis_MPI.py not found; falling back to 3_PROTAC_Analysis.py.", flush=True)
    else:
        raise SystemExit(
            "❌ PROTAC mode selected, but neither 3_PROTAC_Analysis_MPI.py nor 3_PROTAC_Analysis.py was found in the current directory."
        )

    cmd = [sys.executable, protac_script, "--gmx-bin", GMX_BIN]

    if args.headless:
        cmd.extend(["--headless", "--frame-step", "10"])
        print(f"\n🧬 PROTAC mode selected — headless handoff to {os.path.basename(protac_script)} (standard sampling: --frame-step 10)", flush=True)
    else:
        print("\n🧬 PROTAC mode selected — choose analysis depth:", flush=True)
        print("  1. Quick test     (0:500:10)", flush=True)
        print("  2. Standard       (default PROTAC sampling, --frame-step 10)", flush=True)
        print("  3. Full           (--frame-step 1)", flush=True)
        protac_choice = input("Enter choice (1-3) [2]: ").strip() or "2"

        if protac_choice == "1":
            cmd.append("--quick-test")
            print("🧪 Using PROTAC quick-test handoff.", flush=True)
        elif protac_choice == "3":
            cmd.extend(["--frame-step", "1"])
            print("🧬 Using full PROTAC frame-by-frame handoff.", flush=True)
        else:
            cmd.extend(["--frame-step", "10"])
            print("🧬 Using standard PROTAC handoff with frame-step 10.", flush=True)

    print(f"🚀 Launching: {' '.join(cmd)}", flush=True)
    rc = subprocess.call(cmd)
    raise SystemExit(rc)


SIMULATION_TYPE = select_simulation_type(args)
SIM_TYPE_NORM = SIMULATION_TYPE.lower().strip()
PROTEIN_ONLY = SIM_TYPE_NORM in {"just protein", "protein", "biological system", "biological"}
PEPTIDE_LIKE = SIM_TYPE_NORM in {"protein:protein/peptide", "peptide"}
LIGAND_MODE = SIM_TYPE_NORM in {"ligand:protein", "ligand"}
PROTAC_MODE = SIM_TYPE_NORM == "protac"
FORCED_LIGAND_MODE = bool((args.ligand or "").strip()) and not PROTAC_MODE
PARTNER_MODE = LIGAND_MODE or FORCED_LIGAND_MODE
PEPTIDE_INFO = load_peptide_chainmap(".") if PEPTIDE_LIKE else None
COMPONENT_REGISTRY = load_component_registry(".") or {}
SYSTEM_REGISTRY = load_system_registry(".") or {}
REGISTRY_PRIMARY_LIGAND = (COMPONENT_REGISTRY.get("primary_ligand") or None)
REGISTRY_COFACTORS = COMPONENT_REGISTRY.get("cofactors", [])
CONTEXT_COFACTORS = parse_component_list(args.cofactors) if args.cofactors else REGISTRY_COFACTORS
CONTEXT_COFACTORS = [code for code in CONTEXT_COFACTORS if code != REGISTRY_PRIMARY_LIGAND]
BIOMOLECULE_MODE = bool(
    PROTEIN_ONLY and (
        SYSTEM_REGISTRY.get("has_rna")
        or SYSTEM_REGISTRY.get("has_dna")
    )
)
ANALYSIS_POLYMER_GROUP = build_polymer_group_name(SYSTEM_REGISTRY, [], override="auto")

print(f"🧠 Selected analysis type: {SIMULATION_TYPE}")
if COMPONENT_REGISTRY:
    print(f"🧾 Loaded component registry: primary={REGISTRY_PRIMARY_LIGAND} cofactors={REGISTRY_COFACTORS}")
if SYSTEM_REGISTRY:
    print(
        f"🧾 Loaded system registry: class={SYSTEM_REGISTRY.get('system_class')} "
        f"polymer_chains={SYSTEM_REGISTRY.get('chain_types', {})}"
    )
if BIOMOLECULE_MODE:
    print(f"🧬 Biomolecule-aware protein mode enabled via registry group: {ANALYSIS_POLYMER_GROUP}")
    print("🧬 This path currently produces chain-vs-chain biomolecule contacts and distance summaries.")
    print("🧬 Ligand-style per-residue interaction networks are still reserved for ligand/peptide partner analyses.")
if PEPTIDE_LIKE:
    if PEPTIDE_INFO is not None:
        print("🧬 Protein:Protein/Peptide mode detected.")
        print(f"   • chain_id   = {PEPTIDE_INFO.get('chain_id')}")
        print(f"   • name       = {PEPTIDE_INFO.get('name')}")
        print(f"   • group_name = {PEPTIDE_INFO.get('group_name')}")
        print(f"   • atom_range = {PEPTIDE_INFO.get('atom_range')}")
    else:
        print("⚠️ Protein:Protein/Peptide mode selected but .chainmap.peptide.json was not found.")
    if FORCED_LIGAND_MODE:
        print("🧬 Explicit ligand supplied; preserving ligand-centered analysis behavior.")
    else:
        print("🧬 Ligandless multi-chain peptide/protein systems will be evaluated for standalone interface RIN dispatch.")

if PROTAC_MODE:
    dispatch_to_protac_analysis()

NUMBERING_TEMPLATE = args.numbering_template
MIN_CONTACT_FRAC = args.min_contact_frac



def list_numbering_template_files():
    return sorted([
        f for f in os.listdir(".")
        if f.lower().endswith((".pdb", ".cif", ".mmcif")) and os.path.isfile(f)
    ])


def select_numbering_template_interactive():
    template_files = list_numbering_template_files()

    if not template_files:
        print("⚠️ No PDB/CIF/mmCIF files found in current directory.")
        return None

    print("\n🧬 Residue numbering template selection")
    print("Do you want to use a structure file as a residue-numbering template?")
    resp = ask("Use structure template for residue numbering? [y/N]: ", str, default="n")

    if not resp.lower().startswith("y"):
        return None

    print("\n📂 Available template files:")
    for i, f in enumerate(template_files, 1):
        print(f"   {i}) {f}")

    choice = ask(
        "Select a template file by number (or press ENTER to cancel): ",
        str,
        default=""
    )

    if not choice:
        print("⏭️ No template selected — skipping residue remapping.")
        return None

    if choice.isdigit() and 1 <= int(choice) <= len(template_files):
        selected = template_files[int(choice) - 1]
        print(f"✅ Using {selected} as residue numbering template.")
        return selected

    print("⚠️ Invalid selection — skipping residue remapping.")
    return None





# ------------------------------------------------------------
# INTERACTIVE PROMPTS FOR NETWORX ARGUMENTS (if missing)
# ------------------------------------------------------------
import string

def ask(prompt, cast=str, default=None):
    """
    Ask user for input interactively with validation.
    - Strips non-printable characters
    - Retries on invalid input
    - Falls back to default on ENTER or EOF
    """
    while True:
        try:
            raw = input(prompt)

            # Remove non-printable / control characters
            v = "".join(ch for ch in raw if ch in string.printable).strip()

            if v == "":
                return default

            try:
                return cast(v)
            except ValueError:
                print(
                    f"❌ Invalid input: '{raw}'.\n"
                    f"   Please enter a valid {cast.__name__} "
                    f"(digits only, e.g. 0.2)"
                )

        except EOFError:
            return default

# ------------------------------------------------------------
# ALWAYS FORCE NETWORX OUTPUT INTO Analysis_Results/NETWORX
# ------------------------------------------------------------
lig_code = args.ligand or os.path.basename(args.topo).split('.')[0]

NET_OUTDIR = os.path.join("Analysis_Results", "NETWORX")
os.makedirs(NET_OUTDIR, exist_ok=True)

# The base output stem used for all NETWORX panel images (A, B, C)
args.net_out = os.path.join(NET_OUTDIR, lig_code)

if PARTNER_MODE:
    # If no --minfrac provided explicitly → ask
    if args.minfrac is None:
        args.minfrac = ask(
            "📊 Minimum fraction of frames for residue inclusion [0.0–1.0] (default 0.10): ",
            float,
            default=0.10
        )

    # If ellipse radii were not provided → ask interactively
    if args.ellipse_rx is None:
        args.ellipse_rx = ask(
            "⚪ NETWORX ellipse horizontal radius (rx) (press ENTER to use auto layout): ",
            float,
            default=None
        )

    if args.ellipse_ry is None:
        args.ellipse_ry = ask(
            "⚪ NETWORX ellipse vertical radius (ry) (press ENTER to use auto layout): ",
            float,
            default=None
        )

    # If --no-edges was not passed → ask yes/no
    if not args.no_edges:
        v = ask("Hide edges in NETWORX graph? (y/n, default n): ", str, default="n")
        if v.lower().startswith("y"):
            args.no_edges = True
        else:
            args.no_edges = False
else:
    args.minfrac = 0.10 if args.minfrac is None else args.minfrac







PLOT_TYPES = ["H-bond", "Hydrophobic", "Ionic"]
# ------------------------------------------------------------
# Define interaction types and colors for plotting
# ------------------------------------------------------------



# Default color map (safe + readable)
DEFAULT_COLORS = {
    "H-bond": "#1f77b4",        # blue
    "Hydrophobic": "#ff7f0e",   # orange
    "Ionic": "#2ca02c",         # green
}

# Assign colors, fallback to gray if unknown type
COLORS = {
    t: DEFAULT_COLORS.get(t, "#7f7f7f")
    for t in PLOT_TYPES
}

print("🎨 Plot interaction types:", PLOT_TYPES)
print("🎨 Assigned colors:", COLORS)


TOPO_FILE   = args.topo
TRAJ_FILE   = args.traj
POCKET_PDB  = args.pocket_pdb
POCKET_XTC  = args.pocket_xtc
OUTPUT_DIR  = args.outdir
CONTACT_CUTOFF = float(args.contact_cutoff)
POCKET_CUTOFF_ANG = float(args.pocket_cutoff)
POCKET_CUTOFF_NM = POCKET_CUTOFF_ANG / 10.0
HBOND_DISTANCE_CUTOFF = float(args.hbond_distance_cutoff)
HBOND_ANGLE_CUTOFF = float(args.hbond_angle_cutoff)
MIN_CONTACT_FRAC = float(args.min_contact_frac)
MAX_WORKERS = max(1, os.cpu_count() // 2)
print(f"📏 Pocket cutoff: {POCKET_CUTOFF_ANG:.2f} Å")
print(f"📏 Contact cutoff: {CONTACT_CUTOFF:.2f} Å")
print(f"📏 H-bond-like distance cutoff: {HBOND_DISTANCE_CUTOFF:.2f} Å")
print(f"📐 H-bond angle cutoff (documented): {HBOND_ANGLE_CUTOFF:.1f}°")
# ============================================================
# MIAMI COLOR PALETTE (LOCKED)
# ============================================================
MIAMI_ORANGE = "#F47321"   # Protein
MIAMI_GREEN  = "#005030"   # Ligand
MIAMI_BLUE   = "#0033A0"   # Water

# ============================================================


os.makedirs(OUTPUT_DIR, exist_ok=True)
# Create subfolder for distance/contact histograms
DIST_DIR = os.path.join(OUTPUT_DIR, "RESIDUE_Contact_Distance")
os.makedirs(DIST_DIR, exist_ok=True)

PIE_DIR = os.path.join(OUTPUT_DIR, "RESIDUE_Pies")
os.makedirs(PIE_DIR, exist_ok=True)

# ============================================================
# 🧠 OpenMP Thread Management
# ============================================================
if args.threads is None:
    try:
        # Interactive mode: show available cores
        total_cores = multiprocessing.cpu_count()
        print(f"\n🧩 Detected {total_cores} total CPU threads available.")
        user_input = input("⚙️ Enter number of threads to use [press ENTER for auto (66%)]: ").strip()
        if user_input:
            threads = int(user_input)
        else:
            threads = max(1, round(total_cores * 0.66))
    except EOFError:
        # Non-interactive / headless run (e.g., bash pipeline)
        threads = max(1, round(multiprocessing.cpu_count() * 0.66))
else:
    threads = args.threads

# Use either environment or user-provided override
threads_env = os.environ.get("OMP_NUM_THREADS")
if threads_env and threads_env.isdigit():
    print(f"🧠 OMP_NUM_THREADS already set in environment = {threads_env}")
    threads = int(threads_env)
else:
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    os.environ["OMP_PROC_BIND"] = "spread"
    os.environ["OMP_PLACES"] = "cores"

print(f"⚙️ Using {threads} CPU threads for RMSD/RMSF calculations.")

# --- Optional GPU acceleration (CuPy) ---
try:
    import cupy as cp
    GPU_AVAILABLE = True
    print("⚡ CuPy GPU acceleration enabled for distance calculations.")
except ImportError:
    import numpy as cp  # fallback to NumPy API
    GPU_AVAILABLE = False
    print("⚠️ CuPy not found — using CPU (NumPy) for distance calculations.")


# ============================================================
# 🧩 STARTUP: Ligand Handling, Recentered Subset, and Pocket Extraction
# ============================================================
import subprocess, os
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np

# ------------------------- Helper functions -------------------------
def file_exists(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0

def run_cmd(cmd, cwd="."):
    """Run a shell command and stream output live."""
    process = subprocess.Popen(cmd, shell=True, cwd=cwd,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, universal_newlines=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        print(f"⚠️ Command failed (continuing): {cmd}")

# ============================================================
# 🔍 SMART LIGAND DETECTION BLOCK (Auto + User Confirmation)
# ============================================================

def detect_ligand_candidates():
    """Scan directory and topology to find possible ligand candidates."""
    candidates = set()

    # 1) Any file like X_ini.pdb → extract residue code X
    for f in os.listdir("."):
        if f.endswith("_ini.pdb") and len(f.split("_")[0]) <= 4:
            candidates.add(f.split("_")[0].upper())

    # 2) Scan topology for residues NOT standard amino acids or water
    non_standard_res = set()
    try:
        topo = mda.Universe(TOPO_FILE)
        for res in topo.residues:
            if res.resname not in {
                "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS",
                "ILE","LEU","LYS","MET","PHE","PRO","SER","THR",
                "TRP","TYR","VAL",
                "HOH","WAT","SOL"
            }:
                non_standard_res.add(res.resname.strip().upper())
    except Exception:
        pass

    candidates.update(non_standard_res)
    return sorted(list(candidates))


# ------------------------- Ligand Input -------------------------

ligand_code = None
compound_name = None
PARTNER_SELECTION_MDA = None
PARTNER_IS_PEPTIDE = False

if PARTNER_MODE:
    if args.ligand or REGISTRY_PRIMARY_LIGAND:
        ligand_code = (args.ligand or REGISTRY_PRIMARY_LIGAND).strip().upper()
        source = "CLI" if args.ligand else "component registry"
        print(f"🧬 Using ligand from {source}: {ligand_code}")

    else:
        candidates = detect_ligand_candidates()

        if len(candidates) == 1:
            auto = candidates[0]
            print(f"\n🔍 Detected ligand candidate: {auto}")
            resp = input(f"Use {auto} as ligand? [Y/n]: ").strip().lower()
            if resp in ("", "y", "yes"):
                ligand_code = auto
            else:
                ligand_code = input("Enter ligand resname manually: ").strip().upper()

        elif len(candidates) > 1:
            print("\n🔍 Multiple ligand candidates detected:")
            for i, c in enumerate(candidates, 1):
                print(f"   {i}) {c}")

            choice = input("Select a ligand or press ENTER to type manually: ").strip()
            if choice.isdigit() and 1 <= int(choice) <= len(candidates):
                ligand_code = candidates[int(choice)-1]
            else:
                ligand_code = input("Enter ligand resname manually: ").strip().upper()

        else:
            ligand_code = input("❓ Enter ligand resname (e.g., PTC, GDP): ").strip().upper()

        if not ligand_code:
            print("❌ ERROR: No ligand provided and auto-detection failed.")
            exit(1)

    CONTEXT_COMPONENTS = get_all_retained_components(ligand_code, CONTEXT_COFACTORS)
    PARTNER_SELECTION_MDA = mda_partner_selection(ligand_code=ligand_code)
    print(f"🧬 Final ligand selected: {ligand_code}")
    print(f"🧬 Context cofactors retained in trajectory: {', '.join(CONTEXT_COFACTORS) if CONTEXT_COFACTORS else 'None'}")

    print("\n📌 You may now choose a full compound name for labeling graphs/output.")
    print("   • Press ENTER to use the ligand code above (recommended only for short names).")
    print("   • Or type a more descriptive compound name (e.g., Paprotrain, CPD32_Inhibitor).")

    compound_input = input("🏷️  Enter compound display name (or press ENTER to use resname): ").strip()

    if compound_input:
        compound_name = compound_input
        print(f"🏷️ Using custom compound name for all plots + summaries: {compound_name}")
    else:
        compound_name = ligand_code
        print(f"🏷️ Using ligand resname for all plots + summaries: {compound_name}")

    print(f"\n🧬 RESNAME for topology selections: {ligand_code}")
    print(f"🏷️ Compound name for figures/output: {compound_name}\n")

else:
    ligand_code = None
    compound_name = args.compound_name or ("ProteinPeptideInterface" if PEPTIDE_LIKE else "ProteinOnly")
    PARTNER_SELECTION_MDA = None
    CONTEXT_COMPONENTS = []
    print("🧬 No ligand will be used for this analysis mode.")
    print(f"🏷️ Compound/output label for protein-centric analysis: {compound_name}")
# ============================================================
# 🧬 OPTIONAL RESIDUE NUMBERING TEMPLATE SELECTION
# ============================================================

if NUMBERING_TEMPLATE:
    if not os.path.exists(NUMBERING_TEMPLATE):
        print(f"⚠️ Provided numbering template not found: {NUMBERING_TEMPLATE}")
        NUMBERING_TEMPLATE = None
    else:
        print(f"🧬 Using residue numbering template from CLI: {NUMBERING_TEMPLATE}")
else:
    NUMBERING_TEMPLATE = select_numbering_template_interactive()


# ------------------------- Ensure index group -------------------------
def ensure_analysis_group(protein_only=False, ligand_code=None, context_components=None, polymer_group_override=None):
    index_file = "index.ndx"
    merged_group = merged_component_group_name(context_components or [ligand_code] if ligand_code else [], fallback_ligand=ligand_code, protein_only=protein_only)
    wanted_group = polymer_group_override or ("Protein" if protein_only else merged_group)

    if not os.path.exists(index_file):
        raise SystemExit("❌ index.ndx not found. Run MD setup script first.")

    with open(index_file) as f:
        groups = [
            line.strip().strip("[]").strip()
            for line in f
            if line.strip().startswith("[")
        ]

    if wanted_group in groups:
        print(f"✅ Found existing group [{wanted_group}] in index.ndx.")
        return wanted_group

    if protein_only and not polymer_group_override:
        if "Protein" in groups:
            print("✅ Using existing [Protein] group for protein-only mode.")
            return "Protein"
        raise SystemExit("❌ [Protein] group not found in index.ndx.")

    fallback_group = f"Protein_{ligand_code}" if ligand_code else "Protein"
    if fallback_group in groups:
        print(f"✅ Using existing fallback group [{fallback_group}] in index.ndx.")
        return fallback_group

    print(f"⚠️ Group [{wanted_group}] not found — creating it now...")
    proc = subprocess.run(
        gmx_cmd("make_ndx", "-f md_0_1.tpr -o temp.ndx"),
        shell=True,
        capture_output=True,
        text=True,
        input="q\n",
    )
    lig_num = None
    for line in proc.stdout.splitlines():
        if ligand_code and ligand_code in line:
            lig_num = line.strip().split()[0]
            break
    lig_num = lig_num or "13"
    run_cmd(
        f"echo '1 | {lig_num}\nname 21 {wanted_group}\nq' | "
        f"{gmx_cmd('make_ndx', '-f md_0_1.tpr -o index.ndx')}"
    )
    print(f"✅ Created [{wanted_group}] (1 | {lig_num}) in index.ndx.")
    return wanted_group

analysis_group = ensure_analysis_group(
    protein_only=(PROTEIN_ONLY or (PEPTIDE_LIKE and not PARTNER_MODE)),
    ligand_code=ligand_code,
    context_components=CONTEXT_COMPONENTS,
    polymer_group_override=ANALYSIS_POLYMER_GROUP if BIOMOLECULE_MODE else None,
)





import re

def residue_number(res):
    """Extract numeric residue index from labels like ARG329"""
    match = re.search(r"(\d+)", res)
    return int(match.group(1)) if match else -1



from collections import defaultdict

def _parse_resseq_token(token):
    match = re.search(r"-?\d+", str(token))
    return int(match.group(0)) if match else None


def _build_residue_number_map_from_pdb(template_pdb):
    resseqs = []
    last_seen = None

    with open(template_pdb) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue

            if resseq != last_seen:
                resseqs.append(resseq)
                last_seen = resseq

    return resseqs


def _build_residue_number_map_from_cif(template_cif):
    resseqs = []
    in_atom_loop = False
    collecting_headers = False
    atom_headers = []
    last_seen = None

    with open(template_cif) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                if in_atom_loop and atom_headers:
                    break
                continue

            if line == "loop_":
                in_atom_loop = False
                collecting_headers = True
                atom_headers = []
                continue

            if collecting_headers and line.startswith("_"):
                atom_headers.append(line)
                continue

            if collecting_headers:
                collecting_headers = False
                if atom_headers and any(header.startswith("_atom_site.") for header in atom_headers):
                    in_atom_loop = True
                    group_idx = atom_headers.index("_atom_site.group_PDB") if "_atom_site.group_PDB" in atom_headers else None
                    auth_seq_idx = atom_headers.index("_atom_site.auth_seq_id") if "_atom_site.auth_seq_id" in atom_headers else None
                    label_seq_idx = atom_headers.index("_atom_site.label_seq_id") if "_atom_site.label_seq_id" in atom_headers else None
                    auth_chain_idx = atom_headers.index("_atom_site.auth_asym_id") if "_atom_site.auth_asym_id" in atom_headers else None
                    label_chain_idx = atom_headers.index("_atom_site.label_asym_id") if "_atom_site.label_asym_id" in atom_headers else None
                else:
                    in_atom_loop = False

            if not in_atom_loop:
                continue

            if line.startswith("_") or line == "loop_":
                break

            parts = shlex.split(line)
            if len(parts) < len(atom_headers):
                continue

            if group_idx is not None and parts[group_idx].upper() != "ATOM":
                continue

            seq_token = None
            if auth_seq_idx is not None:
                seq_token = parts[auth_seq_idx]
            if (seq_token in {None, ".", "?"}) and label_seq_idx is not None:
                seq_token = parts[label_seq_idx]

            resseq = _parse_resseq_token(seq_token)
            if resseq is None:
                continue

            chain_token = ""
            if auth_chain_idx is not None:
                chain_token = parts[auth_chain_idx]
            if chain_token in {"", ".", "?"} and label_chain_idx is not None:
                chain_token = parts[label_chain_idx]

            current_key = (chain_token, resseq)
            if current_key != last_seen:
                resseqs.append(resseq)
                last_seen = current_key

    return resseqs


def build_residue_number_map(template_path):
    """
    Build an ordered list of polymer residue numbers from a PDB, CIF, or mmCIF template.
    """
    lower = template_path.lower()
    if lower.endswith(".pdb"):
        return _build_residue_number_map_from_pdb(template_path)
    if lower.endswith((".cif", ".mmcif")):
        return _build_residue_number_map_from_cif(template_path)
    raise RuntimeError(f"Unsupported residue-numbering template format: {template_path}")



def apply_template_numbering(pdb_in, pdb_out, resseqs):
    """
    Apply residue numbers sequentially to pdb_in
    using the order defined in resseqs.
    """
    res_idx = -1
    last_resseq = None

    with open(pdb_in) as fin, open(pdb_out, "w") as fout:
        for line in fin:
            # Leave non-coordinate lines intact
            if not line.startswith(("ATOM", "HETATM")):
                fout.write(line)
                continue

            # Do NOT consume template residue numbers for HETATM records.
            # HETATM (ligand/ion) entries should be preserved as-is so they
            # don't shift the protein residue mapping.
            if line.startswith("HETATM"):
                fout.write(line)
                continue

            # Only ATOM records are remapped (protein atoms)
            try:
                cur_resseq = int(line[22:26])
            except ValueError:
                fout.write(line)
                continue

            if cur_resseq != last_resseq:
                res_idx += 1
                if res_idx >= len(resseqs):
                    raise RuntimeError(
                        f"❌ Template has fewer residues ({len(resseqs)}) "
                        f"than trajectory ({res_idx+1})"
                    )
                new_resseq = resseqs[res_idx]
                last_resseq = cur_resseq

            fout.write(
                line[:22] + f"{new_resseq:4d}" + line[26:]
            )





# ============================================================
# STEP A+B — PREPROCESS FINAL TRAJECTORY + POCKET EXTRACTION
# (with correct ATOM → HETATM ligand rewriting)
# ============================================================

phase("STEP A+B — Generate Final_Trajectory.* and binding_pocket_only.*")

lig = ligand_code.upper() if ligand_code else None
group = analysis_group

required = ["index.ndx", "md_0_1.tpr", "md_0_1.xtc"]
for f in required:
    if not os.path.exists(f):
        raise SystemExit(f"❌ Missing required file: {f}")


# ------------------------------------------------------------
# Helper — rewrite ligand to HETATM
# ------------------------------------------------------------
def convert_ligand_to_hetatm(pdb_in, pdb_out, ligand_resname):
    ligand_resname = ligand_resname.upper()
    with open(pdb_in, "r") as fin, open(pdb_out, "w") as fout:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resn = line[17:20].strip().upper()
                if resn == ligand_resname:
                    fout.write("HETATM" + line[6:])
                else:
                    fout.write(line)
            else:
                fout.write(line)
    print(f"🔧 Rewrote ligand {ligand_resname} as HETATM → {pdb_out}")


# ============================================================
# STEP A — Final_Trajectory.xtc / Final_Trajectory.pdb
# ============================================================
if not (file_exists("Final_Trajectory.xtc") and file_exists("Final_Trajectory.pdb")):

    print("🎞️ Creating Final_Trajectory.xtc/pdb ...")

    run_cmd(
        f"printf '{group}\n{group}\n' | "
        f"{gmx_cmd('trjconv', '-s md_0_1.tpr -f md_0_1.xtc')} "
        f"-o Final_Trajectory.xtc -n index.ndx -pbc mol -center"
    )

    run_cmd(
        f"printf '{group}\n' | "
        f"{gmx_cmd('trjconv', '-s md_0_1.tpr -f Final_Trajectory.xtc')} "
        f"-o Final_Trajectory_RAW.pdb -n index.ndx -dump 0"
    )

    if PROTEIN_ONLY or (PEPTIDE_LIKE and not PARTNER_MODE):
        if NUMBERING_TEMPLATE and os.path.exists(NUMBERING_TEMPLATE):
            resmap = build_residue_number_map(NUMBERING_TEMPLATE)
            apply_template_numbering(
                "Final_Trajectory_RAW.pdb",
                "Final_Trajectory.pdb",
                resmap
            )
        else:
            os.rename("Final_Trajectory_RAW.pdb", "Final_Trajectory.pdb")
    else:
        convert_ligand_to_hetatm(
            "Final_Trajectory_RAW.pdb",
            "Final_Trajectory_TMP.pdb",
            lig
        )

        if NUMBERING_TEMPLATE and os.path.exists(NUMBERING_TEMPLATE):
            resmap = build_residue_number_map(NUMBERING_TEMPLATE)
            apply_template_numbering(
                "Final_Trajectory_TMP.pdb",
                "Final_Trajectory.pdb",
                resmap
            )
        else:
            os.rename("Final_Trajectory_TMP.pdb", "Final_Trajectory.pdb")

    if os.path.exists("Final_Trajectory_RAW.pdb"):
        os.remove("Final_Trajectory_RAW.pdb")

    print("✅ Final_Trajectory.xtc/pdb written.\n")
    print("🧬 First 5 residues after renumbering:")
    with open("Final_Trajectory.pdb") as f:
        seen = set()
        for line in f:
            if line.startswith("ATOM"):
                key = line[22:26]
                if key not in seen:
                    print("   ", line[17:26])
                    seen.add(key)
                if len(seen) == 5:
                    break

else:
    print("⏭️ Final_Trajectory.* already exists — skipping.")


# ============================================================
# STEP B — binding_pocket_only.xtc / .pdb
# ============================================================
if not (file_exists("binding_pocket_only.xtc") and file_exists("binding_pocket_only.pdb")):

    if PROTEIN_ONLY or (PEPTIDE_LIKE and not PARTNER_MODE):
        print("🧬 Protein/peptide-centric mode: using full centered protein trajectory as the analysis subset.")
        shutil.copyfile("Final_Trajectory.xtc", "binding_pocket_only.xtc")
        shutil.copyfile("Final_Trajectory.pdb", "binding_pocket_only.pdb")
        print("✅ binding_pocket_only.* created as protein-only working copy.\n")

    else:
        print(f"🔬 Extracting binding pocket using MDTraj with a {POCKET_CUTOFF_ANG:.2f} Å cutoff...")

        traj = md.load("Final_Trajectory.xtc", top="Final_Trajectory.pdb")

        lig_sel = "resname == " + repr(lig)
        lig_atoms = traj.topology.select(lig_sel)
        prot_atoms = traj.topology.select("protein")

        if len(lig_atoms) == 0:
            raise SystemExit(f"❌ Ligand {lig} not found in Final_Trajectory.pdb.")

        neighbors = md.compute_neighbors(
            traj,
            cutoff=POCKET_CUTOFF_NM,
            query_indices=lig_atoms,
            haystack_indices=prot_atoms
        )

        if len(neighbors) == 0:
            raise SystemExit(f"❌ No protein atoms within {POCKET_CUTOFF_ANG:.1f} Å of ligand.")

        pocket_atoms = np.unique(np.concatenate(neighbors))

        pocket_res = sorted(
            {traj.topology.atom(idx).residue.index for idx in pocket_atoms}
        )

        print(f"🎯 Binding pocket contains {len(pocket_res)} residues.")

        sel = traj.topology.select(
            "resid " + " ".join(map(str, pocket_res)) + " or " + ("resname == " + repr(lig))
        )

        pocket = traj.atom_slice(sel)
        pocket.save("binding_pocket_only.xtc")
        pocket[0].save("binding_pocket_only_RAW.pdb")

        convert_ligand_to_hetatm("binding_pocket_only_RAW.pdb",
                                 "binding_pocket_only.pdb", lig)
        os.remove("binding_pocket_only_RAW.pdb")

        print("✅ binding_pocket_only.xtc/pdb written.\n")

else:
    print("⏭️ binding_pocket_only.* already exists — skipping.")

phase("STEP A+B complete")






# ============================================================
# STEP C — Determine Which Chains Interact with the Ligand/Peptide
# ============================================================
phase("STEP C — Identify Binding Chains for Downstream RMSD")

print("\n🔍 Identifying chains for downstream biomolecule analysis...")

u_pocket = mda.Universe("binding_pocket_only.pdb")
protein_residues = u_pocket.select_atoms("protein").residues

atomindex_entries = parse_atomindex_entries("atomIndex.txt")
chain_map = parse_atomindex_file("atomIndex.txt")
chain_entry_meta = {entry["label"]: entry for entry in atomindex_entries}
if not chain_map:
    print("⚠️ atomIndex.txt not found or contained no parseable chain ranges — treating entire protein as one chain")
    chain_map["Protein"] = (protein_residues.atoms[0].index,
                            protein_residues.atoms[-1].index)
    chain_entry_meta["Protein"] = {
        "label": "Protein",
        "start": protein_residues.atoms[0].index,
        "end": protein_residues.atoms[-1].index,
        "current_chain": "Protein",
        "chain_type": "protein",
    }

active_chains = {}

if PROTEIN_ONLY or (PEPTIDE_LIKE and not PARTNER_MODE):
    if BIOMOLECULE_MODE:
        print("🧬 Biomolecule-centric mode: treating all available polymer chains as active.")
    elif PEPTIDE_LIKE:
        print("🧬 Protein-peptide/protein-protein mode: treating all parsed polymer chains as active.")
    else:
        print("🧬 Protein-centric mode: treating all available protein chains as active.")
    active_chains = dict(chain_map)

else:
    lig_atomgroup = u_pocket.select_atoms(PARTNER_SELECTION_MDA or f"resname {ligand_code}")
    if lig_atomgroup.n_atoms == 0:
        print("DEBUG pocket resnames:", sorted(set(u_pocket.residues.resnames)))
        raise SystemExit(f"❌ Could not find ligand '{ligand_code}' in binding_pocket_only.pdb")

    for cname, (start, end) in chain_map.items():
        overlapping = [
            r for r in protein_residues
            if start <= r.atoms[0].index <= end
        ]
        if overlapping:
            active_chains[cname] = (start, end)
            print(f"✅ {cname} participates in ligand binding ({len(overlapping)} pocket residues).")
        else:
            print(f"⚪ {cname} does not interact with ligand.")

if not active_chains:
    print("⚠️ No active protein chains detected!")
else:
    print("📎 Active binding chains:", ", ".join(active_chains.keys()))

polymer_chain_entries = collect_polymer_chain_entries(chain_map, chain_entry_meta, SYSTEM_REGISTRY)
INTERFACE_CHAIN_COUNT = len(polymer_chain_entries)
INTERFACE_RIN_MODE = (not PROTAC_MODE) and (not PARTNER_MODE) and INTERFACE_CHAIN_COUNT >= 2
print(f"🧬 Polymer chains detected from atomIndex.txt: {INTERFACE_CHAIN_COUNT}")
print(f"🧬 Interface RIN mode active: {'yes' if INTERFACE_RIN_MODE else 'no'}")
# ============================================================
# ✅ Ready for downstream RMSD/RMSF and contact analysis
# ============================================================


# ============================================================
# ✅ Ready for downstream RMSD/RMSF and contact analysis
# ============================================================
TOPO_FILE = "Final_Trajectory.pdb"
TRAJ_FILE = "Final_Trajectory.xtc"



# ------------------------- Helpers -------------------------
def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def ensure_interaction_summary_columns(summary, plot_types):
    """Return a summary table with stable interaction-type columns and Total."""
    plot_types = list(plot_types)

    if summary is None or summary.empty:
        normalized = pd.DataFrame(columns=plot_types, dtype=int)
    else:
        normalized = summary.copy().reindex(columns=plot_types, fill_value=0)
        normalized = normalized.astype(int, copy=False)

    normalized["Total"] = normalized.reindex(columns=plot_types, fill_value=0).sum(axis=1)
    return normalized


def infer_element(atom_name: str):
    c = atom_name.strip()[:1]
    return c if c in "ONCHSP" else None

# ============================================================
# Utilities for Checkpointing and RMSD Compatibility
# ============================================================

ANALYSIS_MANIFEST_PATH = os.path.join(OUTPUT_DIR, "analysis_manifest.json")

def checkpoint(path):
    with open(path, "w") as f:
        f.write("done")

def done_checkpoint(path):
    return os.path.exists(path)


def _analysis_input_files():
    return {
        "topology": str(Path(TOPO_FILE).resolve()),
        "trajectory": str(Path(TRAJ_FILE).resolve()),
        "atomindex": str(Path("atomIndex.txt").resolve()),
    }


def _analysis_mode():
    if PROTAC_MODE:
        return "protac_like"
    if INTERFACE_RIN_MODE:
        return "protein_interface"
    if PARTNER_MODE:
        return "ligand_bound"
    return "protein_only"


def _stage_signature(data_parameters=None):
    return compute_signature(
        {
            "inputs": fingerprint_files(_analysis_input_files()),
            "parameters": data_parameters or {},
        }
    )


def _record_stage_manifest(stage_name, numerical_outputs, plot_outputs=None, data_parameters=None, legacy_adopted=False):
    return record_stage_completion(
        ANALYSIS_MANIFEST_PATH,
        stage_name=stage_name,
        analysis_script=Path(__file__).name,
        analysis_source_path=__file__,
        input_files=_analysis_input_files(),
        data_parameters=data_parameters or {},
        plot_parameters={},
        outputs=numerical_outputs,
        plot_outputs=plot_outputs or [],
        system=_analysis_mode(),
        mode=_analysis_mode(),
        legacy_adopted=legacy_adopted,
    )


def _stage_cache_ready(stage_name, checkpoint_name, numerical_outputs, plot_outputs=None, data_parameters=None):
    status = validate_stage_cache(
        manifest_path=ANALYSIS_MANIFEST_PATH,
        stage_name=stage_name,
        required_outputs=numerical_outputs,
        current_data_signature=_stage_signature(data_parameters),
        checkpoint_path=checkpoint_name,
    )
    if status["reusable"]:
        if status["legacy_adopt"]:
            _record_stage_manifest(
                stage_name,
                numerical_outputs=numerical_outputs,
                plot_outputs=plot_outputs or [],
                data_parameters=data_parameters,
                legacy_adopted=True,
            )
        return True

    if os.path.exists(checkpoint_name):
        if status["reason"] == "required_outputs_missing":
            print(f"⚠️ {stage_name}: checkpoint exists but required reusable outputs are missing or empty; rerunning stage.")
        elif status["reason"] == "data_signature_mismatch":
            print(f"⚠️ {stage_name}: checkpoint exists but data-affecting inputs changed; rerunning stage.")
    return False


def _current_chain_names():
    chain_map = parse_atomindex_file("atomIndex.txt")
    if chain_map:
        return list(chain_map.keys())
    if active_chains:
        return list(active_chains.keys())
    return ["Protein"]


def _expected_dssp_chain_labels():
    count = max(1, len(_current_chain_names()))
    labels = []
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for idx in range(count):
        if idx < len(alphabet):
            labels.append(alphabet[idx])
        else:
            labels.append(f"{alphabet[idx // 26 - 1]}{alphabet[idx % 26]}")
    return labels


def _ligand_token():
    return (ligand_code or compound_name or "Ligand").upper()


def _compound_token():
    return compound_name or _ligand_token()


def _stage_d1_outputs():
    return ([os.path.join(OUTPUT_DIR, "Protein_RMSD.csv")], [os.path.join(OUTPUT_DIR, "Protein_RMSD.png")], {})


def _stage_d2_outputs():
    ligand_token = _ligand_token()
    compound_token = _compound_token()
    return (
        [
            os.path.join(OUTPUT_DIR, f"{compound_token}_Ligand_RMSD.csv"),
            os.path.join(OUTPUT_DIR, f"{compound_token}_Ligand_RMSF.csv"),
        ],
        [
            os.path.join(OUTPUT_DIR, f"{ligand_token}_Ligand_RMSD.png"),
            os.path.join(OUTPUT_DIR, f"{ligand_token}_Ligand_RMSF.png"),
        ],
        {"ligand_code": ligand_token},
    )


def _stage_d3_outputs():
    chains = _current_chain_names()
    numerical = []
    plots = []
    for chain in chains:
        numerical.extend(
            [
                os.path.join(OUTPUT_DIR, f"RMSD_{chain}.csv"),
                os.path.join(OUTPUT_DIR, f"RMSF_{chain}.csv"),
                os.path.join(OUTPUT_DIR, f"RMSF_{chain}_per_residue.csv"),
            ]
        )
        plots.extend(
            [
                os.path.join(OUTPUT_DIR, f"RMSD_{chain}.png"),
                os.path.join(OUTPUT_DIR, f"RMSF_{chain}.png"),
                os.path.join(OUTPUT_DIR, f"RMSF_{chain}_per_residue_contacts.png"),
            ]
        )
        if PARTNER_MODE:
            numerical.append(os.path.join(OUTPUT_DIR, f"{chain}_LigandContactPersistence.csv"))
    params = {"chains": chains}
    if PARTNER_MODE:
        params.update({"contact_cutoff_A": float(args.contact_cutoff), "min_contact_frac": float(args.min_contact_frac)})
    return numerical, plots, params


def _stage_d4_outputs():
    ligand_token = _ligand_token()
    return (
        [os.path.join(OUTPUT_DIR, f"{ligand_token}_Complex_RMSD_Overlay.csv")],
        [os.path.join(OUTPUT_DIR, f"{ligand_token}_Complex_RMSD_Overlay.png")],
        {"ligand_code": ligand_token},
    )


def _stage_d4b_outputs():
    if PARTNER_MODE:
        return ([os.path.join(OUTPUT_DIR, "Radius_of_Gyration_Protein_Ligand_Complex.csv")], [os.path.join(OUTPUT_DIR, "Radius_of_Gyration_Overlay.png")], {})
    return ([os.path.join(OUTPUT_DIR, "Radius_of_Gyration_Protein.csv")], [os.path.join(OUTPUT_DIR, "Radius_of_Gyration_Protein.png")], {})


def _stage_d5_outputs():
    return ([os.path.join(OUTPUT_DIR, "FullComplex_RMSD.csv")], [os.path.join(OUTPUT_DIR, "FullComplex_RMSD.png")], {})


def _stage_d6_outputs():
    return (
        [
            os.path.join(OUTPUT_DIR, "DSSP_raw.csv"),
            os.path.join(OUTPUT_DIR, "DSSP_encoded_numeric.csv"),
            os.path.join(OUTPUT_DIR, "SSE_Fractions.csv"),
        ],
        [os.path.join(OUTPUT_DIR, "DSSP_Heatmap.png"), os.path.join(OUTPUT_DIR, "SSE_Fractions.png")],
        {},
    )


def _stage_d6b_outputs():
    numerical = []
    plots = []
    labels = _expected_dssp_chain_labels()
    for label in labels:
        numerical.extend([os.path.join(OUTPUT_DIR, f"DSSP_chain_{label}.csv"), os.path.join(OUTPUT_DIR, f"SSE_chain_{label}_fractions.csv")])
        plots.extend([os.path.join(OUTPUT_DIR, f"DSSP_chain_{label}_heatmap.png"), os.path.join(OUTPUT_DIR, f"SSE_chain_{label}_fractions.png")])
    return numerical, plots, {"chains": labels}


def _stage_d6c_outputs():
    return (
        [os.path.join(OUTPUT_DIR, "Residue_SSE_Persistence.csv")],
        [
            os.path.join(OUTPUT_DIR, "Residue_Helix_Persistence.png"),
            os.path.join(OUTPUT_DIR, "Residue_Sheet_Persistence.png"),
            os.path.join(OUTPUT_DIR, "Residue_Flexibility.png"),
        ],
        {},
    )


def _stage_d6d_outputs():
    ligand_token = _ligand_token()
    return (
        [
            os.path.join(OUTPUT_DIR, f"{ligand_token}_ChangePoint_Summary.csv"),
            os.path.join(OUTPUT_DIR, f"{ligand_token}_SSE_Contact_Coupling_Timeseries.csv"),
        ],
        [os.path.join(OUTPUT_DIR, f"{ligand_token}_SSE_vs_LigandContacts_Dynamics.png")],
        {"contact_cutoff_A": float(args.contact_cutoff), "ligand_code": ligand_token},
    )


def _interface_chain_pairs_for_manifest():
    if args.interface_chain_pairs:
        return [token.strip() for token in str(args.interface_chain_pairs).split(",") if token.strip()]
    labels = []
    indexed_entries = list(enumerate(polymer_chain_entries, start=1))
    for (left_idx, left), (right_idx, right) in combinations(indexed_entries, 2):
        left_label = chain_display_label(left, left_idx)
        right_label = chain_display_label(right, right_idx)
        labels.append(f"{left_label}:{right_label}")
    return labels


def interface_rin_outputs_present(outdir):
    required = [
        os.path.join(outdir, "interface_chain_pair_summary.csv"),
        os.path.join(outdir, "interface_contact_timeseries.csv"),
        os.path.join(outdir, "interface_residue_pair_contacts.csv"),
        os.path.join(outdir, "interface_node_summary.csv"),
        os.path.join(outdir, "interface_frame_index.csv"),
        os.path.join(outdir, "interface_residue_pair_framewise.csv.gz"),
        os.path.join(outdir, "interface_residue_pair_presence.npz"),
    ]
    return all(os.path.exists(path) and os.path.getsize(path) > 0 for path in required)


def dispatch_to_interface_rin():
    interface_script = os.path.join(os.getcwd(), "3C_Interface_RIN.py")
    if not os.path.exists(interface_script):
        raise SystemExit(
            "Interface RIN mode selected, but 3C_Interface_RIN.py was not found in the current directory."
        )

    outdir = os.path.join(args.outdir, "Interface_RIN")
    cmd = [
        sys.executable,
        interface_script,
        "--topo",
        args.topo,
        "--traj",
        args.traj,
        "--atomindex",
        "atomIndex.txt",
        "--outdir",
        outdir,
        "--contact-cutoff",
        str(args.interface_contact_cutoff),
        "--min-contact-frac",
        str(args.interface_min_contact_frac),
        "--frame-step",
        str(args.interface_frame_step),
        "--max-edges",
        str(args.interface_max_edges),
    ]
    if args.interface_chain_pairs:
        cmd.extend(["--chain-pairs", args.interface_chain_pairs])
    if args.headless:
        cmd.append("--headless")

    print(f"🧬 Interface RIN output folder: {outdir}")
    print(f"🚀 Launching interface RIN with command: {' '.join(cmd)}")
    rc = subprocess.call(cmd)
    if rc != 0:
        raise SystemExit(rc)

    if not interface_rin_outputs_present(outdir):
        raise SystemExit("❌ 3C_Interface_RIN.py completed but expected interface outputs were not found.")

    _record_stage_manifest(
        "EINTERFACE",
        numerical_outputs=[
            os.path.join(outdir, "interface_chain_pair_summary.csv"),
            os.path.join(outdir, "interface_contact_timeseries.csv"),
            os.path.join(outdir, "interface_residue_pair_contacts.csv"),
            os.path.join(outdir, "interface_node_summary.csv"),
            os.path.join(outdir, "interface_frame_index.csv"),
            os.path.join(outdir, "interface_residue_pair_framewise.csv.gz"),
            os.path.join(outdir, "interface_residue_pair_presence.npz"),
        ],
        plot_outputs=[
            os.path.join(outdir, "Interface_ChainPair_ContactFractions.png"),
            os.path.join(outdir, "Interface_ContactCount_Timeseries.png"),
            os.path.join(outdir, "Interface_MinDistance_Timeseries.png"),
            os.path.join(outdir, "Interface_ResiduePair_ContactHeatmap.png"),
            os.path.join(outdir, "Interface_RIN_Network.png"),
            os.path.join(outdir, "Interface_ChainLevel_Network.png"),
            os.path.join(outdir, "Interface_TopResiduePairs_Barplot.png"),
        ],
        data_parameters={
            "contact_cutoff_A": float(args.interface_contact_cutoff),
            "frame_step": int(args.interface_frame_step),
            "selected_chain_pairs": _interface_chain_pairs_for_manifest(),
        },
    )
    checkpoint("chk_EINTERFACE.txt")
    print("✅ STEP E-INTERFACE complete.")

def safe_rmsd(traj, ref, frame, atom_indices):
    """Try modern MDTraj rst syntax; fall back when unsupported."""
    try:
        return md.rmsd(traj, ref, frame, atom_indices=atom_indices, parallel=True)
    except TypeError:
        return md.rmsd(traj, ref, frame, atom_indices=atom_indices)

def align_chunk(args):
        traj, bb, idxs = args
        sub = traj.slice(idxs)
        sub.superpose(traj, 0, atom_indices=bb)
        return sub

# ============================================================
# MEMORY-SAFE TRAJECTORY HELPERS
# ============================================================

CHUNK_SIZE = max(1, int(args.chunk_size))

_traj_meta_cache = {
    "ref": None,
    "bb": None,
    "ligand_atoms": None,
    "lig_resname": None,
}

_dssp_numeric_cache = None
_dssp_time_ns_cache = None

def get_traj_metadata():
    """Load only frame 0 and topology-derived selections."""
    global ligand_code

    if _traj_meta_cache["ref"] is not None:
        return (_traj_meta_cache["ref"],
                _traj_meta_cache["bb"],
                _traj_meta_cache["ligand_atoms"],
                _traj_meta_cache["lig_resname"])

    print("\n🌍 Loading reference frame and topology selections…")
    ref = md.load_frame(TRAJ_FILE, 0, top=TOPO_FILE)
    bb = ref.topology.select("backbone")

    if len(bb) == 0:
        raise RuntimeError("❌ No backbone atoms found for alignment.")

    ligand_atoms = np.array([], dtype=int)
    lig_resname = None

    if PARTNER_MODE and ligand_code:
        lig_resname = ligand_code.upper()
        if PEPTIDE_LIKE:
            atoms = mdtraj_partner_indices(ref, peptide_info=PEPTIDE_INFO)
            if len(atoms) == 0:
                raise RuntimeError(
                    "❌ Peptide partner atoms could not be found in the trajectory using .chainmap.peptide.json."
                )
            ligand_atoms = np.array(atoms, dtype=int)
            print(f"🔍 Using peptide partner: {lig_resname} ({len(ligand_atoms)} atoms)")
        else:
            atoms = ref.topology.select("resname == " + repr(lig_resname))
            if len(atoms) == 0:
                raise RuntimeError(
                    f"❌ Ligand '{lig_resname}' not found in trajectory.\n"
                    "Check Final_Trajectory.pdb for residue naming."
                )
            ligand_atoms = np.array(atoms, dtype=int)
            print(f"🔍 Using user-provided ligand: {lig_resname} ({len(ligand_atoms)} atoms)")
    else:
        print("🧬 Protein-centric mode: ligand-specific arrays disabled.")

    _traj_meta_cache.update({
        "ref": ref,
        "bb": bb,
        "ligand_atoms": ligand_atoms,
        "lig_resname": lig_resname,
    })

    return ref, bb, ligand_atoms, lig_resname

def iter_traj_chunks():
    for chunk in md.iterload(TRAJ_FILE, top=TOPO_FILE, chunk=CHUNK_SIZE):
        yield chunk

def iter_aligned_chunks(ref=None, bb=None):
    if ref is None or bb is None:
        ref, bb, _, _ = get_traj_metadata()
    for chunk in iter_traj_chunks():
        chunk.superpose(ref, 0, atom_indices=bb)
        yield chunk

def chunked_rmsd(atom_indices, label="RMSD"):
    ref, bb, _, _ = get_traj_metadata()
    atom_indices = np.asarray(atom_indices, dtype=int)
    times = []
    values = []

    for chunk in tqdm(iter_aligned_chunks(ref, bb), desc=f"📈 {label}", unit="chunk"):
        times.append(np.asarray(chunk.time, dtype=float))
        values.append(safe_rmsd(chunk, ref, 0, atom_indices=atom_indices) * 10.0)

    return np.concatenate(times), np.concatenate(values)

def chunked_multi_rmsd(selection_map, label="RMSD"):
    ref, bb, _, _ = get_traj_metadata()
    times = []
    values = {name: [] for name in selection_map}

    for chunk in tqdm(iter_aligned_chunks(ref, bb), desc=f"📈 {label}", unit="chunk"):
        times.append(np.asarray(chunk.time, dtype=float))
        for name, atom_indices in selection_map.items():
            values[name].append(
                safe_rmsd(chunk, ref, 0, atom_indices=np.asarray(atom_indices, dtype=int)) * 10.0
            )

    return np.concatenate(times), {
        name: np.concatenate(chunks) for name, chunks in values.items()
    }

def chunked_rmsf(atom_indices, label="RMSF"):
    ref, bb, _, _ = get_traj_metadata()
    atom_indices = np.asarray(atom_indices, dtype=int)
    n_atoms = len(atom_indices)
    count = 0
    sum_pos = np.zeros((n_atoms, 3), dtype=np.float64)
    sum_sq = np.zeros((n_atoms, 3), dtype=np.float64)

    for chunk in tqdm(iter_aligned_chunks(ref, bb), desc=f"📈 {label}", unit="chunk"):
        xyz = chunk.xyz[:, atom_indices, :].astype(np.float64, copy=False)
        count += xyz.shape[0]
        sum_pos += xyz.sum(axis=0)
        sum_sq += (xyz * xyz).sum(axis=0)

    mean = sum_pos / max(count, 1)
    variance = np.maximum((sum_sq / max(count, 1)) - (mean * mean), 0.0)
    return np.sqrt(variance.sum(axis=1)) * 10.0, count

def chunked_rg(selection_map):
    times = []
    values = {name: [] for name in selection_map}

    for chunk in tqdm(iter_traj_chunks(), desc="📐 Radius of gyration", unit="chunk"):
        times.append(np.asarray(chunk.time, dtype=float))
        for name, atom_indices in selection_map.items():
            atom_indices = np.asarray(atom_indices, dtype=int)
            sub = chunk.atom_slice(atom_indices) if atom_indices.size != chunk.n_atoms else chunk
            values[name].append(md.compute_rg(sub) * 10.0)

    return np.concatenate(times), {
        name: np.concatenate(chunks) for name, chunks in values.items()
    }

def chunked_contact_counts(ligand_atoms, target_atoms, residue_ids, cutoff_angstrom):
    ligand_atoms = np.asarray(ligand_atoms, dtype=int)
    target_atoms = np.asarray(target_atoms, dtype=int)
    contact_counts = {r: 0 for r in residue_ids}
    n_frames = 0

    for chunk in tqdm(iter_aligned_chunks(), desc="🔍 Ligand contact persistence", unit="chunk"):
        for frame_i in range(chunk.n_frames):
            lig_xyz = chunk.xyz[frame_i, ligand_atoms, :]
            target_xyz = chunk.xyz[frame_i, target_atoms, :]
            d = np.linalg.norm(lig_xyz[:, None, :] - target_xyz[None, :, :], axis=2)
            min_d = d.min(axis=0) * 10.0
            for r, dist in zip(residue_ids, min_d):
                if dist < cutoff_angstrom:
                    contact_counts[r] += 1
        n_frames += chunk.n_frames

    return contact_counts, n_frames

def encode_dssp_numeric(dssp):
    encoded = np.zeros(dssp.shape, dtype=np.int8)
    encoded[np.isin(dssp, ["E", "B"])] = 1
    encoded[np.isin(dssp, ["H", "G", "I"])] = 2
    return encoded

def downsample_frames_for_heatmap(matrix, max_frames=5000):
    if matrix.shape[0] <= max_frames:
        return matrix, 1
    stride = int(np.ceil(matrix.shape[0] / max_frames))
    return matrix[::stride], stride

def compute_dssp_chunked(write_outputs=False):
    global _dssp_numeric_cache, _dssp_time_ns_cache

    raw_path = os.path.join(OUTPUT_DIR, "DSSP_raw.csv")
    encoded_path = os.path.join(OUTPUT_DIR, "DSSP_encoded_numeric.csv")
    if write_outputs:
        for path in (raw_path, encoded_path):
            if os.path.exists(path):
                os.remove(path)

    encoded_chunks = []
    time_chunks = []
    raw_header = True
    encoded_header = True

    for chunk in tqdm(iter_traj_chunks(), desc="📐 DSSP", unit="chunk"):
        dssp = md.compute_dssp(chunk)
        encoded = encode_dssp_numeric(dssp)
        encoded_chunks.append(encoded)
        time_chunks.append(np.asarray(chunk.time / 1000.0, dtype=float))

        if write_outputs:
            pd.DataFrame(dssp).to_csv(
                raw_path, mode="a", header=raw_header, index=False
            )
            pd.DataFrame(encoded).to_csv(
                encoded_path, mode="a", header=encoded_header, index=False
            )
            raw_header = False
            encoded_header = False

    _dssp_numeric_cache = np.vstack(encoded_chunks)
    _dssp_time_ns_cache = np.concatenate(time_chunks)
    return _dssp_numeric_cache, _dssp_time_ns_cache

def get_dssp_numeric_and_time():
    global _dssp_numeric_cache, _dssp_time_ns_cache

    if _dssp_numeric_cache is not None and _dssp_time_ns_cache is not None:
        return _dssp_numeric_cache, _dssp_time_ns_cache

    encoded_path = os.path.join(OUTPUT_DIR, "DSSP_encoded_numeric.csv")
    fractions_path = os.path.join(OUTPUT_DIR, "SSE_Fractions.csv")

    if os.path.exists(encoded_path):
        _dssp_numeric_cache = pd.read_csv(encoded_path).to_numpy(dtype=np.int8)
        if os.path.exists(fractions_path):
            _dssp_time_ns_cache = pd.read_csv(fractions_path)["Time_ns"].to_numpy(dtype=float)
        else:
            _dssp_time_ns_cache = chunked_times_ns()
        return _dssp_numeric_cache, _dssp_time_ns_cache

    return compute_dssp_chunked(write_outputs=True)

def chunked_times_ns():
    times = []
    for chunk in tqdm(iter_traj_chunks(), desc="⏱️ Reading frame times", unit="chunk"):
        times.append(np.asarray(chunk.time / 1000.0, dtype=float))
    return np.concatenate(times)

def dssp_residue_chain_indices(topology, n_dssp_residues):
    all_chains = np.array([res.chain.index for res in topology.residues])
    if all_chains.size == n_dssp_residues:
        return all_chains

    protein_chains = np.array([
        res.chain.index
        for res in topology.residues
        if getattr(res, "is_protein", False)
    ])
    if protein_chains.size == n_dssp_residues:
        return protein_chains

    print(
        "⚠️ DSSP residue count did not exactly match topology residues; "
        "using the overlapping residue range for per-chain SSE plots."
    )
    return all_chains[:n_dssp_residues]

def get_aligned_traj_and_ligand():
    raise RuntimeError(
        "This folder's 3A_AutomateGromacs.py has been patched for chunked analysis. "
        "A remaining code path tried to load the full trajectory at once."
    )
# ============================================================
# STEP D1 — GLOBAL PROTEIN RMSD
# ============================================================

# ------------------ STEP D1 replacement (Protein RMSD) ------------------
_d1_num, _d1_plot, _d1_params = _stage_d1_outputs()
if not _stage_cache_ready("D1", "chk_D1.txt", _d1_num, _d1_plot, _d1_params):
    phase("STEP D1 — Global Protein RMSD")

    ref, bb, ligand_atoms, lig_resname = get_traj_metadata()

    print(f"📈 Computing protein backbone RMSD in chunks of {CHUNK_SIZE} frames…")
    time_ps, prot_rmsd = chunked_rmsd(bb, label="Protein backbone RMSD")

    df = pd.DataFrame({"Time_ps": time_ps, "Protein_RMSD": prot_rmsd})
    df.to_csv(os.path.join(OUTPUT_DIR, "Protein_RMSD.csv"), index=False)

    plt.figure(figsize=(9,5))
    plt.plot(time_ps/1000.0, prot_rmsd, lw=1.8, color=MIAMI_ORANGE)   # use Miami orange for protein
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("Backbone Residue RMSD")
    savefig(os.path.join(OUTPUT_DIR, "Protein_RMSD.png"))

    checkpoint("chk_D1.txt")
    _record_stage_manifest("D1", _d1_num, _d1_plot, _d1_params)
    print("✅ STEP D1 complete.")
else:
    print("⏭️ STEP D1 skipped (checkpoint found).")


# ============================================================
# STEP D2 — LIGAND RMSD & RMSF
# ============================================================

# ------------------ STEP D2 replacement (Ligand RMSD & RMSF) ------------------
if not PARTNER_MODE:
    print("⏭️ STEP D2 skipped (protein-centric mode).")
else:
    _d2_num, _d2_plot, _d2_params = _stage_d2_outputs()
    if _stage_cache_ready("D2", "chk_D2.txt", _d2_num, _d2_plot, _d2_params):
        print("⏭️ STEP D2 skipped (checkpoint found).")
    else:
        phase("STEP D2 — Ligand RMSD & RMSF")

        ref, bb, ligand_atoms, lig_resname = get_traj_metadata()

        print("📈 Ligand RMSD…")
        time_ps, lig_rmsd = chunked_rmsd(ligand_atoms, label="Ligand RMSD")

        print("📈 Ligand RMSF…")
        rmsf, _ = chunked_rmsf(ligand_atoms, label="Ligand RMSF")

        pd.DataFrame({"Time_ps": time_ps, "Ligand_RMSD": lig_rmsd}).to_csv(
            os.path.join(OUTPUT_DIR, f"{compound_name}_Ligand_RMSD.csv"), index=False
        )
        pd.DataFrame({"Atom": ligand_atoms, "RMSF": rmsf}).to_csv(
            os.path.join(OUTPUT_DIR, f"{compound_name}_Ligand_RMSF.csv"), index=False
        )

        plt.figure(figsize=(9,5))
        plt.plot(time_ps/1000.0, lig_rmsd, color=MIAMI_GREEN, lw=1.6)
        plt.title(f"{compound_name} Ligand RMSD")
        plt.xlabel("Time (ns)")
        plt.ylabel("RMSD (Å)")
        savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_Ligand_RMSD.png"))

        plt.figure(figsize=(8,4))
        plt.plot(rmsf, color=MIAMI_GREEN, lw=1.6)
        plt.title(f"{compound_name} Ligand RMSF")
        plt.xlabel("Ligand Atom Index")
        plt.ylabel("RMSF (Å)")
        savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_Ligand_RMSF.png"))

        checkpoint("chk_D2.txt")
        _record_stage_manifest("D2", _d2_num, _d2_plot, _d2_params)
        print("✅ STEP D2 complete.")
# ============================================================
# STEP D3 — PER-CHAIN RMSD / RMSF (ATOM + RESIDUE + PERSISTENCE)
# ============================================================

_d3_num, _d3_plot, _d3_params = _stage_d3_outputs()
if not _stage_cache_ready("D3", "chk_D3.txt", _d3_num, _d3_plot, _d3_params):
    phase("STEP D3 — Per-Chain RMSD / RMSF")

    ref, bb, ligand_atoms, lig_resname = get_traj_metadata()
    n_frames = None

    # --------------------------------------------------------
    # Load chain map
    # --------------------------------------------------------
    atomindex_entries = parse_atomindex_entries("atomIndex.txt")
    chain_map = parse_atomindex_file("atomIndex.txt")
    chain_entry_meta = {entry["label"]: entry for entry in atomindex_entries}

    if not chain_map:
        print("⚠️ atomIndex.txt not found or contained no parseable chain ranges — using full protein as one chain")
        protein_atoms = ref.topology.select("protein")
        if len(protein_atoms) == 0:
            print("❌ No protein atoms found in trajectory.")
            checkpoint("chk_D3.txt")
            print("✅ STEP D3 complete.")
        else:
            chain_map["Protein"] = (int(protein_atoms.min()), int(protein_atoms.max()))

    # --------------------------------------------------------
    # Loop over chains
    # --------------------------------------------------------
    for chain, (start, end) in chain_map.items():

        # -------------------------
        # Clamp invalid ranges
        # -------------------------
        orig_start, orig_end = start, end

        start = max(0, start)
        end = min(end, ref.n_atoms - 1)

        if orig_start != start or orig_end != end:
            print(
                f"⚠️ {chain}: atomIndex range {orig_start}-{orig_end} "
                f"adjusted to valid range {start}-{end}"
            )

        if start > end:
            print(f"⚠️ {chain}: invalid atom range after clamping {start}-{end}")
            continue

        sel = np.arange(start, end + 1, dtype=int)

        if sel.size == 0:
            print(f"⚠️ {chain}: empty atom selection after clamping")
            continue

        # =========================
        # RMSD (Å)
        # =========================
        time_ps, rmsd = chunked_rmsd(sel, label=f"{chain} RMSD")

        pd.DataFrame({
            "Time_ps": time_ps,
            "RMSD_Å": rmsd
        }).to_csv(
            os.path.join(OUTPUT_DIR, f"RMSD_{chain}.csv"),
            index=False
        )

        plt.figure(figsize=(9, 5))
        plt.plot(time_ps / 1000.0, rmsd, color=MIAMI_ORANGE, lw=1.8)
        plt.title(f"{chain} RMSD")
        plt.xlabel("Time (ns)")
        plt.ylabel("RMSD (Å)")
        savefig(os.path.join(OUTPUT_DIR, f"RMSD_{chain}.png"))
        plt.close()

        # =========================
        # ATOM-LEVEL RMSF (Å)
        # =========================
        rmsf_atom, n_frames = chunked_rmsf(sel, label=f"{chain} atom RMSF")

        pd.DataFrame({
            "Atom_Index": sel,
            "RMSF_Å": rmsf_atom
        }).to_csv(
            os.path.join(OUTPUT_DIR, f"RMSF_{chain}.csv"),
            index=False
        )

        plt.figure(figsize=(9, 4))
        plt.plot(rmsf_atom, color=MIAMI_ORANGE, lw=1.6)
        plt.title(f"{chain} RMSF (per atom)")
        plt.xlabel("Atom Index")
        plt.ylabel("RMSF (Å)")
        savefig(os.path.join(OUTPUT_DIR, f"RMSF_{chain}.png"))
        plt.close()

        # =========================
        # RESIDUE-LEVEL RMSF (anchor atoms)
        # =========================
        entry_meta = chain_entry_meta.get(chain, {})
        current_chain_id = entry_meta.get("current_chain")
        chain_type = entry_meta.get("chain_type")
        if not chain_type and current_chain_id and SYSTEM_REGISTRY.get("chain_types"):
            chain_type = SYSTEM_REGISTRY["chain_types"].get(current_chain_id)
        if not chain_type:
            chain_type = "mixed_or_unknown"

        anchor_indices, res_ids, anchor_name = select_chain_anchor_atoms(
            ref, start, end, chain_type=chain_type
        )

        if anchor_indices.size == 0:
            print(
                f"⚠️ {chain}: no residue anchor atoms found for "
                f"{polymer_type_label(chain_type)} chain"
            )
            continue

        rmsf_res, n_frames = chunked_rmsf(anchor_indices, label=f"{chain} residue RMSF")

        pd.DataFrame({
            "Residue": res_ids,
            "RMSF_Å": rmsf_res
        }).to_csv(
            os.path.join(OUTPUT_DIR, f"RMSF_{chain}_per_residue.csv"),
            index=False
        )

        # =========================
        # CONTACT PERSISTENCE
        # Only for ligand-bearing runs
        # =========================
        persistent_residues = set()

        has_ligand_for_contacts = (
            ligand_atoms is not None and
            len(ligand_atoms) > 0 and
            (not globals().get("PROTEIN_CENTRIC", False))
        )

        if has_ligand_for_contacts:
            INTERACTION_CUTOFF = float(args.contact_cutoff)
            contact_counts, n_frames = chunked_contact_counts(
                ligand_atoms, anchor_indices, res_ids, INTERACTION_CUTOFF
            )

            contact_frac = {
                r: contact_counts[r] / n_frames
                for r in res_ids
            }

            persistent_residues = {
                r for r, frac in contact_frac.items()
                if frac >= MIN_CONTACT_FRAC
            }

            pd.DataFrame({
                "Residue": res_ids,
                "Contact_Fraction": [contact_frac[r] for r in res_ids]
            }).to_csv(
                os.path.join(OUTPUT_DIR, f"{chain}_LigandContactPersistence.csv"),
                index=False
            )

        # =========================
        # BUILD RMSF LINE WITH TRUE GAPS
        # =========================
        res_plot = []
        rmsf_plot = []

        for i, r in enumerate(res_ids):
            if i > 0 and r - res_ids[i - 1] > 1:
                res_plot.append(np.nan)
                rmsf_plot.append(np.nan)
            res_plot.append(r)
            rmsf_plot.append(rmsf_res[i])

        res_plot = np.array(res_plot, dtype=float)
        rmsf_plot = np.array(rmsf_plot, dtype=float)

        from matplotlib.lines import Line2D

        # =========================
        # RESIDUE RMSF PLOT
        # =========================
        plt.figure(figsize=(10, 4))

        # RMSF curve (broken at gaps)
        plt.plot(res_plot, rmsf_plot, color=MIAMI_ORANGE, lw=1.8)

        # Persistent ligand-contact markers only when ligand exists
        if persistent_residues:
            for r, rmsf_val in zip(res_ids, rmsf_res):
                if r in persistent_residues:
                    plt.vlines(
                        x=r,
                        ymin=0,
                        ymax=rmsf_val,
                        color=MIAMI_GREEN,
                        linestyles="dashed",
                        alpha=0.6,
                        lw=1.2
                    )

            legend_elements = [
                Line2D(
                    [0], [0],
                    color=MIAMI_GREEN,
                    lw=1.5,
                    linestyle="dashed",
                    label=f"Ligand contact ≥ {MIN_CONTACT_FRAC:.0%} frames"
                )
            ]

            plt.legend(
                handles=legend_elements,
                loc="upper right",
                frameon=False
            )

        plt.title(f"{chain} RMSF per Residue")
        plt.xlabel("Residue Number")
        plt.ylabel("RMSF (Å)")
        if anchor_name:
            plt.suptitle(f"Anchor atom: {anchor_name}", y=1.02, fontsize=10)
        plt.tight_layout()

        savefig(
            os.path.join(
                OUTPUT_DIR,
                f"RMSF_{chain}_per_residue_contacts.png"
            )
        )
        plt.close()

    checkpoint("chk_D3.txt")
    _record_stage_manifest("D3", _d3_num, _d3_plot, _d3_params)
    print("✅ STEP D3 complete.")

else:
    print("⏭️ STEP D3 skipped (checkpoint found).")


# ============================================================
# STEP D4 — BOUND COMPLEX RMSD (Dual-Axis Protein + Ligand)
# ============================================================

# ------------------ STEP D4 replacement (Bound Complex RMSD) ------------------
if not PARTNER_MODE:
    print("⏭️ STEP D4 skipped (protein-centric mode).")
elif not _stage_cache_ready("D4", "chk_D4.txt", *_stage_d4_outputs()):
    phase("STEP D4 — Bound Complex RMSD (Dual-Axis)")

    ref, bb, ligand_atoms, lig_resname = get_traj_metadata()
    time_ps, rmsd_map = chunked_multi_rmsd(
        {"protein": bb, "ligand": ligand_atoms},
        label="Protein/ligand RMSD"
    )
    prot_rmsd = rmsd_map["protein"]
    lig_rmsd = rmsd_map["ligand"]

    time_ns = time_ps / 1000.0

    df = pd.DataFrame({
        "Time_ns": time_ns,
        "Protein_RMSD": prot_rmsd,
        "Ligand_RMSD": lig_rmsd
    })

    df.to_csv(os.path.join(OUTPUT_DIR,
                           f"{lig_resname}_Complex_RMSD_Overlay.csv"),
                           index=False)

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(time_ns, prot_rmsd, color=MIAMI_ORANGE, label="Protein RMSD")
    ax1.set_ylabel("Protein RMSD (Å)", color=MIAMI_ORANGE)
    ax1.tick_params(axis="y", labelcolor=MIAMI_ORANGE)

    ax2 = ax1.twinx()
    ax2.plot(time_ns, lig_rmsd, color=MIAMI_GREEN, label="Ligand RMSD")
    ax2.set_ylabel("Ligand RMSD (Å)", color=MIAMI_GREEN)
    ax2.tick_params(axis="y", labelcolor=MIAMI_GREEN)

    chain_label = ", ".join(active_chains.keys()) if active_chains else "Protein"
    ax1.set_title(f"RMSD of {compound_name} bound with {chain_label}")
    ax1.set_xlabel("Time (ns)")

    fig.tight_layout()
    savefig(os.path.join(OUTPUT_DIR,
                         f"{lig_resname}_Complex_RMSD_Overlay.png"))

    checkpoint("chk_D4.txt")
    _d4_num, _d4_plot, _d4_params = _stage_d4_outputs()
    _record_stage_manifest("D4", _d4_num, _d4_plot, _d4_params)
    print("✅ STEP D4 complete.")
else:
    print("⏭️ STEP D4 skipped (checkpoint found).")






# ============================================================
# STEP D4B — Radius of Gyration (Protein / Ligand / Complex)
# ============================================================

_d4b_num, _d4b_plot, _d4b_params = _stage_d4b_outputs()
if not _stage_cache_ready("D4B", "chk_D4B.txt", _d4b_num, _d4b_plot, _d4b_params):
    phase("STEP D4B — Radius of Gyration")

    ref, bb, ligand_atoms, lig_resname = get_traj_metadata()

    print("📐 Computing radius of gyration…")

    rg_selections = {"protein": bb}
    if PARTNER_MODE:
        rg_selections["ligand"] = ligand_atoms
        rg_selections["complex"] = np.arange(ref.n_atoms, dtype=int)

    time_ps, rg_map = chunked_rg(rg_selections)
    time_ns = time_ps / 1000.0
    rg_protein = rg_map["protein"]

    if not PARTNER_MODE:
        df = pd.DataFrame({
            "Time_ns": time_ns,
            "Rg_Protein_Å": rg_protein
        })
        df.to_csv(
            os.path.join(OUTPUT_DIR, "Radius_of_Gyration_Protein.csv"),
            index=False
        )

        plt.figure(figsize=(10,6))
        plt.plot(time_ns, rg_protein, lw=2.0, label="Protein Rg", color="#1E88E5")
        plt.xlabel("Time (ns)")
        plt.ylabel("Radius of Gyration (Å)")
        plt.title("Radius of Gyration — Protein")
        plt.legend(frameon=False)
        plt.grid(alpha=0.25)

        savefig(os.path.join(OUTPUT_DIR, "Radius_of_Gyration_Protein.png"))
    else:
        rg_ligand  = rg_map["ligand"]
        rg_complex = rg_map["complex"]

        df = pd.DataFrame({
            "Time_ns": time_ns,
            "Rg_Protein_Å": rg_protein,
            "Rg_Ligand_Å": rg_ligand,
            "Rg_Complex_Å": rg_complex
        })

        df.to_csv(
            os.path.join(OUTPUT_DIR, "Radius_of_Gyration_Protein_Ligand_Complex.csv"),
            index=False
        )

        plt.figure(figsize=(10,6))
        plt.plot(time_ns, rg_protein, lw=2.0, label="Protein Rg", color="#1E88E5")
        plt.plot(time_ns, rg_ligand,  lw=2.0, label="Ligand Rg",  color="#D32F2F")
        plt.plot(time_ns, rg_complex, lw=2.0, label="Complex Rg", color="#2E7D32")

        plt.xlabel("Time (ns)")
        plt.ylabel("Radius of Gyration (Å)")
        plt.title("Radius of Gyration — Protein, Ligand, and Complex")
        plt.legend(frameon=False)
        plt.grid(alpha=0.25)

        savefig(os.path.join(OUTPUT_DIR, "Radius_of_Gyration_Overlay.png"))

    checkpoint("chk_D4B.txt")
    _record_stage_manifest("D4B", _d4b_num, _d4b_plot, _d4b_params)
    print("✅ STEP D4B complete.")

else:
    print("⏭️ STEP D4B skipped (checkpoint found).")







# ============================================================
# STEP D5 — FULL SYSTEM RMSD
# ============================================================

# ------------------ STEP D5 replacement (Full System RMSD) ------------------
_d5_num, _d5_plot, _d5_params = _stage_d5_outputs()
if not _stage_cache_ready("D5", "chk_D5.txt", _d5_num, _d5_plot, _d5_params):
    phase("STEP D5 — Full Complex RMSD")

    ref, bb, ligand_atoms, lig_resname = get_traj_metadata()

    time_ps, full_rmsd = chunked_rmsd(
        np.arange(ref.n_atoms, dtype=int),
        label="Full complex RMSD"
    )
    pd.DataFrame({"Time_ps": time_ps, "RMSD_Å": full_rmsd}).to_csv(
        os.path.join(OUTPUT_DIR, "FullComplex_RMSD.csv"), index=False
    )

    plt.figure(figsize=(9,5))
    plt.plot(time_ps/1000., full_rmsd, color=MIAMI_ORANGE, lw=1.6)  # full-system (protein-majority) in orange
    plt.title("Backbone RMSD of ALL Chains in Complex")
    plt.ylabel("RMSD (Å)")
    plt.xlabel("Time (ns)")

    savefig(os.path.join(OUTPUT_DIR, "FullComplex_RMSD.png"))

    checkpoint("chk_D5.txt")
    _record_stage_manifest("D5", _d5_num, _d5_plot, _d5_params)
    print("✅ STEP D5 complete.")
else:
    print("⏭️ STEP D5 skipped (checkpoint found).")


# ============================================================
# STEP D6 — Secondary Structure Evolution (DSSP)
# ============================================================

_d6_num, _d6_plot, _d6_params = _stage_d6_outputs()
if not _stage_cache_ready("D6", "chk_D6.txt", _d6_num, _d6_plot, _d6_params):
    phase("STEP D6 — Secondary Structure Evolution (DSSP)")

    print(f"📐 Computing DSSP in chunks of {CHUNK_SIZE} frames...")
    int_dssp, time_ns = compute_dssp_chunked(write_outputs=True)

    # ============================================================
    # 2) HEATMAP — SSE Timeline (Residue × Time)
    # ============================================================

    from matplotlib.colors import ListedColormap, BoundaryNorm

    print("🎨 Plotting DSSP heatmap...")

    # Coil = yellow, Sheet = blue, Helix = red
    cmap = ListedColormap(["#FFD700", "#1E88E5", "#D32F2F"])
    bounds = [-0.5, 0.5, 1.5, 2.5]          # 0,1,2 centered
    norm = BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(figsize=(16, 6))

    heatmap_dssp, heatmap_stride = downsample_frames_for_heatmap(int_dssp)
    if heatmap_stride > 1:
        print(f"ℹ️ DSSP heatmap downsampled every {heatmap_stride} frames for plotting; CSV outputs keep all frames.")

    im = ax.imshow(
        heatmap_dssp.T,
        aspect="auto",
        cmap=cmap,
        norm=norm,
        origin="upper"  # keep residue 0 at the top
    )

    ax.set_xlabel("Frame")
    ax.set_ylabel("Residue Index")
    ax.set_title("Secondary Structure Evolution (Coil / Sheet / Helix)")

    # Colorbar with fixed ticks/labels
    cbar = fig.colorbar(im, ax=ax, boundaries=bounds, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(["Coil", "Sheet", "Helix"])
    cbar.set_label("SSE", rotation=270, labelpad=15)

    fig.tight_layout()
    savefig(os.path.join(OUTPUT_DIR, "DSSP_Heatmap.png"))

    # ============================================================
    # 3) FRACTIONS — Helix / Sheet / Coil Over Time
    # ============================================================

    print("📊 Computing SSE fractions per frame...")

    helix_frac = (int_dssp == 2).mean(axis=1)
    sheet_frac = (int_dssp == 1).mean(axis=1)
    coil_frac  = (int_dssp == 0).mean(axis=1)

    df_sse = pd.DataFrame({
        "Time_ns": time_ns,
        "Helix": helix_frac,
        "Sheet": sheet_frac,
        "Coil":  coil_frac,
    })

    df_sse.to_csv(os.path.join(OUTPUT_DIR, "SSE_Fractions.csv"), index=False)

    print("🎨 Plotting SSE fraction timeline...")
    plt.figure(figsize=(12, 6))
    plt.plot(df_sse["Time_ns"], df_sse["Helix"], label="Helix", lw=2, color="#D32F2F")
    plt.plot(df_sse["Time_ns"], df_sse["Sheet"], label="Sheet", lw=2, color="#1E88E5")
    plt.plot(df_sse["Time_ns"], df_sse["Coil"],  label="Coil",  lw=2, color="#FFD700")

    plt.xlabel("Time (ns)")
    plt.ylabel("Fraction of Protein")
    plt.title("Secondary Structure Fractions Over Time")
    plt.legend()

    savefig(os.path.join(OUTPUT_DIR, "SSE_Fractions.png"))

    checkpoint("chk_D6.txt")
    _record_stage_manifest("D6", _d6_num, _d6_plot, _d6_params)
    print("✅ STEP D6 complete — DSSP analysis finished.\n")

else:
    print("⏭️ STEP D6 skipped (checkpoint found).")



# ============================================================
# STEP D6B — Per-Chain Secondary Structure Evolution (DSSP)
# ============================================================

import string

def chain_index_to_letter(idx):
    letters = string.ascii_uppercase
    if idx < len(letters):
        return letters[idx]
    else:
        return f"{letters[idx // 26 - 1]}{letters[idx % 26]}"



_d6b_num, _d6b_plot, _d6b_params = _stage_d6b_outputs()
if not _stage_cache_ready("D6B", "chk_D6B.txt", _d6b_num, _d6b_plot, _d6b_params):
    phase("STEP D6B — Per-Chain Secondary Structure Evolution")

    ref, bb, ligand_atoms, lig_resname = get_traj_metadata()

    print("📐 Loading DSSP numeric data (reusing D6 results if available)...")
    int_dssp_all, time_ns = get_dssp_numeric_and_time()

    # Map each residue to its chain index
    residue_chains = dssp_residue_chain_indices(ref.topology, int_dssp_all.shape[1])
    unique_chains = np.unique(residue_chains)

    from matplotlib.colors import ListedColormap, BoundaryNorm

    # Shared color mapping
    cmap = ListedColormap(["#FFD700", "#1E88E5", "#D32F2F"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = BoundaryNorm(bounds, cmap.N)

    for chain_id in unique_chains:
        chain_label = chain_index_to_letter(chain_id)
        print(f"\n🔍 Processing chain {chain_label} ...")


        mask = residue_chains == chain_id
        int_sse = int_dssp_all[:, mask]          # shape: (frames, residues_in_chain)

        # Save numeric matrix
        pd.DataFrame(int_sse).to_csv(
            os.path.join(OUTPUT_DIR, f"DSSP_chain_{chain_label}.csv"),
            index=False
        )


        # ---------------------------
        # Heatmap per chain
        # ---------------------------
        fig, ax = plt.subplots(figsize=(16, 6))
        heatmap_sse, heatmap_stride = downsample_frames_for_heatmap(int_sse)
        if heatmap_stride > 1:
            print(f"ℹ️ Chain {chain_label} heatmap downsampled every {heatmap_stride} frames for plotting.")
        im = ax.imshow(
            heatmap_sse.T,
            aspect="auto",
            cmap=cmap,
            norm=norm,
            origin="upper"
        )

        ax.set_ylabel("Residue Index (Chain-local)")
        ax.set_xlabel("Frame")
        ax.set_title(f"Chain {chain_label} — Secondary Structure Evolution")


        cbar = fig.colorbar(im, ax=ax, boundaries=bounds, ticks=[0, 1, 2])
        cbar.ax.set_yticklabels(["Coil", "Sheet", "Helix"])
        cbar.set_label("SSE", rotation=270, labelpad=15)

        fig.tight_layout()
        savefig(os.path.join(
            OUTPUT_DIR, f"DSSP_chain_{chain_label}_heatmap.png"
        ))


        # ---------------------------
        # Chain-level fractions
        # ---------------------------
        helix_frac = (int_sse == 2).mean(axis=1)
        sheet_frac = (int_sse == 1).mean(axis=1)
        coil_frac = (int_sse == 0).mean(axis=1)

        df = pd.DataFrame({
            "Time_ns": time_ns,
            "Helix": helix_frac,
            "Sheet": sheet_frac,
            "Coil": coil_frac
        })

        df.to_csv(
            os.path.join(OUTPUT_DIR, f"SSE_chain_{chain_label}_fractions.csv"),
            index=False
        )


        plt.figure(figsize=(10,5))
        plt.plot(df["Time_ns"], df["Helix"], label="Helix", color="#D32F2F")
        plt.plot(df["Time_ns"], df["Sheet"], label="Sheet", color="#1E88E5")
        plt.plot(df["Time_ns"], df["Coil"],  label="Coil",  color="#FFD700")
        plt.title(f"Chain {chain_label} — Secondary Structure Fractions Over Time")
        plt.xlabel("Time (ns)")
        plt.ylabel("Fraction")
        plt.legend()

        savefig(os.path.join(
            OUTPUT_DIR, f"SSE_chain_{chain_label}_fractions.png"
        ))


    checkpoint("chk_D6B.txt")
    _record_stage_manifest("D6B", _d6b_num, _d6b_plot, _d6b_params)
    print("✅ STEP D6B complete.\n")

else:
    print("⏭️ STEP D6B skipped (checkpoint found).")










# ============================================================
# STEP D6C — Residue-Level SSE Persistence Ranking
# ============================================================

_d6c_num, _d6c_plot, _d6c_params = _stage_d6c_outputs()
if not _stage_cache_ready("D6C", "chk_D6C.txt", _d6c_num, _d6c_plot, _d6c_params):
    phase("STEP D6C — Residue-Level SSE Persistence Ranking")

    int_dssp, time_ns = get_dssp_numeric_and_time()

    helix_persist = (int_dssp == 2).mean(axis=0)
    sheet_persist = (int_dssp == 1).mean(axis=0)
    coil_persist = (int_dssp == 0).mean(axis=0)

    df = pd.DataFrame({
        "Residue_Index": np.arange(int_dssp.shape[1]),
        "Helix_Persistence": helix_persist,
        "Sheet_Persistence": sheet_persist,
        "Coil_Persistence": coil_persist
    })

    df.to_csv(os.path.join(OUTPUT_DIR, "Residue_SSE_Persistence.csv"), index=False)

    # Plot ranking (Helix)
    plt.figure(figsize=(14,6))
    plt.bar(df["Residue_Index"], df["Helix_Persistence"], color="red")
    plt.title("Residue Helix Persistence")
    plt.xlabel("Residue Index")
    plt.ylabel("Fraction of Time in Helix")
    savefig(os.path.join(OUTPUT_DIR, "Residue_Helix_Persistence.png"))

    # Plot ranking (Sheet)
    plt.figure(figsize=(14,6))
    plt.bar(df["Residue_Index"], df["Sheet_Persistence"], color="blue")
    plt.title(" Residue Sheet Persistence")
    plt.xlabel("Residue Index")
    plt.ylabel("Fraction of Time in Sheet")
    savefig(os.path.join(OUTPUT_DIR, "Residue_Sheet_Persistence.png"))

    # Flexibility = 1 - max(structured)
    df["Flexibility"] = 1 - df[[
        "Helix_Persistence",
        "Sheet_Persistence"
    ]].max(axis=1)

    plt.figure(figsize=(14,6))
    plt.bar(df["Residue_Index"], df["Flexibility"], color="green")
    plt.title(" Residue Flexibility Index (1 = Highly Flexible)")
    plt.xlabel("Residue Index")
    plt.ylabel("Flexibility")
    savefig(os.path.join(OUTPUT_DIR, "Residue_Flexibility.png"))

    checkpoint("chk_D6C.txt")
    _record_stage_manifest("D6C", _d6c_num, _d6c_plot, _d6c_params)
    print("✅ STEP D6C complete.\n")

else:
    print("⏭️ STEP D6C skipped (checkpoint found).")







phase("STEP E-INTERFACE — Protein-Protein / Protein-Peptide Interface RIN Analysis")
if INTERFACE_RIN_MODE and not args.skip_interface_rin:
    interface_outdir = os.path.join(args.outdir, "Interface_RIN")
    _e_num = [
        os.path.join(interface_outdir, "interface_chain_pair_summary.csv"),
        os.path.join(interface_outdir, "interface_contact_timeseries.csv"),
        os.path.join(interface_outdir, "interface_residue_pair_contacts.csv"),
        os.path.join(interface_outdir, "interface_node_summary.csv"),
        os.path.join(interface_outdir, "interface_frame_index.csv"),
        os.path.join(interface_outdir, "interface_residue_pair_framewise.csv.gz"),
        os.path.join(interface_outdir, "interface_residue_pair_presence.npz"),
    ]
    _e_plot = [
        os.path.join(interface_outdir, "Interface_ChainPair_ContactFractions.png"),
        os.path.join(interface_outdir, "Interface_ContactCount_Timeseries.png"),
        os.path.join(interface_outdir, "Interface_MinDistance_Timeseries.png"),
        os.path.join(interface_outdir, "Interface_ResiduePair_ContactHeatmap.png"),
        os.path.join(interface_outdir, "Interface_RIN_Network.png"),
        os.path.join(interface_outdir, "Interface_ChainLevel_Network.png"),
        os.path.join(interface_outdir, "Interface_TopResiduePairs_Barplot.png"),
    ]
    _e_params = {
        "contact_cutoff_A": float(args.interface_contact_cutoff),
        "frame_step": int(args.interface_frame_step),
        "selected_chain_pairs": _interface_chain_pairs_for_manifest(),
    }
    if _stage_cache_ready("EINTERFACE", "chk_EINTERFACE.txt", _e_num, _e_plot, _e_params):
        print("⏭️ STEP E-INTERFACE skipped (checkpoint found and reusable data cache is valid).")
    else:
        if interface_rin_outputs_present(interface_outdir):
            print("⚠️ chk_EINTERFACE.txt exists but interface outputs are missing; rerunning 3C_Interface_RIN.py.")
            dispatch_to_interface_rin()
        else:
            if args.headless:
                dispatch_to_interface_rin()
            else:
                resp = input(
                    f"Detected a ligandless multi-chain system with {INTERFACE_CHAIN_COUNT} chains. "
                    "Run protein-protein / protein-peptide interface RIN analysis? [Y/n]: "
                ).strip().lower()
                if resp in {"", "y", "yes"}:
                    dispatch_to_interface_rin()
                else:
                    print("⏭️ STEP E-INTERFACE skipped by user choice.")
elif INTERFACE_RIN_MODE and args.skip_interface_rin:
    print("⏭️ STEP E-INTERFACE skipped due to --skip-interface-rin.")
elif not PARTNER_MODE:
    print("⏭️ Skipping ligand-contact/NETWORX analysis because this is apo single-chain or protein-centric mode.")
else:
    # ============================================================
    # STEP E — FULL BINDING-POCKET ANALYSIS WITH CSV EXPORTS
    # ============================================================
    _, _, _, lig_resname = get_traj_metadata()
    chain_label = ", ".join(active_chains.keys()) if active_chains else "Protein"
    phase("STEP E — Binding Pocket Analysis")

    u = mda.Universe("binding_pocket_only.pdb", "binding_pocket_only.xtc")
    pocket_ligand  = u.select_atoms(PARTNER_SELECTION_MDA or f"resname {lig_resname}")
    if PEPTIDE_LIKE and pocket_ligand.n_atoms > 0:
        pocket_protein = u.select_atoms(f"protein and not ({PARTNER_SELECTION_MDA})")
    else:
        pocket_protein = u.select_atoms("protein")
    total_frames = len(u.trajectory)

    CONTACT_CUTOFF = float(args.contact_cutoff)
    MIN_CONTACT_FRAC = float(args.min_contact_frac)

    # Storage tables
    interaction_data = []          # (Time_ps, Residue)
    framewise_interaction_rows = []  # (Time_ns, Residue, Type)

    # ============================================================
    # STEP E1 — CONTACT DETECTION (RESIDUE LEVEL)
    # ============================================================

    phase("STEP E1 — Pocket Contact Detection")

    interaction_data = []

    try:
        if pocket_ligand.n_atoms == 0:
            print(f"⚠️ No ligand atoms found for resname {lig_resname}. Skipping STEP E1/E1A/E1B.")
        elif pocket_protein.n_atoms == 0:
            print("⚠️ No protein atoms found in binding pocket selection. Skipping STEP E1/E1A/E1B.")
        else:
            for ts in tqdm(u.trajectory, desc="🔍 Detecting contacts"):
                time_ps = ts.time

                # IMPORTANT:
                # In peptide mode, residue atom indices are universe-global indices,
                # while a distance matrix built against pocket_protein.positions is
                # pocket-local. Using global indices against that matrix causes the
                # exact out-of-bounds failure seen for index 3197 with size 3197.
                #
                # To keep ligand mode and peptide mode both safe, compute distances
                # residue-by-residue against res.atoms.positions instead of slicing a
                # pocket-wide protein distance matrix.
                ligand_positions = pocket_ligand.positions

                for res in pocket_protein.residues:
                    res_atoms = res.atoms
                    if res_atoms.n_atoms == 0:
                        continue

                    dists = distances.distance_array(
                        ligand_positions,
                        res_atoms.positions
                    )

                    if dists.size == 0:
                        continue

                    ligand_atom_idx, protein_atom_idx = np.unravel_index(
                        np.argmin(dists),
                        dists.shape
                    )

                    min_dist = dists[ligand_atom_idx, protein_atom_idx]

                    if min_dist < CONTACT_CUTOFF:
                        ligand_atom = pocket_ligand.atoms[ligand_atom_idx]
                        protein_atom = res_atoms[protein_atom_idx]

                        interaction_data.append((
                            time_ps,
                            f"{res.resname}{res.resid}",
                            ligand_atom.name,
                            ligand_atom.index,
                            protein_atom.name,
                            protein_atom.index,
                            float(min_dist)
                        ))

    except Exception as e:
        print(f"⚠️ STEP E1 contact detection failed — continuing pipeline.")
        print(f"   ↳ {type(e).__name__}: {e}")


    # ============================================================
    # DATAFRAME CONSTRUCTION
    # ============================================================

    if len(interaction_data) == 0:
        print("⚠️ No ligand–protein contacts detected under current cutoff.")
        interaction_df = pd.DataFrame(columns=[
            "Time_ps", "Residue",
            "LigAtomName", "LigAtomIndex",
            "ProtAtomName", "ProtAtomIndex",
            "Distance", "Time_ns", "ResidNum"
        ])
    else:
        interaction_df = pd.DataFrame(
            interaction_data,
            columns=[
                "Time_ps", "Residue",
                "LigAtomName", "LigAtomIndex",
                "ProtAtomName", "ProtAtomIndex",
                "Distance"
            ]
        )

        interaction_df["Time_ns"] = interaction_df["Time_ps"] / 1000.0

        # Extract numeric residue number for sorting
        interaction_df["ResidNum"] = (
            interaction_df["Residue"]
            .str.extract(r"(\d+)")
            .astype(int)
        )

    # Save full framewise contacts whether empty or not
    interaction_df.to_csv(
        os.path.join(OUTPUT_DIR, "AllContacts_Framewise.csv"),
        index=False
    )

    print(f"📊 interaction_df rows: {len(interaction_df)}")
    print(f"📊 unique contacting residues: {interaction_df['Residue'].nunique() if not interaction_df.empty else 0}")


    # ============================================================
    # FILTER BY CONTACT PERSISTENCE
    # ============================================================

    if interaction_df.empty:
        keep_res = []
        filtered_df = interaction_df.copy()
    else:
        counts = interaction_df["Residue"].value_counts()
        keep_res = counts[counts >= total_frames * MIN_CONTACT_FRAC].index
        filtered_df = interaction_df[interaction_df["Residue"].isin(keep_res)].copy()

    filtered_df.to_csv(
        os.path.join(OUTPUT_DIR, "FilteredContacts_Framewise.csv"),
        index=False
    )

    print(f"📊 residues passing persistence filter ({MIN_CONTACT_FRAC:.1%}): {len(keep_res)}")
    print(f"📊 filtered_df rows: {len(filtered_df)}")


    # ============================================================
    # STEP E1A — CONTACT HEATMAP (NUMERICALLY SORTED)
    # ============================================================

    phase("STEP E1A — Contact Heatmap")

    try:
        if filtered_df.empty:
            print("⚠️ No persistent contacts found — skipping STEP E1A.")
        else:
            pivot = filtered_df.pivot_table(
                index=["ResidNum", "Residue"],
                columns="Time_ns",
                aggfunc="size",
                fill_value=0
            )

            if pivot.empty or pivot.shape[0] == 0 or pivot.shape[1] == 0:
                print("⚠️ Empty pivot table — skipping STEP E1A.")
            else:
                # Sort residues numerically
                pivot = pivot.sort_index(level="ResidNum")

                # Drop numeric index level (keep clean labels)
                pivot.index = pivot.index.get_level_values("Residue")

                plt.figure(figsize=(16, 9))
                sns.heatmap(pivot, cmap="coolwarm")
                plt.title(f"{compound_name}–{chain_label} Contact Map (Binding Pocket Only)")
                plt.xlabel("Time (ns)")
                plt.ylabel("Residue")

                savefig(
                    os.path.join(
                        OUTPUT_DIR,
                        f"{lig_resname}_contact_map_residue_time.png"
                    )
                )

                # Optional: save pivot used for plotting
                pivot.to_csv(
                    os.path.join(
                        OUTPUT_DIR,
                        f"{lig_resname}_contact_map_residue_time.csv"
                    )
                )

                print("✅ STEP E1A complete.")

    except Exception as e:
        print("⚠️ STEP E1A failed — continuing pipeline.")
        print(f"   ↳ {type(e).__name__}: {e}")


    # ============================================================
    # STEP E1B — CONTACT FREQUENCY BARPLOT (MATCHING ORDER)
    # ============================================================

    phase("STEP E1B — Contact Frequency Barplot")

    try:
        if filtered_df.empty:
            print("⚠️ No persistent contacts found — skipping STEP E1B.")
        else:
            freq = (
                filtered_df
                .drop_duplicates(["Residue", "ResidNum", "Time_ns"])
                .groupby(["ResidNum", "Residue"])
                .size()
                .sort_index(level="ResidNum")
            )

            if freq.empty:
                print("⚠️ Frequency table is empty — skipping STEP E1B.")
            else:
                # Clean index for plotting
                freq.index = freq.index.get_level_values("Residue")

                # Convert counts → percent of total frames
                freq_pct = (freq.astype(float) / float(total_frames)) * 100.0

                # Save numeric table
                freq_pct.to_csv(
                    os.path.join(
                        OUTPUT_DIR,
                        f"{lig_resname}_contact_frequency_filtered.csv"
                    ),
                    header=["PercentFrames"]
                )

                plt.figure(figsize=(8, 12))
                ax = freq_pct.plot(kind="barh", color=MIAMI_ORANGE)
                ax.set_title(f"{compound_name}-{chain_label} Contact Frequency")
                ax.set_xlabel("Percent of frames (%)")
                ax.set_xlim(0, 100)

                # Annotate bars
                for p in ax.patches:
                    width = p.get_width()
                    if np.isfinite(width):
                        ax.text(
                            width + 0.8,
                            p.get_y() + p.get_height() / 2,
                            f"{width:.1f}%",
                            va="center",
                            fontsize=8
                        )

                savefig(
                    os.path.join(
                        OUTPUT_DIR,
                        f"{lig_resname}_contact_frequency_filtered.png"
                    )
                )

                print("✅ STEP E1B complete.")

    except Exception as e:
        print("⚠️ STEP E1B failed — continuing pipeline.")
        print(f"   ↳ {type(e).__name__}: {e}")

    # ============================================================
    # STEP E2 — INTERACTION TYPE CLASSIFICATION (SAFE SINGLE THREAD)
    # ============================================================

    phase("STEP E2 — Classifying Interaction Types")

    def classify_frame(ts):
        time_ns = ts.time / 1000.0
        ligand_positions = pocket_ligand.positions
        rows = []

        for res in pocket_protein.residues:
            rid = f"{res.resname}{res.resid}"
            res_atoms = res.atoms
            if res_atoms.n_atoms == 0:
                continue

            # Same index-safety fix as STEP E1: use residue-local coordinates
            # instead of slicing a pocket-wide matrix with universe-global indices.
            dists = distances.distance_array(
                ligand_positions,
                res_atoms.positions
            )
            if dists.size == 0:
                continue

            lig_i, prot_j = np.unravel_index(np.argmin(dists), dists.shape)
            min_dist = dists[lig_i, prot_j]

            lig_atom = pocket_ligand.atoms[lig_i]
            prot_atom = res_atoms[prot_j]

            itype = None

            # H-bond detection
            if infer_element(prot_atom.name) in ("O","N") and infer_element(lig_atom.name) in ("O","N"):
                if min_dist < HBOND_DISTANCE_CUTOFF:
                    itype = "H-bond"

            # Hydrophobic
            if itype is None and res.resname in {"PHE","LEU","ILE","VAL","MET","TRP","TYR","ALA"} and min_dist < 4.5:
                itype = "Hydrophobic"

            # Ionic
            if itype is None and res.resname in {"ASP","GLU","ARG","LYS","HIS"} and min_dist < 6.0:
                itype = "Ionic"

            if itype:
                rows.append((
                    time_ns, rid, itype,
                    lig_atom.name, lig_atom.index,
                    prot_atom.name, prot_atom.index,
                    float(min_dist)
                ))

        return rows



    # Execute (safe, interruptible)
    for ts in tqdm(u.trajectory, desc="🔬 Classifying interactions"):
        framewise_interaction_rows.extend(classify_frame(ts))


    # Convert & save
    int_df = pd.DataFrame(
        framewise_interaction_rows,
        columns=[
            "Time_ns", "Residue", "Type",
            "LigAtomName", "LigAtomIndex",
            "ProtAtomName", "ProtAtomIndex",
            "Distance"
        ]
    )

    int_df.to_csv(os.path.join(OUTPUT_DIR, "InteractionTypes_Framewise.csv"), index=False)

    if int_df.empty:
        print("⚠️ No ligand:protein interactions were classified in STEP E2. Creating an empty interaction summary and skipping data-dependent interaction plots.")
        summary = ensure_interaction_summary_columns(pd.DataFrame(), PLOT_TYPES)
    else:
        summary = int_df.groupby(["Residue","Type"]).size().unstack(fill_value=0)
        summary = ensure_interaction_summary_columns(summary, PLOT_TYPES)

    print("📊 Interaction summary columns:", list(summary.columns))
    print("📊 Interaction type counts:", int_df["Type"].value_counts().to_dict() if not int_df.empty else {})
    summary.to_csv(os.path.join(OUTPUT_DIR, "InteractionTypes_Summary.csv"))

    stacked_frac = summary[PLOT_TYPES].div(total_frames).fillna(0)

    import re

    def residue_sort_key(res):
        # Extract trailing residue number (e.g. ARG329 → 329)
        m = re.search(r"(\d+)$", res)
        return int(m.group(1)) if m else 0

    # Sort stacked_frac by residue number
    stacked_frac = stacked_frac.loc[
        sorted(stacked_frac.index, key=residue_sort_key)
    ]

    # ============================================================
    # STEP E2A — Grouped bar plot (comparison view)
    # ============================================================

    if stacked_frac.empty:
        print("⚠️ No interaction summary rows available for STEP E2A grouped bar plot — skipping figure.")
    else:
        fig, ax = plt.subplots(figsize=(16, 6), dpi=150)

        stacked_frac.plot(
            kind="bar",
            ax=ax,
            stacked=False,
            width=0.75,
            color=[COLORS[c] for c in PLOT_TYPES],
            edgecolor="black",
            linewidth=0.4,
        )

        ax.set_title(
            f"Normalized {chain_label}–{compound_name} Interaction Frequencies",
            fontsize=16,
            pad=20
        )
        ax.set_ylabel("Fraction of Frames", fontsize=12)
        ax.set_xlabel("Residue", fontsize=12)

        # Comparison emphasis
        ax.yaxis.grid(True, linestyle="--", alpha=0.4)
        ax.set_axisbelow(True)

        ax.tick_params(axis="x", rotation=90, labelsize=8)
        ax.tick_params(axis="y", labelsize=10)

        # Legend ABOVE axes, right side of title (figure-level)
        fig.legend(
            handles=ax.get_legend_handles_labels()[0],
            labels=ax.get_legend_handles_labels()[1],
            title="Interaction Type",
            fontsize=11,
            title_fontsize=12,
            ncol=len(PLOT_TYPES),
            loc="upper right",
            bbox_to_anchor=(0.98, 1.01),
            frameon=False
        )

        # Remove axis legend
        ax.get_legend().remove()

        ax.margins(x=0.01)
        plt.tight_layout(rect=[0, 0, 1, 0.88])

        savefig(os.path.join(OUTPUT_DIR, "interaction_grouped_normalized.png"))
        plt.close()


    # ============================================================
    # STEP E2B — Stacked bar plot (composition view)
    # ============================================================

    if stacked_frac.empty:
        print("⚠️ No interaction summary rows available for STEP E2B stacked bar plot — skipping figure.")
    else:
        fig, ax = plt.subplots(figsize=(16, 6), dpi=150)

        stacked_frac.plot(
            kind="bar",
            ax=ax,
            stacked=True,
            width=0.75,
            color=[COLORS[c] for c in PLOT_TYPES],
            edgecolor="white",
            linewidth=0.6,
        )

        ax.set_title(
            f"Normalized {chain_label}–{compound_name} Interaction Composition",
            fontsize=16,
            pad=20
        )
        ax.set_ylabel("Fraction of Frames", fontsize=12)
        ax.set_xlabel("Residue", fontsize=12)

        # Composition emphasis
        ax.yaxis.grid(False)

        ax.tick_params(axis="x", rotation=90, labelsize=8)
        ax.tick_params(axis="y", labelsize=10)

        # Legend ABOVE axes, right side of title (figure-level)
        fig.legend(
            handles=ax.get_legend_handles_labels()[0],
            labels=ax.get_legend_handles_labels()[1],
            title="Interaction Type",
            fontsize=11,
            title_fontsize=12,
            ncol=len(PLOT_TYPES),
            loc="upper right",
            bbox_to_anchor=(0.98, 1.01),
            frameon=False
        )

        # Remove axis legend
        ax.get_legend().remove()

        ax.margins(x=0.01)
        plt.tight_layout(rect=[0, 0, 1, 0.88])

        savefig(os.path.join(OUTPUT_DIR, "interaction_stacked_normalized.png"))
        plt.close()



    # ============================================================
    # STEP E2.5 — Call NETWORX (NEWNETWORX_noLP.py)
    # ============================================================
    phase("STEP E3 —Running Advanced NETWORX Renderer")

    NETWORX_SCRIPT = "3B_NETWORX.py"   # or 3B_NETWORX.py if preferred

    if not os.path.exists(NETWORX_SCRIPT):
        print(f"❌ NETWORX script {NETWORX_SCRIPT} not found — skipping NETWORX rendering.")
    else:
        output_name = args.net_out or f"{lig_resname}"
        minfrac = args.minfrac           # same threshold as analyzer
        ellipse_rx = args.ellipse_rx
        ellipse_ry = args.ellipse_ry
    
        # Build CLI call string
        cmd = f"python {NETWORX_SCRIPT} -l {lig_resname} -f {minfrac} -o {output_name}"
    
        if ellipse_rx:
            cmd += f" --ellipse-rx {ellipse_rx}"
        if ellipse_ry:
            cmd += f" --ellipse-ry {ellipse_ry}"
        if args.no_edges:
            cmd += " --no-edges"

        print("➡️ Running NETWORX with command:")
        print("   ", cmd)

        run_cmd(cmd)
        print("✅ NETWORX rendering complete.")


    # ============================================================
    # STEP E3 — RUN ADVANCED NETWORX RENDERER
    # ============================================================
    phase("STEP E3 — Running Advanced NETWORX Renderer")

    NETWORX_SCRIPT = "3B_NETWORX.py"

    if not os.path.exists(NETWORX_SCRIPT):
        print(f"❌ NETWORX script {NETWORX_SCRIPT} not found — skipping NETWORX rendering.")
    else:
        # Protein-only jobs should not try to render ligand/peptide NETWORX
        if not PARTNER_MODE:
            print("⏭️ Skipping NETWORX renderer for protein-only mode.")
        else:
            output_name = args.net_out or f"{lig_resname or ligand_code or compound_name or 'NETWORX'}"
            minfrac = args.minfrac if args.minfrac is not None else MIN_CONTACT_FRAC

            # Decide which NETWORX mode to send
            if PEPTIDE_LIKE:
                net_mode = "peptide"
                net_label = (ligand_code or "PEPTIDE").strip()
            else:
                net_mode = "ligand"
                net_label = (lig_resname or ligand_code or compound_name or "LIG").strip()

            cmd_parts = [
                "python",
                NETWORX_SCRIPT,
                "-l", str(net_label),
                "-f", str(minfrac),
                "-o", str(output_name),
                "--mode", net_mode,
            ]

            if args.ellipse_rx is not None:
                cmd_parts.extend(["--ellipse-rx", str(args.ellipse_rx)])

            if args.ellipse_ry is not None:
                cmd_parts.extend(["--ellipse-ry", str(args.ellipse_ry)])

            if args.no_edges:
                cmd_parts.append("--no-edges")

            cmd = " ".join(cmd_parts)

            print("➡️ Running NETWORX with command:")
            print("   ", cmd)

            run_cmd(cmd)
            print("✅ NETWORX rendering complete.")



    # ============================================================
    # STEP E4 — INTERACTION TIMELINE HEATMAP (Priority Encoded)
    # ============================================================

    phase("STEP E4 — Interaction Timeline Heatmap")

    # -------------------------------
    # Priority encoding
    # -------------------------------
    PRIORITY_MAP = {
        "H-bond": 3,
        "Ionic": 2,
        "Hydrophobic": 1
    }

    if int_df.empty:
        print("⚠️ No interaction rows available for STEP E4 — skipping interaction timeline heatmap.")
    else:
        # Encode interaction strength numerically
        int_df["InteractionCode"] = int_df["Type"].map(PRIORITY_MAP).fillna(0)

        # ------------------------------------------------------------
        # OPTIONAL (HIGHLY RECOMMENDED): bin time to reduce sparsity
        # ------------------------------------------------------------
        TIME_BIN_NS = 2.0  # adjust if needed
        int_df["Time_bin"] = (
            (int_df["Time_ns"] / TIME_BIN_NS)
            .round(0)
            * TIME_BIN_NS
        )

        # ------------------------------------------------------------
        # Pivot using MAX priority per residue per time bin
        # ------------------------------------------------------------
        timeline_matrix = (
            int_df
            .pivot_table(
                index="Residue",
                columns="Time_bin",
                values="InteractionCode",
                aggfunc="max",
                fill_value=0
            )
            .sort_index(
                key=lambda x: x.map(residue_number),
                ascending=False   # 🔥 highest residue number on top
            )
        )

        # ------------------------------------------------------------
        # Plot
        # ------------------------------------------------------------
        plt.figure(figsize=(20, 10))

        sns.heatmap(
            timeline_matrix,
            cmap="viridis",
            cbar_kws={
                "label": "Interaction Strength (3=H-bond, 2=Ionic, 1=Hydrophobic)"
            },
            linewidths=0.05
        )

        plt.title(f"{compound_name} Interaction Timeline Heatmap")
        plt.xlabel("Time (ns)")
        plt.ylabel("Residue")
        plt.tight_layout()

        savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_interaction_timeline.png"))

        print("✅ STEP E4 complete — interaction timeline heatmap created.")


    # ============================================================
    # STEP E4B — INTERACTION EVENT TIMELINE (Figure 24)
    # ============================================================

    phase("STEP E4B — Interaction Event Timeline (Figure 25)")

    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    if int_df.empty:
        print("⚠️ No interaction rows available for STEP E4B — skipping interaction event timeline.")
    else:
        # Sort data
        df = int_df.sort_values(["Residue", "Time_ns"])

        # Residue ordering
        residues = sorted(
            df["Residue"].unique(),
            key=residue_number,
            reverse=False   # 🔥 highest number at top
        )

        res_to_y = {res: i for i, res in enumerate(residues)}

        fig, ax = plt.subplots(figsize=(22, 12))
        BAR_HEIGHT = 0.8

        # Color map
        TYPE_COLORS = {
            "H-bond": "#FDE725",
            "Ionic": "#5EC962",
            "Hydrophobic": "#31688E"
        }

        # Gap threshold (ns) to define event break
        GAP_NS = 1.5  # ns; maximum allowed gap between frames for a continuous event

        # Build all contiguous residue/type events before drawing.  Calling
        # ``broken_barh`` once per event creates tens of thousands of separate
        # Matplotlib collections for long trajectories and becomes very slow.
        from matplotlib.collections import PolyCollection

        event_start = (
            df["Residue"].ne(df["Residue"].shift())
            | df["Type"].ne(df["Type"].shift())
            | df["Time_ns"].sub(df["Time_ns"].shift()).gt(GAP_NS)
        )
        events = (
            df.assign(_event_id=event_start.cumsum())
            .groupby("_event_id", sort=False, as_index=False)
            .agg(
                Residue=("Residue", "first"),
                Type=("Type", "first"),
                Start=("Time_ns", "first"),
                End=("Time_ns", "last"),
            )
        )

        event_polygons = []
        event_colors = []
        half_height = BAR_HEIGHT / 2
        for residue, itype, start_time, end_time in events.itertuples(
            index=False, name=None
        ):
            y = res_to_y[residue]
            event_polygons.append([
                (start_time, y - half_height),
                (end_time, y - half_height),
                (end_time, y + half_height),
                (start_time, y + half_height),
            ])
            event_colors.append(TYPE_COLORS.get(itype, "gray"))

        ax.add_collection(
            PolyCollection(event_polygons, facecolors=event_colors, edgecolors="none")
        )
        print(f"  ✓ Rendering {len(events):,} interaction events as one collection.")

        # Axis formatting
        ax.set_yticks(range(len(residues)))
        ax.set_yticklabels(residues)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Residue")
        ax.set_title(f"{compound_name} Interaction Event Timeline")

        ax.set_ylim(-1, len(residues))
        ax.set_xlim(df["Time_ns"].min(), df["Time_ns"].max())

        # Legend
        legend_elements = [
            Patch(facecolor=c, label=t) for t, c in TYPE_COLORS.items()
        ]
        ax.legend(
            handles=legend_elements,
            title="Interaction Type",
            loc="upper right"
        )

        plt.tight_layout()
        savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_interaction_event_timeline.png"))

        print("✅ STEP E4B complete — Figure 25 interaction event timeline created.")



    # ============================================================
    # STEP E5 — Per-Residue Interaction Composition Pie Charts
    # ============================================================

    phase("STEP E5 — Residue Pie Charts")
    print("\n📊 Generating pie charts for per-residue interaction composition...")

    if summary.empty or summary["Total"].sum() == 0:
        print("⚠️ No interaction summary data available for STEP E5 — skipping residue pie charts.")
    else:
        top_residues = summary.sort_values("Total", ascending=False).head(20).index
        # top_residues = summary.sort_values("Total", ascending=False).index

        for res in top_residues:
            print(f"   • Creating pie chart for {res} ...")

            vals = summary.loc[res].reindex(PLOT_TYPES, fill_value=0)
            vals = vals[vals > 0]

            if len(vals) == 0:
                print(f"     ⚠️  Skipping {res} — no nonzero interactions.")
                continue

            plt.figure(figsize=(4,4))
            plt.pie(
                vals,
                labels=vals.index,
                autopct="%1.1f%%",
                startangle=90,
                pctdistance=0.75
            )
            plt.title(f"{res} Interaction Composition")
            plt.tight_layout()

            savefig(os.path.join(PIE_DIR, f"pie_{res}.png"))

    print("✅ STEP E5 complete — Pie charts saved.\n")


    # ============================================================
    # STEP E6 — Interaction Distance Distributions
    # ============================================================

    phase("STEP E6 — Interaction Distance Distributions")
    print("\n📈 Generating contact distance histograms per residue...")

    for res in filtered_df["Residue"].unique():
        d = filtered_df.loc[filtered_df["Residue"] == res, "Distance"]

        plt.figure(figsize=(5,4))
        sns.histplot(d, bins=30, kde=True)

        plt.xlabel("Distance (Å)")
        plt.ylabel("Count")
        plt.title(f"{res} Contact Distance Distribution")

        savefig(os.path.join(DIST_DIR, f"dist_{res}.png"))
        plt.close()

    print("✅ STEP E6 complete — Distance histograms saved.\n")



    # ============================================================
    # STEP E7 — Residue–Ligand Contact Persistence Distribution
    # ============================================================

    phase("STEP E7 — Interaction Degree Distribution")
    print("\n📊 Generating residue–ligand contact persistence distribution...")
    print("\n🔍 int_df columns:")
    print(int_df.columns.tolist())

    # --- Any-contact persistence ---
    contact_persistence = (
        int_df
        .groupby("Residue")["Time_ns"]
        .nunique()
    )

    contact_persistence = contact_persistence[
        contact_persistence >= total_frames * MIN_CONTACT_FRAC
    ]

    # --- Hydrogen-bond persistence ---
    hbond_df = int_df[int_df["Type"].str.lower().isin(["hbond", "hydrogen_bond"])]


    hbond_persistence = (
        hbond_df
        .groupby("Residue")["Time_ns"]
        .nunique()
    )

    # Apply same persistence threshold
    hbond_persistence = hbond_persistence[
        hbond_persistence >= total_frames * MIN_CONTACT_FRAC
    ]

    plt.figure(figsize=(6, 4))

    sns.histplot(
        contact_persistence,
        bins=15,
        color="steelblue",
        alpha=0.6,
        label="Any ligand contact"
    )

    sns.histplot(
        hbond_persistence,
        bins=15,
        color="darkorange",
        alpha=0.6,
        label="Hydrogen bond contact"
    )

    plt.xlabel("Number of Simulation Frames with Contact")
    plt.ylabel("Number of Residues")
    plt.title(f"{compound_name} Residue–Ligand Contact Persistence")

    plt.legend(frameon=False)
    savefig(os.path.join(DIST_DIR, "interaction_degree_distribution.png"))
    plt.close()

    print("✅ STEP E7 complete — Contact persistence distribution saved.\n")



    # ============================================================
    # STEP E7B — Composite Distance Histogram
    # ============================================================

    phase("STEP E7B — Composite Distance Histogram")
    print("\n📊 Generating composite histogram of ALL interaction distances...")

    if int_df.empty:
        print("⚠️ No interaction rows available for STEP E7B — skipping composite distance histogram.")
    else:
        plt.figure(figsize=(10, 7))

        for itype, color in zip(
            PLOT_TYPES,
            ["#4DB6AC", "#64B5F6", "#FFD54F"]
        ):
            d = int_df.loc[int_df["Type"] == itype, "Distance"]
            sns.histplot(d, bins=40, kde=True, label=itype, color=color, alpha=0.6)

        plt.xlabel("Interaction Distance (Å)")
        plt.ylabel("Count")
        plt.title(f"{compound_name} – Composite Interaction Distance Distributions")
        plt.legend()

        savefig(os.path.join(OUTPUT_DIR, "Composite_Distance_Distributions.png"))
        plt.close()

    print("✅ STEP E7B complete — Composite distance histogram saved.\n")


    # ============================================================
    # STEP E8 — Distance Violin Plot by Interaction Type
    # ============================================================

    phase("STEP E8 — Violin Plot by Interaction Type")
    print("\n🎻 Generating violin plot for interaction type distance distributions...")

    violin_df = int_df[int_df["Type"].isin(PLOT_TYPES)] if not int_df.empty else int_df
    if violin_df.empty:
        print("⚠️ No interaction rows available for STEP E8 — skipping violin plot.")
    else:
        plt.figure(figsize=(9,6))
        sns.violinplot(
            data=violin_df,
            x="Type",
            y="Distance",
            palette="Set2"
        )

        plt.title(f"{compound_name} — Distance Distribution by Interaction Type")
        plt.ylabel("Distance (Å)")

        savefig(os.path.join(OUTPUT_DIR, "interaction_type_violin.png"))
        plt.close()

    print("✅ STEP E8 complete — Violin plot saved.\n")

    # ============================================================
    # STEP E9 — Interaction Heatmap (Residue × Type)
    # ============================================================

    phase("STEP E9 — Interaction Heatmap")
    print("\n🔥 Generating interaction frequency heatmap...")

    # Normalize counts → fraction of frames
    heatmap_df = summary.reindex(columns=PLOT_TYPES, fill_value=0).div(total_frames)

    if heatmap_df.empty:
        print("⚠️ No interaction summary rows available for STEP E9 — skipping interaction heatmap.")
    else:
        # Order residues by overall interaction frequency (NOT a Total column)
        order = heatmap_df.sum(axis=1).sort_values(ascending=False).index
        heatmap_df = heatmap_df.loc[order]

        plt.figure(figsize=(10, 8))
        sns.heatmap(
            heatmap_df,
            cmap="magma",
            vmin=0,
            vmax=1,
            cbar_kws={
                "label": "Interaction Frequency\n(fraction of simulation frames)"
            }
        )

        plt.title(f"{compound_name} Interaction Type vs Residue Heatmap")
        plt.xlabel("Interaction Type")
        plt.ylabel("Residue")

        savefig(os.path.join(OUTPUT_DIR, "interaction_heatmap_normalized.png"))
        plt.close()

    print("✅ STEP E9 complete — Normalized interaction heatmap saved.\n")






    # ============================================================
    # STEP E10 — Residue Binding Importance Ranking
    # ============================================================

    phase("STEP E10 — Binding Importance Ranking")
    print("\n🏆 Computing and plotting binding importance ranking...")

    importance = (
        3 * summary.reindex(columns=PLOT_TYPES, fill_value=0)["H-bond"] +
        2 * summary.reindex(columns=PLOT_TYPES, fill_value=0)["Ionic"] +
        1 * summary.reindex(columns=PLOT_TYPES, fill_value=0)["Hydrophobic"]
    ).sort_values(ascending=False)

    if importance.empty:
        print("⚠️ No interaction summary rows available for STEP E10 — skipping binding importance chart.")
    else:
        plt.figure(figsize=(10,6))
        importance.plot(kind="bar", color="slateblue")
        plt.title("Residue Binding Importance Score")
        plt.ylabel("Weighted Interaction Score")
        plt.xlabel("Residue")

        savefig(os.path.join(OUTPUT_DIR, "binding_importance_ranking.png"))

    print("✅ STEP E10 complete — Binding importance chart saved.\n")




    # ============================================================
    # STEP E11 — PERSISTENCE VS DISTANCE SCATTER (FIXED)
    # ============================================================

    phase("STEP E11 — Persistence vs Distance Scatter")

    if int_df.empty or summary.empty:
        print("⚠️ No interaction rows available for STEP E11 — skipping persistence vs distance scatter.")
    else:
        # --- Compute per-residue mean interaction distance ---
        mean_dist_series = int_df.groupby("Residue")["Distance"].mean()

        # --- Compute persistence fraction for the same residues ---
        persistence_series = summary["Total"] / total_frames

        # --- Align indices so x and y match ---
        common_residues = mean_dist_series.index.intersection(persistence_series.index)

        mean_dist = mean_dist_series.loc[common_residues]
        frac = persistence_series.loc[common_residues]

        # --- Make sure lengths match ---
        print("Residues used:", len(common_residues))

        plt.figure(figsize=(10, 6))
        plt.scatter(mean_dist, frac, s=120, edgecolor="k")

        for r in common_residues:
            plt.text(mean_dist[r], frac[r], r, fontsize=8)

        plt.xlabel("Mean Minimum Distance (Å)")
        plt.ylabel("Persistence Fraction")
        plt.title(f"{compound_name} — Persistence vs Distance")

        savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_persistence_vs_distance.png"))

    print("✅ STEP E11 complete.")





    # ============================================================
    # STEP D6D — FULL SSE ↔ Ligand Contact Dynamics Analysis (SAFE)
    # Includes:
    #  - Stacked SSE + Contacts plots
    #  - Pearson & Spearman correlations
    #  - Rolling correlations aligned to trajectory
    #  - Memory-safe change-point detection
    #  - Segment statistics export
    # ============================================================

    if not PARTNER_MODE:
        print("⏭️ STEP D6D skipped (protein-centric mode).")

    elif not _stage_cache_ready("D6D", "chk_D6D.txt", *_stage_d6d_outputs()):
        phase("STEP D6D — FULL SSE–Ligand Contact Coupling Analysis")

        ref, bb, ligand_atoms, lig_resname = get_traj_metadata()

        import ruptures as rpt
        from scipy.stats import pearsonr, spearmanr

        # ------------------------------------------------------------
        # 1) LOAD DSSP FROM STEP D6 OUTPUT IF AVAILABLE
        #    (avoid recomputing DSSP on large trajectories)
        # ------------------------------------------------------------
        dssp_csv = os.path.join(OUTPUT_DIR, "DSSP_raw.csv")

        helix_codes = {"H", "G", "I"}
        sheet_codes = {"E", "B"}
        coil_codes  = {"S", "T", "C", " "}

        print("📐 Loading DSSP data for SSE/contact coupling analysis...")

        if os.path.exists(dssp_csv):
            dssp_df = pd.read_csv(dssp_csv)
            dssp = dssp_df.astype(str).fillna(" ").values
            print(f"✅ Reused DSSP from STEP D6: {dssp.shape[0]} frames × {dssp.shape[1]} residues")
        else:
            print("⚠️ DSSP_raw.csv not found; recomputing DSSP now...")
            compute_dssp_chunked(write_outputs=True)
            dssp_df = pd.read_csv(dssp_csv)
            dssp = dssp_df.astype(str).fillna(" ").values
            print(f"✅ DSSP recomputed and saved: {dssp.shape[0]} frames × {dssp.shape[1]} residues")

        helix_frac, sheet_frac, coil_frac = [], [], []

        for row in dssp:
            total = len(row) if len(row) > 0 else 1
            helix_frac.append(sum(c in helix_codes for c in row) / total)
            sheet_frac.append(sum(c in sheet_codes for c in row) / total)
            coil_frac.append(sum(c in coil_codes for c in row) / total)

        helix_frac = np.asarray(helix_frac, dtype=float)
        sheet_frac = np.asarray(sheet_frac, dtype=float)
        coil_frac  = np.asarray(coil_frac, dtype=float)
        if os.path.exists(os.path.join(OUTPUT_DIR, "SSE_Fractions.csv")):
            time_ns = pd.read_csv(os.path.join(OUTPUT_DIR, "SSE_Fractions.csv"))["Time_ns"].to_numpy(dtype=float)
        else:
            time_ns = chunked_times_ns()

        # ------------------------------------------------------------
        # 2) LOAD CONTACTS FROM STEP E2 / E1 OUTPUT
        # ------------------------------------------------------------
        print("🔗 Loading ligand-contact counts over time...")

        if "int_df" in globals() and isinstance(int_df, pd.DataFrame) and not int_df.empty:
            contact_source_df = int_df.copy()
        else:
            fallback_csvs = [
                os.path.join(OUTPUT_DIR, "All_Interaction_Types_Framewise.csv"),
                os.path.join(OUTPUT_DIR, "AllContacts_Framewise.csv"),
                os.path.join(OUTPUT_DIR, "FilteredContacts_Framewise.csv"),
            ]

            contact_source_df = None
            for csv_path in fallback_csvs:
                if os.path.exists(csv_path):
                    tmp = pd.read_csv(csv_path)
                    if not tmp.empty and ("Time_ns" in tmp.columns or "Time_ps" in tmp.columns):
                        contact_source_df = tmp
                        print(f"✅ Loaded contacts from: {os.path.basename(csv_path)}")
                        break

            if contact_source_df is None:
                raise RuntimeError(
                    "❌ No interaction dataframe available for STEP D6D. "
                    "Run STEP E1/E2 first so contact counts exist."
                )

        if "Time_ns" not in contact_source_df.columns:
            if "Time_ps" in contact_source_df.columns:
                contact_source_df["Time_ns"] = contact_source_df["Time_ps"] / 1000.0
            else:
                raise RuntimeError("❌ Contact dataframe has neither Time_ns nor Time_ps.")

        # Snap contact times to trajectory times to avoid float mismatch issues
        rounded_contact_time = np.round(contact_source_df["Time_ns"].astype(float).values, 6)
        rounded_traj_time    = np.round(time_ns, 6)

        contact_counts = (
            pd.Series(rounded_contact_time)
            .value_counts()
            .sort_index()
            .reindex(rounded_traj_time, fill_value=0)
        )

        contact_counts.index = time_ns

        # Smooth contacts for interpretability
        smooth_window = 51
        if len(contact_counts) < smooth_window:
            smooth_window = max(5, len(contact_counts) // 5 * 2 + 1)

        contacts_smooth = (
            contact_counts
            .rolling(window=smooth_window, center=True, min_periods=1)
            .mean()
        )

        # ------------------------------------------------------------
        # 3) CORRELATION ANALYSIS
        # ------------------------------------------------------------
        pear_corr, pear_p = pearsonr(contacts_smooth.values, sheet_frac)
        spear_corr, spear_p = spearmanr(contacts_smooth.values, sheet_frac)

        print("\n📈 SSE ↔ Contact Correlations")
        print(f"   Pearson  (contacts vs sheet): {pear_corr:.3f}  [p={pear_p:.3e}]")
        print(f"   Spearman (contacts vs sheet): {spear_corr:.3f}  [p={spear_p:.3e}]")

        # ------------------------------------------------------------
        # 4) ROLLING CORRELATION
        # ------------------------------------------------------------
        contacts_series = pd.Series(contacts_smooth.values, index=time_ns)
        sheet_series    = pd.Series(sheet_frac, index=time_ns)

        roll_window_frames = min(500, max(25, len(time_ns) // 100))
        if roll_window_frames % 2 == 0:
            roll_window_frames += 1

        rolling_corr = (
            contacts_series
            .rolling(roll_window_frames, center=True, min_periods=max(5, roll_window_frames // 5))
            .corr(sheet_series)
        )

        rolling_corr = rolling_corr.reindex(time_ns)
        rolling_corr = rolling_corr.bfill().ffill()

        # ------------------------------------------------------------
        # 5) MEMORY-SAFE CHANGE-POINT DETECTION
        # ------------------------------------------------------------
        print("\n🔍 Detecting change-points in ligand contacts (memory-safe mode)...")

        signal = contacts_smooth.fillna(0).astype(np.float32).values
        n_frames = len(signal)

        # Downsample aggressively for very long trajectories
        if n_frames > 40000:
            ds_factor = 20
        elif n_frames > 20000:
            ds_factor = 10
        elif n_frames > 10000:
            ds_factor = 5
        elif n_frames > 5000:
            ds_factor = 2
        else:
            ds_factor = 1

        if ds_factor > 1:
            signal_ds = signal[::ds_factor]
            time_ds = time_ns[::ds_factor]
        else:
            signal_ds = signal
            time_ds = time_ns

        # reshape for ruptures
        signal_ds_2d = signal_ds.reshape(-1, 1)

        # choose a sane minimum segment size
        min_size = max(10, len(signal_ds) // 100)

        try:
            algo = rpt.Pelt(model="l2", min_size=min_size, jump=1).fit(signal_ds_2d)
            change_points_ds = algo.predict(pen=3)
        except Exception as e:
            print(f"⚠️ Change-point detection failed: {type(e).__name__}: {e}")
            change_points_ds = [len(signal_ds)]

        # Remove final endpoint if present
        if change_points_ds and change_points_ds[-1] == len(signal_ds):
            change_points_ds = change_points_ds[:-1]

        # Map back to original frame indices
        change_points = sorted(set(int(cp * ds_factor) for cp in change_points_ds if cp * ds_factor < n_frames))

        print(f"   Downsample factor used: {ds_factor}")
        print(f"   Change-points detected at frames: {change_points}")

        # ------------------------------------------------------------
        # 6) CHANGE-POINT SUMMARY CSV
        # ------------------------------------------------------------
        cp_out = os.path.join(OUTPUT_DIR, f"{lig_resname}_ChangePoint_Summary.csv")
        with open(cp_out, "w") as f:
            f.write("Segment,StartFrame,EndFrame,StartTime_ns,EndTime_ns,MeanContacts,MeanHelix,MeanSheet,MeanCoil\n")

            boundaries = change_points + [n_frames]
            start = 0

            for i, cp in enumerate(boundaries):
                if cp <= start:
                    continue

                seg_contacts = float(np.mean(signal[start:cp]))
                seg_helix    = float(np.mean(helix_frac[start:cp]))
                seg_sheet    = float(np.mean(sheet_frac[start:cp]))
                seg_coil     = float(np.mean(coil_frac[start:cp]))

                f.write(
                    f"{i},{start},{cp},{time_ns[start]:.6f},{time_ns[cp-1]:.6f},"
                    f"{seg_contacts:.6f},{seg_helix:.6f},{seg_sheet:.6f},{seg_coil:.6f}\n"
                )
                start = cp

        print(f"📄 Change-point summary saved: {cp_out}")

        # ------------------------------------------------------------
        # 7) PLOTTING — SSE TOP, CONTACTS MIDDLE, ROLLING CORR BOTTOM
        # ------------------------------------------------------------
        fig, axes = plt.subplots(
            3, 1,
            figsize=(16, 14),
            sharex=True,
            gridspec_kw={"height_ratios": [2.0, 1.2, 1.0]}
        )

        ax1, ax2, ax3 = axes

        # --- TOP: SSE fractions
        ax1.plot(time_ns, helix_frac, lw=2, color="blue",   label="Helix")
        ax1.plot(time_ns, sheet_frac, lw=2, color="orange", label="Sheet")
        ax1.plot(time_ns, coil_frac,  lw=2, color="green",  label="Coil")
        ax1.set_ylabel("Secondary Structure Fraction")
        ax1.set_title(
            f"{compound_name} — SSE vs Ligand Contacts\n"
            f"Pearson={pear_corr:.3f}   Spearman={spear_corr:.3f}"
        )
        ax1.legend(loc="upper right", frameon=False)

        for cp in change_points:
            ax1.axvline(time_ns[cp], color="black", linestyle="--", alpha=0.35, lw=1)

        # --- MIDDLE: raw + smoothed contacts
        ax2.plot(time_ns, contact_counts.values, lw=0.8, alpha=0.35, color="gray", label="Raw contacts")
        ax2.plot(time_ns, contacts_smooth.values, lw=2.0, color="crimson", label="Smoothed contacts")
        ax2.set_ylabel("Ligand Contacts / Frame")
        ax2.legend(loc="upper right", frameon=False)

        for cp in change_points:
            ax2.axvline(time_ns[cp], color="black", linestyle="--", alpha=0.35, lw=1)

        # --- BOTTOM: rolling correlation
        ax3.plot(time_ns, rolling_corr.values, lw=2.0, color="purple", label="Rolling corr (contacts vs sheet)")
        ax3.axhline(0, color="black", linestyle=":", lw=1)
        ax3.set_ylabel("Rolling Corr")
        ax3.set_xlabel("Time (ns)")
        ax3.legend(loc="upper right", frameon=False)

        for cp in change_points:
            ax3.axvline(time_ns[cp], color="black", linestyle="--", alpha=0.35, lw=1)

        savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_SSE_vs_LigandContacts_Dynamics.png"))

        # ------------------------------------------------------------
        # 8) EXPORT FULL TIME-SERIES TABLE
        # ------------------------------------------------------------
        out_df = pd.DataFrame({
            "Time_ns": time_ns,
            "Helix_Fraction": helix_frac,
            "Sheet_Fraction": sheet_frac,
            "Coil_Fraction": coil_frac,
            "Raw_Contacts": contact_counts.values,
            "Smoothed_Contacts": contacts_smooth.values,
            "Rolling_Corr_Contacts_vs_Sheet": rolling_corr.values
        })

        out_csv = os.path.join(OUTPUT_DIR, f"{lig_resname}_SSE_Contact_Coupling_Timeseries.csv")
        out_df.to_csv(out_csv, index=False)
        print(f"📄 Time-series coupling table saved: {out_csv}")

        checkpoint("chk_D6D.txt")
        _d6d_num, _d6d_plot, _d6d_params = _stage_d6d_outputs()
        _record_stage_manifest("D6D", _d6d_num, _d6d_plot, _d6d_params)
        print("✅ STEP D6D complete.\n")

    else:
        print("⏭️ STEP D6D skipped (checkpoint found).")







if "u" in globals():
    del u
gc.collect()
print("\n✅ Pocket analysis complete — all memory released.")
print("✅ All analyses complete. Outputs saved to:", os.path.abspath(OUTPUT_DIR))

print("\n🧾 Pipeline Summary:")
print("✅ Step A: Final_Trajectory.xtc/pdb — OK")
print("✅ Step B: binding_pocket_only.xtc/pdb — OK")
print("✅ Step C: RMSD/RMSF — OK")
print("✅ Step D: Contact + Network Analysis — OK" if PARTNER_MODE else "✅ Step D: Protein-centric analysis path completed — OK")
print(f"📂 Outputs saved to {os.path.abspath(OUTPUT_DIR)}")
