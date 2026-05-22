#!/usr/bin/env python3
"""
Generic PROTAC / ternary-complex MD analysis.

This script is designed for PROTAC systems prepared by the current
automation workflow where `atomIndex.txt` is the source of truth for
protein atom ranges and component ordering:

  - first protein entry   -> Ligase
  - later protein entries -> Target_A, Target_B, Target_C, ...

The ligand / PROTAC identity is detected dynamically when possible and
can always be overridden with `--ligand`.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import shutil
import shlex
import sys
import time
from collections import defaultdict
from contextlib import suppress
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

PLOTTING_IMPORT_ERROR = None
try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
except Exception as exc:  # pragma: no cover - environment-dependent
    matplotlib = None
    plt = None
    sns = None
    PLOTTING_IMPORT_ERROR = exc

NETWORKX_IMPORT_ERROR = None
try:
    import networkx as nx
except Exception as exc:  # pragma: no cover - environment-dependent
    nx = None
    NETWORKX_IMPORT_ERROR = exc

import numpy as np
PANDAS_IMPORT_ERROR = None
try:
    import pandas as pd
except Exception as exc:  # pragma: no cover - environment-dependent
    pd = None
    PANDAS_IMPORT_ERROR = exc

MDANALYSIS_IMPORT_ERROR = None
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
    from MDAnalysis.lib.distances import distance_array
except Exception as exc:  # pragma: no cover - environment-dependent
    mda = None
    align = None
    distance_array = None
    MDANALYSIS_IMPORT_ERROR = exc

with suppress(ImportError):
    from tqdm.auto import tqdm  # type: ignore

if "tqdm" not in globals():
    tqdm = None

with suppress(ImportError):
    from adjustText import adjust_text  # type: ignore

if "adjust_text" not in globals():
    adjust_text = None


STANDARD_AMINO_ACIDS = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

WATER_NAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P", "SPC", "SPCE"}
ION_NAMES = {
    "NA",
    "CL",
    "K",
    "MG",
    "CA",
    "ZN",
    "MN",
    "FE",
    "CU",
    "CO",
    "NI",
    "CD",
}
SOLVENT_LIKE_NAMES = {
    "DMS",
    "DMSO",
    "PEG",
    "EDO",
    "GOL",
    "EOH",
    "IPA",
    "ACT",
    "ACN",
    "SO4",
    "PO4",
    "MES",
    "HEP",
    "TRS",
}

RESNAME_EXCLUDE = STANDARD_AMINO_ACIDS | WATER_NAMES | ION_NAMES | SOLVENT_LIKE_NAMES

HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "ILE", "LEU", "MET", "PRO"}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}
POLAR_RESIDUES = {"SER", "THR", "ASN", "GLN", "CYS"}
POSITIVE_RESIDUES = {"ARG", "LYS", "HIS"}
NEGATIVE_RESIDUES = {"ASP", "GLU"}

FILE_KEYS = [
    "mol2",
    "str",
    "itp",
    "prm",
    "ini_pdb",
    "pose_match_pdb",
    "from_pdb_pdb",
    "pose_h_pdb",
    "gro",
    "posre",
    "index",
]

TOPOLOGY_IGNORE_STEMS = {
    "forcefield",
    "posre",
    "ions",
    "tip3p",
    "spce",
    "watermodels",
    "atomtypes",
}

SYSTEM_GRO_IGNORE = {
    "npt",
    "nvt",
    "em",
    "complex",
    "newbox",
    "solv",
    "solv_ions",
    "protein_processed",
    "final_trajectory",
    "binding_pocket_only",
    "md_0_1",
    "ions",
}

DESCRIPTION_KEYS = {
    "protein_rmsd": "protac_rmsd",
    "protein_rmsf": "protac_rmsf",
    "protein_rmsf_contacts": "protac_rmsf_contacts",
    "component_rmsf_contacts": "protac_component_rmsf_contacts",
    "ligand_rmsd": "protac_rmsd",
    "ligand_rmsf": "protac_rmsf",
    "radius_of_gyration": "protac_radius_of_gyration",
    "contacts_timeseries": "protac_contacts_timeseries",
    "contact_heatmap": "protac_contact_heatmap",
    "interaction_network": "protac_interaction_network",
    "interaction_composition": "protac_interaction_composition",
    "pca": "protac_pca",
    "water_bridges": "protac_water_bridges",
    "water_bridge_frequency": "protac_water_bridge_frequency",
    "water_bridge_network": "protac_water_bridge_network",
}

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
    return None


def gmx_cmd(subcommand, args=""):
    if GMX_BIN is None:
        raise RuntimeError("GROMACS executable has not been resolved.")
    base = f"{shlex.quote(GMX_BIN)} {subcommand}"
    return f"{base} {args}".strip()


def gmx_argv(subcommand, args=None):
    if GMX_BIN is None:
        raise RuntimeError("GROMACS executable has not been resolved.")
    return [GMX_BIN, subcommand] + list(args or [])


def print_gmx_preflight():
    if GMX_BIN is None:
        log("⚠️ No GROMACS executable resolved; MDAnalysis fallback will be used when possible.")
        return

    log(f"🧬 Using GROMACS executable: {GMX_BIN}")
    try:
        result = subprocess.run([GMX_BIN, "--version"], text=True, capture_output=True)
    except Exception as exc:
        raise RuntimeError(f"Failed to execute '{GMX_BIN} --version': {exc}") from exc

    if result.returncode != 0:
        stderr = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(f"GROMACS executable is not callable: {GMX_BIN}\n{stderr}")

    version_line = next((line.strip() for line in result.stdout.splitlines() if line.strip()), None)
    if version_line:
        log(f"🧪 GROMACS version: {version_line}")


def log(message: str = "") -> None:
    print(message, flush=True)


def ensure_plotting_available() -> None:
    if PLOTTING_IMPORT_ERROR is not None or plt is None or sns is None:
        raise RuntimeError(f"Plotting dependencies are unavailable: {PLOTTING_IMPORT_ERROR}")


def ensure_pandas_available() -> None:
    if PANDAS_IMPORT_ERROR is not None or pd is None:
        raise RuntimeError(
            f"Pandas is unavailable: {PANDAS_IMPORT_ERROR}. "
            "Activate the MDAnalysis Conda environment (see environment_mdanalysis.yml) before running this script."
        )


def ensure_mdanalysis_available() -> None:
    if MDANALYSIS_IMPORT_ERROR is not None or mda is None or align is None or distance_array is None:
        raise RuntimeError(
            f"MDAnalysis is unavailable: {MDANALYSIS_IMPORT_ERROR}. "
            "Activate the MDAnalysis Conda environment (see environment_mdanalysis.yml) before running this script."
        )


def ensure_networkx_available() -> None:
    if NETWORKX_IMPORT_ERROR is not None or nx is None:
        raise RuntimeError(f"NetworkX is unavailable: {NETWORKX_IMPORT_ERROR}")


class RuntimeProfiler:
    def __init__(self, output_path: Path):
        self.output_path = output_path
        self.records: List[dict] = []

    def record(self, stage: str, elapsed_s: float, status: str = "OK", details: str = "") -> None:
        self.records.append(
            {
                "stage": stage,
                "elapsed_s": elapsed_s,
                "status": status,
                "details": details,
            }
        )
        self.write()

    def write(self) -> None:
        total_elapsed = 0.0
        lines = ["PROTAC analysis runtime profile", ""]
        for record in self.records:
            total_elapsed += record["elapsed_s"]
            suffix = f" | {record['details']}" if record["details"] else ""
            lines.append(
                f"{record['status']:>7} | {record['elapsed_s']:9.2f} s | {record['stage']}{suffix}"
            )
        lines.append("")
        lines.append(f"Total recorded runtime: {total_elapsed:.2f} s")
        self.output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class FigureRegistry:
    def __init__(self, qc_dir: Path):
        self.qc_dir = qc_dir
        self.records: List[dict] = []

    def add(
        self,
        figure_path: Path,
        figure_title: str,
        figure_category: str,
        ligand_code: str,
        protein_components: Sequence[str],
        description_key: str,
        notes_hint: str,
        include_in_pdf: bool,
        sort_order: int,
    ) -> None:
        self.records.append(
            {
                "figure_path": str(Path(figure_path).resolve()),
                "figure_title": figure_title,
                "figure_category": figure_category,
                "ligand_code": ligand_code,
                "protein_components": ",".join(protein_components),
                "description_key": description_key,
                "notes_hint": notes_hint,
                "include_in_pdf": bool(include_in_pdf),
                "sort_order": int(sort_order),
            }
        )

    def write(self) -> Tuple[Path, Path]:
        csv_path = self.qc_dir / "protac_figure_manifest.csv"
        json_path = self.qc_dir / "protac_figure_manifest.json"
        df = pd.DataFrame(self.records)
        if df.empty:
            df = pd.DataFrame(
                columns=[
                    "figure_path",
                    "figure_title",
                    "figure_category",
                    "ligand_code",
                    "protein_components",
                    "description_key",
                    "notes_hint",
                    "include_in_pdf",
                    "sort_order",
                ]
            )
        df = df.sort_values(by=["sort_order", "figure_title"], kind="stable")
        df.to_csv(csv_path, index=False)
        json_path.write_text(df.to_json(orient="records", indent=2), encoding="utf-8")
        return csv_path, json_path


def execute_stage(profiler: RuntimeProfiler, stage: str, func, *args, **kwargs):
    log(f"▶ {stage}...")
    started = time.perf_counter()
    try:
        result = func(*args, **kwargs)
    except Exception as exc:
        elapsed = time.perf_counter() - started
        profiler.record(stage, elapsed, status="FAILED", details=str(exc))
        log(f"❌ {stage} failed after {elapsed:.2f} s: {exc}")
        raise
    elapsed = time.perf_counter() - started
    profiler.record(stage, elapsed, status="OK")
    log(f"✅ {stage} completed in {elapsed:.2f} s")
    return result


def iter_with_progress(items: Sequence[int], label: str) -> Iterable[Tuple[int, int]]:
    total = len(items)
    if total == 0:
        return

    if tqdm is not None:
        iterator = tqdm(items, total=total, desc=label, unit="frame", leave=False)
        for index, item in enumerate(iterator, start=1):
            yield index, item
        return

    every = max(1, total // 10)
    log(f"   ⏳ {label}: 0/{total}")
    for index, item in enumerate(items, start=1):
        if index == 1 or index % every == 0 or index == total:
            log(f"   ⏳ {label}: {index}/{total}")
        yield index, item


def stream_subprocess(command: Sequence[str], input_text: Optional[str] = None, cwd: Optional[Path] = None) -> int:
    process = subprocess.Popen(
        list(command),
        cwd=str(cwd) if cwd else None,
        stdin=subprocess.PIPE if input_text is not None else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        universal_newlines=True,
        bufsize=1,
    )
    if input_text is not None and process.stdin is not None:
        process.stdin.write(input_text)
        process.stdin.close()
    if process.stdout is not None:
        for line in process.stdout:
            log(line.rstrip())
    return process.wait()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def file_exists(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def copy_if_needed(source: Path, target: Path, overwrite: bool = False) -> str:
    if not file_exists(source):
        return "missing_source"
    if file_exists(target) and not overwrite:
        return "reused"
    ensure_dir(target.parent)
    shutil.copy2(source, target)
    return "created" if overwrite or not target.exists() else "created"


def default_pdb_reference() -> Optional[str]:
    candidate = Path("Complex.pdb")
    return str(candidate) if candidate.exists() else None


def normalize_candidate_code(value: str) -> Optional[str]:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", value.strip())
    if 2 <= len(cleaned) <= 5:
        return cleaned.upper()
    return None


def parse_name_list(value: str) -> List[str]:
    return [item.strip().upper() for item in value.split(",") if item.strip()]


def role_name_for_target(index: int) -> str:
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if index < len(letters):
        return f"Target_{letters[index]}"
    return f"Target_{index + 1}"


def parse_range_metadata(value: Optional[str]) -> Tuple[Optional[int], Optional[int]]:
    if not value:
        return None, None
    match = re.match(r"^\s*(-?\d+)\s*-\s*(-?\d+)\s*$", value)
    if not match:
        return None, None
    return int(match.group(1)), int(match.group(2))


def apply_runtime_mode_defaults(args: argparse.Namespace) -> argparse.Namespace:
    if args.contacts_only:
        args.skip_pca = True
        args.skip_network = True
        args.skip_contact_map = True
        args.skip_rmsf = True

    if args.basic_only:
        args.skip_pca = True
        args.skip_network = True
        args.skip_contact_map = True
        args.skip_rmsf = True

    if args.quick_test:
        if args.start_frame is None:
            args.start_frame = 0
        if args.end_frame is None:
            args.end_frame = 500
        if args.frame_step is None:
            args.frame_step = 10

    if args.frame_step is None:
        args.frame_step = 10

    return args


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze a PROTAC ternary-complex trajectory using atomIndex.txt-defined ligase/target roles.",
        epilog=(
            "Examples:\n"
            "  python 3_PROTAC_Analysis.py --dry-run --headless\n"
            "  python 3_PROTAC_Analysis.py --quick-test --headless\n"
            "  python 3_PROTAC_Analysis.py --basic-only --frame-step 10 --headless\n"
            "  python 3_PROTAC_Analysis.py --frame-step 1 --headless"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-l", "--ligand", default=None, help="Ligand / PROTAC residue name. Optional when auto-detection is confident.")
    parser.add_argument("--gmx-bin", default="auto", help="GROMACS executable to use: auto, gmx, gmx_mpi, gmx-mpi, or an explicit path.")
    parser.add_argument("--topo", default="npt.gro", help="Topology/coordinate file for MDAnalysis (default: npt.gro).")
    parser.add_argument("--traj", default="md_0_1.xtc", help="Trajectory file for MDAnalysis (default: md_0_1.xtc).")
    parser.add_argument("--pdb-reference", default=None, help="Optional reference PDB. Defaults to Complex.pdb when present.")
    parser.add_argument("--atomindex", default="atomIndex.txt", help="atomIndex file with ordered protein atom ranges (default: atomIndex.txt).")
    parser.add_argument("--outdir", default="Analysis_Results/PROTAC", help="Output directory root (default: Analysis_Results/PROTAC).")
    parser.add_argument("--distance-cutoff", type=float, default=4.0, help="Protein-ligand contact cutoff in angstroms (default: 4.0).")
    parser.add_argument("--headless", action="store_true", help="Disable interactive prompts and fail clearly when detection is ambiguous.")
    parser.add_argument("--start-frame", type=int, default=None, help="First frame index to analyze for trajectory-scan analyses.")
    parser.add_argument("--end-frame", type=int, default=None, help="Last frame index (exclusive) to analyze for trajectory-scan analyses.")
    parser.add_argument("--frame-step", type=int, default=None, help="Frame stride for expensive analyses (default for PROTAC analysis: 10).")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs, run ligand detection, QC generation, and exit before expensive analysis.")
    parser.add_argument("--quick-test", action="store_true", help="Smoke-test mode: defaults to start=0, end=500, step=10 unless frame flags were explicitly provided.")
    parser.add_argument("--skip-pca", action="store_true", help="Skip PCA analysis.")
    parser.add_argument("--skip-network", action="store_true", help="Skip ligand interaction network generation.")
    parser.add_argument("--skip-contact-map", action="store_true", help="Skip ligand contact heatmap generation.")
    parser.add_argument("--skip-rmsf", action="store_true", help="Skip RMSF outputs.")
    parser.add_argument("--contacts-only", action="store_true", help="Run only ligand contact analysis and QC contact summaries.")
    parser.add_argument("--basic-only", action="store_true", help="Run RMSD, radius of gyration, ligand contacts, and QC only; skip PCA, network, contact-map, and RMSF outputs.")
    parser.add_argument("--residue-numbering", choices=["source", "current"], default="source", help="Residue numbering style for display labels (default: source).")
    parser.add_argument("--network-min-contact-fraction", type=float, default=0.05, help="Minimum contact fraction for residues shown in the network figure (default: 0.05).")
    parser.add_argument("--network-max-residues", type=int, default=40, help="Maximum number of residue nodes shown in the network figure (default: 40).")
    parser.add_argument("--contact-map-min-contact-fraction", type=float, default=0.02, help="Minimum contact fraction for residues shown in the contact heatmap (default: 0.02).")
    parser.add_argument("--contact-map-max-residues", type=int, default=80, help="Maximum number of residues shown in the contact heatmap (default: 80).")
    parser.add_argument("--skip-preprocess", action="store_true", help="Skip writing PROTAC subset structures and binding-pocket preprocessing outputs.")
    parser.add_argument("--preprocess-only", action="store_true", help="Generate subset structures, groups, and binding-pocket outputs, then exit before RMSD/contact plots.")
    parser.add_argument("--pocket-cutoff", type=float, default=5.0, help="Pocket residue cutoff from the PROTAC in angstroms (default: 5.0).")
    parser.add_argument("--pocket-mode", choices=["first_frame", "any_sampled_frame", "persistent"], default="first_frame", help="Pocket selection mode (default: first_frame).")
    parser.add_argument("--pocket-min-fraction", type=float, default=0.05, help="Minimum contact fraction for persistent pocket mode (default: 0.05).")
    parser.add_argument("--interaction-mode", choices=["residue_chemistry", "geometry"], default="residue_chemistry", help="Broad interaction classification mode (default: residue_chemistry).")
    parser.add_argument("--plot-allprotein", action="store_true", help="Include the full-protein aggregate trace in user-facing contact plots.")
    parser.add_argument("--write-subset-trajectories", dest="write_subset_trajectories", action="store_true", help="Write subset XTC trajectories in addition to subset PDB structures.")
    parser.add_argument("--no-write-subset-trajectories", dest="write_subset_trajectories", action="store_false", help="Do not write subset XTC trajectories.")
    parser.add_argument("--write-final-trajectory", dest="write_final_trajectory", action="store_true", help="Write centered Final_Trajectory compatibility files (default: enabled).")
    parser.add_argument("--no-write-final-trajectory", dest="write_final_trajectory", action="store_false", help="Disable Final_Trajectory compatibility-file generation.")
    parser.add_argument("--use-final-trajectory", dest="use_final_trajectory", action="store_true", help="Reload the centered Final_Trajectory for downstream analysis (default: enabled).")
    parser.add_argument("--no-use-final-trajectory", dest="use_final_trajectory", action="store_false", help="Keep using the original topology/trajectory after preprocessing.")
    parser.add_argument("--overwrite-final-trajectory", action="store_true", help="Overwrite existing Final_Trajectory files and remapping outputs.")
    parser.add_argument("--overwrite-structures", action="store_true", help="Overwrite existing subset structure outputs such as binding_pocket_only aliases.")
    parser.add_argument("--final-trajectory-name", default="Final_Trajectory", help="Base name for compatibility final-trajectory files (default: Final_Trajectory).")
    parser.add_argument("--rmsf-contact-threshold", type=float, default=0.05, help="Persistent-contact threshold used for RMSF contact overlays (default: 0.05).")
    parser.add_argument("--rmsf-label-top-contacts", type=int, default=10, help="Maximum number of contacting residues to annotate on RMSF contact figures (default: 10).")
    parser.add_argument("--rmsf-contact-marker-mode", choices=["vlines", "scatter", "both"], default="both", help="Marker style for RMSF contact overlays (default: both).")
    parser.add_argument("--analyze-water-bridges", action="store_true", help="Run optional ligand-water-protein bridge analysis on the original raw trajectory.")
    parser.add_argument("--water-bridge-cutoff", type=float, default=3.5, help="Cutoff in angstroms for ligand-water and water-protein bridge detection (default: 3.5).")
    parser.add_argument("--water-pocket-cutoff", type=float, default=6.0, help="Cutoff in angstroms for optional water-pocket trajectory selection (default: 6.0).")
    parser.add_argument("--write-water-pocket-trajectory", action="store_true", help="Write a raw-trajectory water-pocket subset containing ligand, nearby protein residues, and relevant waters.")
    parser.add_argument("--water-names", default="SOL,WAT,HOH,TIP3,TIP3P,SPC,SPCE", help="Comma-separated water residue names for water-bridge analysis.")
    parser.set_defaults(
        write_subset_trajectories=False,
        write_final_trajectory=True,
        use_final_trajectory=True,
    )
    return apply_runtime_mode_defaults(parser.parse_args())


def parse_atomindex_entries(atomindex_path: str) -> List[dict]:
    path = Path(atomindex_path)
    if not path.exists():
        raise FileNotFoundError(
            f"atomIndex file not found: {path}. This script requires atomIndex.txt from the current preparation workflow."
        )

    entries: List[dict] = []
    invalid_lines: List[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue

            match = re.match(r"^(.*?)\s+(\d+)\s*-\s*(\d+)(?:\s*#(.*))?$", stripped)
            if not match:
                invalid_lines.append(f"line {line_number}: {raw.rstrip()}")
                continue

            name = match.group(1).strip()
            start_1based = int(match.group(2))
            end_1based = int(match.group(3))
            metadata_comment = (match.group(4) or "").strip()

            if not name or start_1based <= 0 or end_1based < start_1based:
                invalid_lines.append(f"line {line_number}: {raw.rstrip()}")
                continue

            metadata = dict(re.findall(r"([A-Za-z_]+)=([^\s#]+)", metadata_comment))
            current_res_start, current_res_end = parse_range_metadata(metadata.get("current_residues"))
            source_res_start, source_res_end = parse_range_metadata(metadata.get("source_residues"))

            source_rescount = metadata.get("source_rescount")
            try:
                source_rescount_i = int(source_rescount) if source_rescount is not None else None
            except ValueError:
                source_rescount_i = None

            entries.append(
                {
                    "name": name,
                    "display_name": name,
                    "start_1based": start_1based,
                    "end_1based": end_1based,
                    "start0": start_1based - 1,
                    "end0": end_1based - 1,
                    "current_chain": metadata.get("current_chain"),
                    "current_res_start": current_res_start,
                    "current_res_end": current_res_end,
                    "source_chain": metadata.get("source_chain"),
                    "source_res_start": source_res_start,
                    "source_res_end": source_res_end,
                    "source_rescount": source_rescount_i,
                    "metadata_comment": metadata_comment,
                }
            )

    if invalid_lines:
        joined = "\n".join(invalid_lines[:10])
        raise ValueError(
            "Failed to parse atomIndex.txt. Expected lines like 'Name 1-1975 # current_chain=A ...'.\n"
            f"Unparsable entries:\n{joined}"
        )

    if not entries:
        raise ValueError(
            f"No atom ranges were parsed from {path}. Please regenerate atomIndex.txt with the current workflow."
        )

    return entries


def assign_protac_roles(entries: Sequence[dict]) -> List[dict]:
    assigned: List[dict] = []
    for index, entry in enumerate(entries):
        role = "Ligase" if index == 0 else role_name_for_target(index - 1)
        enriched = dict(entry)
        enriched["role"] = role
        enriched["display_name"] = entry.get("name") or role
        enriched["component_order"] = index
        assigned.append(enriched)
    return assigned


def component_entries(entries: Sequence[dict]) -> List[dict]:
    return list(entries)


def target_entries(entries: Sequence[dict]) -> List[dict]:
    return [entry for entry in entries if entry["role"].startswith("Target_")]


def has_multiple_targets(entries: Sequence[dict]) -> bool:
    return len(target_entries(entries)) > 1


def plot_display_name_for_entry(entry: dict) -> str:
    return entry.get("display_name") or entry.get("name") or entry.get("role") or "Protein"


def atomgroup_for_entry(universe: mda.Universe, entry: dict) -> mda.AtomGroup:
    start0 = int(entry["start0"])
    end0 = int(entry["end0"])
    if start0 < 0 or end0 >= universe.atoms.n_atoms:
        raise ValueError(
            f"atomIndex range {entry['name']} {entry['start_1based']}-{entry['end_1based']} falls outside topology atom count ({universe.atoms.n_atoms})."
        )
    return universe.atoms[start0 : end0 + 1]


def overlap_count(indices: np.ndarray, start0: int, end0: int) -> int:
    return int(np.count_nonzero((indices >= start0) & (indices <= end0)))


def entry_for_residue(residue: mda.core.groups.Residue, entries: Sequence[dict]) -> Optional[dict]:
    indices = residue.atoms.indices
    best_entry = None
    best_overlap = 0
    for entry in entries:
        count = overlap_count(indices, int(entry["start0"]), int(entry["end0"]))
        if count > best_overlap:
            best_overlap = count
            best_entry = entry
    return best_entry


def residue_role_for_residue(residue: mda.core.groups.Residue, entries: Sequence[dict]) -> Tuple[Optional[str], Optional[str]]:
    entry = entry_for_residue(residue, entries)
    if not entry:
        return None, None
    return entry["role"], entry["name"]


def residue_first_atom_key(residue: mda.core.groups.Residue) -> int:
    return int(residue.atoms.indices[0])


def get_display_residue_id(residue: mda.core.groups.Residue, entry: dict) -> int:
    current_resid = int(residue.resid)
    current_start = entry.get("current_res_start")
    current_end = entry.get("current_res_end")
    source_start = entry.get("source_res_start")
    source_end = entry.get("source_res_end")

    if None not in (current_start, current_end, source_start, source_end):
        if current_start <= current_resid <= current_end:
            mapped = source_start + (current_resid - current_start)
            if source_start <= mapped <= source_end:
                return int(mapped)
    return current_resid


def residue_label(
    residue: mda.core.groups.Residue,
    entry: dict,
    numbering: str = "source",
) -> str:
    protein_name = plot_display_name_for_entry(entry)
    current_chain = entry.get("current_chain")
    current_resid = int(residue.resid)
    source_resid = get_display_residue_id(residue, entry)
    residue_name = residue.resname

    if numbering == "current":
        chain_part = f"{current_chain}:" if current_chain else ""
        return f"{protein_name}:{chain_part}{residue_name}{current_resid}_current"
    return f"{protein_name}:{residue_name}{source_resid}"


def select_ligand(universe: mda.Universe, ligand_code: str) -> mda.AtomGroup:
    code = ligand_code.strip().upper()
    indices: List[int] = []
    for residue in universe.residues:
        if residue.resname.strip().upper() == code:
            indices.extend(residue.atoms.indices.tolist())
    if not indices:
        return universe.atoms[np.array([], dtype=int)]
    unique_indices = np.array(sorted(set(indices)), dtype=int)
    return universe.atoms[unique_indices]


def validate_protac_inputs(universe: mda.Universe, ligand: mda.AtomGroup, entries: Sequence[dict]) -> None:
    if universe.atoms.n_atoms == 0:
        raise ValueError("Loaded universe has zero atoms.")
    if not entries:
        raise ValueError("No atomIndex entries were provided.")
    for entry in entries:
        atomgroup_for_entry(universe, entry)
    if ligand.n_atoms == 0:
        raise ValueError("Selected ligand / PROTAC has zero atoms in the loaded topology.")


def parse_topol_top(topol_path: Path) -> dict:
    result = {"include_candidates": set(), "molecule_candidates": set()}
    if not topol_path.exists():
        return result

    in_molecules = False
    with topol_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue

            if line.lower().startswith("#include"):
                match = re.search(r'"([^"]+)"', raw)
                if match:
                    include_name = Path(match.group(1)).name
                    stem = Path(include_name).stem
                    normalized = normalize_candidate_code(stem)
                    lower_stem = stem.lower()
                    if include_name.lower().endswith((".itp", ".prm")) and normalized:
                        if lower_stem not in TOPOLOGY_IGNORE_STEMS and not lower_stem.startswith("posre"):
                            result["include_candidates"].add(normalized)
                continue

            if line.startswith("["):
                in_molecules = line.lower().startswith("[ molecules ]")
                continue

            if in_molecules:
                tokens = line.split()
                if len(tokens) >= 2:
                    normalized = normalize_candidate_code(tokens[0])
                    if normalized and normalized not in RESNAME_EXCLUDE:
                        result["molecule_candidates"].add(normalized)
    return result


def discover_ligand_files(run_dir: Path) -> Dict[str, dict]:
    candidates: Dict[str, dict] = {}

    def get_candidate(code: str) -> dict:
        if code not in candidates:
            candidates[code] = {
                "score": 0,
                "sources": [],
                "files": {key: None for key in FILE_KEYS},
                "warnings": [],
                "errors": [],
                "residue_atom_count": 0,
                "residue_count": 0,
            }
        return candidates[code]

    for path in run_dir.iterdir():
        if not path.is_file():
            continue

        name_lower = path.name.lower()
        entry_code: Optional[str] = None
        file_key: Optional[str] = None

        if name_lower.endswith(".cgenff.mol2"):
            entry_code = normalize_candidate_code(path.name[: -len(".cgenff.mol2")])
            file_key = "mol2"
        elif name_lower.endswith(".str"):
            entry_code = normalize_candidate_code(path.stem)
            file_key = "str"
        elif name_lower.endswith(".itp"):
            if name_lower.startswith("posre_"):
                entry_code = normalize_candidate_code(path.stem[6:])
                file_key = "posre"
            else:
                entry_code = normalize_candidate_code(path.stem)
                file_key = "itp"
        elif name_lower.endswith(".prm"):
            entry_code = normalize_candidate_code(path.stem)
            file_key = "prm"
        elif name_lower.endswith(".ndx"):
            if name_lower.startswith("index_"):
                entry_code = normalize_candidate_code(path.stem[6:])
                file_key = "index"
        elif name_lower.endswith("_ini.pdb"):
            entry_code = normalize_candidate_code(path.name[: -len("_ini.pdb")])
            file_key = "ini_pdb"
        elif name_lower.endswith("_pose_match.pdb"):
            entry_code = normalize_candidate_code(path.name[: -len("_pose_match.pdb")])
            file_key = "pose_match_pdb"
        elif name_lower.endswith("_from_pdb.pdb"):
            entry_code = normalize_candidate_code(path.name[: -len("_from_pdb.pdb")])
            file_key = "from_pdb_pdb"
        elif name_lower.endswith("_pose_h.pdb"):
            entry_code = normalize_candidate_code(path.name[: -len("_pose_h.pdb")])
            file_key = "pose_h_pdb"
        elif name_lower.endswith(".gro"):
            if path.stem.lower() not in SYSTEM_GRO_IGNORE:
                entry_code = normalize_candidate_code(path.stem)
                file_key = "gro"

        if entry_code and file_key:
            get_candidate(entry_code)["files"][file_key] = str(path.resolve())

    for code, data in candidates.items():
        files = data["files"]
        if all(files[key] for key in ["mol2", "str", "itp", "prm", "ini_pdb", "pose_match_pdb", "gro"]):
            data["score"] += 4
            data["sources"].append("cgenff_file_set")
        if files["gro"]:
            data["score"] += 2
            data["sources"].append("gro_residue")
        if files["ini_pdb"] or files["pose_match_pdb"] or files["from_pdb_pdb"] or files["pose_h_pdb"]:
            data["score"] += 2
            data["sources"].append("cgenff_file_set")
        if files["posre"] or files["index"]:
            data["score"] += 1
            data["sources"].append("index_ndx")

    return candidates


def scan_structure_for_ligands(
    topo_path: Optional[str],
    traj_path: Optional[str],
) -> Tuple[Dict[str, dict], List[str], Optional[str]]:
    residue_info: Dict[str, dict] = {}
    warnings: List[str] = []
    source_label: Optional[str] = None

    if not topo_path or not Path(topo_path).exists():
        return residue_info, warnings, source_label

    try:
        if traj_path and Path(traj_path).exists():
            universe = mda.Universe(topo_path, traj_path)
            source_label = "trajectory_residue"
        else:
            universe = mda.Universe(topo_path)
            source_label = "gro_residue"
    except Exception as exc:
        warnings.append(f"Could not inspect topology/trajectory for ligand residues: {exc}")
        return residue_info, warnings, source_label

    counts = defaultdict(lambda: {"atoms": 0, "residues": 0})
    for residue in universe.residues:
        resname = residue.resname.strip().upper()
        if resname in RESNAME_EXCLUDE:
            continue
        normalized = normalize_candidate_code(resname)
        if not normalized:
            continue
        counts[normalized]["atoms"] += residue.atoms.n_atoms
        counts[normalized]["residues"] += 1

    for code, info in counts.items():
        residue_info[code] = {
            "score": 2,
            "sources": [source_label] if source_label else [],
            "residue_atom_count": info["atoms"],
            "residue_count": info["residues"],
        }

    if len(residue_info) > 1:
        ranked = sorted(residue_info.items(), key=lambda item: (-item[1]["residue_atom_count"], item[0]))
        summary = ", ".join(
            f"{code}({data['residue_count']} res / {data['residue_atom_count']} atoms)"
            for code, data in ranked
        )
        warnings.append(f"Multiple ligand-like residues detected in structure: {summary}")

    return residue_info, warnings, source_label


def merge_candidate_evidence(base: Dict[str, dict], evidence: Dict[str, dict]) -> Dict[str, dict]:
    merged = {key: value for key, value in base.items()}
    for code, data in evidence.items():
        if code not in merged:
            merged[code] = {
                "score": 0,
                "sources": [],
                "files": {key: None for key in FILE_KEYS},
                "warnings": [],
                "errors": [],
                "residue_atom_count": 0,
                "residue_count": 0,
            }
        merged[code]["score"] += data.get("score", 0)
        merged[code]["sources"].extend(data.get("sources", []))
        merged[code]["residue_atom_count"] = max(merged[code]["residue_atom_count"], data.get("residue_atom_count", 0))
        merged[code]["residue_count"] = max(merged[code]["residue_count"], data.get("residue_count", 0))
    return merged


def source_string(sources: Sequence[str]) -> str:
    ordered = []
    for source in ["cli", "cgenff_file_set", "topol_top", "gro_residue", "trajectory_residue", "index_ndx", "manual"]:
        if source in sources and source not in ordered:
            ordered.append(source)
    for source in sources:
        if source not in ordered:
            ordered.append(source)
    return "|".join(ordered)


def score_confidence(best_score: int, lead: int, has_cli: bool) -> str:
    if has_cli or best_score >= 9 or (best_score >= 7 and lead >= 2):
        return "high"
    if best_score >= 5 and lead >= 1:
        return "medium"
    return "low"


def prompt_user_for_ligand(candidates: Sequence[Tuple[str, dict]]) -> str:
    log("\nLigand auto-detection is ambiguous. Candidate summary:")
    for index, (code, data) in enumerate(candidates, start=1):
        log(
            f"  {index}. {code} | score={data['score']} | residues={data.get('residue_count', 0)} | atoms={data.get('residue_atom_count', 0)} | sources={source_string(data.get('sources', [])) or 'manual'}"
        )
    while True:
        choice = input("Select a ligand by number or type a residue name: ").strip()
        if choice.isdigit():
            number = int(choice)
            if 1 <= number <= len(candidates):
                return candidates[number - 1][0]
        normalized = normalize_candidate_code(choice)
        if normalized:
            return normalized
        log("Please enter a valid choice.")


def detect_protac_ligand_setup(
    run_dir: str = ".",
    topo_path: Optional[str] = None,
    traj_path: Optional[str] = None,
    cli_ligand: Optional[str] = None,
    headless: bool = False,
) -> dict:
    run_path = Path(run_dir).resolve()
    candidates = discover_ligand_files(run_path)
    topol_info = parse_topol_top(run_path / "topol.top")

    for code in topol_info["include_candidates"] | topol_info["molecule_candidates"]:
        if code not in candidates:
            candidates[code] = {
                "score": 0,
                "sources": [],
                "files": {key: None for key in FILE_KEYS},
                "warnings": [],
                "errors": [],
                "residue_atom_count": 0,
                "residue_count": 0,
            }
        candidates[code]["score"] += 3
        candidates[code]["sources"].append("topol_top")

    structure_evidence, warnings, structure_source = scan_structure_for_ligands(topo_path, traj_path)
    candidates = merge_candidate_evidence(candidates, structure_evidence)

    cli_code = normalize_candidate_code(cli_ligand) if cli_ligand else None
    if cli_code:
        if cli_code not in candidates:
            candidates[cli_code] = {
                "score": 0,
                "sources": [],
                "files": {key: None for key in FILE_KEYS},
                "warnings": [],
                "errors": [],
                "residue_atom_count": 0,
                "residue_count": 0,
            }
        candidates[cli_code]["score"] += 5
        candidates[cli_code]["sources"].append("cli")

    ranked = sorted(
        candidates.items(),
        key=lambda item: (-item[1]["score"], -item[1].get("residue_atom_count", 0), -item[1].get("residue_count", 0), item[0]),
    )

    selected_code: Optional[str] = cli_code
    errors: List[str] = []

    if not selected_code:
        if not ranked:
            if headless:
                errors.append("Could not auto-detect a PROTAC ligand from the current directory. Pass --ligand <RESNAME> in headless mode.")
            else:
                manual = input("Ligand could not be auto-detected. Enter ligand residue name: ").strip()
                selected_code = normalize_candidate_code(manual)
                if not selected_code:
                    errors.append("No valid ligand residue name was provided.")
        else:
            best_code, best_data = ranked[0]
            if len(ranked) > 1 and best_data["score"] == ranked[1][1]["score"]:
                if headless:
                    summary = "; ".join(
                        f"{code}: score={data['score']} sources={source_string(data['sources']) or 'manual'}"
                        for code, data in ranked[:5]
                    )
                    errors.append(f"Ligand detection is ambiguous in headless mode. Top candidates: {summary}. Pass --ligand explicitly.")
                else:
                    selected_code = prompt_user_for_ligand(ranked[:5])
            elif best_data["score"] < 3:
                if headless:
                    errors.append("Ligand auto-detection confidence is too low for headless mode. Pass --ligand explicitly.")
                else:
                    selected_code = prompt_user_for_ligand(ranked[:5])
            else:
                selected_code = best_code

    selected = candidates.get(selected_code) if selected_code else None
    selected_files = {key: None for key in FILE_KEYS}
    confidence = "low"
    source = "manual" if selected_code else ""

    if selected:
        selected_files.update(selected["files"])
        best_score = selected.get("score", 0)
        second_score = ranked[1][1]["score"] if len(ranked) > 1 and ranked[0][0] == selected_code else 0
        confidence = score_confidence(best_score, best_score - second_score, bool(cli_code))
        source = source_string(selected.get("sources", []))

    residue_candidates = {code: data for code, data in structure_evidence.items() if data.get("residue_atom_count", 0) > 0}
    selected_warnings = list(warnings)
    if selected_code and residue_candidates and selected_code not in residue_candidates:
        selected_warnings.append(
            f"Selected ligand code {selected_code} was not observed as a residue in the loaded structure. Detected non-protein residues: "
            + ", ".join(f"{code}({data['residue_atom_count']} atoms)" for code, data in sorted(residue_candidates.items()))
        )

    if cli_code and structure_source and cli_code not in residue_candidates:
        errors.append(f"--ligand {cli_code} was provided, but that residue name was not found in the loaded topology/trajectory.")

    if selected_code and len(residue_candidates) > 1:
        others = [code for code in residue_candidates if code != selected_code]
        if others:
            selected_warnings.append(f"Additional ligand-like residues were detected: {', '.join(sorted(others))}")

    return {
        "ligand_code": selected_code,
        "confidence": confidence,
        "source": source,
        "files": selected_files,
        "warnings": selected_warnings,
        "errors": errors,
        "topology_residue_atoms": selected.get("residue_atom_count", 0) if selected else 0,
        "topology_residue_count": selected.get("residue_count", 0) if selected else 0,
        "candidate_scores": [
            {
                "ligand_code": code,
                "score": data["score"],
                "sources": source_string(data.get("sources", [])),
                "residue_atom_count": data.get("residue_atom_count", 0),
                "residue_count": data.get("residue_count", 0),
            }
            for code, data in ranked
        ],
    }


def write_ligand_preflight_reports(ligand_setup: dict, qc_dir: Path) -> Tuple[Path, Path]:
    ensure_dir(qc_dir)
    json_path = qc_dir / "protac_ligand_preflight.json"
    txt_path = qc_dir / "protac_ligand_preflight.txt"
    json_path.write_text(json.dumps(ligand_setup, indent=2), encoding="utf-8")

    lines = [
        "PROTAC ligand preflight",
        f"Ligand selected: {ligand_setup.get('ligand_code') or 'None'}",
        f"Confidence: {ligand_setup.get('confidence', 'unknown')}",
        f"Source: {ligand_setup.get('source', '') or 'manual'}",
        "Found files:",
    ]
    for key in FILE_KEYS:
        value = ligand_setup.get("files", {}).get(key)
        marker = "✓" if value else "-"
        lines.append(f"  {marker} {key}: {value or 'missing'}")
    lines.append(f"Topology residue atoms: {ligand_setup.get('topology_residue_atoms', 0)}")
    lines.append(f"Topology residue count: {ligand_setup.get('topology_residue_count', 0)}")
    warnings = ligand_setup.get("warnings", [])
    errors = ligand_setup.get("errors", [])
    lines.append("Warnings:")
    lines.extend([f"  - {warning}" for warning in warnings] or ["  none"])
    lines.append("Errors:")
    lines.extend([f"  - {error}" for error in errors] or ["  none"])
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, txt_path


def print_preflight_block(ligand_setup: dict) -> None:
    source_display = (ligand_setup.get("source") or "manual").replace("|", " + ")
    log("\n" + "=" * 72)
    log("🧬 PROTAC ligand preflight")
    log(f"Ligand selected: {ligand_setup.get('ligand_code') or 'None'}")
    log(f"Confidence: {ligand_setup.get('confidence', 'unknown')}")
    log(f"Source: {source_display}")
    log("Found files:")
    any_file = False
    for key in FILE_KEYS:
        value = ligand_setup.get("files", {}).get(key)
        if value:
            any_file = True
            log(f"  ✓ {Path(value).name}")
    if not any_file:
        log("  none")
    log(f"Topology residue atoms: {ligand_setup.get('topology_residue_atoms', 0)}")
    warnings = ligand_setup.get("warnings", [])
    if warnings:
        log("Warnings:")
        for warning in warnings:
            log(f"  - {warning}")
    else:
        log("Warnings: none")
    log("=" * 72)


def resolve_frame_window(
    n_frames: int,
    start_frame: Optional[int],
    end_frame: Optional[int],
    frame_step: int,
) -> Tuple[int, int, int, List[int]]:
    start = 0 if start_frame is None else max(0, start_frame)
    end = n_frames if end_frame is None else min(n_frames, end_frame)
    step = max(1, frame_step)
    if start >= end:
        raise ValueError(f"Invalid frame window: start={start}, end={end}, total_frames={n_frames}")
    return start, end, step, list(range(start, end, step))


def describe_frame_selection(start: int, end: int, step: int, sampled_count: int) -> str:
    return f"Analyzing frames {start}:{end}:{step} ({sampled_count} sampled frames)"


def superpose_positions(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    mobile_centered = mobile - mobile.mean(axis=0)
    reference_centered = reference - reference.mean(axis=0)
    if len(mobile) < 3:
        return mobile_centered
    rotation, _ = align.rotation_matrix(mobile_centered, reference_centered)
    return np.dot(mobile_centered, rotation)


def gather_aligned_positions(
    universe: mda.Universe,
    atomgroup: mda.AtomGroup,
    frame_indices: Sequence[int],
    label: str,
) -> np.ndarray:
    if atomgroup.n_atoms == 0:
        return np.empty((0, 0, 3))
    log(f"   ⏳ {label}: {len(frame_indices)} frames, {atomgroup.n_atoms} atoms")
    universe.trajectory[frame_indices[0]]
    reference = atomgroup.positions.copy()
    aligned_frames = []
    for _, frame in iter_with_progress(frame_indices, label):
        universe.trajectory[frame]
        aligned_frames.append(superpose_positions(atomgroup.positions.copy(), reference))
    return np.asarray(aligned_frames)


def rmsd_from_aligned_positions(aligned_positions: np.ndarray) -> np.ndarray:
    if aligned_positions.size == 0:
        return np.array([])
    reference = aligned_positions[0]
    diffs = aligned_positions - reference
    return np.sqrt(np.mean(np.sum(diffs * diffs, axis=2), axis=1))


def rmsf_from_aligned_positions(aligned_positions: np.ndarray) -> np.ndarray:
    if aligned_positions.size == 0:
        return np.array([])
    mean_positions = aligned_positions.mean(axis=0)
    diffs = aligned_positions - mean_positions
    return np.sqrt(np.mean(np.sum(diffs * diffs, axis=2), axis=0))


def classify_interaction_type(residue: mda.core.groups.Residue, min_distance: float, mode: str = "residue_chemistry") -> str:
    resname = residue.resname.strip().upper()
    labels: List[str] = []
    if resname in HYDROPHOBIC_RESIDUES and min_distance <= 4.5:
        labels.append("Hydrophobic")
    if resname in AROMATIC_RESIDUES and min_distance <= 4.8:
        labels.append("Aromatic")
    if resname in POSITIVE_RESIDUES | NEGATIVE_RESIDUES and min_distance <= 4.0:
        labels.append("Ionic")
    if mode == "geometry" and resname in POLAR_RESIDUES | POSITIVE_RESIDUES | NEGATIVE_RESIDUES and min_distance <= 3.3:
        labels.append("H-bond-like")
    elif resname in POLAR_RESIDUES | POSITIVE_RESIDUES | NEGATIVE_RESIDUES and min_distance <= 3.8:
        labels.append("Polar")
    return ";".join(labels) if labels else "Contact"


def build_protein_residue_records(
    protein: mda.AtomGroup,
    entries: Sequence[dict],
    numbering: str,
) -> List[dict]:
    local_lookup = {global_index: local_index for local_index, global_index in enumerate(protein.indices)}
    records: List[dict] = []
    for residue in protein.residues:
        entry = entry_for_residue(residue, entries)
        if not entry:
            continue
        protein_local_indices = np.array(
            [local_lookup[index] for index in residue.atoms.indices if index in local_lookup],
            dtype=int,
        )
        if protein_local_indices.size == 0:
            continue

        current_resid = int(residue.resid)
        source_resid = get_display_residue_id(residue, entry)
        current_chain = entry.get("current_chain")
        source_chain = entry.get("source_chain") or current_chain
        label_current = residue_label(residue, entry, numbering="current")
        label_source = residue_label(residue, entry, numbering="source")
        label_display = residue_label(residue, entry, numbering=numbering)

        records.append(
            {
                "residue": residue,
                "entry": entry,
                "role": entry["role"],
                "protein_name": entry["name"],
                "display_name": plot_display_name_for_entry(entry),
                "current_chain": current_chain,
                "current_resid": current_resid,
                "source_chain": source_chain,
                "source_resid": source_resid,
                "residue_name": residue.resname,
                "residue_label_current": label_current,
                "residue_label_source": label_source,
                "residue_label_display": label_display,
                "protein_local_indices": protein_local_indices,
                "component_order": entry["component_order"],
                "residue_sort_key": source_resid if numbering == "source" else current_resid,
                "residue_key": residue_first_atom_key(residue),
            }
        )
    return records


def write_system_roles(entries: Sequence[dict], qc_dir: Path) -> Path:
    rows = []
    for entry in entries:
        rows.append(
            {
                "Name": entry["name"],
                "DisplayName": entry["display_name"],
                "Role": entry["role"],
                "Start_1based": entry["start_1based"],
                "End_1based": entry["end_1based"],
                "Start0": entry["start0"],
                "End0": entry["end0"],
                "CurrentChain": entry.get("current_chain"),
                "CurrentResidStart": entry.get("current_res_start"),
                "CurrentResidEnd": entry.get("current_res_end"),
                "SourceChain": entry.get("source_chain"),
                "SourceResidStart": entry.get("source_res_start"),
                "SourceResidEnd": entry.get("source_res_end"),
                "SourceResidCount": entry.get("source_rescount"),
            }
        )
    path = qc_dir / "protac_system_roles.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def write_ligand_detection_summary(ligand: mda.AtomGroup, ligand_code: str, qc_dir: Path) -> Path:
    lines = [
        f"Ligand resname: {ligand_code}",
        f"Atom count: {ligand.n_atoms}",
        f"Residue count: {ligand.residues.n_residues}",
    ]
    path = qc_dir / "protac_ligand_detection.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_labeling_audit(entries: Sequence[dict], qc_dir: Path) -> Path:
    rows = []
    for entry in entries:
        example_current = entry.get("current_res_start") or ""
        example_source = entry.get("source_res_start") or entry.get("current_res_start") or ""
        example_label = ""
        if example_source:
            example_label = f"{plot_display_name_for_entry(entry)}:RES{example_source}"
        rows.append(
            {
                "Role": entry["role"],
                "ProteinName": entry["name"],
                "DisplayName": plot_display_name_for_entry(entry),
                "AtomRange": f"{entry['start_1based']}-{entry['end_1based']}",
                "CurrentChain": entry.get("current_chain"),
                "CurrentResidueRange": (
                    f"{entry.get('current_res_start')}-{entry.get('current_res_end')}"
                    if entry.get("current_res_start") is not None and entry.get("current_res_end") is not None
                    else ""
                ),
                "SourceChain": entry.get("source_chain"),
                "SourceResidueRange": (
                    f"{entry.get('source_res_start')}-{entry.get('source_res_end')}"
                    if entry.get("source_res_start") is not None and entry.get("source_res_end") is not None
                    else ""
                ),
                "ExampleCurrentResid": example_current,
                "ExampleSourceResid": example_source,
                "ExampleDisplayLabel": example_label,
            }
        )
    path = qc_dir / "protac_labeling_audit.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def print_startup_summary(
    ligand_code: str,
    entries: Sequence[dict],
    topo_path: str,
    traj_path: str,
    frame_count: int,
) -> None:
    ligase = entries[0]
    targets = entries[1:]
    target_display = ", ".join(
        f"{entry['role']}={entry['name']}({entry['start_1based']}-{entry['end_1based']})"
        for entry in targets
    )
    target_display = target_display or "none"
    log("\nStartup summary")
    log(f"Ligand: {ligand_code}")
    log(f"Ligase: {ligase['name']} ({ligase['start_1based']}-{ligase['end_1based']})")
    log(f"Targets: {target_display}")
    log(f"Topology: {topo_path}")
    log(f"Trajectory: {traj_path}")
    log(f"Frames: {frame_count}")


def component_color_map(entries: Sequence[dict]) -> Dict[str, str]:
    palette = [
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#ff7f0e",
        "#8c564b",
        "#9467bd",
        "#17becf",
    ]
    colors = {"Full complex": "#7f7f7f"}
    for index, entry in enumerate(entries):
        colors[plot_display_name_for_entry(entry)] = palette[index % len(palette)]
    if has_multiple_targets(entries):
        colors["Combined target side"] = "#555555"
    return colors


def slugify(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("_") or "component"


def combine_atomgroups(groups: Sequence[mda.AtomGroup]) -> mda.AtomGroup:
    if not groups:
        raise ValueError("No atom groups were provided for combination.")
    universe = groups[0].universe
    indices = np.array(sorted(set(np.concatenate([group.indices for group in groups]).tolist())), dtype=int)
    return universe.atoms[indices]


def protein_atomgroup_from_entries(universe: mda.Universe, entries: Sequence[dict]) -> mda.AtomGroup:
    return combine_atomgroups([atomgroup_for_entry(universe, entry) for entry in entries])


def convert_ligand_records_to_hetatm(source_path: Path, target_path: Path, ligand_code: str) -> None:
    ligand_code = ligand_code.strip().upper()
    with source_path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()

    converted = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            resname = line[17:20].strip().upper()
            if resname == ligand_code:
                converted.append("HETATM" + line[6:])
                continue
        converted.append(line)

    target_path.write_text("".join(converted), encoding="utf-8")


def write_atomgroup_pdb(
    universe: mda.Universe,
    atomgroup: mda.AtomGroup,
    output_path: Path,
    ligand_code: str,
    frame: int,
) -> None:
    ensure_dir(output_path.parent)
    temp_path = output_path.with_suffix(".raw.pdb")
    universe.trajectory[frame]
    atomgroup.write(str(temp_path))
    convert_ligand_records_to_hetatm(temp_path, output_path, ligand_code)
    temp_path.unlink(missing_ok=True)


def write_atomgroup_xtc(
    universe: mda.Universe,
    atomgroup: mda.AtomGroup,
    output_path: Path,
    frame_indices: Sequence[int],
) -> None:
    ensure_dir(output_path.parent)
    with mda.Writer(str(output_path), atomgroup.n_atoms) as writer:
        for _, frame in iter_with_progress(frame_indices, f"Writing {output_path.name}"):
            universe.trajectory[frame]
            writer.write(atomgroup)


def residue_mapping_rows_for_output(output_name: str, residue_records: Sequence[dict]) -> List[dict]:
    rows = []
    for output_resid, record in enumerate(
        sorted(
            residue_records,
            key=lambda row: (row["component_order"], row["source_resid"], row["current_resid"], row["residue_name"]),
        ),
        start=1,
    ):
        rows.append(
            {
                "OutputFile": output_name,
                "Role": record["role"],
                "ProteinName": record["protein_name"],
                "OutputResid": output_resid,
                "CurrentResid": record["current_resid"],
                "SourceResid": record["source_resid"],
                "ResidueName": record["residue_name"],
                "ResidueLabel_Display": record["residue_label_display"],
            }
        )
    return rows


def write_ndx(groups: Dict[str, np.ndarray], output_path: Path) -> Path:
    ensure_dir(output_path.parent)
    with output_path.open("w", encoding="utf-8") as handle:
        for name, indices in groups.items():
            handle.write(f"[ {name} ]\n")
            one_based = [str(int(index) + 1) for index in np.array(indices, dtype=int)]
            for start in range(0, len(one_based), 15):
                handle.write(" ".join(one_based[start : start + 15]) + "\n")
            handle.write("\n")
    return output_path


def entry_for_atom_index(atom_index: int, entries: Sequence[dict]) -> Optional[dict]:
    for entry in entries:
        if int(entry["start0"]) <= atom_index <= int(entry["end0"]):
            return entry
    return None


def build_atom_mapping_for_subset(
    atomgroup: mda.AtomGroup,
    entries: Sequence[dict],
    ligand_code: str,
) -> pd.DataFrame:
    rows = []
    ligand_code = ligand_code.strip().upper()
    for output_index0, atom in enumerate(atomgroup.atoms):
        original_index0 = int(atom.index)
        entry = entry_for_atom_index(original_index0, entries)
        residue = atom.residue
        is_ligand = residue.resname.strip().upper() == ligand_code or entry is None
        role = entry["role"] if entry else "Ligand"
        protein_name = entry["name"] if entry else ligand_code
        display_name = plot_display_name_for_entry(entry) if entry else ligand_code
        current_chain = entry.get("current_chain") if entry else ""
        source_chain = entry.get("source_chain") if entry else current_chain
        current_resid = int(residue.resid)
        source_resid = get_display_residue_id(residue, entry) if entry else current_resid
        rows.append(
            {
                "OriginalAtomIndex0": original_index0,
                "OriginalAtomIndex1": original_index0 + 1,
                "OutputAtomIndex0": output_index0,
                "OutputAtomIndex1": output_index0 + 1,
                "Role": role,
                "ProteinName": protein_name,
                "DisplayName": display_name,
                "CurrentChain": current_chain,
                "SourceChain": source_chain,
                "AtomName": atom.name,
                "ResidueName": residue.resname,
                "CurrentResid": current_resid,
                "SourceResid": source_resid,
                "IsLigand": bool(is_ligand),
            }
        )
    return pd.DataFrame(rows)


def remap_entries_for_subset(entries: Sequence[dict], mapping_df: pd.DataFrame) -> List[dict]:
    remapped = []
    protein_rows = mapping_df[~mapping_df["IsLigand"]].copy()
    for entry in entries:
        mask = (
            (protein_rows["OriginalAtomIndex0"] >= int(entry["start0"]))
            & (protein_rows["OriginalAtomIndex0"] <= int(entry["end0"]))
        )
        subset = protein_rows[mask].sort_values(by="OutputAtomIndex0", kind="stable")
        if subset.empty:
            raise ValueError(
                f"Component {entry['name']} lost all atoms during Final_Trajectory remapping. "
                "Check the centered subset generation and atomIndex.txt alignment."
            )
        enriched = dict(entry)
        enriched["original_start0"] = entry["start0"]
        enriched["original_end0"] = entry["end0"]
        enriched["original_start_1based"] = entry["start_1based"]
        enriched["original_end_1based"] = entry["end_1based"]
        enriched["start0"] = int(subset["OutputAtomIndex0"].min())
        enriched["end0"] = int(subset["OutputAtomIndex0"].max())
        enriched["start_1based"] = int(subset["OutputAtomIndex1"].min())
        enriched["end_1based"] = int(subset["OutputAtomIndex1"].max())
        remapped.append(enriched)
    return remapped


def write_final_mapping_audit(entries: Sequence[dict], remapped_entries: Sequence[dict], qc_dir: Path) -> Path:
    rows = []
    for raw_entry, mapped_entry in zip(entries, remapped_entries):
        rows.append(
            {
                "Role": raw_entry["role"],
                "ProteinName": raw_entry["name"],
                "OriginalAtomRange": f"{raw_entry['start_1based']}-{raw_entry['end_1based']}",
                "OutputAtomRange": f"{mapped_entry['start_1based']}-{mapped_entry['end_1based']}",
                "OriginalAtomCount": int(raw_entry["end0"]) - int(raw_entry["start0"]) + 1,
                "OutputAtomCount": int(mapped_entry["end0"]) - int(mapped_entry["start0"]) + 1,
                "CurrentResidueRange": (
                    f"{raw_entry.get('current_res_start')}-{raw_entry.get('current_res_end')}"
                    if raw_entry.get("current_res_start") is not None and raw_entry.get("current_res_end") is not None
                    else ""
                ),
                "SourceResidueRange": (
                    f"{raw_entry.get('source_res_start')}-{raw_entry.get('source_res_end')}"
                    if raw_entry.get("source_res_start") is not None and raw_entry.get("source_res_end") is not None
                    else ""
                ),
            }
        )
    path = qc_dir / "protac_final_trajectory_mapping_audit.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def compute_com_summary(universe: mda.Universe, atomgroup: mda.AtomGroup, frame_indices: Sequence[int], label: str) -> dict:
    if atomgroup.n_atoms == 0:
        return {"Label": label, "Frames": 0, "XRange": 0.0, "YRange": 0.0, "ZRange": 0.0, "Displacement": 0.0}
    coords = []
    for _, frame in iter_with_progress(frame_indices, f"COM drift ({label})"):
        universe.trajectory[frame]
        coords.append(atomgroup.center_of_mass())
    array = np.asarray(coords, dtype=float)
    ranges = np.ptp(array, axis=0) if len(array) else np.zeros(3, dtype=float)
    displacement = float(np.linalg.norm(array[-1] - array[0])) if len(array) > 1 else 0.0
    return {
        "Label": label,
        "Frames": len(array),
        "XRange": float(ranges[0]),
        "YRange": float(ranges[1]),
        "ZRange": float(ranges[2]),
        "Displacement": displacement,
    }


def write_visualization_helpers(
    structures_dir: Path,
    final_pdb: Path,
    final_xtc: Optional[Path],
    pocket_pdb: Optional[Path],
    ligand_code: str,
    entries: Sequence[dict],
) -> Tuple[Path, Path]:
    pml_path = structures_dir / "load_protac_visualization.pml"
    readme_path = structures_dir / "README_visualization.txt"
    lines = [
        f'load "{final_pdb.name}", final_complex',
    ]
    if final_xtc and file_exists(final_xtc):
        lines.append(f'try: load_traj "{final_xtc.name}", final_complex')
    if pocket_pdb and file_exists(pocket_pdb):
        lines.append(f'load "{pocket_pdb.name}", binding_pocket')
    lines.extend(
        [
            f"select protac_ligand, resn {ligand_code}",
            "hide everything",
            "show cartoon, final_complex",
            "show sticks, protac_ligand",
            "color tv_orange, protac_ligand",
        ]
    )
    palette = ["marine", "forest", "slate", "teal", "olive", "deepblue"]
    for idx, entry in enumerate(entries):
        current_chain = entry.get("current_chain")
        selector = f"chain {current_chain}" if current_chain else f"resi {entry.get('current_res_start', '')}-{entry.get('current_res_end', '')}"
        lines.append(f"select component_{slugify(entry['name'])}, ({selector}) and final_complex")
        lines.append(f"color {palette[idx % len(palette)]}, component_{slugify(entry['name'])}")
    if pocket_pdb and file_exists(pocket_pdb):
        lines.extend(
            [
                "show sticks, binding_pocket and polymer.protein",
                "color yellow, binding_pocket and polymer.protein",
            ]
        )
    pml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    readme_path.write_text(
        "\n".join(
            [
                "PROTAC visualization helper files",
                "",
                f"{final_pdb.name} / {final_xtc.name if final_xtc else 'Final_Trajectory.xtc'}: centered no-water protein+PROTAC trajectory for visualization.",
                "binding_pocket_only.pdb / binding_pocket_only.xtc: reduced PROTAC binding pocket subset when generated.",
                "protac_structure_residue_mapping.csv: maps output subset residue numbering back to current/source numbering metadata.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return pml_path, readme_path


def generate_centered_final_trajectory(
    topo_path: str,
    traj_path: str,
    universe: mda.Universe,
    ligand: mda.AtomGroup,
    entries: Sequence[dict],
    structures_dir: Path,
    qc_dir: Path,
    ligand_code: str,
    frame_indices: Sequence[int],
    final_trajectory_name: str = "Final_Trajectory",
    write_xtc: bool = True,
    overwrite: bool = False,
    prefer_gromacs: bool = True,
) -> dict:
    ensure_dir(structures_dir)
    final_name = final_trajectory_name.strip() or "Final_Trajectory"
    cwd_final_pdb = Path.cwd() / f"{final_name}.pdb"
    cwd_final_xtc = Path.cwd() / f"{final_name}.xtc"
    structures_final_pdb = structures_dir / f"{final_name}.pdb"
    structures_final_xtc = structures_dir / f"{final_name}.xtc"
    protac_complex_pdb = structures_dir / "protac_complex_nowater.pdb"
    protac_complex_xtc = structures_dir / "protac_complex_nowater.xtc"
    mapping_csv = structures_dir / "final_trajectory_atom_mapping.csv"

    need_xtc = write_xtc
    have_reusable = file_exists(cwd_final_pdb) and (file_exists(cwd_final_xtc) if need_xtc else True) and file_exists(mapping_csv)
    protein_group = protein_atomgroup_from_entries(universe, entries)
    complex_group = combine_atomgroups([protein_group, ligand])

    groups = build_protac_structure_groups(universe, ligand, entries, None)
    ndx_path = write_ndx(groups, structures_dir / "protac_analysis_groups.ndx")

    warnings: List[str] = []
    status = "reused" if have_reusable and not overwrite else "created"
    method = "reused"

    if not have_reusable or overwrite:
        method = "mdanalysis"
        gmx_available = GMX_BIN is not None
        tpr_candidates = [Path("md_0_1.tpr"), Path(traj_path).with_suffix(".tpr"), Path(topo_path).with_suffix(".tpr")]
        tpr_path = next((candidate for candidate in tpr_candidates if candidate.exists()), None)
        raw_pdb = structures_dir / f"{final_name}_RAW.pdb"

        if prefer_gromacs and gmx_available and tpr_path is not None and write_xtc:
            log("   ⏳ Attempting GROMACS-centered Final_Trajectory generation")
            method = "gromacs"
            rc_xtc = stream_subprocess(
                gmx_argv(
                    "trjconv",
                    [
                        "-s",
                        str(tpr_path),
                        "-f",
                        str(Path(traj_path).resolve()),
                        "-o",
                        str(cwd_final_xtc),
                        "-n",
                        str(ndx_path),
                        "-pbc",
                        "mol",
                        "-center",
                    ],
                ),
                input_text="PROTAC_Complex_NoWater\nPROTAC_Complex_NoWater\n",
            )
            rc_pdb = stream_subprocess(
                gmx_argv(
                    "trjconv",
                    [
                        "-s",
                        str(tpr_path),
                        "-f",
                        str(cwd_final_xtc),
                        "-o",
                        str(raw_pdb),
                        "-n",
                        str(ndx_path),
                        "-dump",
                        "0",
                    ],
                ),
                input_text="PROTAC_Complex_NoWater\n",
            ) if rc_xtc == 0 else 1

            if rc_xtc == 0 and rc_pdb == 0 and file_exists(cwd_final_xtc) and file_exists(raw_pdb):
                convert_ligand_records_to_hetatm(raw_pdb, cwd_final_pdb, ligand_code)
                raw_pdb.unlink(missing_ok=True)
            else:
                warnings.append("GROMACS Final_Trajectory centering failed; falling back to MDAnalysis subset writing.")
                method = "mdanalysis"
                cwd_final_xtc.unlink(missing_ok=True)
                raw_pdb.unlink(missing_ok=True)

        if method == "mdanalysis":
            log("   ⏳ Writing Final_Trajectory via MDAnalysis fallback")
            write_atomgroup_pdb(universe, complex_group, cwd_final_pdb, ligand_code, frame_indices[0])
            if write_xtc:
                try:
                    write_atomgroup_xtc(universe, complex_group, cwd_final_xtc, frame_indices)
                except Exception as exc:
                    warnings.append(f"MDAnalysis XTC writing failed for Final_Trajectory; PDB was still written: {exc}")
                    cwd_final_xtc.unlink(missing_ok=True)

        copy_if_needed(cwd_final_pdb, structures_final_pdb, overwrite=True)
        if file_exists(cwd_final_xtc):
            copy_if_needed(cwd_final_xtc, structures_final_xtc, overwrite=True)
        copy_if_needed(cwd_final_pdb, protac_complex_pdb, overwrite=True)
        if file_exists(cwd_final_xtc):
            copy_if_needed(cwd_final_xtc, protac_complex_xtc, overwrite=True)

        mapping_df = build_atom_mapping_for_subset(complex_group, entries, ligand_code)
        mapping_df.to_csv(mapping_csv, index=False)
    else:
        copy_if_needed(cwd_final_pdb, structures_final_pdb, overwrite=False)
        if need_xtc and file_exists(cwd_final_xtc):
            copy_if_needed(cwd_final_xtc, structures_final_xtc, overwrite=False)
        copy_if_needed(cwd_final_pdb, protac_complex_pdb, overwrite=False)
        if need_xtc and file_exists(cwd_final_xtc):
            copy_if_needed(cwd_final_xtc, protac_complex_xtc, overwrite=False)
        mapping_df = pd.read_csv(mapping_csv)

    remapped_entries = remap_entries_for_subset(entries, mapping_df)
    mapping_audit = write_final_mapping_audit(entries, remapped_entries, qc_dir)

    final_universe = mda.Universe(str(cwd_final_pdb), str(cwd_final_xtc)) if file_exists(cwd_final_xtc) else mda.Universe(str(cwd_final_pdb))
    final_ligand = select_ligand(final_universe, ligand_code)
    final_protein = protein_atomgroup_from_entries(final_universe, remapped_entries)
    if final_ligand.n_atoms == 0:
        raise ValueError("Final_Trajectory validation failed: ligand atoms were not found after centered subset generation.")
    for entry in remapped_entries:
        if atomgroup_for_entry(final_universe, entry).n_atoms == 0:
            raise ValueError(f"Final_Trajectory validation failed: component {entry['name']} has zero atoms after remapping.")

    final_resnames = {res.resname.strip().upper() for res in final_universe.residues}
    forbidden = sorted(final_resnames & (WATER_NAMES | ION_NAMES | SOLVENT_LIKE_NAMES))
    if forbidden:
        warnings.append(f"Final_Trajectory still contains excluded solvent/ion residue names: {', '.join(forbidden)}")

    raw_com = compute_com_summary(universe, complex_group, frame_indices, "RawComplex")
    final_frame_indices = list(range(min(len(frame_indices), final_universe.trajectory.n_frames)))
    final_com = compute_com_summary(final_universe, combine_atomgroups([final_protein, final_ligand]), final_frame_indices, "FinalTrajectory")
    centering_qc_path = qc_dir / "protac_centering_qc.csv"
    pd.DataFrame([raw_com, final_com]).to_csv(centering_qc_path, index=False)

    structure_rows = [
        {
            "OutputName": final_name,
            "PDBPath": str(structures_final_pdb.resolve()),
            "XTCCreated": file_exists(structures_final_xtc),
            "XTCPath": str(structures_final_xtc.resolve()) if file_exists(structures_final_xtc) else "",
            "ResidueCount": final_protein.residues.n_residues,
            "Status": status,
            "Method": method,
        }
    ]

    return {
        "status": status,
        "method": method,
        "warnings": warnings,
        "cwd_final_pdb": cwd_final_pdb,
        "cwd_final_xtc": cwd_final_xtc,
        "structures_final_pdb": structures_final_pdb,
        "structures_final_xtc": structures_final_xtc,
        "protac_complex_pdb": protac_complex_pdb,
        "protac_complex_xtc": protac_complex_xtc,
        "mapping_csv": mapping_csv,
        "mapping_audit": mapping_audit,
        "centering_qc": centering_qc_path,
        "ndx_path": ndx_path,
        "remapped_entries": remapped_entries,
        "structure_rows": structure_rows,
    }


def compute_residue_distance_summary(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    ligand: mda.AtomGroup,
    residue_records: Sequence[dict],
    frame_indices: Sequence[int],
    cutoff: float,
    mode: str,
    persistent_fraction: float,
) -> pd.DataFrame:
    residue_local_indices = [record["protein_local_indices"] for record in residue_records]

    if mode == "first_frame":
        iter_frames = [frame_indices[0]]
    else:
        iter_frames = list(frame_indices)

    hit_counts = np.zeros(len(residue_records), dtype=int)
    min_observed = np.full(len(residue_records), np.inf, dtype=float)
    total_frames = max(1, len(iter_frames))

    for _, frame in iter_with_progress(iter_frames, f"Pocket scan ({mode})"):
        universe.trajectory[frame]
        dist_matrix = distance_array(ligand.positions, protein.positions)
        residue_min = np.array(
            [float(np.min(dist_matrix[:, local_indices])) for local_indices in residue_local_indices],
            dtype=float,
        )
        min_observed = np.minimum(min_observed, residue_min)
        hit_counts += (residue_min <= cutoff).astype(int)

    fractions = hit_counts / float(total_frames)
    rows = []
    for index, record in enumerate(residue_records):
        rows.append(
            {
                "Role": record["role"],
                "ProteinName": record["protein_name"],
                "DisplayName": record["display_name"],
                "CurrentChain": record["current_chain"],
                "CurrentResid": record["current_resid"],
                "SourceChain": record["source_chain"],
                "SourceResid": record["source_resid"],
                "ResidueName": record["residue_name"],
                "ResidueLabel_Current": record["residue_label_current"],
                "ResidueLabel_Source": record["residue_label_source"],
                "ResidueLabel_Display": record["residue_label_display"],
                "MinDistanceObserved": float(min_observed[index]),
                "ContactFraction": float(fractions[index]),
                "ContactCount": int(hit_counts[index]),
                "FramesObserved": total_frames,
                "ComponentOrder": record["component_order"],
                "ResidueKey": record["residue_key"],
            }
        )

    df = pd.DataFrame(rows)
    if mode == "persistent":
        df = df[df["ContactFraction"] >= persistent_fraction].copy()
    else:
        df = df[df["ContactCount"] > 0].copy()
    return order_residue_summary_for_display(df)


def component_specs(entries: Sequence[dict], plot_allprotein: bool) -> List[dict]:
    specs = []
    ligase = entries[0]
    targets = target_entries(entries)
    multiple_targets = len(targets) > 1

    specs.append(
        {
            "component": "Ligase",
            "role": ligase["role"],
            "protein_name": ligase["name"],
            "display_name": plot_display_name_for_entry(ligase),
            "entries": [ligase],
            "user_plot": True,
        }
    )

    if targets:
        target_side_name = " + ".join(plot_display_name_for_entry(entry) for entry in targets)
        specs.append(
            {
                "component": "TargetSide",
                "role": "TargetSide",
                "protein_name": "TargetSide",
                "display_name": "Target side" if multiple_targets else plot_display_name_for_entry(targets[0]),
                "entries": list(targets),
                "user_plot": multiple_targets,
            }
        )

    for entry in targets:
        specs.append(
            {
                "component": entry["name"],
                "role": entry["role"],
                "protein_name": entry["name"],
                "display_name": plot_display_name_for_entry(entry),
                "entries": [entry],
                "user_plot": True,
            }
        )

    specs.append(
        {
            "component": "AllProtein",
            "role": "All_Protein",
            "protein_name": "All_Protein",
            "display_name": "Full complex",
            "entries": list(entries),
            "user_plot": plot_allprotein,
        }
    )
    return specs


def build_protac_structure_groups(
    universe: mda.Universe,
    ligand: mda.AtomGroup,
    entries: Sequence[dict],
    binding_pocket_group: Optional[mda.AtomGroup],
) -> Dict[str, np.ndarray]:
    groups: Dict[str, np.ndarray] = {}
    protein_group = protein_atomgroup_from_entries(universe, entries)
    groups["PROTAC_Ligand"] = ligand.indices
    groups["PROTAC_AllProtein"] = protein_group.indices
    groups["PROTAC_Complex_NoWater"] = combine_atomgroups([protein_group, ligand]).indices

    ligase_group = atomgroup_for_entry(universe, entries[0])
    groups["PROTAC_Ligase"] = ligase_group.indices

    target_groups = [atomgroup_for_entry(universe, entry) for entry in target_entries(entries)]
    if target_groups:
        groups["PROTAC_TargetSide"] = combine_atomgroups(target_groups).indices

    for entry in entries:
        pname = slugify(plot_display_name_for_entry(entry))
        pgroup = atomgroup_for_entry(universe, entry)
        groups[f"PROTAC_{pname}"] = pgroup.indices
        groups[f"PROTAC_{pname}_Complex"] = combine_atomgroups([pgroup, ligand]).indices

    if binding_pocket_group is not None and binding_pocket_group.n_atoms > 0:
        groups["PROTAC_BindingPocket"] = binding_pocket_group.indices
    return groups


def write_preprocessing_summary(
    structures_dir: Path,
    structure_rows: Sequence[dict],
    pocket_df: pd.DataFrame,
    ligand_code: str,
    pocket_mode: str,
    pocket_cutoff: float,
) -> Tuple[Path, Path, Path]:
    structure_csv = structures_dir.parent / "QC" / "protac_structure_outputs.csv"
    binding_csv = structures_dir.parent / "QC" / "protac_binding_pocket_summary.csv"
    summary_txt = structures_dir.parent / "QC" / "protac_preprocessing_summary.txt"

    pd.DataFrame(structure_rows).to_csv(structure_csv, index=False)
    pocket_df.to_csv(binding_csv, index=False)

    lines = [
        "PROTAC preprocessing summary",
        f"Ligand: {ligand_code}",
        f"Pocket mode: {pocket_mode}",
        f"Pocket cutoff (A): {pocket_cutoff}",
        f"Structure outputs: {len(structure_rows)}",
        f"Pocket residues: {len(pocket_df)}",
    ]
    summary_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return structure_csv, binding_csv, summary_txt


def run_preprocessing(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    ligand: mda.AtomGroup,
    entries: Sequence[dict],
    residue_records: Sequence[dict],
    frame_indices: Sequence[int],
    structures_dir: Path,
    qc_dir: Path,
    ligand_code: str,
    write_subset_trajectories: bool,
    pocket_cutoff: float,
    pocket_mode: str,
    pocket_min_fraction: float,
    overwrite_structures: bool = False,
    final_trajectory_name: str = "Final_Trajectory",
    initial_structure_rows: Optional[Sequence[dict]] = None,
) -> dict:
    ensure_dir(structures_dir)
    protein_group = protein
    complex_nowater = combine_atomgroups([protein_group, ligand])
    ligase_group = combine_atomgroups([atomgroup_for_entry(universe, entries[0]), ligand])

    target_only_entries = target_entries(entries)
    structure_rows = list(initial_structure_rows or [])
    mapping_rows = []

    def record_output(name: str, pdb_path: Path, xtc_path: Optional[Path], residue_subset: Sequence[dict]) -> None:
        structure_rows.append(
            {
                "OutputName": name,
                "PDBPath": str(pdb_path.resolve()),
                "XTCCreated": bool(xtc_path and xtc_path.exists()),
                "XTCPath": str(xtc_path.resolve()) if xtc_path and xtc_path.exists() else "",
                "ResidueCount": len(residue_subset),
            }
        )
        mapping_rows.extend(residue_mapping_rows_for_output(name, residue_subset))

    log("   ⏳ Determining binding-pocket residues")
    pocket_df = compute_residue_distance_summary(
        universe,
        protein_group,
        ligand,
        residue_records,
        frame_indices,
        cutoff=pocket_cutoff,
        mode=pocket_mode,
        persistent_fraction=pocket_min_fraction,
    )
    pocket_residue_keys = set(pocket_df["ResidueKey"].tolist()) if not pocket_df.empty else set()
    binding_pocket_records = [record for record in residue_records if record["residue_key"] in pocket_residue_keys]

    binding_pocket_group = None
    if binding_pocket_records:
        pocket_atomgroups = [record["residue"].atoms for record in binding_pocket_records]
        binding_pocket_group = combine_atomgroups([combine_atomgroups(pocket_atomgroups), ligand])

    subset_definitions = [
        (
            "protac_complex_nowater",
            complex_nowater,
            list(residue_records),
        ),
        (
            "protein_only_nowater",
            protein_group,
            list(residue_records),
        ),
        (
            "ligase_protac",
            ligase_group,
            [record for record in residue_records if record["role"] == entries[0]["role"]],
        ),
    ]

    if target_only_entries:
        target_side_group = combine_atomgroups([protein_atomgroup_from_entries(universe, target_only_entries), ligand])
        subset_definitions.append(
            (
                "targetside_protac",
                target_side_group,
                [record for record in residue_records if record["role"].startswith("Target_")],
            )
        )

    for entry in target_only_entries:
        pname = slugify(plot_display_name_for_entry(entry))
        subset_definitions.append(
            (
                f"target_{pname}_protac",
                combine_atomgroups([atomgroup_for_entry(universe, entry), ligand]),
                [record for record in residue_records if record["role"] == entry["role"]],
            )
        )

    if binding_pocket_group is not None and binding_pocket_group.n_atoms > 0:
        subset_definitions.append(
            (
                "binding_pocket_only",
                binding_pocket_group,
                list(binding_pocket_records),
            )
        )

    for output_name, atomgroup, residue_subset in subset_definitions:
        log(f"   ⏳ Writing subset structure: {output_name}")
        pdb_path = structures_dir / f"{output_name}.pdb"
        xtc_path = structures_dir / f"{output_name}.xtc"
        if overwrite_structures or not file_exists(pdb_path):
            write_atomgroup_pdb(universe, atomgroup, pdb_path, ligand_code, frame_indices[0])
        if write_subset_trajectories or output_name == "binding_pocket_only":
            if overwrite_structures or not file_exists(xtc_path):
                write_atomgroup_xtc(universe, atomgroup, xtc_path, frame_indices)
        record_output(output_name, pdb_path, xtc_path if xtc_path.exists() else None, residue_subset)

    pocket_residue_csv = structures_dir / "binding_pocket_residues.csv"
    if not pocket_df.empty:
        pocket_export = pocket_df.drop(columns=["ComponentOrder", "ResidueKey"], errors="ignore")
        pocket_export.to_csv(pocket_residue_csv, index=False)
    else:
        pd.DataFrame(columns=["Role", "ProteinName", "CurrentChain", "CurrentResid", "SourceChain", "SourceResid", "ResidueName", "ResidueLabel_Display", "MinDistanceObserved", "ContactFraction"]).to_csv(pocket_residue_csv, index=False)

    mapping_csv = structures_dir / "protac_structure_residue_mapping.csv"
    pd.DataFrame(mapping_rows).to_csv(mapping_csv, index=False)

    groups = build_protac_structure_groups(universe, ligand, entries, binding_pocket_group)
    ndx_path = write_ndx(groups, structures_dir / "protac_analysis_groups.ndx")

    binding_pocket_pdb = structures_dir / "binding_pocket_only.pdb"
    binding_pocket_xtc = structures_dir / "binding_pocket_only.xtc"
    final_pdb = structures_dir / f"{final_trajectory_name}.pdb"
    final_xtc = structures_dir / f"{final_trajectory_name}.xtc"
    if file_exists(binding_pocket_pdb):
        copy_if_needed(binding_pocket_pdb, Path.cwd() / "binding_pocket_only.pdb", overwrite=overwrite_structures)
    if file_exists(binding_pocket_xtc):
        copy_if_needed(binding_pocket_xtc, Path.cwd() / "binding_pocket_only.xtc", overwrite=overwrite_structures)
    pml_path, vis_readme = write_visualization_helpers(
        structures_dir,
        final_pdb if file_exists(final_pdb) else (structures_dir / "protac_complex_nowater.pdb"),
        final_xtc if file_exists(final_xtc) else None,
        binding_pocket_pdb if file_exists(binding_pocket_pdb) else None,
        ligand_code,
        entries,
    )

    structure_csv, pocket_summary_csv, preprocessing_txt = write_preprocessing_summary(
        structures_dir,
        structure_rows,
        pocket_df.drop(columns=["ComponentOrder", "ResidueKey"], errors="ignore"),
        ligand_code,
        pocket_mode,
        pocket_cutoff,
    )

    return {
        "binding_pocket_df": pocket_df.drop(columns=["ComponentOrder", "ResidueKey"], errors="ignore"),
        "binding_pocket_records": binding_pocket_records,
        "binding_pocket_group": binding_pocket_group,
        "ndx_path": ndx_path,
        "structure_csv": structure_csv,
        "pocket_summary_csv": pocket_summary_csv,
        "preprocessing_txt": preprocessing_txt,
        "mapping_csv": mapping_csv,
        "pocket_residue_csv": pocket_residue_csv,
        "visualization_pml": pml_path,
        "visualization_readme": vis_readme,
    }


def plot_title_contacts(ligand_code: str, entries: Sequence[dict]) -> str:
    names = [plot_display_name_for_entry(entry) for entry in entries]
    if len(names) == 2:
        return f"{ligand_code} contacts over time: {names[0]} vs {names[1]}"
    return f"{ligand_code} contacts over time by protein component"


def make_alias(source: Path, target: Path) -> None:
    if not source.exists():
        return
    ensure_dir(target.parent)
    shutil.copy2(source, target)


def register_standard_figure(
    registry: FigureRegistry,
    path: Path,
    title: str,
    category: str,
    ligand_code: str,
    entries: Sequence[dict],
    notes_hint: str,
    sort_order: int,
    include_in_pdf: bool = True,
) -> None:
    registry.add(
        figure_path=path,
        figure_title=title,
        figure_category=category,
        ligand_code=ligand_code,
        protein_components=[plot_display_name_for_entry(entry) for entry in entries],
        description_key=DESCRIPTION_KEYS[category],
        notes_hint=notes_hint,
        include_in_pdf=include_in_pdf,
        sort_order=sort_order,
    )


def merge_rmsf_contact_annotations(
    rmsf_df: pd.DataFrame,
    contact_frequency_df: Optional[pd.DataFrame],
    threshold: float,
) -> pd.DataFrame:
    if rmsf_df.empty:
        return rmsf_df
    merged = rmsf_df.copy()
    defaults = {
        "ContactFraction": 0.0,
        "ContactPercent": 0.0,
        "IsPersistentContact": False,
        "DominantInteractionType": "",
        "MeanMinDistance": math.nan,
        "MinDistanceObserved": math.nan,
    }
    if contact_frequency_df is None or contact_frequency_df.empty:
        for key, value in defaults.items():
            merged[key] = value
        return merged

    contact_subset = contact_frequency_df[
        [
            "Component",
            "CurrentResid",
            "SourceResid",
            "ResidueName",
            "ResidueLabel_Display",
            "ContactFraction",
            "ContactPercent",
            "DominantInteractionType",
            "MeanMinDistance",
            "MinDistanceObserved",
        ]
    ].copy()
    contact_subset["ComponentPreference"] = contact_subset["Component"].map(
        lambda value: 2 if value in {"AllProtein", "TargetSide"} else 0
    )
    contact_subset = contact_subset.sort_values(
        by=["ComponentPreference", "ContactFraction"],
        ascending=[True, False],
        kind="stable",
    ).drop_duplicates(
        subset=["CurrentResid", "SourceResid", "ResidueName", "ResidueLabel_Display"],
        keep="first",
    )
    contact_subset = contact_subset.drop(columns=["ComponentPreference"], errors="ignore")
    contact_subset["SourceResid"] = pd.to_numeric(contact_subset["SourceResid"], errors="coerce")
    contact_subset["CurrentResid"] = pd.to_numeric(contact_subset["CurrentResid"], errors="coerce")
    merged["SourceResid"] = pd.to_numeric(merged["SourceResid"], errors="coerce")
    merged["CurrentResid"] = pd.to_numeric(merged["CurrentResid"], errors="coerce")
    merged = merged.merge(
        contact_subset,
        how="left",
        on=["CurrentResid", "SourceResid", "ResidueName", "ResidueLabel_Display"],
        suffixes=("", "_contact"),
    )
    for key, value in defaults.items():
        if key in merged:
            merged[key] = merged[key].fillna(value)
        else:
            merged[key] = value
    merged["IsPersistentContact"] = merged["ContactFraction"] >= threshold
    return merged


def plot_rmsf_with_contact_overlay(
    data: pd.DataFrame,
    output_path: Path,
    title: str,
    colors: Dict[str, str],
    threshold: float,
    marker_mode: str,
    label_top_n: int,
    labeled_rows: List[dict],
    text_box: Optional[str] = None,
) -> None:
    ensure_plotting_available()
    plt.figure(figsize=(11, 5))
    texts = []
    for display_name, subset in data.groupby("DisplayName", sort=False):
        subset = subset.sort_values(by=["SourceResid", "CurrentResid"], kind="stable")
        x_values = subset["SourceResid"].fillna(subset["CurrentResid"])
        plt.plot(
            x_values,
            subset["ResidueRMSF"],
            label=display_name,
            linewidth=1.8,
            color=colors.get(display_name, "#444444"),
        )

        contacts = subset[subset["ContactFraction"] >= threshold].copy()
        if contacts.empty:
            continue
        marker_sizes = 40 + 200 * contacts["ContactFraction"].clip(0.0, 1.0)
        if marker_mode in {"scatter", "both"}:
            plt.scatter(
                contacts["SourceResid"].fillna(contacts["CurrentResid"]),
                contacts["ResidueRMSF"],
                s=marker_sizes,
                color=colors.get(display_name, "#444444"),
                alpha=0.65,
                edgecolors="black",
                linewidths=0.4,
                zorder=3,
            )
        if marker_mode in {"vlines", "both"}:
            for _, row in contacts.iterrows():
                x_value = row["SourceResid"] if pd.notna(row["SourceResid"]) else row["CurrentResid"]
                plt.vlines(
                    x=x_value,
                    ymin=0,
                    ymax=row["ResidueRMSF"],
                    colors=colors.get(display_name, "#444444"),
                    linestyles="dashed",
                    linewidth=1.0,
                    alpha=min(0.85, max(0.2, float(row["ContactFraction"]))),
                    zorder=1,
                )

        for _, row in contacts.nlargest(label_top_n, "ContactFraction").iterrows():
            x_value = row["SourceResid"] if pd.notna(row["SourceResid"]) else row["CurrentResid"]
            y_value = row["ResidueRMSF"]
            text = plt.text(x_value, y_value + 0.08, row["ResidueLabel_Display"], fontsize=8, color=colors.get(display_name, "#444444"))
            texts.append(text)
            labeled_rows.append(
                {
                    "Figure": output_path.name,
                    "DisplayName": display_name,
                    "ResidueLabel_Display": row["ResidueLabel_Display"],
                    "SourceResid": row["SourceResid"],
                    "CurrentResid": row["CurrentResid"],
                    "ResidueRMSF": row["ResidueRMSF"],
                    "ContactFraction": row["ContactFraction"],
                    "ContactPercent": row["ContactPercent"],
                    "DominantInteractionType": row["DominantInteractionType"],
                }
            )

    if adjust_text is not None and texts:
        adjust_text(texts, only_move={"points": "y", "texts": "y"})

    if text_box:
        plt.gca().text(
            1.02,
            0.98,
            text_box,
            transform=plt.gca().transAxes,
            fontsize=8,
            va="top",
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "#aaaaaa"},
        )
    plt.xlabel("Residue number")
    plt.ylabel("RMSF (Å)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def run_rmsd_and_rmsf(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    ligand: mda.AtomGroup,
    entries: Sequence[dict],
    residue_records: Sequence[dict],
    frame_indices: Sequence[int],
    rms_dir: Path,
    ligand_code: str,
    registry: FigureRegistry,
    include_rmsf: bool,
    contact_frequency_df: Optional[pd.DataFrame],
    rmsf_contact_threshold: float,
    rmsf_label_top_contacts: int,
    rmsf_contact_marker_mode: str,
) -> None:
    ensure_plotting_available()
    ensure_dir(rms_dir)
    colors = component_color_map(entries)

    log("   ⏳ Preparing RMSD/RMSF frame times...")
    times = []
    for _, frame in iter_with_progress(frame_indices, "Collecting frame times"):
        universe.trajectory[frame]
        times.append(float(universe.trajectory.time))

    rows_rmsd = []
    rows_rmsf = []
    residue_record_by_key = {record["residue_key"]: record for record in residue_records}

    protein_ca = protein.select_atoms("name CA")
    if protein_ca.n_atoms == 0:
        protein_ca = protein

    protein_positions = gather_aligned_positions(universe, protein_ca, frame_indices, "Aligning full-complex coordinates for RMSD/RMSF")
    protein_rmsd = rmsd_from_aligned_positions(protein_positions)
    for frame, time_ps, value in zip(frame_indices, times, protein_rmsd):
        rows_rmsd.append(
            {
                "Frame": frame,
                "Time_ps": time_ps,
                "Role": "All_Protein",
                "ProteinName": "All_Protein",
                "DisplayName": "Full complex",
                "RMSD": float(value),
            }
        )

    if include_rmsf:
        protein_atom_rmsf = rmsf_from_aligned_positions(protein_positions)
        cursor = 0
        for residue in protein_ca.residues:
            count = residue.atoms.n_atoms
            residue_rmsf = float(np.mean(protein_atom_rmsf[cursor : cursor + count])) if count else math.nan
            cursor += count
            record = residue_record_by_key.get(residue_first_atom_key(residue))
            if record:
                rows_rmsf.append(
                    {
                        "Role": "All_Protein",
                        "ProteinName": "All_Protein",
                        "DisplayName": "Full complex",
                        "CurrentChain": record["current_chain"],
                        "CurrentResid": record["current_resid"],
                        "SourceChain": record["source_chain"],
                        "SourceResid": record["source_resid"],
                        "ResidueName": record["residue_name"],
                        "ResidueLabel_Current": record["residue_label_current"],
                        "ResidueLabel_Source": record["residue_label_source"],
                        "ResidueLabel_Display": record["residue_label_display"],
                        "ResidueRMSF": residue_rmsf,
                    }
                )

    rmsd_series = {"Full complex": protein_rmsd}

    for entry in entries:
        entry_group = atomgroup_for_entry(universe, entry).select_atoms("protein and name CA")
        if entry_group.n_atoms == 0:
            entry_group = atomgroup_for_entry(universe, entry).select_atoms("protein")
        if entry_group.n_atoms == 0:
            continue
        entry_name = plot_display_name_for_entry(entry)
        entry_positions = gather_aligned_positions(universe, entry_group, frame_indices, f"Aligning {entry_name} coordinates for RMSD/RMSF")
        entry_rmsd = rmsd_from_aligned_positions(entry_positions)
        rmsd_series[entry_name] = entry_rmsd
        for frame, time_ps, value in zip(frame_indices, times, entry_rmsd):
            rows_rmsd.append(
                {
                    "Frame": frame,
                    "Time_ps": time_ps,
                    "Role": entry["role"],
                    "ProteinName": entry["name"],
                    "DisplayName": entry_name,
                    "RMSD": float(value),
                }
            )

        if include_rmsf:
            entry_atom_rmsf = rmsf_from_aligned_positions(entry_positions)
            cursor = 0
            for residue in entry_group.residues:
                count = residue.atoms.n_atoms
                residue_rmsf = float(np.mean(entry_atom_rmsf[cursor : cursor + count])) if count else math.nan
                cursor += count
                record = residue_record_by_key.get(residue_first_atom_key(residue))
                if record:
                    rows_rmsf.append(
                        {
                            "Role": entry["role"],
                            "ProteinName": entry["name"],
                            "DisplayName": entry_name,
                            "CurrentChain": record["current_chain"],
                            "CurrentResid": record["current_resid"],
                            "SourceChain": record["source_chain"],
                            "SourceResid": record["source_resid"],
                            "ResidueName": record["residue_name"],
                            "ResidueLabel_Current": record["residue_label_current"],
                            "ResidueLabel_Source": record["residue_label_source"],
                            "ResidueLabel_Display": record["residue_label_display"],
                            "ResidueRMSF": residue_rmsf,
                        }
                    )

    ligand_positions = gather_aligned_positions(universe, ligand, frame_indices, "Aligning ligand coordinates for RMSD/RMSF")
    ligand_rmsd = rmsd_from_aligned_positions(ligand_positions)
    ligand_rmsd_df = pd.DataFrame(
        {
            "Frame": frame_indices,
            "Time_ps": times,
            "Role": "Ligand",
            "ProteinName": "Ligand",
            "DisplayName": ligand_code,
            "RMSD": ligand_rmsd,
        }
    )
    ligand_rmsd_df.to_csv(rms_dir / "ligand_rmsd.csv", index=False)

    rmsd_df = pd.DataFrame(rows_rmsd).sort_values(by=["Frame", "DisplayName"], kind="stable")
    rmsd_df.to_csv(rms_dir / "protein_rmsd_by_component.csv", index=False)

    plt.figure(figsize=(10, 5))
    for label, values in rmsd_series.items():
        alpha = 0.45 if label == "Full complex" else 1.0
        linestyle = "--" if label == "Full complex" else "-"
        plt.plot(frame_indices, values, label=label, color=colors.get(label, "#444444"), linewidth=2, alpha=alpha, linestyle=linestyle)
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title("Protein RMSD by protein component")
    plt.legend()
    plt.tight_layout()
    rmsd_plot_path = rms_dir / "protein_rmsd_by_component.png"
    plt.savefig(rmsd_plot_path, dpi=300)
    plt.close()

    make_alias(rmsd_plot_path, rms_dir / "Protein_RMSD.png")
    register_standard_figure(
        registry,
        rmsd_plot_path,
        "Protein RMSD by protein component",
        "protein_rmsd",
        ligand_code,
        entries,
        "Full complex plus component-resolved protein RMSD.",
        10,
    )

    plt.figure(figsize=(9, 5))
    plt.plot(frame_indices, ligand_rmsd, color="#9467bd", linewidth=2)
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title(f"{ligand_code} RMSD")
    plt.tight_layout()
    ligand_rmsd_path = rms_dir / "ligand_rmsd.png"
    plt.savefig(ligand_rmsd_path, dpi=300)
    plt.close()
    register_standard_figure(
        registry,
        ligand_rmsd_path,
        f"{ligand_code} RMSD",
        "ligand_rmsd",
        ligand_code,
        entries,
        "Ligand / PROTAC RMSD over time.",
        15,
    )

    if include_rmsf:
        ligand_rmsf = rmsf_from_aligned_positions(ligand_positions)
        pd.DataFrame(
            {
                "AtomIndex": np.arange(ligand.n_atoms),
                "AtomName": ligand.names,
                "RMSF": ligand_rmsf,
            }
        ).to_csv(rms_dir / "ligand_rmsf.csv", index=False)

        rmsf_df = pd.DataFrame(rows_rmsf).sort_values(
            by=["DisplayName", "SourceResid", "CurrentResid", "ResidueName"],
            kind="stable",
        )
        rmsf_df = merge_rmsf_contact_annotations(rmsf_df, contact_frequency_df, rmsf_contact_threshold)
        rmsf_df.to_csv(rms_dir / "protein_rmsf_by_component.csv", index=False)

        plt.figure(figsize=(10, 5))
        for display_name, subset in rmsf_df[rmsf_df["DisplayName"] != "Full complex"].groupby("DisplayName", sort=False):
            x_values = subset["SourceResid"].fillna(subset["CurrentResid"])
            plt.plot(x_values, subset["ResidueRMSF"], label=display_name, linewidth=1.8, color=colors.get(display_name, "#444444"))
        plt.xlabel("Residue number")
        plt.ylabel("RMSF (Å)")
        plt.title("Protein RMSF by protein component")
        plt.legend()
        plt.tight_layout()
        protein_rmsf_path = rms_dir / "protein_rmsf_by_component.png"
        plt.savefig(protein_rmsf_path, dpi=300)
        plt.close()
        register_standard_figure(
            registry,
            protein_rmsf_path,
            "Protein RMSF by protein component",
            "protein_rmsf",
            ligand_code,
            entries,
            "Component-resolved residue flexibility using display numbering.",
            20,
        )

        labeled_rows: List[dict] = []
        overlay_path = rms_dir / "protein_rmsf_by_component_contacts.png"
        overlay_subset = rmsf_df[rmsf_df["DisplayName"] != "Full complex"].copy()
        plot_rmsf_with_contact_overlay(
            overlay_subset,
            overlay_path,
            "Protein RMSF by component with PROTAC-contact overlays",
            colors,
            rmsf_contact_threshold,
            rmsf_contact_marker_mode,
            rmsf_label_top_contacts,
            labeled_rows,
        )
        register_standard_figure(
            registry,
            overlay_path,
            "Protein RMSF by component with PROTAC-contact overlays",
            "protein_rmsf_contacts",
            ligand_code,
            entries,
            "RMSF with persistent PROTAC-contact overlays for component residues.",
            22,
        )

        for entry in entries:
            display_name = plot_display_name_for_entry(entry)
            entry_subset = rmsf_df[rmsf_df["DisplayName"] == display_name].copy()
            if entry_subset.empty:
                continue
            top_contacts = entry_subset[entry_subset["ContactFraction"] >= rmsf_contact_threshold].nlargest(
                min(rmsf_label_top_contacts, max(1, len(entry_subset))),
                "ContactFraction",
            )
            text_box = (
                f"Persistent contact threshold: {rmsf_contact_threshold:.0%}\n"
                + "Top contact residues: "
                + (", ".join(top_contacts["ResidueLabel_Display"].tolist()) if not top_contacts.empty else "none")
            )
            component_path = rms_dir / f"protein_rmsf_contacts_{slugify(display_name)}.png"
            plot_rmsf_with_contact_overlay(
                entry_subset,
                component_path,
                f"{display_name} RMSF with PROTAC-contact overlays",
                colors,
                rmsf_contact_threshold,
                rmsf_contact_marker_mode,
                rmsf_label_top_contacts,
                labeled_rows,
                text_box=text_box,
            )
            register_standard_figure(
                registry,
                component_path,
                f"{display_name} RMSF with PROTAC-contact overlays",
                "component_rmsf_contacts",
                ligand_code,
                [entry],
                "Component RMSF with labeled persistent PROTAC-contact residues.",
                23 + int(entry.get("component_order", 0)),
                include_in_pdf=False,
            )

        labeled_df = pd.DataFrame(
            labeled_rows,
            columns=[
                "Figure",
                "DisplayName",
                "ResidueLabel_Display",
                "SourceResid",
                "CurrentResid",
                "ResidueRMSF",
                "ContactFraction",
                "ContactPercent",
                "DominantInteractionType",
            ],
        ).drop_duplicates(subset=["Figure", "ResidueLabel_Display"], keep="first")
        labeled_df.to_csv(rms_dir / "rmsf_contact_labeled_residues.csv", index=False)

        plt.figure(figsize=(9, 5))
        plt.plot(np.arange(ligand.n_atoms), ligand_rmsf, color="#9467bd", linewidth=1.8)
        plt.xlabel("Ligand atom index")
        plt.ylabel("RMSF (Å)")
        plt.title(f"{ligand_code} RMSF")
        plt.tight_layout()
        ligand_rmsf_path = rms_dir / "ligand_rmsf.png"
        plt.savefig(ligand_rmsf_path, dpi=300)
        plt.close()
        register_standard_figure(
            registry,
            ligand_rmsf_path,
            f"{ligand_code} RMSF",
            "ligand_rmsf",
            ligand_code,
            entries,
            "Ligand / PROTAC per-atom RMSF.",
            25,
        )


def run_radius_of_gyration(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    entries: Sequence[dict],
    frame_indices: Sequence[int],
    rms_dir: Path,
    geometry_dir: Path,
    ligand_code: str,
    registry: FigureRegistry,
) -> None:
    ensure_plotting_available()
    ensure_dir(rms_dir)
    rows = []
    for _, frame in iter_with_progress(frame_indices, "Computing radius of gyration"):
        universe.trajectory[frame]
        rows.append(
            {
                "Frame": frame,
                "Time_ps": float(universe.trajectory.time),
                "Role": "All_Protein",
                "ProteinName": "All_Protein",
                "DisplayName": "Full complex",
                "RadiusOfGyration": float(protein.radius_of_gyration()),
            }
        )
        for entry in entries:
            entry_group = atomgroup_for_entry(universe, entry).select_atoms("protein")
            if entry_group.n_atoms == 0:
                continue
            rows.append(
                {
                    "Frame": frame,
                    "Time_ps": float(universe.trajectory.time),
                    "Role": entry["role"],
                    "ProteinName": entry["name"],
                    "DisplayName": plot_display_name_for_entry(entry),
                    "RadiusOfGyration": float(entry_group.radius_of_gyration()),
                }
            )

    df = pd.DataFrame(rows).sort_values(by=["Frame", "DisplayName"], kind="stable")
    df.to_csv(rms_dir / "radius_of_gyration_by_component.csv", index=False)

    colors = component_color_map(entries)
    plt.figure(figsize=(10, 5))
    for display_name, subset in df.groupby("DisplayName", sort=False):
        alpha = 0.45 if display_name == "Full complex" else 1.0
        linestyle = "--" if display_name == "Full complex" else "-"
        plt.plot(subset["Frame"], subset["RadiusOfGyration"], label=display_name, linewidth=2, alpha=alpha, linestyle=linestyle, color=colors.get(display_name, "#444444"))
    plt.xlabel("Frame")
    plt.ylabel("Radius of gyration (Å)")
    plt.title("Radius of gyration by protein component")
    plt.legend()
    plt.tight_layout()
    rg_plot_path = rms_dir / "radius_of_gyration_by_component.png"
    plt.savefig(rg_plot_path, dpi=300)
    plt.close()

    make_alias(rg_plot_path, geometry_dir / "Radius_of_Gyration_Overlay.png")
    register_standard_figure(
        registry,
        rg_plot_path,
        "Radius of gyration by protein component",
        "radius_of_gyration",
        ligand_code,
        entries,
        "Full complex and per-component compactness traces.",
        30,
    )


def order_residue_summary_for_display(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    return df.sort_values(
        by=["ComponentOrder", "SourceResid", "CurrentResid", "ResidueName", "ResidueLabel_Display"],
        kind="stable",
    )


def write_contact_csvs(
    detail_df: pd.DataFrame,
    timeseries_df: pd.DataFrame,
    residue_summary_df: pd.DataFrame,
    contacts_dir: Path,
    qc_dir: Path,
) -> None:
    detail_df.to_csv(contacts_dir / "ligand_contacts_detailed.csv", index=False)
    timeseries_df.to_csv(contacts_dir / "ligand_contacts_by_component_timeseries.csv", index=False)
    residue_summary_df.to_csv(contacts_dir / "ligand_contact_residue_summary.csv", index=False)

    role_summary = (
        timeseries_df.groupby(["Role", "ProteinName", "DisplayName"], dropna=False)
        .agg(
            FramesWithContact=("HasContact", "sum"),
            ContactFraction=("HasContact", "mean"),
            ContactPercent=("HasContact", lambda s: float(np.mean(s)) * 100.0),
            MeanContactResidues=("ContactResidues", "mean"),
            MeanMinDistance=("MinDistance", "mean"),
        )
        .reset_index()
    )
    role_summary.to_csv(qc_dir / "protac_contact_summary.csv", index=False)


def run_contact_analysis(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    ligand: mda.AtomGroup,
    ligand_code: str,
    residue_records: Sequence[dict],
    entries: Sequence[dict],
    frame_indices: Sequence[int],
    distance_cutoff: float,
    contacts_dir: Path,
    networks_dir: Path,
    geometry_dir: Path,
    qc_dir: Path,
    registry: FigureRegistry,
    numbering: str,
    generate_contact_map: bool,
    generate_network: bool,
    interaction_mode: str,
    plot_allprotein: bool,
    network_min_contact_fraction: float,
    network_max_residues: int,
    contact_map_min_contact_fraction: float,
    contact_map_max_residues: int,
) -> None:
    ensure_plotting_available()
    if generate_network:
        ensure_networkx_available()
    ensure_dir(contacts_dir)
    ensure_dir(networks_dir)
    ensure_dir(geometry_dir)
    log(f"   ⏳ Preparing contact analysis structures for {len(frame_indices)} sampled frames")

    colors = component_color_map(entries)
    multiple_targets = has_multiple_targets(entries)
    specs = component_specs(entries, plot_allprotein=plot_allprotein)
    protein_local_lookup = {atom_index: local_index for local_index, atom_index in enumerate(protein.indices)}

    component_groups = []
    for spec in specs:
        entry_groups = [atomgroup_for_entry(universe, entry).select_atoms("protein") for entry in spec["entries"]]
        entry_groups = [group for group in entry_groups if group.n_atoms > 0]
        if not entry_groups:
            continue
        group = combine_atomgroups(entry_groups)
        local_indices = np.array(
            [protein_local_lookup[index] for index in group.indices if index in protein_local_lookup],
            dtype=int,
        )
        if local_indices.size == 0:
            continue
        enriched = dict(spec)
        enriched["local_indices"] = local_indices
        component_groups.append(enriched)

    residue_local_indices = [record["protein_local_indices"] for record in residue_records]
    residue_objects = [record["residue"] for record in residue_records]

    component_masks = {}
    for spec in component_groups:
        allowed_roles = {entry["role"] for entry in spec["entries"]}
        if spec["component"] == "AllProtein":
            component_masks[spec["component"]] = np.ones(len(residue_records), dtype=bool)
        else:
            component_masks[spec["component"]] = np.array(
                [record["role"] in allowed_roles for record in residue_records],
                dtype=bool,
            )

    detail_rows: List[dict] = []
    timeseries_rows: List[dict] = []
    residue_distance_history: Dict[str, Dict[str, List[float]]] = {
        spec["component"]: defaultdict(list) for spec in component_groups
    }

    for _, frame in iter_with_progress(frame_indices, "Computing ligand contacts"):
        universe.trajectory[frame]
        dist_matrix = distance_array(ligand.positions, protein.positions)
        residue_min_distances = np.array(
            [float(np.min(dist_matrix[:, local_indices])) for local_indices in residue_local_indices],
            dtype=float,
        )
        residue_contact_mask = residue_min_distances <= distance_cutoff

        for spec in component_groups:
            component_name = spec["component"]
            submatrix = dist_matrix[:, spec["local_indices"]]
            min_distance = float(np.min(submatrix))
            contact_atom_pairs = int(np.count_nonzero(submatrix <= distance_cutoff))
            residue_mask = residue_contact_mask & component_masks[component_name]
            residue_count = int(np.count_nonzero(residue_mask))
            timeseries_rows.append(
                {
                    "Frame": frame,
                    "Time_ps": float(universe.trajectory.time),
                    "Component": component_name,
                    "Role": spec["role"],
                    "ProteinName": spec["protein_name"],
                    "DisplayName": spec["display_name"],
                    "CurrentChain": "",
                    "CurrentResid": "",
                    "SourceChain": "",
                    "SourceResid": "",
                    "ResidueName": "",
                    "ResidueLabel_Current": "",
                    "ResidueLabel_Source": "",
                    "ResidueLabel_Display": "",
                    "MinDistance": min_distance,
                    "ContactBinary": int(contact_atom_pairs > 0),
                    "InteractionType": "Contact",
                    "ContactResidues": residue_count,
                    "ContactAtomPairs": contact_atom_pairs,
                }
            )

        for residue_index in np.where(residue_contact_mask)[0]:
            record = residue_records[residue_index]
            residue = residue_objects[residue_index]
            min_distance = float(residue_min_distances[residue_index])
            interaction_type = classify_interaction_type(residue, min_distance, mode=interaction_mode)

            for spec in component_groups:
                if not component_masks[spec["component"]][residue_index]:
                    continue
                residue_distance_history[spec["component"]][record["residue_label_display"]].append(min_distance)
                detail_rows.append(
                    {
                        "Frame": frame,
                        "Time_ps": float(universe.trajectory.time),
                        "Component": spec["component"],
                        "Role": record["role"],
                        "ProteinName": record["protein_name"],
                        "DisplayName": record["display_name"],
                        "CurrentChain": record["current_chain"],
                        "CurrentResid": record["current_resid"],
                        "SourceChain": record["source_chain"],
                        "SourceResid": record["source_resid"],
                        "ResidueName": record["residue_name"],
                        "ResidueLabel_Current": record["residue_label_current"],
                        "ResidueLabel_Source": record["residue_label_source"],
                        "ResidueLabel_Display": record["residue_label_display"],
                        "MinDistance": min_distance,
                        "ContactBinary": 1,
                        "InteractionType": interaction_type,
                        "ComponentOrder": record["component_order"],
                        "ContactFraction": np.nan,
                        "ContactPercent": np.nan,
                    }
                )

    detail_df = pd.DataFrame(detail_rows)
    timeseries_df = pd.DataFrame(timeseries_rows)
    if detail_df.empty:
        raise ValueError("No ligand-protein contacts were detected. Check the ligand residue name or distance cutoff.")

    total_sampled_frames = len(frame_indices)

    residue_frequency = (
        detail_df.groupby(
            [
                "Component",
                "Role",
                "ProteinName",
                "DisplayName",
                "CurrentChain",
                "CurrentResid",
                "SourceChain",
                "SourceResid",
                "ResidueName",
                "ResidueLabel_Current",
                "ResidueLabel_Source",
                "ResidueLabel_Display",
                "ComponentOrder",
            ],
            dropna=False,
        )
        .agg(
            ContactCount=("Frame", "nunique"),
            MeanMinDistance=("MinDistance", "mean"),
            MedianMinDistance=("MinDistance", "median"),
            MinDistanceObserved=("MinDistance", "min"),
            DominantInteractionType=("InteractionType", lambda vals: pd.Series(vals).mode().iloc[0] if len(pd.Series(vals).mode()) else "Contact"),
        )
        .reset_index()
    )
    residue_frequency["TotalSampledFrames"] = total_sampled_frames
    residue_frequency["ContactFraction"] = residue_frequency["ContactCount"] / total_sampled_frames
    residue_frequency["ContactPercent"] = residue_frequency["ContactFraction"] * 100.0
    residue_frequency = residue_frequency.sort_values(
        by=["Component", "ComponentOrder", "SourceResid", "CurrentResid", "ResidueLabel_Display"],
        kind="stable",
    )

    detail_df["ContactFraction"] = np.nan
    detail_df["ContactPercent"] = np.nan
    timeseries_df["ContactFraction"] = timeseries_df["ContactBinary"]
    timeseries_df["ContactPercent"] = timeseries_df["ContactBinary"] * 100.0

    protac_long_path = contacts_dir / "protac_contacts_long.csv"
    detail_df.to_csv(protac_long_path, index=False)
    frequency_path = contacts_dir / "protac_contact_frequency_by_residue.csv"
    residue_frequency.to_csv(frequency_path, index=False)

    component_summary_rows = []
    for spec in component_groups:
        component_name = spec["component"]
        ts_subset = timeseries_df[timeseries_df["Component"] == component_name]
        rf_subset = residue_frequency[residue_frequency["Component"] == component_name]
        top_residues = ", ".join(rf_subset.nlargest(5, "ContactFraction")["ResidueLabel_Display"].tolist())
        component_summary_rows.append(
            {
                "Component": component_name,
                "ProteinName": spec["protein_name"],
                "DisplayName": spec["display_name"],
                "TotalContacts": int(ts_subset["ContactAtomPairs"].sum()) if not ts_subset.empty else 0,
                "MeanContactsPerFrame": float(ts_subset["ContactResidues"].mean()) if not ts_subset.empty else 0.0,
                "MedianContactsPerFrame": float(ts_subset["ContactResidues"].median()) if not ts_subset.empty else 0.0,
                "PersistentResidueCount": int(np.count_nonzero(rf_subset["ContactFraction"] >= max(contact_map_min_contact_fraction, 0.1))) if not rf_subset.empty else 0,
                "TopResidues": top_residues,
            }
        )

    component_summary_df = pd.DataFrame(component_summary_rows)
    component_summary_path = qc_dir / "protac_contact_summary_by_component.csv"
    component_summary_df.to_csv(component_summary_path, index=False)

    interaction_summary_path = qc_dir / "protac_interaction_analysis_summary.csv"
    interaction_summary_df = component_summary_df.copy()
    interaction_summary_df["InteractionMode"] = interaction_mode
    interaction_summary_df["DistanceCutoff"] = distance_cutoff
    interaction_summary_df.to_csv(interaction_summary_path, index=False)

    log("   ⏳ Rendering main contact timeseries plot")
    plt.figure(figsize=(10, 5))
    plotted_any = False
    for spec in component_groups:
        if not spec["user_plot"]:
            continue
        subset = timeseries_df[timeseries_df["Component"] == spec["component"]]
        if subset.empty:
            continue
        plotted_any = True
        plt.plot(
            subset["Frame"],
            subset["ContactResidues"],
            label=spec["display_name"],
            color=colors.get(spec["display_name"], "#444444"),
            linewidth=2,
        )
    if not plotted_any:
        for spec in component_groups:
            subset = timeseries_df[timeseries_df["Component"] == spec["component"]]
            if subset.empty:
                continue
            plt.plot(subset["Frame"], subset["ContactResidues"], label=spec["display_name"], linewidth=2)
    plt.xlabel("Frame")
    plt.ylabel("Residues contacting ligand")
    plt.title(plot_title_contacts(ligand_code, entries))
    plt.legend()
    plt.tight_layout()
    main_timeseries_path = contacts_dir / "protac_contacts_by_component_timeseries.png"
    plt.savefig(main_timeseries_path, dpi=300)
    plt.close()
    register_standard_figure(
        registry,
        main_timeseries_path,
        plot_title_contacts(ligand_code, entries),
        "contacts_timeseries",
        ligand_code,
        entries,
        "Component-resolved PROTAC contact counts over time.",
        40,
    )

    top_frequency = residue_frequency[residue_frequency["Component"].isin([spec["component"] for spec in component_groups if spec["user_plot"]])].copy()
    top_frequency = top_frequency.nlargest(30, "ContactFraction").sort_values(
        by=["ContactFraction", "ComponentOrder", "SourceResid"],
        ascending=[False, True, True],
        kind="stable",
    )
    log("   ⏳ Rendering contact-frequency plot")
    plt.figure(figsize=(12, max(5, len(top_frequency) * 0.25)))
    bar_colors = [colors.get(row["DisplayName"], "#444444") for _, row in top_frequency.iterrows()]
    plt.barh(top_frequency["ResidueLabel_Display"], top_frequency["ContactPercent"], color=bar_colors)
    plt.xlabel("Contact percent")
    plt.ylabel("Residue")
    plt.title(f"{ligand_code} contact frequency by residue")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    frequency_plot_path = contacts_dir / "protac_contact_frequency_by_residue.png"
    plt.savefig(frequency_plot_path, dpi=300)
    plt.close()
    register_standard_figure(
        registry,
        frequency_plot_path,
        f"{ligand_code} contact frequency by residue",
        "contacts_timeseries",
        ligand_code,
        entries,
        "Top residues ranked by PROTAC contact persistence.",
        45,
    )

    if generate_contact_map:
        log("   ⏳ Rendering contact heatmaps")
        heatmap_specs = [("AllComponents", "protac_contact_heatmap_all_components.png", None)]
        for spec in component_groups:
            if spec["component"] in {"AllProtein"}:
                continue
            label = spec["component"] if spec["component"] != "TargetSide" else "TargetSide"
            heatmap_specs.append((spec["component"], f"protac_contact_heatmap_{slugify(label)}.png", spec["component"]))

        for key, filename, component_name in heatmap_specs:
            if component_name is None:
                subset = residue_frequency[residue_frequency["Component"].isin([spec["component"] for spec in component_groups if spec["user_plot"] or spec["component"] == "TargetSide"])].copy()
            else:
                subset = residue_frequency[residue_frequency["Component"] == component_name].copy()
            subset = subset[subset["ContactFraction"] >= contact_map_min_contact_fraction].copy()
            subset = subset.sort_values(by=["ComponentOrder", "SourceResid", "CurrentResid"], kind="stable")
            if len(subset) > contact_map_max_residues:
                subset = subset.nlargest(contact_map_max_residues, "ContactFraction").sort_values(
                    by=["ComponentOrder", "SourceResid", "CurrentResid"],
                    kind="stable",
                )
            if subset.empty:
                continue
            ordered_labels = subset["ResidueLabel_Display"].tolist()
            detail_subset = detail_df if component_name is None else detail_df[detail_df["Component"] == component_name]
            heatmap_df = (
                detail_subset[detail_subset["ResidueLabel_Display"].isin(ordered_labels)]
                .assign(ContactValue=1)
                .pivot_table(index="ResidueLabel_Display", columns="Frame", values="ContactValue", aggfunc="max", fill_value=0)
                .reindex(ordered_labels)
            )
            plt.figure(figsize=(12, max(5, len(heatmap_df) * 0.22)))
            sns.heatmap(heatmap_df, cmap="mako", cbar=True)
            plt.xlabel("Frame")
            plt.ylabel("Residue")
            title = f"{ligand_code} contact heatmap" if component_name is None else f"{ligand_code} contact heatmap: {subset['DisplayName'].iloc[0]}"
            plt.title(title)
            plt.tight_layout()
            output_path = contacts_dir / filename if component_name is None else contacts_dir / slugify(component_name) / filename
            ensure_dir(output_path.parent)
            plt.savefig(output_path, dpi=300)
            plt.close()
            if component_name is None:
                make_alias(output_path, contacts_dir / f"{ligand_code}_contact_map_residue_time.png")
                register_standard_figure(
                    registry,
                    output_path,
                    title,
                    "contact_heatmap",
                    ligand_code,
                    entries,
                    "Filtered contact heatmap ordered by component and biological residue numbering.",
                    50,
                )
    else:
        log("   ⏭️ Skipping contact heatmap by request")

    component_labels = [spec["display_name"] for spec in component_groups if spec["component"] not in {"AllProtein", "TargetSide"} or spec["user_plot"]]
    composition_df = (
        detail_df[detail_df["DisplayName"].isin(component_labels)]
        .groupby(["DisplayName", "InteractionType"], dropna=False)
        .size()
        .reset_index(name="Count")
    )
    composition_df.to_csv(geometry_dir / "ligand_interaction_type_summary.csv", index=False)
    composition_matrix = composition_df.pivot_table(index="DisplayName", columns="InteractionType", values="Count", fill_value=0).reindex(component_labels).fillna(0)
    composition_fraction = composition_matrix.div(composition_matrix.sum(axis=1), axis=0).fillna(0.0)
    composition_fraction.to_csv(geometry_dir / "protac_interaction_composition_by_component.csv")

    log("   ⏳ Rendering interaction-composition plot")
    ax = composition_fraction.plot(kind="bar", stacked=True, figsize=(10, 5), colormap="tab20")
    ax.set_ylabel("Fraction of residue-contact events")
    ax.set_xlabel("Protein component")
    ax.set_title(f"{ligand_code} interaction composition by component")
    ax.legend(title="Interaction type", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    composition_path = geometry_dir / "protac_interaction_composition_by_component.png"
    plt.savefig(composition_path, dpi=300)
    plt.close()
    register_standard_figure(
        registry,
        composition_path,
        f"{ligand_code} interaction composition by component",
        "interaction_composition",
        ligand_code,
        entries,
        "Broad chemistry-derived interaction classes across protein components.",
        60,
    )

    min_distance_candidates = residue_frequency[residue_frequency["Component"].isin([spec["component"] for spec in component_groups if spec["user_plot"]])].copy()
    min_distance_candidates = min_distance_candidates.nlargest(12, "ContactFraction")
    log("   ⏳ Rendering min-distance plot for top residues")
    plt.figure(figsize=(10, 5))
    for _, row in min_distance_candidates.iterrows():
        distances = residue_distance_history[row["Component"]].get(row["ResidueLabel_Display"], [])
        if not distances:
            continue
        plt.plot(
            range(len(distances)),
            distances,
            label=row["ResidueLabel_Display"],
            linewidth=1.3,
            color=colors.get(row["DisplayName"], "#444444"),
        )
    plt.xlabel("Contact event index")
    plt.ylabel("Min distance (Å)")
    plt.title(f"{ligand_code} minimum distances for top contacting residues")
    plt.legend(fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    min_distance_path = contacts_dir / "protac_min_distance_top_residues.png"
    plt.savefig(min_distance_path, dpi=300)
    plt.close()
    register_standard_figure(
        registry,
        min_distance_path,
        f"{ligand_code} minimum distances for top residues",
        "contacts_timeseries",
        ligand_code,
        entries,
        "Distance traces for top persistent residues.",
        65,
    )

    if generate_network:
        log("   ⏳ Building interaction networks")
        network_specs = [("AllComponents", None, "protac_interaction_network_all_components.png")]
        for spec in component_groups:
            if spec["component"] == "AllProtein":
                continue
            network_specs.append((spec["display_name"], spec["component"], f"protac_interaction_network_{slugify(spec['component'])}.png"))

        for display_label, component_name, filename in network_specs:
            subset = residue_frequency if component_name is None else residue_frequency[residue_frequency["Component"] == component_name]
            subset = subset[subset["ContactFraction"] >= network_min_contact_fraction].copy()
            subset = subset.nlargest(network_max_residues, "ContactFraction").sort_values(
                by=["ContactFraction", "ComponentOrder", "SourceResid"],
                ascending=[False, True, True],
                kind="stable",
            )
            if subset.empty:
                continue
            graph = nx.Graph()
            graph.add_node(ligand_code, role="Ligand", display_name=ligand_code)
            for _, row in subset.iterrows():
                node = row["ResidueLabel_Display"]
                graph.add_node(node, role=row["Role"], display_name=row["DisplayName"])
                graph.add_edge(ligand_code, node, weight=float(row["ContactPercent"]), interaction=row["DominantInteractionType"])

            positions = nx.spring_layout(graph, seed=42, weight="weight")
            positions[ligand_code] = np.array([0.0, 0.0])
            plt.figure(figsize=(12, 10))
            edge_widths = [max(1.0, graph[u][v]["weight"] / 12.0) for u, v in graph.edges()]
            node_sizes = [2600 if node == ligand_code else 400 + 15 * graph.degree(node) for node in graph.nodes()]
            node_colors = [
                "#9467bd" if node == ligand_code else colors.get(graph.nodes[node].get("display_name", ""), "#444444")
                for node in graph.nodes()
            ]
            nx.draw_networkx_edges(graph, positions, alpha=0.4, width=edge_widths, edge_color="#999999")
            nx.draw_networkx_nodes(graph, positions, node_color=node_colors, node_size=node_sizes, edgecolors="black")
            nx.draw_networkx_labels(graph, positions, font_size=8)
            title = f"{ligand_code} interaction network" if component_name is None else f"{ligand_code} interaction network: {display_label}"
            plt.title(title)
            plt.axis("off")
            plt.tight_layout()
            output_path = networks_dir / filename if component_name is None else networks_dir / slugify(component_name) / filename
            ensure_dir(output_path.parent)
            plt.savefig(output_path, dpi=300)
            plt.close()
            if component_name is None:
                register_standard_figure(
                    registry,
                    output_path,
                    title,
                    "interaction_network",
                    ligand_code,
                    entries,
                    "Filtered interaction network for top persistent residues.",
                    70,
                )
    else:
        log("   ⏭️ Skipping network generation by request")

    overview_candidates = residue_frequency[residue_frequency["Component"].isin([spec["component"] for spec in component_groups if spec["user_plot"]])].nlargest(20, "ContactFraction")
    if not overview_candidates.empty:
        log("   ⏳ Rendering binding-pocket overview plot")
        plt.figure(figsize=(10, 5))
        scatter = plt.scatter(
            overview_candidates["ContactPercent"],
            overview_candidates["MinDistanceObserved"],
            c=[colors.get(name, "#444444") for name in overview_candidates["DisplayName"]],
            s=60,
        )
        for _, row in overview_candidates.iterrows():
            plt.text(row["ContactPercent"], row["MinDistanceObserved"], row["ResidueLabel_Display"], fontsize=8)
        plt.xlabel("Contact percent")
        plt.ylabel("Minimum observed distance (Å)")
        plt.title(f"{ligand_code} binding-pocket overview")
        plt.tight_layout()
        overview_path = contacts_dir / "protac_binding_pocket_overview.png"
        plt.savefig(overview_path, dpi=300)
        plt.close()
        register_standard_figure(
            registry,
            overview_path,
            f"{ligand_code} binding-pocket overview",
            "contacts_timeseries",
            ligand_code,
            entries,
            "Overview of persistent pocket residues by contact percent and closest observed approach.",
            67,
            include_in_pdf=False,
        )

    for spec in component_groups:
        component_name = spec["component"]
        component_contacts_dir = contacts_dir / slugify(component_name)
        ensure_dir(component_contacts_dir)
        detail_df[detail_df["Component"] == component_name].to_csv(component_contacts_dir / "contacts_long.csv", index=False)
        residue_frequency[residue_frequency["Component"] == component_name].to_csv(component_contacts_dir / "contact_frequency_by_residue.csv", index=False)
        timeseries_df[timeseries_df["Component"] == component_name].to_csv(component_contacts_dir / "timeseries.csv", index=False)
        persistent_df = residue_frequency[
            (residue_frequency["Component"] == component_name) & (residue_frequency["ContactFraction"] >= max(contact_map_min_contact_fraction, 0.1))
        ]
        persistent_df.to_csv(component_contacts_dir / "persistent_binding_residues.csv", index=False)

    return {
        "detail_df": detail_df,
        "timeseries_df": timeseries_df,
        "residue_frequency_df": residue_frequency,
        "component_summary_df": component_summary_df,
        "component_summary_path": component_summary_path,
        "interaction_summary_path": interaction_summary_path,
    }


def find_water_oxygen_indices(universe: mda.Universe, water_names: Sequence[str]) -> np.ndarray:
    water_set = {name.strip().upper() for name in water_names if name.strip()}
    oxygen_indices: List[int] = []
    for residue in universe.residues:
        if residue.resname.strip().upper() not in water_set:
            continue
        oxygen_atoms = [atom.index for atom in residue.atoms if atom.name.strip().upper().startswith("O")]
        if not oxygen_atoms and residue.atoms.n_atoms:
            oxygen_atoms = [int(residue.atoms.indices[0])]
        oxygen_indices.extend(oxygen_atoms[:1])
    return np.array(sorted(set(oxygen_indices)), dtype=int)


def write_water_bridge_summary_txt(
    output_path: Path,
    water_names: Sequence[str],
    water_count: int,
    sampled_frames: int,
    cutoff: float,
    total_events: int,
    top_rows: pd.DataFrame,
) -> Path:
    lines = [
        "PROTAC water-bridge summary",
        f"Water names searched: {','.join(water_names)}",
        f"Water molecules in raw topology: {water_count}",
        f"Sampled frames: {sampled_frames}",
        f"Water-bridge cutoff (A): {cutoff}",
        f"Total bridge events: {total_events}",
        "Top water-bridged residues:",
    ]
    if top_rows.empty:
        lines.append("  none")
    else:
        for _, row in top_rows.iterrows():
            lines.append(
                f"  - {row['ResidueLabel_Display']}: {row['BridgePercent']:.1f}% ({int(row['BridgeCount'])} frames)"
            )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def run_water_bridge_analysis(
    raw_universe: mda.Universe,
    raw_protein: mda.AtomGroup,
    raw_ligand: mda.AtomGroup,
    raw_entries: Sequence[dict],
    raw_residue_records: Sequence[dict],
    frame_indices: Sequence[int],
    ligand_code: str,
    water_dir: Path,
    structures_dir: Path,
    qc_dir: Path,
    registry: FigureRegistry,
    water_names: Sequence[str],
    bridge_cutoff: float,
    water_pocket_cutoff: float,
    write_water_pocket_trajectory: bool,
    overwrite_structures: bool,
) -> dict:
    ensure_plotting_available()
    ensure_dir(water_dir)
    ensure_dir(structures_dir)
    colors = component_color_map(raw_entries)
    water_oxygen_indices = find_water_oxygen_indices(raw_universe, water_names)
    if water_oxygen_indices.size == 0:
        raise ValueError(
            f"No water oxygen atoms were found using water names: {', '.join(water_names)}. "
            "Adjust --water-names for this trajectory."
        )

    water_oxygens = raw_universe.atoms[water_oxygen_indices]
    protein_local_lookup = {atom_index: local_index for local_index, atom_index in enumerate(raw_protein.indices)}
    residue_local_indices = [record["protein_local_indices"] for record in raw_residue_records]
    component_specs_all = component_specs(raw_entries, plot_allprotein=False)
    component_masks = {}
    for spec in component_specs_all:
        allowed_roles = {entry["role"] for entry in spec["entries"]}
        if spec["component"] == "AllProtein":
            component_masks[spec["component"]] = np.ones(len(raw_residue_records), dtype=bool)
        else:
            component_masks[spec["component"]] = np.array(
                [record["role"] in allowed_roles for record in raw_residue_records],
                dtype=bool,
            )

    local_atom_to_residue = np.full(raw_protein.n_atoms, -1, dtype=int)
    for residue_index, local_indices in enumerate(residue_local_indices):
        local_atom_to_residue[local_indices] = residue_index

    detail_rows: List[dict] = []
    timeseries_rows: List[dict] = []
    bridge_residue_keys = set()
    pocket_residue_keys = set()
    water_residue_ids = set()
    pocket_water_residue_ids = set()

    for _, frame in iter_with_progress(frame_indices, "Computing water bridges"):
        raw_universe.trajectory[frame]
        ligand_water_dist = distance_array(raw_ligand.positions, water_oxygens.positions)
        ligand_water_min = ligand_water_dist.min(axis=0)
        close_water_indices = np.where(ligand_water_min <= bridge_cutoff)[0]
        pocket_water_indices = np.where(ligand_water_min <= water_pocket_cutoff)[0]
        for water_index in pocket_water_indices:
            pocket_water_residue_ids.add(int(water_oxygens[water_index].residue.ix))
        frame_component_counts = defaultdict(set)
        if close_water_indices.size == 0:
            for spec in component_specs_all:
                if spec["component"] == "AllProtein":
                    continue
                timeseries_rows.append(
                    {
                        "Frame": frame,
                        "Time_ps": float(raw_universe.trajectory.time),
                        "Component": spec["component"],
                        "ProteinName": spec["protein_name"],
                        "DisplayName": spec["display_name"],
                        "WaterBridgeResidues": 0,
                    }
                )
            continue

        close_waters = water_oxygens[close_water_indices]
        dist_wp = distance_array(close_waters.positions, raw_protein.positions)
        ligand_pocket_dist = distance_array(raw_ligand.positions, raw_protein.positions)
        residue_ligand_min = np.array(
            [float(np.min(ligand_pocket_dist[:, local_indices])) for local_indices in residue_local_indices],
            dtype=float,
        )
        for residue_index, distance_value in enumerate(residue_ligand_min):
            if distance_value <= water_pocket_cutoff:
                pocket_residue_keys.add(raw_residue_records[residue_index]["residue_key"])

        for water_offset, water_atom in enumerate(close_waters):
            contact_atom_locals = np.where(dist_wp[water_offset] <= bridge_cutoff)[0]
            if contact_atom_locals.size == 0:
                continue
            water_residue_ids.add(int(water_atom.residue.ix))
            contacted_residue_indices = np.unique(local_atom_to_residue[contact_atom_locals])
            contacted_residue_indices = contacted_residue_indices[contacted_residue_indices >= 0]
            ligand_water_distance = float(ligand_water_min[close_water_indices[water_offset]])
            for residue_index in contacted_residue_indices:
                record = raw_residue_records[int(residue_index)]
                bridge_residue_keys.add(record["residue_key"])
                frame_component_counts[record["display_name"]].add(record["residue_key"])
                if record["role"].startswith("Target_"):
                    frame_component_counts["TargetSide"].add(record["residue_key"])
                water_protein_distance = float(np.min(dist_wp[water_offset, residue_local_indices[int(residue_index)]]))
                detail_rows.append(
                    {
                        "Frame": frame,
                        "Time_ps": float(raw_universe.trajectory.time),
                        "WaterResid": int(water_atom.residue.resid),
                        "WaterResname": water_atom.residue.resname,
                        "Role": record["role"],
                        "ProteinName": record["protein_name"],
                        "DisplayName": record["display_name"],
                        "CurrentChain": record["current_chain"],
                        "CurrentResid": record["current_resid"],
                        "SourceChain": record["source_chain"],
                        "SourceResid": record["source_resid"],
                        "ResidueName": record["residue_name"],
                        "ResidueLabel_Current": record["residue_label_current"],
                        "ResidueLabel_Source": record["residue_label_source"],
                        "ResidueLabel_Display": record["residue_label_display"],
                        "LigandWaterDistance": ligand_water_distance,
                        "WaterProteinDistance": water_protein_distance,
                        "BridgeBinary": 1,
                        "ComponentOrder": record["component_order"],
                    }
                )

        for spec in component_specs_all:
            if spec["component"] == "AllProtein":
                continue
            if spec["component"] == "TargetSide":
                subset_count = len(frame_component_counts.get("TargetSide", set()))
            else:
                subset_count = len(frame_component_counts.get(spec["display_name"], set()))
            timeseries_rows.append(
                {
                    "Frame": frame,
                    "Time_ps": float(raw_universe.trajectory.time),
                    "Component": spec["component"],
                    "ProteinName": spec["protein_name"],
                    "DisplayName": spec["display_name"],
                    "WaterBridgeResidues": subset_count,
                }
            )

    detail_df = pd.DataFrame(detail_rows)
    timeseries_df = pd.DataFrame(timeseries_rows)
    if detail_df.empty:
        detail_df = pd.DataFrame(
            columns=[
                "Frame",
                "Time_ps",
                "WaterResid",
                "WaterResname",
                "Role",
                "ProteinName",
                "DisplayName",
                "CurrentChain",
                "CurrentResid",
                "SourceChain",
                "SourceResid",
                "ResidueName",
                "ResidueLabel_Current",
                "ResidueLabel_Source",
                "ResidueLabel_Display",
                "LigandWaterDistance",
                "WaterProteinDistance",
                "BridgeBinary",
            ]
        )

    long_csv = water_dir / "protac_water_bridges_long.csv"
    detail_df.to_csv(long_csv, index=False)

    total_sampled_frames = len(frame_indices)
    if detail_df.empty:
        summary_df = pd.DataFrame(
            columns=[
                "Role",
                "ProteinName",
                "ResidueLabel_Display",
                "SourceResid",
                "CurrentResid",
                "BridgeCount",
                "TotalSampledFrames",
                "BridgeFraction",
                "BridgePercent",
                "MeanLigandWaterDistance",
                "MeanWaterProteinDistance",
                "MinLigandWaterDistance",
                "MinWaterProteinDistance",
                "DisplayName",
                "ComponentOrder",
            ]
        )
    else:
        summary_df = (
            detail_df.groupby(
                [
                    "Role",
                    "ProteinName",
                    "DisplayName",
                    "CurrentResid",
                    "SourceResid",
                    "ResidueName",
                    "ResidueLabel_Current",
                    "ResidueLabel_Source",
                    "ResidueLabel_Display",
                    "ComponentOrder",
                ],
                dropna=False,
            )
            .agg(
                BridgeCount=("Frame", "nunique"),
                MeanLigandWaterDistance=("LigandWaterDistance", "mean"),
                MeanWaterProteinDistance=("WaterProteinDistance", "mean"),
                MinLigandWaterDistance=("LigandWaterDistance", "min"),
                MinWaterProteinDistance=("WaterProteinDistance", "min"),
            )
            .reset_index()
        )
        summary_df["TotalSampledFrames"] = total_sampled_frames
        summary_df["BridgeFraction"] = summary_df["BridgeCount"] / float(max(1, total_sampled_frames))
        summary_df["BridgePercent"] = summary_df["BridgeFraction"] * 100.0
        summary_df = summary_df.sort_values(
            by=["ComponentOrder", "SourceResid", "CurrentResid", "ResidueLabel_Display"],
            kind="stable",
        )
    frequency_csv = water_dir / "protac_water_bridge_frequency_by_residue.csv"
    summary_df.to_csv(frequency_csv, index=False)

    summary_component_rows = []
    for spec in component_specs_all:
        if spec["component"] == "AllProtein":
            continue
        if spec["component"] == "TargetSide":
            allowed_roles = {entry["role"] for entry in spec["entries"]}
            subset = summary_df[summary_df["Role"].isin(allowed_roles)].copy()
        else:
            subset = summary_df[summary_df["DisplayName"] == spec["display_name"]].copy()
        ts_subset = timeseries_df[timeseries_df["Component"] == spec["component"]].copy()
        top_residues = ", ".join(subset.nlargest(5, "BridgeFraction")["ResidueLabel_Display"].tolist()) if not subset.empty else ""
        summary_component_rows.append(
            {
                "Component": spec["component"],
                "ProteinName": spec["protein_name"],
                "TotalWaterBridgeEvents": int(subset["BridgeCount"].sum()) if not subset.empty else 0,
                "MeanWaterBridgesPerFrame": float(ts_subset["WaterBridgeResidues"].mean()) if not ts_subset.empty else 0.0,
                "PersistentWaterBridgeResidueCount": int(np.count_nonzero(subset["BridgeFraction"] >= 0.05)) if not subset.empty else 0,
                "TopWaterBridgeResidues": top_residues,
            }
        )
    component_summary_df = pd.DataFrame(summary_component_rows)
    component_summary_csv = water_dir / "protac_water_bridge_summary_by_component.csv"
    component_summary_df.to_csv(component_summary_csv, index=False)

    log("   ⏳ Rendering water-bridge timeseries")
    plt.figure(figsize=(10, 5))
    for spec in component_specs_all:
        if spec["component"] == "AllProtein" or not spec["user_plot"]:
            continue
        subset = timeseries_df[timeseries_df["Component"] == spec["component"]]
        if subset.empty:
            continue
        plt.plot(
            subset["Frame"],
            subset["WaterBridgeResidues"],
            label=spec["display_name"],
            linewidth=2,
            color=colors.get(spec["display_name"], "#444444"),
        )
    plt.xlabel("Frame")
    plt.ylabel("Water-bridged residues")
    plt.title(f"{ligand_code} water bridges over time")
    plt.legend()
    plt.tight_layout()
    timeseries_path = water_dir / "protac_water_bridges_over_time.png"
    plt.savefig(timeseries_path, dpi=300)
    plt.close()
    registry.add(
        figure_path=timeseries_path,
        figure_title=f"{ligand_code} water bridges over time",
        figure_category="water_bridges",
        ligand_code=ligand_code,
        protein_components=[plot_display_name_for_entry(entry) for entry in raw_entries],
        description_key=DESCRIPTION_KEYS["water_bridges"],
        notes_hint="Time-resolved ligand-water-protein bridge counts across protein components.",
        include_in_pdf=True,
        sort_order=72,
    )

    plot_summary = summary_df[summary_df["BridgeFraction"] > 0].copy()
    plot_summary = plot_summary.nlargest(30, "BridgeFraction")
    if not plot_summary.empty:
        log("   ⏳ Rendering water-bridge frequency plot")
        plt.figure(figsize=(12, max(5, len(plot_summary) * 0.25)))
        bar_colors = [colors.get(row["DisplayName"], "#444444") for _, row in plot_summary.iterrows()]
        plt.barh(plot_summary["ResidueLabel_Display"], plot_summary["BridgePercent"], color=bar_colors)
        plt.xlabel("Bridge percent")
        plt.ylabel("Residue")
        plt.title(f"{ligand_code} water-bridge frequency by residue")
        plt.gca().invert_yaxis()
        plt.tight_layout()
        frequency_path = water_dir / "protac_water_bridge_frequency_by_residue.png"
        plt.savefig(frequency_path, dpi=300)
        plt.close()
        registry.add(
            figure_path=frequency_path,
            figure_title=f"{ligand_code} water-bridge frequency by residue",
            figure_category="water_bridge_frequency",
            ligand_code=ligand_code,
            protein_components=[plot_display_name_for_entry(entry) for entry in raw_entries],
            description_key=DESCRIPTION_KEYS["water_bridge_frequency"],
            notes_hint="Water-bridge persistence by biological residue label.",
            include_in_pdf=True,
            sort_order=73,
        )

        heatmap_labels = plot_summary.sort_values(by=["ComponentOrder", "SourceResid", "CurrentResid"], kind="stable")["ResidueLabel_Display"].tolist()
        heatmap_df = (
            detail_df[detail_df["ResidueLabel_Display"].isin(heatmap_labels)]
            .assign(BridgeValue=1)
            .pivot_table(index="ResidueLabel_Display", columns="Frame", values="BridgeValue", aggfunc="max", fill_value=0)
            .reindex(heatmap_labels)
        )
        log("   ⏳ Rendering water-bridge heatmap")
        plt.figure(figsize=(12, max(5, len(heatmap_df) * 0.22)))
        sns.heatmap(heatmap_df, cmap="crest", cbar=True)
        plt.xlabel("Frame")
        plt.ylabel("Residue")
        plt.title(f"{ligand_code} water-bridge heatmap")
        plt.tight_layout()
        heatmap_path = water_dir / "protac_water_bridge_heatmap.png"
        plt.savefig(heatmap_path, dpi=300)
        plt.close()
        registry.add(
            figure_path=heatmap_path,
            figure_title=f"{ligand_code} water-bridge heatmap",
            figure_category="water_bridges",
            ligand_code=ligand_code,
            protein_components=[plot_display_name_for_entry(entry) for entry in raw_entries],
            description_key=DESCRIPTION_KEYS["water_bridges"],
            notes_hint="Residue-vs-frame map for ligand-water-protein bridge events.",
            include_in_pdf=True,
            sort_order=74,
        )

        if nx is not None:
            log("   ⏳ Rendering water-bridge network")
            net_subset = plot_summary[plot_summary["BridgeFraction"] >= 0.02].nlargest(30, "BridgeFraction")
            graph = nx.Graph()
            graph.add_node(ligand_code, role="Ligand")
            graph.add_node("Water", role="Water")
            graph.add_edge(ligand_code, "Water", weight=100.0)
            for _, row in net_subset.iterrows():
                node = row["ResidueLabel_Display"]
                graph.add_node(node, role=row["Role"], display_name=row["DisplayName"])
                graph.add_edge("Water", node, weight=float(row["BridgePercent"]))
            positions = nx.spring_layout(graph, seed=42, weight="weight")
            plt.figure(figsize=(12, 10))
            node_colors = []
            for node in graph.nodes():
                if node == ligand_code:
                    node_colors.append("#9467bd")
                elif node == "Water":
                    node_colors.append("#4fa3ff")
                else:
                    node_colors.append(colors.get(graph.nodes[node].get("display_name", ""), "#444444"))
            edge_widths = [max(1.0, graph[u][v]["weight"] / 15.0) for u, v in graph.edges()]
            nx.draw_networkx_edges(graph, positions, alpha=0.4, width=edge_widths, edge_color="#999999")
            nx.draw_networkx_nodes(graph, positions, node_size=800, node_color=node_colors, edgecolors="black")
            nx.draw_networkx_labels(graph, positions, font_size=8)
            plt.title(f"{ligand_code} water-bridge network")
            plt.axis("off")
            plt.tight_layout()
            network_path = water_dir / "protac_water_bridge_network.png"
            plt.savefig(network_path, dpi=300)
            plt.close()
            registry.add(
                figure_path=network_path,
                figure_title=f"{ligand_code} water-bridge network",
                figure_category="water_bridge_network",
                ligand_code=ligand_code,
                protein_components=[plot_display_name_for_entry(entry) for entry in raw_entries],
                description_key=DESCRIPTION_KEYS["water_bridge_network"],
                notes_hint="Network view of residues linked to the PROTAC through bridging waters.",
                include_in_pdf=True,
                sort_order=75,
            )
        else:
            log("   ⏭️ Skipping water-bridge network because NetworkX is unavailable")

    pocket_outputs = {}
    if write_water_pocket_trajectory:
        selected_residue_records = [
            record for record in raw_residue_records if record["residue_key"] in (pocket_residue_keys | bridge_residue_keys)
        ]
        selected_water_residue_ids = pocket_water_residue_ids | water_residue_ids
        if selected_residue_records:
            selected_groups = [record["residue"].atoms for record in selected_residue_records]
            if selected_water_residue_ids:
                water_atom_indices = []
                for residue in raw_universe.residues:
                    if int(residue.ix) in selected_water_residue_ids:
                        water_atom_indices.extend(residue.atoms.indices.tolist())
                if water_atom_indices:
                    selected_groups.append(raw_universe.atoms[np.array(sorted(set(water_atom_indices)), dtype=int)])
            selected_groups.append(raw_ligand)
            pocket_group = combine_atomgroups(selected_groups)
            pocket_pdb = structures_dir / "water_bridge_pocket.pdb"
            pocket_xtc = structures_dir / "water_bridge_pocket.xtc"
            if overwrite_structures or not file_exists(pocket_pdb):
                write_atomgroup_pdb(raw_universe, pocket_group, pocket_pdb, ligand_code, frame_indices[0])
            if overwrite_structures or not file_exists(pocket_xtc):
                write_atomgroup_xtc(raw_universe, pocket_group, pocket_xtc, frame_indices)
            pocket_outputs = {"pdb": pocket_pdb, "xtc": pocket_xtc}

    summary_txt = write_water_bridge_summary_txt(
        qc_dir / "protac_water_bridge_summary.txt",
        water_names,
        water_oxygens.n_atoms,
        total_sampled_frames,
        bridge_cutoff,
        int(detail_df["BridgeBinary"].sum()) if not detail_df.empty else 0,
        summary_df.nlargest(10, "BridgeFraction") if not summary_df.empty else summary_df,
    )

    return {
        "detail_df": detail_df,
        "summary_df": summary_df,
        "timeseries_df": timeseries_df,
        "component_summary_df": component_summary_df,
        "long_csv": long_csv,
        "frequency_csv": frequency_csv,
        "component_summary_csv": component_summary_csv,
        "summary_txt": summary_txt,
        "pocket_outputs": pocket_outputs,
    }


def run_pca(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    frame_indices: Sequence[int],
    qc_dir: Path,
    ligand_code: str,
    entries: Sequence[dict],
    registry: FigureRegistry,
) -> None:
    ensure_plotting_available()
    protein_ca = protein.select_atoms("name CA")
    if protein_ca.n_atoms == 0:
        protein_ca = protein

    log(f"   ⏳ Preparing PCA with {len(frame_indices)} frames and {protein_ca.n_atoms} atoms")
    aligned_positions = gather_aligned_positions(universe, protein_ca, frame_indices, "Aligning protein coordinates for PCA")
    if aligned_positions.size == 0:
        return

    matrix = aligned_positions.reshape(aligned_positions.shape[0], -1)
    matrix = matrix - matrix.mean(axis=0)
    _, _, vt = np.linalg.svd(matrix, full_matrices=False)
    projected = np.dot(matrix, vt[:2].T)
    pd.DataFrame({"Frame": frame_indices, "PC1": projected[:, 0], "PC2": projected[:, 1]}).to_csv(qc_dir / "protein_pca.csv", index=False)

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(projected[:, 0], projected[:, 1], c=frame_indices, cmap="viridis", s=30)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Protein PCA projection")
    plt.colorbar(scatter, label="Frame")
    plt.tight_layout()
    pca_path = qc_dir / "protein_pca.png"
    plt.savefig(pca_path, dpi=300)
    plt.close()
    register_standard_figure(
        registry,
        pca_path,
        "Protein PCA projection",
        "pca",
        ligand_code,
        entries,
        "Principal component projection of protein motions.",
        80,
    )


def write_empty_manifest(registry: FigureRegistry) -> Tuple[Path, Path]:
    return registry.write()


def main() -> int:
    args = parse_args()
    global GMX_BIN
    try:
        GMX_BIN = resolve_gmx_binary(args.gmx_bin)
        print_gmx_preflight()
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1
    args.pdb_reference = args.pdb_reference or default_pdb_reference()

    try:
        ensure_pandas_available()
        ensure_mdanalysis_available()
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1

    outdir = ensure_dir(Path(args.outdir).resolve())
    rms_dir = ensure_dir(outdir / "RMSD_RMSF")
    contacts_dir = ensure_dir(outdir / "Contacts")
    networks_dir = ensure_dir(outdir / "Networks")
    geometry_dir = ensure_dir(outdir / "Geometry")
    structures_dir = ensure_dir(outdir / "Structures")
    water_dir = ensure_dir(outdir / "Water_Bridges")
    qc_dir = ensure_dir(outdir / "QC")

    runtime_profiler = RuntimeProfiler(qc_dir / "protac_runtime_profile.txt")
    figure_registry = FigureRegistry(qc_dir)

    log("🧪 Starting generic PROTAC analysis")
    log(
        f"Command-line mode: dry_run={args.dry_run}, quick_test={args.quick_test}, contacts_only={args.contacts_only}, basic_only={args.basic_only}"
    )
    log(f"Frame controls (requested): start={args.start_frame}, end={args.end_frame}, step={args.frame_step}")

    ligand_setup = execute_stage(
        runtime_profiler,
        "Ligand preflight / auto-detection",
        detect_protac_ligand_setup,
        run_dir=".",
        topo_path=args.topo,
        traj_path=args.traj,
        cli_ligand=args.ligand,
        headless=args.headless,
    )
    execute_stage(runtime_profiler, "Write ligand preflight QC files", write_ligand_preflight_reports, ligand_setup, qc_dir)
    print_preflight_block(ligand_setup)

    if ligand_setup["errors"]:
        for error in ligand_setup["errors"]:
            print(f"ERROR: {error}", file=sys.stderr, flush=True)
        return 1
    if not ligand_setup["ligand_code"]:
        print("ERROR: Ligand code could not be determined. Pass --ligand <RESNAME>.", file=sys.stderr, flush=True)
        return 1

    raw_entries = execute_stage(runtime_profiler, "Parse atomIndex entries", parse_atomindex_entries, args.atomindex)
    entries = execute_stage(runtime_profiler, "Assign PROTAC roles from atomIndex order", assign_protac_roles, raw_entries)

    try:
        universe = execute_stage(runtime_profiler, "Load topology and trajectory", mda.Universe, args.topo, args.traj)
    except Exception as exc:
        print(f"ERROR: Could not load topology/trajectory: {exc}", file=sys.stderr, flush=True)
        return 1

    ligand = execute_stage(runtime_profiler, "Select ligand atoms", select_ligand, universe, ligand_setup["ligand_code"])
    protein = execute_stage(runtime_profiler, "Select protein atoms from atomIndex ranges", protein_atomgroup_from_entries, universe, entries)
    raw_universe = universe
    raw_entries = [dict(entry) for entry in entries]
    raw_ligand = ligand
    raw_protein = protein

    try:
        execute_stage(runtime_profiler, "Validate PROTAC inputs and role ranges", validate_protac_inputs, universe, ligand, entries)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1

    start, end, step, frame_indices = resolve_frame_window(
        universe.trajectory.n_frames,
        args.start_frame,
        args.end_frame,
        args.frame_step,
    )

    execute_stage(runtime_profiler, "Write role-assignment QC", write_system_roles, entries, qc_dir)
    execute_stage(runtime_profiler, "Write labeling audit QC", write_labeling_audit, entries, qc_dir)
    execute_stage(runtime_profiler, "Write ligand detection QC", write_ligand_detection_summary, ligand, ligand_setup["ligand_code"], qc_dir)

    execute_stage(runtime_profiler, "Initialize figure manifest", write_empty_manifest, figure_registry)

    if args.dry_run:
        print_startup_summary(ligand_setup["ligand_code"], entries, args.topo, args.traj, universe.trajectory.n_frames)
        log(describe_frame_selection(start, end, step, len(frame_indices)))
        log("Dry-run complete. QC outputs were written; no full analysis was executed.")
        return 0

    final_outputs = None
    analysis_topo = args.topo
    analysis_traj = args.traj

    if args.skip_preprocess:
        log("Skipping preprocessing by request.")
        if args.use_final_trajectory:
            log("Centered Final_Trajectory reuse is skipped because --skip-preprocess was requested.")
    elif args.write_final_trajectory:
        final_outputs = execute_stage(
            runtime_profiler,
            "Generate centered Final_Trajectory",
            generate_centered_final_trajectory,
            args.topo,
            args.traj,
            universe,
            ligand,
            entries,
            structures_dir,
            qc_dir,
            ligand_setup["ligand_code"],
            frame_indices,
            args.final_trajectory_name,
            True,
            args.overwrite_final_trajectory,
            True,
        )
        for warning in final_outputs.get("warnings", []):
            log(f"⚠️ {warning}")
        log(
            f"Final trajectory: {final_outputs['status']} via {final_outputs['method']}\n"
            f"Final topology: {final_outputs['cwd_final_pdb'].name}\n"
            f"Final trajectory file: {final_outputs['cwd_final_xtc'].name if file_exists(final_outputs['cwd_final_xtc']) else 'not written'}"
        )
        if args.use_final_trajectory and file_exists(final_outputs["cwd_final_pdb"]) and file_exists(final_outputs["cwd_final_xtc"]):
            reloaded_universe = execute_stage(
                runtime_profiler,
                "Reload centered Final_Trajectory universe",
                mda.Universe,
                str(final_outputs["cwd_final_pdb"]),
                str(final_outputs["cwd_final_xtc"]),
            )
            universe = reloaded_universe
            entries = final_outputs["remapped_entries"]
            ligand = execute_stage(runtime_profiler, "Reselect ligand atoms in Final_Trajectory", select_ligand, universe, ligand_setup["ligand_code"])
            protein = execute_stage(runtime_profiler, "Remap protein atoms in Final_Trajectory", protein_atomgroup_from_entries, universe, entries)
            execute_stage(runtime_profiler, "Validate remapped Final_Trajectory inputs", validate_protac_inputs, universe, ligand, entries)
            analysis_topo = final_outputs["cwd_final_pdb"].name
            analysis_traj = final_outputs["cwd_final_xtc"].name if file_exists(final_outputs["cwd_final_xtc"]) else analysis_traj
        else:
            if args.use_final_trajectory and not file_exists(final_outputs["cwd_final_xtc"]):
                log("⚠️ Final_Trajectory.xtc was not available, so the original trajectory will remain the downstream analysis input.")
            else:
                log("Final_Trajectory was generated but the original trajectory will remain the downstream analysis input.")
    else:
        log("Final trajectory generation disabled by request.")

    print_startup_summary(ligand_setup["ligand_code"], entries, analysis_topo, analysis_traj, universe.trajectory.n_frames)
    log(describe_frame_selection(start, end, step, len(frame_indices)))
    log(
        "Analysis input: "
        + (
            f"{analysis_topo} + {analysis_traj}"
            if analysis_traj
            else analysis_topo
        )
    )

    residue_records = execute_stage(
        runtime_profiler,
        "Build residue-to-role mapping",
        build_protein_residue_records,
        protein,
        entries,
        args.residue_numbering,
    )
    if not residue_records:
        print(
            "ERROR: No protein residues could be mapped onto atomIndex roles. Check that atomIndex.txt matches the analyzed topology.",
            file=sys.stderr,
            flush=True,
        )
        return 1

    preprocessing_outputs = None
    if not args.skip_preprocess:
        preprocessing_outputs = execute_stage(
            runtime_profiler,
            "Preprocessing / subset structure generation",
            run_preprocessing,
            universe,
            protein,
            ligand,
            entries,
            residue_records,
            frame_indices,
            structures_dir,
            qc_dir,
            ligand_setup["ligand_code"],
            args.write_subset_trajectories,
            args.pocket_cutoff,
            args.pocket_mode,
            args.pocket_min_fraction,
            args.overwrite_structures or args.overwrite_final_trajectory,
            args.final_trajectory_name,
            final_outputs.get("structure_rows", []) if final_outputs else [],
        )

    if args.preprocess_only:
        execute_stage(runtime_profiler, "Write figure manifest", figure_registry.write)
        execute_stage(runtime_profiler, "Write runtime profile", runtime_profiler.write)
        log("Preprocess-only run complete.")
        return 0

    contact_outputs = execute_stage(
        runtime_profiler,
        "Ligand contacts" + ("" if args.skip_contact_map else " / contact map") + ("" if args.skip_network else " / network"),
        run_contact_analysis,
        universe,
        protein,
        ligand,
        ligand_setup["ligand_code"],
        residue_records,
        entries,
        frame_indices,
        args.distance_cutoff,
        contacts_dir,
        networks_dir,
        geometry_dir,
        qc_dir,
        figure_registry,
        args.residue_numbering,
        not args.skip_contact_map,
        not args.skip_network,
        args.interaction_mode,
        args.plot_allprotein,
        args.network_min_contact_fraction,
        args.network_max_residues,
        args.contact_map_min_contact_fraction,
        args.contact_map_max_residues,
    )

    if args.analyze_water_bridges:
        raw_residue_records = execute_stage(
            runtime_profiler,
            "Build raw residue-to-role mapping for water bridges",
            build_protein_residue_records,
            raw_protein,
            raw_entries,
            args.residue_numbering,
        )
        execute_stage(
            runtime_profiler,
            "Water-bridge analysis",
            run_water_bridge_analysis,
            raw_universe,
            raw_protein,
            raw_ligand,
            raw_entries,
            raw_residue_records,
            frame_indices,
            ligand_setup["ligand_code"],
            water_dir,
            structures_dir,
            qc_dir,
            figure_registry,
            parse_name_list(args.water_names),
            args.water_bridge_cutoff,
            args.water_pocket_cutoff,
            args.write_water_pocket_trajectory,
            args.overwrite_structures or args.overwrite_final_trajectory,
        )
    else:
        log("Skipping water-bridge analysis (use --analyze-water-bridges to enable).")

    if args.contacts_only:
        log("Contacts-only mode enabled: skipping RMSD, RMSF, radius of gyration, and PCA.")
    else:
        execute_stage(
            runtime_profiler,
            "RMSD" + ("" if args.skip_rmsf else " / RMSF"),
            run_rmsd_and_rmsf,
            universe,
            protein,
            ligand,
            entries,
            residue_records,
            frame_indices,
            rms_dir,
            ligand_setup["ligand_code"],
            figure_registry,
            not args.skip_rmsf,
            contact_outputs.get("residue_frequency_df"),
            args.rmsf_contact_threshold,
            args.rmsf_label_top_contacts,
            args.rmsf_contact_marker_mode,
        )
        execute_stage(
            runtime_profiler,
            "Radius of gyration",
            run_radius_of_gyration,
            universe,
            protein,
            entries,
            frame_indices,
            rms_dir,
            geometry_dir,
            ligand_setup["ligand_code"],
            figure_registry,
        )

    if not args.contacts_only and not args.skip_pca:
        execute_stage(
            runtime_profiler,
            "PCA",
            run_pca,
            universe,
            protein,
            frame_indices,
            qc_dir,
            ligand_setup["ligand_code"],
            entries,
            figure_registry,
        )
    elif args.skip_pca:
        log("Skipping PCA by request.")

    execute_stage(runtime_profiler, "Write figure manifest", figure_registry.write)
    execute_stage(runtime_profiler, "Write runtime profile", runtime_profiler.write)
    log(f"Analysis complete. Outputs written under {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
