#!/usr/bin/env python3
"""
Generic PROTAC / ternary-complex MD analysis.

This script is designed for PROTAC systems prepared by the current
automation workflow where `atomIndex.txt` is the source of truth for
protein atom ranges:

  - first protein entry  -> Ligase
  - later protein entries -> Target_A, Target_B, Target_C, ...

The ligand / PROTAC identity is detected dynamically when possible and
can always be overridden with `--ligand`.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import MDAnalysis as mda

from MDAnalysis.analysis import align
from MDAnalysis.lib.distances import distance_array


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze a PROTAC ternary-complex trajectory using atomIndex.txt-defined ligase/target roles."
    )
    parser.add_argument(
        "-l",
        "--ligand",
        default=None,
        help="Ligand / PROTAC residue name. Optional when auto-detection is confident.",
    )
    parser.add_argument(
        "--topo",
        default="npt.gro",
        help="Topology/coordinate file for MDAnalysis (default: npt.gro).",
    )
    parser.add_argument(
        "--traj",
        default="md_0_1.xtc",
        help="Trajectory file for MDAnalysis (default: md_0_1.xtc).",
    )
    parser.add_argument(
        "--pdb-reference",
        default=None,
        help="Optional reference PDB. Defaults to Complex.pdb when present.",
    )
    parser.add_argument(
        "--atomindex",
        default="atomIndex.txt",
        help="atomIndex file with ordered protein atom ranges (default: atomIndex.txt).",
    )
    parser.add_argument(
        "--outdir",
        default="Analysis_Results/PROTAC",
        help="Output directory root (default: Analysis_Results/PROTAC).",
    )
    parser.add_argument(
        "--distance-cutoff",
        type=float,
        default=4.0,
        help="Protein-ligand contact cutoff in angstroms (default: 4.0).",
    )
    parser.add_argument(
        "--headless",
        action="store_true",
        help="Disable interactive prompts and fail clearly when detection is ambiguous.",
    )
    parser.add_argument(
        "--start-frame",
        type=int,
        default=None,
        help="First frame index to analyze for trajectory-scan analyses.",
    )
    parser.add_argument(
        "--end-frame",
        type=int,
        default=None,
        help="Last frame index (exclusive) to analyze for trajectory-scan analyses.",
    )
    parser.add_argument(
        "--frame-step",
        type=int,
        default=1,
        help="Frame stride for expensive analyses (default: 1).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs, run ligand preflight, write QC, and exit without full analysis.",
    )
    return parser.parse_args()


def default_pdb_reference() -> Optional[str]:
    candidate = Path("Complex.pdb")
    return str(candidate) if candidate.exists() else None


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def normalize_candidate_code(value: str) -> Optional[str]:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", value.strip())
    if 2 <= len(cleaned) <= 5:
        return cleaned.upper()
    return None


def role_name_for_target(index: int) -> str:
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if index < len(letters):
        return f"Target_{letters[index]}"
    return f"Target_{index + 1}"


def parse_atomindex_entries(atomindex_path: str) -> List[dict]:
    path = Path(atomindex_path)
    if not path.exists():
        raise FileNotFoundError(
            f"atomIndex file not found: {path}. "
            "This script requires atomIndex.txt from the current preparation workflow."
        )

    entries: List[dict] = []
    invalid_lines: List[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue

            content = stripped.split("#", 1)[0].strip()
            if not content:
                continue

            match = re.match(r"^(.*?)\s+(\d+)\s*-\s*(\d+)$", content)
            if not match:
                invalid_lines.append(f"line {line_number}: {raw.rstrip()}")
                continue

            name = match.group(1).strip()
            start_1based = int(match.group(2))
            end_1based = int(match.group(3))

            if not name or start_1based <= 0 or end_1based < start_1based:
                invalid_lines.append(f"line {line_number}: {raw.rstrip()}")
                continue

            entries.append(
                {
                    "name": name,
                    "start_1based": start_1based,
                    "end_1based": end_1based,
                    "start0": start_1based - 1,
                    "end0": end_1based - 1,
                }
            )

    if invalid_lines:
        joined = "\n".join(invalid_lines[:10])
        raise ValueError(
            "Failed to parse atomIndex.txt. Expected lines like 'Name 1-1975'.\n"
            f"Unparsable entries:\n{joined}"
        )

    if not entries:
        raise ValueError(
            f"No atom ranges were parsed from {path}. "
            "Please regenerate atomIndex.txt with the current workflow."
        )

    return entries


def assign_protac_roles(entries: Sequence[dict]) -> List[dict]:
    assigned: List[dict] = []
    for index, entry in enumerate(entries):
        role = "Ligase" if index == 0 else role_name_for_target(index - 1)
        new_entry = dict(entry)
        new_entry["role"] = role
        assigned.append(new_entry)
    return assigned


def atomgroup_for_entry(universe: mda.Universe, entry: dict) -> mda.AtomGroup:
    start0 = int(entry["start0"])
    end0 = int(entry["end0"])
    if start0 < 0 or end0 >= universe.atoms.n_atoms:
        raise ValueError(
            f"atomIndex range {entry['name']} {entry['start_1based']}-{entry['end_1based']} "
            f"falls outside topology atom count ({universe.atoms.n_atoms})."
        )
    return universe.atoms[start0 : end0 + 1]


def overlap_count(indices: np.ndarray, start0: int, end0: int) -> int:
    return int(np.count_nonzero((indices >= start0) & (indices <= end0)))


def residue_role_for_residue(residue: mda.core.groups.Residue, entries: Sequence[dict]) -> Tuple[Optional[str], Optional[str]]:
    indices = residue.atoms.indices
    best_role = None
    best_name = None
    best_overlap = 0

    for entry in entries:
        count = overlap_count(indices, int(entry["start0"]), int(entry["end0"]))
        if count > best_overlap:
            best_overlap = count
            best_role = entry["role"]
            best_name = entry["name"]

    return best_role, best_name


def residue_label(residue: mda.core.groups.Residue, role: str) -> str:
    return f"{role}:{residue.resname}{residue.resid}"


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
    result = {
        "include_candidates": set(),
        "molecule_candidates": set(),
        "include_lines": [],
    }
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
                    include_name = os.path.basename(match.group(1))
                    stem = Path(include_name).stem
                    normalized = normalize_candidate_code(stem)
                    lower_stem = stem.lower()
                    result["include_lines"].append(include_name)

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
                    molname = tokens[0]
                    normalized = normalize_candidate_code(molname)
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
        ranked = sorted(
            residue_info.items(),
            key=lambda item: (-item[1]["residue_atom_count"], item[0]),
        )
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
        merged[code]["residue_atom_count"] = max(
            merged[code]["residue_atom_count"], data.get("residue_atom_count", 0)
        )
        merged[code]["residue_count"] = max(
            merged[code]["residue_count"], data.get("residue_count", 0)
        )
        for key, value in data.get("files", {}).items():
            if value and not merged[code]["files"].get(key):
                merged[code]["files"][key] = value
        merged[code]["warnings"].extend(data.get("warnings", []))
        merged[code]["errors"].extend(data.get("errors", []))
    return merged


def source_string(sources: Sequence[str]) -> str:
    ordered = []
    for source in [
        "cli",
        "cgenff_file_set",
        "topol_top",
        "gro_residue",
        "trajectory_residue",
        "index_ndx",
        "manual",
    ]:
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
    print("\nLigand auto-detection is ambiguous. Candidate summary:")
    for index, (code, data) in enumerate(candidates, start=1):
        print(
            f"  {index}. {code} | score={data['score']} | "
            f"residues={data.get('residue_count', 0)} | atoms={data.get('residue_atom_count', 0)} | "
            f"sources={source_string(data.get('sources', [])) or 'manual'}"
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
        print("Please enter a valid choice.")


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
    for code in topol_info["include_candidates"]:
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

    for code in topol_info["molecule_candidates"]:
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
        key=lambda item: (
            -item[1]["score"],
            -item[1].get("residue_atom_count", 0),
            -item[1].get("residue_count", 0),
            item[0],
        ),
    )

    selected_code: Optional[str] = cli_code
    errors: List[str] = []

    if not selected_code:
        if not ranked:
            if not headless:
                manual = input("Ligand could not be auto-detected. Enter ligand residue name: ").strip()
                selected_code = normalize_candidate_code(manual)
                if selected_code:
                    candidates[selected_code] = {
                        "score": 0,
                        "sources": ["manual"],
                        "files": {key: None for key in FILE_KEYS},
                        "warnings": [],
                        "errors": [],
                        "residue_atom_count": 0,
                        "residue_count": 0,
                    }
                else:
                    errors.append("No valid ligand residue name was provided.")
            else:
                errors.append(
                    "Could not auto-detect a PROTAC ligand from the current directory. "
                    "Pass --ligand <RESNAME> in headless mode."
                )
        else:
            best_code, best_data = ranked[0]
            lead = best_data["score"] - ranked[1][1]["score"] if len(ranked) > 1 else best_data["score"]
            if len(ranked) > 1 and best_data["score"] == ranked[1][1]["score"]:
                if headless:
                    summary = "; ".join(
                        f"{code}: score={data['score']} sources={source_string(data['sources']) or 'manual'}"
                        for code, data in ranked[:5]
                    )
                    errors.append(
                        "Ligand detection is ambiguous in headless mode. "
                        f"Top candidates: {summary}. Pass --ligand explicitly."
                    )
                else:
                    selected_code = prompt_user_for_ligand(ranked[:5])
                    candidates[selected_code]["sources"].append("manual")
            elif best_data["score"] < 3:
                if headless:
                    errors.append(
                        "Ligand auto-detection confidence is too low for headless mode. "
                        "Pass --ligand explicitly."
                    )
                else:
                    selected_code = prompt_user_for_ligand(ranked[:5])
                    candidates[selected_code]["sources"].append("manual")
            else:
                selected_code = best_code

    selected = candidates.get(selected_code, None) if selected_code else None

    if selected_code and selected is None:
        selected = {
            "score": 0,
            "sources": ["manual"],
            "files": {key: None for key in FILE_KEYS},
            "warnings": [],
            "errors": [],
            "residue_atom_count": 0,
            "residue_count": 0,
        }

    selected_files = {key: None for key in FILE_KEYS}
    confidence = "low"
    source = "manual" if selected_code else ""

    if selected:
        selected_files.update(selected["files"])
        best_score = selected.get("score", 0)
        second_score = ranked[1][1]["score"] if len(ranked) > 1 and ranked[0][0] == selected_code else ranked[0][1]["score"] if ranked and ranked[0][0] != selected_code else 0
        lead = best_score - second_score if ranked else best_score
        confidence = score_confidence(best_score, lead, bool(cli_code))
        source = source_string(selected.get("sources", []))

    selected_warnings = list(warnings)
    if selected:
        selected_warnings.extend(selected.get("warnings", []))

    residue_candidates = {
        code: data for code, data in structure_evidence.items() if data.get("residue_atom_count", 0) > 0
    }

    if selected_code and residue_candidates and selected_code not in residue_candidates:
        residue_summary = ", ".join(
            f"{code}({data['residue_atom_count']} atoms)"
            for code, data in sorted(residue_candidates.items())
        )
        selected_warnings.append(
            f"Selected ligand code {selected_code} was not observed as a residue in the loaded structure. "
            f"Detected non-protein residues: {residue_summary}"
        )

    if selected_code and (run_path / "topol.top").exists():
        include_match = selected_files.get("itp") or selected_files.get("prm")
        if not include_match and selected_code not in topol_info["include_candidates"]:
            selected_warnings.append(
                f"topol.top did not show a clear include for ligand code {selected_code}."
            )

    if cli_code and structure_source and cli_code not in residue_candidates:
        errors.append(
            f"--ligand {cli_code} was provided, but that residue name was not found in the loaded topology/trajectory."
        )

    if selected_code and len(residue_candidates) > 1:
        other_codes = [code for code in residue_candidates if code != selected_code]
        if other_codes:
            selected_warnings.append(
                f"Additional ligand-like residues were detected: {', '.join(sorted(other_codes))}"
            )

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

    lines.append(
        f"Topology residue atoms: {ligand_setup.get('topology_residue_atoms', 0)}"
    )
    lines.append(
        f"Topology residue count: {ligand_setup.get('topology_residue_count', 0)}"
    )

    warnings = ligand_setup.get("warnings", [])
    errors = ligand_setup.get("errors", [])
    lines.append("Warnings:")
    if warnings:
        lines.extend([f"  - {warning}" for warning in warnings])
    else:
        lines.append("  none")
    lines.append("Errors:")
    if errors:
        lines.extend([f"  - {error}" for error in errors])
    else:
        lines.append("  none")

    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, txt_path


def print_preflight_block(ligand_setup: dict) -> None:
    source_display = (ligand_setup.get("source") or "manual").replace("|", " + ")
    print("\n" + "=" * 72)
    print("🧬 PROTAC ligand preflight")
    print(f"Ligand selected: {ligand_setup.get('ligand_code') or 'None'}")
    print(f"Confidence: {ligand_setup.get('confidence', 'unknown')}")
    print(f"Source: {source_display}")
    print("Found files:")
    for key in FILE_KEYS:
        value = ligand_setup.get("files", {}).get(key)
        if value:
            print(f"  ✓ {Path(value).name}")
    if not any(ligand_setup.get("files", {}).values()):
        print("  none")
    print(f"Topology residue atoms: {ligand_setup.get('topology_residue_atoms', 0)}")
    warnings = ligand_setup.get("warnings", [])
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("Warnings: none")
    print("=" * 72)


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
        raise ValueError(
            f"Invalid frame window: start={start}, end={end}, total_frames={n_frames}"
        )
    frame_indices = list(range(start, end, step))
    return start, end, step, frame_indices


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
) -> np.ndarray:
    if atomgroup.n_atoms == 0:
        return np.empty((0, 0, 3))

    universe.trajectory[frame_indices[0]]
    reference = atomgroup.positions.copy()
    aligned_frames = []

    for frame in frame_indices:
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


def atom_rmsf_to_residue_rmsf(atomgroup: mda.AtomGroup, atom_rmsf: np.ndarray) -> pd.DataFrame:
    rows = []
    cursor = 0
    for residue in atomgroup.residues:
        count = residue.atoms.n_atoms
        residue_values = atom_rmsf[cursor : cursor + count]
        cursor += count
        rows.append(
            {
                "ResidueName": residue.resname,
                "ResidueID": residue.resid,
                "ResidueRMSF": float(np.mean(residue_values)) if len(residue_values) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def classify_interaction_type(residue: mda.core.groups.Residue, min_distance: float) -> str:
    resname = residue.resname.strip().upper()
    labels: List[str] = []

    if resname in HYDROPHOBIC_RESIDUES and min_distance <= 4.5:
        labels.append("Hydrophobic")
    if resname in AROMATIC_RESIDUES and min_distance <= 4.8:
        labels.append("Aromatic")
    if resname in POSITIVE_RESIDUES | NEGATIVE_RESIDUES and min_distance <= 4.0:
        labels.append("Ionic")
    if resname in POLAR_RESIDUES | POSITIVE_RESIDUES | NEGATIVE_RESIDUES and min_distance <= 3.5:
        labels.append("Hydrophilic")

    if not labels:
        return "Contact"
    return ";".join(labels)


def build_protein_residue_records(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    entries: Sequence[dict],
) -> List[dict]:
    local_lookup = {global_index: local_index for local_index, global_index in enumerate(protein.indices)}
    records: List[dict] = []

    for residue in protein.residues:
        role, protein_name = residue_role_for_residue(residue, entries)
        if not role:
            continue
        local_indices = np.array(
            [local_lookup[index] for index in residue.atoms.indices if index in local_lookup],
            dtype=int,
        )
        if local_indices.size == 0:
            continue
        records.append(
            {
                "residue": residue,
                "role": role,
                "protein_name": protein_name,
                "residue_label": residue_label(residue, role),
                "protein_local_indices": local_indices,
            }
        )
    return records


def write_system_roles(entries: Sequence[dict], qc_dir: Path) -> Path:
    rows = []
    for entry in entries:
        rows.append(
            {
                "Name": entry["name"],
                "Role": entry["role"],
                "Start_1based": entry["start_1based"],
                "End_1based": entry["end_1based"],
                "Start0": entry["start0"],
                "End0": entry["end0"],
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
    if not target_display:
        target_display = "none"

    print("\nStartup summary")
    print(f"Ligand: {ligand_code}")
    print(
        f"Ligase: {ligase['name']} ({ligase['start_1based']}-{ligase['end_1based']})"
    )
    print(f"Targets: {target_display}")
    print(f"Topology: {topo_path}")
    print(f"Trajectory: {traj_path}")
    print(f"Frames: {frame_count}")


def role_color(role: str) -> str:
    if role == "Ligase":
        return "#1f77b4"
    return "#d62728" if role == "Target_A" else "#2ca02c"


def plot_series(
    x: Sequence[float],
    y_series: Dict[str, Sequence[float]],
    xlabel: str,
    ylabel: str,
    title: str,
    output_path: Path,
) -> None:
    plt.figure(figsize=(9, 5))
    for label, values in y_series.items():
        plt.plot(x, values, label=label, linewidth=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def run_rmsd_and_rmsf(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    ligand: mda.AtomGroup,
    entries: Sequence[dict],
    frame_indices: Sequence[int],
    rms_dir: Path,
) -> None:
    ensure_dir(rms_dir)
    rows_rmsd = []
    rows_rmsf = []
    times = []
    for frame in frame_indices:
        universe.trajectory[frame]
        times.append(float(universe.trajectory.time))

    protein_ca = protein.select_atoms("name CA")
    if protein_ca.n_atoms == 0:
        protein_ca = protein

    protein_positions = gather_aligned_positions(universe, protein_ca, frame_indices)
    protein_rmsd = rmsd_from_aligned_positions(protein_positions)
    protein_atom_rmsf = rmsf_from_aligned_positions(protein_positions)
    protein_residue_rmsf = atom_rmsf_to_residue_rmsf(protein_ca, protein_atom_rmsf)
    protein_residue_rmsf["Role"] = "All_Protein"
    protein_residue_rmsf["ProteinName"] = "All_Protein"
    rows_rmsf.append(protein_residue_rmsf)

    for frame, time_ps, value in zip(frame_indices, times, protein_rmsd):
        rows_rmsd.append(
            {
                "Frame": frame,
                "Time_ps": time_ps,
                "Role": "All_Protein",
                "ProteinName": "All_Protein",
                "RMSD": float(value),
            }
        )

    entry_rmsd_series = {"All_Protein": protein_rmsd}

    for entry in entries:
        entry_group = atomgroup_for_entry(universe, entry).select_atoms("protein and name CA")
        if entry_group.n_atoms == 0:
            entry_group = atomgroup_for_entry(universe, entry).select_atoms("protein")
        if entry_group.n_atoms == 0:
            continue

        entry_positions = gather_aligned_positions(universe, entry_group, frame_indices)
        entry_rmsd = rmsd_from_aligned_positions(entry_positions)
        entry_rmsf = rmsf_from_aligned_positions(entry_positions)
        entry_residue_rmsf = atom_rmsf_to_residue_rmsf(entry_group, entry_rmsf)
        entry_residue_rmsf["Role"] = entry["role"]
        entry_residue_rmsf["ProteinName"] = entry["name"]
        rows_rmsf.append(entry_residue_rmsf)
        entry_rmsd_series[entry["role"]] = entry_rmsd

        for frame, time_ps, value in zip(frame_indices, times, entry_rmsd):
            rows_rmsd.append(
                {
                    "Frame": frame,
                    "Time_ps": time_ps,
                    "Role": entry["role"],
                    "ProteinName": entry["name"],
                    "RMSD": float(value),
                }
            )

    ligand_positions = gather_aligned_positions(universe, ligand, frame_indices)
    ligand_rmsd = rmsd_from_aligned_positions(ligand_positions)
    ligand_rmsf = rmsf_from_aligned_positions(ligand_positions)
    ligand_df = pd.DataFrame(
        {
            "AtomIndex": np.arange(ligand.n_atoms),
            "AtomName": ligand.names,
            "RMSF": ligand_rmsf,
        }
    )
    ligand_df.to_csv(rms_dir / "ligand_rmsf.csv", index=False)

    ligand_rmsd_df = pd.DataFrame(
        {
            "Frame": frame_indices,
            "Time_ps": times,
            "Role": "Ligand",
            "ProteinName": "Ligand",
            "RMSD": ligand_rmsd,
        }
    )
    ligand_rmsd_df.to_csv(rms_dir / "ligand_rmsd.csv", index=False)

    rmsd_df = pd.DataFrame(rows_rmsd)
    rmsf_df = pd.concat(rows_rmsf, ignore_index=True)
    rmsd_df.to_csv(rms_dir / "protein_rmsd_by_role.csv", index=False)
    rmsf_df.to_csv(rms_dir / "protein_rmsf_by_role.csv", index=False)

    plot_series(
        frame_indices,
        entry_rmsd_series,
        xlabel="Frame",
        ylabel="RMSD (Å)",
        title="Protein RMSD by PROTAC Role",
        output_path=rms_dir / "protein_rmsd_by_role.png",
    )

    plt.figure(figsize=(9, 5))
    for role, subset in rmsf_df.groupby("Role"):
        plt.plot(subset["ResidueID"], subset["ResidueRMSF"], label=role, linewidth=1.6)
    plt.xlabel("Residue ID")
    plt.ylabel("RMSF (Å)")
    plt.title("Protein RMSF by PROTAC Role")
    plt.legend()
    plt.tight_layout()
    plt.savefig(rms_dir / "protein_rmsf_by_role.png", dpi=300)
    plt.close()

    plt.figure(figsize=(9, 5))
    plt.plot(frame_indices, ligand_rmsd, color="#9467bd", linewidth=2)
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title("Ligand RMSD")
    plt.tight_layout()
    plt.savefig(rms_dir / "ligand_rmsd.png", dpi=300)
    plt.close()

    plt.figure(figsize=(9, 5))
    plt.plot(np.arange(ligand.n_atoms), ligand_rmsf, color="#9467bd", linewidth=1.8)
    plt.xlabel("Ligand atom index")
    plt.ylabel("RMSF (Å)")
    plt.title("Ligand RMSF")
    plt.tight_layout()
    plt.savefig(rms_dir / "ligand_rmsf.png", dpi=300)
    plt.close()


def run_radius_of_gyration(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    entries: Sequence[dict],
    frame_indices: Sequence[int],
    rms_dir: Path,
) -> None:
    rows = []
    for frame in frame_indices:
        universe.trajectory[frame]
        rows.append(
            {
                "Frame": frame,
                "Time_ps": float(universe.trajectory.time),
                "Role": "All_Protein",
                "ProteinName": "All_Protein",
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
                    "RadiusOfGyration": float(entry_group.radius_of_gyration()),
                }
            )

    df = pd.DataFrame(rows)
    df.to_csv(rms_dir / "radius_of_gyration_by_role.csv", index=False)

    plt.figure(figsize=(9, 5))
    for role, subset in df.groupby("Role"):
        plt.plot(subset["Frame"], subset["RadiusOfGyration"], label=role, linewidth=1.8)
    plt.xlabel("Frame")
    plt.ylabel("Radius of gyration (Å)")
    plt.title("Protein compactness by PROTAC role")
    plt.legend()
    plt.tight_layout()
    plt.savefig(rms_dir / "radius_of_gyration_by_role.png", dpi=300)
    plt.close()


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
) -> None:
    ensure_dir(contacts_dir)
    ensure_dir(networks_dir)
    ensure_dir(geometry_dir)

    role_groups = {}
    target_indices = []
    protein_local_lookup = {atom_index: local_index for local_index, atom_index in enumerate(protein.indices)}
    for entry in entries:
        group = atomgroup_for_entry(universe, entry).select_atoms("protein")
        if group.n_atoms == 0:
            continue
        role_groups[entry["role"]] = {
            "protein_name": entry["name"],
            "local_indices": np.array(
                [protein_local_lookup[atom_index] for atom_index in group.indices if atom_index in protein_local_lookup],
                dtype=int,
            ),
        }
        if entry["role"].startswith("Target_"):
            target_indices.extend(group.indices.tolist())

    if target_indices:
        target_local = np.array(
            [protein_local_lookup[atom_index] for atom_index in sorted(set(target_indices)) if atom_index in protein_local_lookup],
            dtype=int,
        )
        role_groups["Target_Combined"] = {
            "protein_name": "Target_Combined",
            "local_indices": target_local,
        }

    role_groups["All_Protein"] = {
        "protein_name": "All_Protein",
        "local_indices": np.arange(protein.n_atoms, dtype=int),
    }

    detail_rows = []
    role_timeseries_rows = []
    per_role_frame_contact_counts = defaultdict(list)

    for frame in frame_indices:
        universe.trajectory[frame]
        dist_matrix = distance_array(ligand.positions, protein.positions)

        for role, group_info in role_groups.items():
            local_indices = group_info["local_indices"]
            if local_indices.size == 0:
                continue
            submatrix = dist_matrix[:, local_indices]
            min_distance = float(np.min(submatrix))
            contact_atom_pairs = int(np.count_nonzero(submatrix <= distance_cutoff))

            if role in {"All_Protein", "Target_Combined"}:
                residue_count = 0
                for record in residue_records:
                    if role == "Target_Combined" and not record["role"].startswith("Target_"):
                        continue
                    res_min = float(np.min(dist_matrix[:, record["protein_local_indices"]]))
                    if res_min <= distance_cutoff:
                        residue_count += 1
            else:
                residue_count = 0
                for record in residue_records:
                    if record["role"] != role:
                        continue
                    res_min = float(np.min(dist_matrix[:, record["protein_local_indices"]]))
                    if res_min <= distance_cutoff:
                        residue_count += 1

            role_timeseries_rows.append(
                {
                    "Frame": frame,
                    "Time_ps": float(universe.trajectory.time),
                    "Role": role,
                    "ProteinName": group_info["protein_name"],
                    "ResidueName": "",
                    "ResidueID": "",
                    "ResidueLabel": "",
                    "InteractionType": "Contact",
                    "MinDistance": min_distance,
                    "ContactResidues": residue_count,
                    "ContactAtomPairs": contact_atom_pairs,
                }
            )
            per_role_frame_contact_counts[role].append(int(contact_atom_pairs > 0))

        for record in residue_records:
            residue = record["residue"]
            min_distance = float(np.min(dist_matrix[:, record["protein_local_indices"]]))
            if min_distance > distance_cutoff:
                continue
            interaction_type = classify_interaction_type(residue, min_distance)
            detail_rows.append(
                {
                    "Frame": frame,
                    "Time_ps": float(universe.trajectory.time),
                    "Role": record["role"],
                    "ProteinName": record["protein_name"],
                    "ResidueName": residue.resname,
                    "ResidueID": residue.resid,
                    "ResidueLabel": record["residue_label"],
                    "InteractionType": interaction_type,
                    "MinDistance": min_distance,
                    "ContactFraction": np.nan,
                    "ContactPercent": np.nan,
                }
            )

    detail_df = pd.DataFrame(detail_rows)
    timeseries_df = pd.DataFrame(role_timeseries_rows)
    detail_path = contacts_dir / "ligand_contacts_detailed.csv"
    timeseries_path = contacts_dir / "ligand_contacts_by_role_timeseries.csv"
    detail_df.to_csv(detail_path, index=False)
    timeseries_df.to_csv(timeseries_path, index=False)

    if detail_df.empty:
        raise ValueError(
            "No ligand-protein contacts were detected. "
            "Check the ligand residue name or distance cutoff."
        )

    residue_summary = (
        detail_df.groupby(
            ["Role", "ProteinName", "ResidueName", "ResidueID", "ResidueLabel", "InteractionType"],
            dropna=False,
        )
        .agg(ContactFrames=("Frame", "nunique"), MinDistance=("MinDistance", "min"))
        .reset_index()
    )
    total_frames = len(frame_indices)
    residue_summary["ContactFraction"] = residue_summary["ContactFrames"] / total_frames
    residue_summary["ContactPercent"] = residue_summary["ContactFraction"] * 100.0
    residue_summary.to_csv(contacts_dir / "ligand_contact_residue_summary.csv", index=False)

    role_summary_rows = []
    for role, values in per_role_frame_contact_counts.items():
        frames_with_contact = int(np.sum(values))
        fraction = frames_with_contact / total_frames if total_frames else 0.0
        role_subset = timeseries_df[timeseries_df["Role"] == role]
        role_summary_rows.append(
            {
                "Role": role,
                "FramesWithContact": frames_with_contact,
                "ContactFraction": fraction,
                "ContactPercent": fraction * 100.0,
                "MeanContactResidues": float(role_subset["ContactResidues"].mean()) if not role_subset.empty else 0.0,
                "MeanMinDistance": float(role_subset["MinDistance"].mean()) if not role_subset.empty else math.nan,
            }
        )

    role_summary_df = pd.DataFrame(role_summary_rows)
    role_summary_df.to_csv(qc_dir / "protac_contact_summary.csv", index=False)

    plt.figure(figsize=(10, 5))
    for role, subset in timeseries_df.groupby("Role"):
        plt.plot(subset["Frame"], subset["ContactResidues"], label=role, linewidth=1.8)
    plt.xlabel("Frame")
    plt.ylabel("Residues contacting ligand")
    plt.title(f"{ligand_code} contacts over time by PROTAC role")
    plt.legend()
    plt.tight_layout()
    plt.savefig(contacts_dir / "ligand_contacts_over_time_by_role.png", dpi=300)
    plt.close()

    heatmap_df = (
        detail_df.assign(ContactValue=1)
        .pivot_table(
            index="ResidueLabel",
            columns="Frame",
            values="ContactValue",
            aggfunc="max",
            fill_value=0,
        )
        .sort_index()
    )
    plt.figure(figsize=(12, max(5, len(heatmap_df) * 0.2)))
    sns.heatmap(heatmap_df, cmap="mako", cbar=True)
    plt.xlabel("Frame")
    plt.ylabel("Residue")
    plt.title(f"{ligand_code} contact map by residue and frame")
    plt.tight_layout()
    plt.savefig(contacts_dir / "ligand_contact_map.png", dpi=300)
    plt.close()

    interaction_summary = (
        detail_df.groupby(["Role", "InteractionType"])
        .size()
        .reset_index(name="Count")
    )
    interaction_summary.to_csv(geometry_dir / "ligand_interaction_type_summary.csv", index=False)

    interaction_plot = (
        interaction_summary.pivot_table(
            index="Role",
            columns="InteractionType",
            values="Count",
            fill_value=0,
        )
        .sort_index()
    )
    interaction_plot = interaction_plot.div(interaction_plot.sum(axis=1), axis=0).fillna(0.0)
    interaction_plot.to_csv(geometry_dir / "ligand_interaction_type_fraction_by_role.csv")
    interaction_plot.plot(kind="bar", stacked=True, figsize=(10, 5), colormap="tab20")
    plt.ylabel("Fraction of residue-contact events")
    plt.title(f"{ligand_code} interaction type composition by PROTAC role")
    plt.tight_layout()
    plt.savefig(geometry_dir / "ligand_interaction_type_fraction_by_role.png", dpi=300)
    plt.close()

    network_df = (
        detail_df.groupby(
            ["Role", "ProteinName", "ResidueName", "ResidueID", "ResidueLabel"],
            dropna=False,
        )
        .agg(
            ContactFrames=("Frame", "nunique"),
            MinDistance=("MinDistance", "min"),
            InteractionType=("InteractionType", lambda values: ";".join(sorted(set(values)))),
        )
        .reset_index()
    )
    network_df["ContactFraction"] = network_df["ContactFrames"] / total_frames
    network_df["ContactPercent"] = network_df["ContactFraction"] * 100.0
    network_df["EdgeWeight"] = network_df["ContactPercent"]
    network_df.to_csv(networks_dir / "ligand_interaction_network_edges.csv", index=False)

    graph = nx.Graph()
    ligand_node = ligand_code
    graph.add_node(ligand_node, role="Ligand")

    for _, row in network_df.iterrows():
        node = row["ResidueLabel"]
        graph.add_node(node, role=row["Role"])
        graph.add_edge(
            ligand_node,
            node,
            weight=float(row["ContactPercent"]),
            interaction=row["InteractionType"],
        )

    positions = nx.spring_layout(graph, seed=42, weight="weight")
    positions[ligand_node] = np.array([0.0, 0.0])

    plt.figure(figsize=(12, 10))
    edge_widths = [max(1.0, graph[u][v]["weight"] / 10.0) for u, v in graph.edges()]
    node_colors = [
        "#9467bd" if node == ligand_node else role_color(graph.nodes[node].get("role", "Target_A"))
        for node in graph.nodes()
    ]
    node_sizes = [
        2500 if node == ligand_node else 300 + 20 * graph.degree(node)
        for node in graph.nodes()
    ]
    nx.draw_networkx_edges(graph, positions, alpha=0.5, width=edge_widths, edge_color="#999999")
    nx.draw_networkx_nodes(graph, positions, node_color=node_colors, node_size=node_sizes, edgecolors="black")
    nx.draw_networkx_labels(graph, positions, font_size=8)
    plt.title(f"{ligand_code} interaction network by PROTAC role")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(networks_dir / "ligand_interaction_network.png", dpi=300)
    plt.close()


def run_pca(
    universe: mda.Universe,
    protein: mda.AtomGroup,
    frame_indices: Sequence[int],
    qc_dir: Path,
) -> None:
    protein_ca = protein.select_atoms("name CA")
    if protein_ca.n_atoms == 0:
        protein_ca = protein

    aligned_positions = gather_aligned_positions(universe, protein_ca, frame_indices)
    if aligned_positions.size == 0:
        return

    matrix = aligned_positions.reshape(aligned_positions.shape[0], -1)
    matrix = matrix - matrix.mean(axis=0)
    _, _, vt = np.linalg.svd(matrix, full_matrices=False)
    components = vt[:2].T
    projected = np.dot(matrix, components)

    pca_df = pd.DataFrame(
        {
            "Frame": frame_indices,
            "PC1": projected[:, 0],
            "PC2": projected[:, 1],
        }
    )
    pca_df.to_csv(qc_dir / "protein_pca.csv", index=False)

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(projected[:, 0], projected[:, 1], c=frame_indices, cmap="viridis", s=30)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Protein PCA projection")
    plt.colorbar(scatter, label="Frame")
    plt.tight_layout()
    plt.savefig(qc_dir / "protein_pca.png", dpi=300)
    plt.close()


def main() -> int:
    args = parse_args()
    args.pdb_reference = args.pdb_reference or default_pdb_reference()

    outdir = ensure_dir(Path(args.outdir).resolve())
    rms_dir = ensure_dir(outdir / "RMSD_RMSF")
    contacts_dir = ensure_dir(outdir / "Contacts")
    networks_dir = ensure_dir(outdir / "Networks")
    geometry_dir = ensure_dir(outdir / "Geometry")
    qc_dir = ensure_dir(outdir / "QC")
    _ = geometry_dir

    ligand_setup = detect_protac_ligand_setup(
        run_dir=".",
        topo_path=args.topo,
        traj_path=args.traj,
        cli_ligand=args.ligand,
        headless=args.headless,
    )
    write_ligand_preflight_reports(ligand_setup, qc_dir)
    print_preflight_block(ligand_setup)

    if ligand_setup["errors"]:
        for error in ligand_setup["errors"]:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    if not ligand_setup["ligand_code"]:
        print(
            "ERROR: Ligand code could not be determined. Pass --ligand <RESNAME>.",
            file=sys.stderr,
        )
        return 1

    entries = assign_protac_roles(parse_atomindex_entries(args.atomindex))

    try:
        universe = mda.Universe(args.topo, args.traj)
    except Exception as exc:
        print(f"ERROR: Could not load topology/trajectory: {exc}", file=sys.stderr)
        return 1

    ligand = select_ligand(universe, ligand_setup["ligand_code"])
    protein = universe.select_atoms("protein")

    try:
        validate_protac_inputs(universe, ligand, entries)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    start, end, step, frame_indices = resolve_frame_window(
        universe.trajectory.n_frames,
        args.start_frame,
        args.end_frame,
        args.frame_step,
    )
    print_startup_summary(
        ligand_setup["ligand_code"],
        entries,
        args.topo,
        args.traj,
        universe.trajectory.n_frames,
    )
    print(f"Analyzing frames {start}:{end}:{step} ({len(frame_indices)} sampled frames)")

    write_system_roles(entries, qc_dir)
    write_ligand_detection_summary(ligand, ligand_setup["ligand_code"], qc_dir)

    if args.dry_run:
        print("Dry-run complete. QC outputs were written; no full analysis was executed.")
        return 0

    residue_records = build_protein_residue_records(universe, protein, entries)
    if not residue_records:
        print(
            "ERROR: No protein residues could be mapped onto atomIndex roles. "
            "Check that atomIndex.txt matches the analyzed topology.",
            file=sys.stderr,
        )
        return 1

    run_rmsd_and_rmsf(universe, protein, ligand, entries, frame_indices, rms_dir)
    run_radius_of_gyration(universe, protein, entries, frame_indices, rms_dir)
    run_contact_analysis(
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
    )
    run_pca(universe, protein, frame_indices, qc_dir)

    print(f"Analysis complete. Outputs written under {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
