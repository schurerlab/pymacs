#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import math
import os
import re
import sys
import time
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cbook as mpl_cbook
import numpy as np
import pandas as pd
import seaborn as sns

from pymacs_analysis_cache import (
    compute_signature,
    describe_output_files,
    fingerprint_files,
    load_npz_dict,
    read_json_manifest,
    record_stage_completion,
    save_npz_atomic,
    validate_stage_cache,
    write_json_atomic,
)
from pymacs_progress import ProgressReporter
from pymacs_interaction_chemistry import (
    EVENT_COLUMNS,
    aggregate_event_occupancy,
    classify_protein_residue_pair_events,
)

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    tqdm = None

try:
    import networkx as nx
except Exception:  # pragma: no cover
    nx = None

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance, distance_array
except Exception:  # pragma: no cover
    mda = None
    capped_distance = None
    distance_array = None


if not hasattr(mpl_cbook, "iterable"):
    def _mpl_iterable_compat(obj):
        if isinstance(obj, (str, bytes)):
            return False
        try:
            iter(obj)
            return True
        except TypeError:
            return False

    mpl_cbook.iterable = _mpl_iterable_compat


HYDROPHOBIC = {"ALA", "VAL", "ILE", "LEU", "MET", "PRO"}
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
POSITIVE = {"ARG", "LYS", "HIS"}
NEGATIVE = {"ASP", "GLU"}
POLAR = {"SER", "THR", "ASN", "GLN", "CYS"}
AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def build_parser():
    parser = argparse.ArgumentParser(
        description="Standalone protein-protein / protein-peptide interface residue interaction network analysis."
    )
    parser.add_argument("--topo", default="Final_Trajectory.pdb")
    parser.add_argument("--traj", default="Final_Trajectory.xtc")
    parser.add_argument("--atomindex", default="atomIndex.txt")
    parser.add_argument("--outdir", default="Analysis_Results/Interface_RIN")
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    parser.add_argument("--min-contact-frac", type=float, default=0.10)
    parser.add_argument("--frame-step", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=1000,
                        help="Frames per resumable interface scan shard (default: 1000).")
    parser.add_argument("--max-edges", type=int, default=100)
    parser.add_argument("--chain-pairs", default=None)
    parser.add_argument("--peptide-chains", default="", help="Comma-separated current chain IDs explicitly designated as peptides.")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--headless", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def parse_atomindex_entries(path="atomIndex.txt"):
    entries = []
    atomindex_path = Path(path)
    if not atomindex_path.exists():
        return entries

    with atomindex_path.open() as handle:
        for raw in handle:
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

            match = re.match(r"^(\S+)\s+(\d+)\s*-\s*(\d+)$", line)
            if not match:
                print(f"⚠️ Skipping unparsable atomIndex line: {raw.rstrip()}")
                continue

            meta = {}
            for key, value in re.findall(r"(\w+)=([^\s#]+)", comment):
                meta[key] = value

            entries.append(
                {
                    "label": match.group(1),
                    "start": int(match.group(2)) - 1,
                    "end": int(match.group(3)) - 1,
                    "current_chain": meta.get("current_chain"),
                    "chain_type": (meta.get("chain_type") or "mixed_or_unknown").strip().lower(),
                    "current_residues": meta.get("current_residues"),
                    "source_chain": meta.get("source_chain"),
                    "source_residues": meta.get("source_residues"),
                    "source_rescount": meta.get("source_rescount"),
                }
            )
    return entries


def normalize_chain_type(chain_type):
    return (chain_type or "mixed_or_unknown").strip().lower()


def is_polymer_entry(entry):
    chain_type = normalize_chain_type(entry.get("chain_type"))
    if chain_type in {"", "mixed_or_unknown"}:
        return True
    return chain_type not in {"ligand", "small_molecule", "cofactor", "ion", "water", "solvent", "lipid"}


def chain_display_label(entry, index):
    """Keep the submitted component name visible beside its simulation chain."""
    name = entry.get("label") or entry.get("source_chain") or f"Chain_{index}"
    chain = entry.get("current_chain") or entry.get("source_chain") or "?"
    role = "peptide" if normalize_chain_type(entry.get("chain_type")) == "peptide" else "protein"
    return f"{name} ({role}, chain {chain})"


def apply_peptide_chain_hints(entries, peptide_chains):
    """Mark known peptide chains without guessing from their length."""
    requested = {item.strip() for item in str(peptide_chains or "").split(",") if item.strip()}
    for entry in entries:
        if entry.get("current_chain") in requested:
            entry["chain_type"] = "peptide"
    return entries


def detect_polymer_entries(atomindex_entries):
    polymer_entries = []
    for idx, entry in enumerate(atomindex_entries, start=1):
        if not is_polymer_entry(entry):
            continue
        polymer_entry = dict(entry)
        polymer_entry["display_label"] = chain_display_label(entry, idx)
        polymer_entry["chain_type"] = normalize_chain_type(polymer_entry.get("chain_type"))
        polymer_entries.append(polymer_entry)
    return polymer_entries


def load_universe(topo_path, traj_path):
    print(f"🧬 Topology loaded: {topo_path}")
    print(f"🧬 Trajectory loaded: {traj_path}")
    return mda.Universe(str(topo_path), str(traj_path))


def atomgroup_for_entry(universe, entry):
    n_atoms = len(universe.atoms)
    start = max(0, int(entry["start"]))
    end = min(n_atoms - 1, int(entry["end"]))
    if end < start:
        return universe.atoms[[]]
    return universe.atoms[start:end + 1]


def atom_is_hydrogen(atom):
    try:
        element = (atom.element or "").strip().upper()
        if element == "H":
            return True
    except Exception:
        pass
    name = (getattr(atom, "name", "") or "").strip().upper()
    if re.match(r"^\d*H", name):
        return True
    return False


def heavy_atomgroup(atomgroup):
    heavy_indices = [atom.index for atom in atomgroup if not atom_is_hydrogen(atom)]
    if heavy_indices:
        return atomgroup.universe.atoms[heavy_indices]
    return atomgroup


def residue_number(residue):
    for attr in ("resnum", "resid"):
        value = getattr(residue, attr, None)
        if value is not None:
            try:
                return int(value)
            except Exception:
                pass
    return int(residue.ix) + 1


def residue_label(chain_label, residue):
    return f"{chain_label}:{residue.resname}{residue_number(residue)}"


def classify_residue_pair(residue_a, residue_b, frame=0, time_ns=0.0):
    """Compatibility wrapper returning evidence-backed interface labels."""
    events = classify_protein_residue_pair_events(residue_a, residue_b, frame, time_ns)
    validated = sorted({event["InteractionClass"] for event in events if event["ValidationStatus"] == "validated"})
    return ";".join(validated) if validated else "Other Contact"


def build_chain_data(universe, entry):
    atomgroup = atomgroup_for_entry(universe, entry)
    heavy_group = heavy_atomgroup(atomgroup)
    residues = []
    atom_to_residue = np.full(heavy_group.n_atoms, -1, dtype=int)

    for residue in heavy_group.residues:
        residue_atom_positions = [i for i, atom in enumerate(heavy_group) if atom.residue.ix == residue.ix]
        if not residue_atom_positions:
            continue
        residue_index = len(residues)
        for pos in residue_atom_positions:
            atom_to_residue[pos] = residue_index
        residues.append(
            {
                "residue": residue,
                "label": residue_label(entry["display_label"], residue),
                "resname": residue.resname,
                "resnum": residue_number(residue),
            }
        )

    return {
        "entry": entry,
        "atomgroup": atomgroup,
        "heavy_group": heavy_group,
        "residues": residues,
        "atom_to_residue": atom_to_residue,
    }


def parse_requested_chain_pairs(chain_pairs_arg, chain_entries):
    if not chain_pairs_arg:
        return [(left, right) for left, right in itertools.combinations(chain_entries, 2)]

    alias_map = {}
    for entry in chain_entries:
        aliases = {
            entry["display_label"],
            entry.get("label"),
            entry.get("current_chain"),
            entry.get("source_chain"),
        }
        for alias in aliases:
            if alias:
                alias_map[alias.strip()] = entry

    selected = []
    seen = set()
    for token in chain_pairs_arg.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            print(f"⚠️ Skipping invalid chain pair token: {token}")
            continue
        left_key, right_key = [part.strip() for part in token.split(":", 1)]
        left = alias_map.get(left_key)
        right = alias_map.get(right_key)
        if left is None or right is None:
            print(f"⚠️ Skipping unknown chain pair token: {token}")
            continue
        if left["display_label"] == right["display_label"]:
            print(f"⚠️ Skipping self-pair token: {token}")
            continue
        pair_key = tuple(sorted([left["display_label"], right["display_label"]]))
        if pair_key in seen:
            continue
        seen.add(pair_key)
        selected.append((left, right))
    return selected


def iter_frame_indices(n_frames, frame_step):
    return list(range(0, n_frames, max(1, int(frame_step))))


def interface_numerical_output_paths(outdir: Path) -> list[Path]:
    return [
        outdir / "interface_chain_pair_summary.csv",
        outdir / "interface_contact_timeseries.csv",
        outdir / "interface_residue_pair_contacts.csv",
        outdir / "interface_node_summary.csv",
        outdir / "interface_frame_index.csv",
        outdir / "interface_residue_pair_framewise.csv.gz",
        outdir / "interface_residue_pair_presence.npz",
        outdir / "ValidatedInterfaceInteractionEvents.csv",
        outdir / "ValidatedInterfaceInteractionOccupancy.csv",
    ]


def interface_plot_output_paths(outdir: Path) -> list[Path]:
    return [
        outdir / "Interface_ChainPair_ContactFractions.png",
        outdir / "Interface_ContactCount_Timeseries.png",
        outdir / "Interface_MinDistance_Timeseries.png",
        outdir / "Interface_ResiduePair_ContactHeatmap.png",
        outdir / "Interface_RIN_Network.png",
        outdir / "Interface_ChainLevel_Network.png",
        outdir / "Interface_TopResiduePairs_Barplot.png",
        outdir / "Interface_RIN_Explorer.html",
    ]


def build_interface_frame_index(frame_indices: Sequence[int], analyzed_times_ns: Sequence[float]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Frame": list(frame_indices),
            "Time_ns": [float(value) for value in analyzed_times_ns],
        }
    )


def make_pair_id(chain_a_label, chain_b_label, residue_a_label, residue_b_label) -> str:
    return f"{chain_a_label}__{chain_b_label}__{residue_a_label}__{residue_b_label}"


def build_presence_payload(
    residue_df: pd.DataFrame,
    pair_binary_series: Dict[str, Sequence[int]],
    frame_indices: Sequence[int],
    analyzed_times_ns: Sequence[float],
) -> Dict[str, np.ndarray]:
    frames = len(frame_indices)
    rows = []
    pair_ids = []
    chain_a_labels = []
    chain_b_labels = []
    residue_a_labels = []
    residue_b_labels = []

    for _, row in residue_df.iterrows():
        series_key = str(row.get("_series_key", ""))
        binary_series = np.asarray(pair_binary_series.get(series_key, [0] * frames), dtype=np.uint8)
        if binary_series.shape[0] != frames:
            raise ValueError(
                f"Residue-pair series length mismatch for {series_key}: "
                f"{binary_series.shape[0]} vs expected {frames}"
            )
        rows.append(binary_series)
        pair_ids.append(str(row["Pair_ID"]))
        chain_a_labels.append(str(row["ChainA_Label"]))
        chain_b_labels.append(str(row["ChainB_Label"]))
        residue_a_labels.append(str(row["ResidueA"]))
        residue_b_labels.append(str(row["ResidueB"]))

    contact_matrix = np.vstack(rows).astype(np.uint8, copy=False) if rows else np.zeros((0, frames), dtype=np.uint8)
    return {
        "contact_matrix": contact_matrix,
        "frame_indices": np.asarray(frame_indices, dtype=np.int64),
        "times_ns": np.asarray(analyzed_times_ns, dtype=np.float64),
        "pair_ids": np.asarray(pair_ids),
        "chain_a_labels": np.asarray(chain_a_labels),
        "chain_b_labels": np.asarray(chain_b_labels),
        "residue_a_labels": np.asarray(residue_a_labels),
        "residue_b_labels": np.asarray(residue_b_labels),
    }


def current_data_signature(input_files: Dict[str, str], data_parameters: Dict[str, Any]) -> str:
    return compute_signature(
        {
            "inputs": fingerprint_files(input_files),
            "parameters": data_parameters,
        }
    )


def interface_plot_parameters(args) -> Dict[str, Any]:
    return {
        "min_contact_frac": float(args.min_contact_frac),
        "max_edges": int(args.max_edges),
    }


def load_cached_interface_outputs(outdir: Path) -> Dict[str, Any]:
    residue_df = pd.read_csv(outdir / "interface_residue_pair_contacts.csv")
    if "Pair_ID" not in residue_df.columns:
        residue_df["Pair_ID"] = [
            make_pair_id(
                row["ChainA_Label"],
                row["ChainB_Label"],
                row["ResidueA"],
                row["ResidueB"],
            )
            for _, row in residue_df.iterrows()
        ]
    return {
        "summary_df": pd.read_csv(outdir / "interface_chain_pair_summary.csv"),
        "timeseries_df": pd.read_csv(outdir / "interface_contact_timeseries.csv"),
        "residue_df": residue_df,
        "node_df": pd.read_csv(outdir / "interface_node_summary.csv"),
        "frame_index_df": pd.read_csv(outdir / "interface_frame_index.csv"),
        "presence_payload": load_npz_dict(outdir / "interface_residue_pair_presence.npz"),
    }


def write_interface_outputs(
    outdir: Path,
    summary_df: pd.DataFrame,
    timeseries_df: pd.DataFrame,
    residue_df: pd.DataFrame,
    node_df: pd.DataFrame,
    frame_index_df: pd.DataFrame,
    framewise_df: pd.DataFrame,
    presence_payload: Dict[str, np.ndarray],
    runtime_lines: Sequence[str],
) -> None:
    summary_df.to_csv(outdir / "interface_chain_pair_summary.csv", index=False)
    timeseries_df.to_csv(outdir / "interface_contact_timeseries.csv", index=False)
    residue_export_df = residue_df.drop(columns=["_pair_key", "_pair_id", "_series_key"], errors="ignore")
    residue_export_df.to_csv(outdir / "interface_residue_pair_contacts.csv", index=False)
    node_df.to_csv(outdir / "interface_node_summary.csv", index=False)
    frame_index_df.to_csv(outdir / "interface_frame_index.csv", index=False)
    framewise_df.to_csv(outdir / "interface_residue_pair_framewise.csv.gz", index=False, compression="gzip")
    save_npz_atomic(outdir / "interface_residue_pair_presence.npz", **presence_payload)
    (outdir / "interface_runtime_summary.txt").write_text("\n".join(runtime_lines) + "\n")


# The work directory deliberately contains only complete shard files.  A shard
# is written to a temporary name and renamed atomically, so an interrupted job
# can safely resume without treating a partial CSV as analysis evidence.
FRAMEWISE_COLUMNS = ["Pair_ID", "Frame", "Time_ns", "ChainA_Label", "ChainB_Label", "ResidueA", "ResidueB",
                     "ResidueA_Name", "ResidueB_Name", "ResidueA_Index", "ResidueB_Index", "MinDistance_A", "InteractionClass"]
TIMESERIES_COLUMNS = ["Time_ns", "Frame", "ChainA_Label", "ChainB_Label", "ContactCount", "MinDistance_A", "InContact"]


def _atomic_csv(path, rows, columns):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    pd.DataFrame(rows, columns=columns).to_csv(tmp, index=False)
    os.replace(tmp, path)


def _append_csv_shards(shards, destination, columns, gzip_output=False):
    opener = gzip.open if gzip_output else open
    mode = "wt"
    with opener(destination, mode, newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=columns)
        writer.writeheader()
        for shard in shards:
            with open(shard, newline="", encoding="utf-8") as source:
                for row in csv.DictReader(source):
                    writer.writerow({name: row.get(name, "") for name in columns})


def _work_signature(data_signature, args):
    return compute_signature({"data_signature": data_signature, "chunk_size": int(args.chunk_size)})


def scan_interface_shards(universe, chain_pairs, frame_indices, cutoff_ang, workdir, state):
    """Scan bounded frame ranges, retaining only one shard's Python rows."""
    workdir = Path(workdir)
    prepared = []
    for left_entry, right_entry in chain_pairs:
        left, right = build_chain_data(universe, left_entry), build_chain_data(universe, right_entry)
        if left["heavy_group"].n_atoms and right["heavy_group"].n_atoms and left["residues"] and right["residues"]:
            prepared.append((left_entry, right_entry, left, right))
    size = max(1, int(state["chunk_size"]))
    completed = set(state.get("completed_chunks", []))
    for chunk_no, start in enumerate(range(0, len(frame_indices), size)):
        stop = min(len(frame_indices), start + size)
        token = f"{chunk_no:06d}"
        expected = [workdir / f"framewise_{token}.csv", workdir / f"timeseries_{token}.csv", workdir / f"events_{token}.csv"]
        if token in completed and all(path.exists() for path in expected):
            continue
        frame_rows, time_rows, event_rows = [], [], []
        for global_pos, frame_index in enumerate(frame_indices[start:stop], start):
            universe.trajectory[frame_index]
            time_ns = float(universe.trajectory.time) / 1000.0
            for left_entry, right_entry, left, right in prepared:
                pair_label = f"{left_entry['display_label']}__{right_entry['display_label']}"
                contacts, distances = capped_distance(left["heavy_group"].positions, right["heavy_group"].positions,
                                                      max_cutoff=cutoff_ang, return_distances=True)
                minima = {}
                for (i, j), dist in zip(contacts, distances):
                    a, b = left["atom_to_residue"][i], right["atom_to_residue"][j]
                    if a >= 0 and b >= 0 and ((a, b) not in minima or dist < minima[(a, b)]):
                        minima[(a, b)] = float(dist)
                minimum = float(np.min(distances)) if len(distances) else (float(np.min(distance_array(left["heavy_group"].positions, right["heavy_group"].positions))) if left["heavy_group"].n_atoms and right["heavy_group"].n_atoms else math.nan)
                for (a, b), dist in minima.items():
                    ra, rb = left["residues"][a], right["residues"][b]
                    pair_id = make_pair_id(left_entry["display_label"], right_entry["display_label"], ra["label"], rb["label"])
                    cls = classify_residue_pair(ra["residue"], rb["residue"], frame_index, time_ns)
                    frame_rows.append({"Pair_ID": pair_id, "Frame": frame_index, "Time_ns": time_ns,
                                       "ChainA_Label": left_entry["display_label"], "ChainB_Label": right_entry["display_label"],
                                       "ResidueA": ra["label"], "ResidueB": rb["label"], "ResidueA_Name": ra["resname"], "ResidueB_Name": rb["resname"],
                                       "ResidueA_Index": ra["resnum"], "ResidueB_Index": rb["resnum"], "MinDistance_A": dist, "InteractionClass": cls})
                    for event in classify_protein_residue_pair_events(ra["residue"], rb["residue"], frame_index, time_ns):
                        event = dict(event); event.update({"Partner": ra["label"], "Component": pair_label, "InterfaceResidueA": ra["label"], "InterfaceResidueB": rb["label"]})
                        event_rows.append(event)
                time_rows.append({"Time_ns": time_ns, "Frame": frame_index, "ChainA_Label": left_entry["display_label"], "ChainB_Label": right_entry["display_label"], "ContactCount": len(minima), "MinDistance_A": minimum, "InContact": int(bool(minima))})
        _atomic_csv(expected[0], frame_rows, FRAMEWISE_COLUMNS)
        _atomic_csv(expected[1], time_rows, TIMESERIES_COLUMNS)
        _atomic_csv(expected[2], event_rows, EVENT_COLUMNS)
        state.setdefault("completed_chunks", []).append(token)
        write_json_atomic(workdir / "progress.json", state)


def reduce_interface_shards(outdir, workdir, frame_indices, analyzed_times_ns, chain_pairs, cutoff_ang):
    """Reduce shard CSVs into public artifacts without retaining raw rows."""
    workdir, outdir = Path(workdir), Path(outdir)
    frame_shards = sorted(workdir.glob("framewise_*.csv")); time_shards = sorted(workdir.glob("timeseries_*.csv")); event_shards = sorted(workdir.glob("events_*.csv"))
    _append_csv_shards(frame_shards, outdir / "interface_residue_pair_framewise.csv.gz", FRAMEWISE_COLUMNS, gzip_output=True)
    _append_csv_shards(time_shards, outdir / "interface_contact_timeseries.csv", TIMESERIES_COLUMNS)
    _append_csv_shards(event_shards, outdir / "ValidatedInterfaceInteractionEvents.csv", EVENT_COLUMNS)
    stats, metadata = {}, {}
    for shard in frame_shards:
        for row in pd.read_csv(shard).to_dict("records"):
            key = row["Pair_ID"]; stats.setdefault(key, [0, 0.0, float("inf")]); stats[key][0] += 1; stats[key][1] += float(row["MinDistance_A"]); stats[key][2] = min(stats[key][2], float(row["MinDistance_A"])); metadata[key] = row
    residue_rows = []
    for key, values in stats.items():
        row = metadata[key]; n = values[0]
        residue_rows.append({"Pair_ID": key, "ChainA_Label": row["ChainA_Label"], "ChainB_Label": row["ChainB_Label"], "ResidueA": row["ResidueA"], "ResidueB": row["ResidueB"], "ResidueA_Name": row["ResidueA_Name"], "ResidueB_Name": row["ResidueB_Name"], "ResidueA_Index": row["ResidueA_Index"], "ResidueB_Index": row["ResidueB_Index"], "ContactFrames": n, "Frames_Analyzed": len(frame_indices), "ContactFraction": n / max(1, len(frame_indices)), "MeanMinDistance_A": values[1] / n, "MinObservedDistance_A": values[2], "InteractionClass": "Evidence in interface_residue_pair_framewise.csv"})
    residue_df = pd.DataFrame(residue_rows).sort_values(["ContactFraction", "ContactFrames"], ascending=False) if residue_rows else pd.DataFrame(columns=["Pair_ID"])
    residue_df.to_csv(outdir / "interface_residue_pair_contacts.csv", index=False)
    pair_to_row = {str(row.Pair_ID): i for i, row in residue_df.iterrows()}
    matrix = np.zeros((len(residue_df), len(frame_indices)), dtype=np.uint8); frame_to_col = {int(frame): i for i, frame in enumerate(frame_indices)}
    for shard in frame_shards:
        for row in pd.read_csv(shard, usecols=["Pair_ID", "Frame"]).itertuples(index=False): matrix[pair_to_row[str(row.Pair_ID)], frame_to_col[int(row.Frame)]] = 1
    save_npz_atomic(outdir / "interface_residue_pair_presence.npz", contact_matrix=matrix, frame_indices=np.asarray(frame_indices, dtype=np.int64), times_ns=np.asarray(analyzed_times_ns, dtype=float), pair_ids=residue_df.get("Pair_ID", pd.Series(dtype=str)).astype(str).to_numpy(), chain_a_labels=residue_df.get("ChainA_Label", pd.Series(dtype=str)).astype(str).to_numpy(), chain_b_labels=residue_df.get("ChainB_Label", pd.Series(dtype=str)).astype(str).to_numpy(), residue_a_labels=residue_df.get("ResidueA", pd.Series(dtype=str)).astype(str).to_numpy(), residue_b_labels=residue_df.get("ResidueB", pd.Series(dtype=str)).astype(str).to_numpy())
    frame_index_df = build_interface_frame_index(frame_indices, analyzed_times_ns); frame_index_df.to_csv(outdir / "interface_frame_index.csv", index=False)
    timeseries_df = pd.read_csv(outdir / "interface_contact_timeseries.csv")
    summary = timeseries_df.groupby(["ChainA_Label", "ChainB_Label"], as_index=False).agg(Frames_Analyzed=("Frame", "size"), Frames_In_Contact=("InContact", "sum"), ContactFraction=("InContact", "mean"), MeanContactsPerFrame=("ContactCount", "mean"), MaxContactsInFrame=("ContactCount", "max"), MeanMinDistance_A=("MinDistance_A", "mean"), MinObservedDistance_A=("MinDistance_A", "min"))
    types = {entry["display_label"]: entry.get("chain_type") for pair in chain_pairs for entry in pair}; summary["ChainA_Type"] = summary["ChainA_Label"].map(types); summary["ChainB_Type"] = summary["ChainB_Label"].map(types); summary["ContactCutoff_A"] = float(cutoff_ang); summary = summary.sort_values(["ContactFraction", "MeanContactsPerFrame"], ascending=False); summary.to_csv(outdir / "interface_chain_pair_summary.csv", index=False)
    node_df = build_node_summary(residue_df); node_df.to_csv(outdir / "interface_node_summary.csv", index=False)
    events = pd.read_csv(outdir / "ValidatedInterfaceInteractionEvents.csv")
    aggregate_event_occupancy(events[events.get("ValidationStatus", pd.Series(dtype=str)) == "validated"], len(frame_indices)).to_csv(outdir / "ValidatedInterfaceInteractionOccupancy.csv", index=False)
    return summary, timeseries_df, residue_df, node_df, frame_index_df, {"contact_matrix": matrix, "frame_indices": np.asarray(frame_indices), "times_ns": np.asarray(analyzed_times_ns)}, events


def analyze_chain_pair_contacts(universe, left_entry, right_entry, frame_indices, cutoff_ang):
    left = build_chain_data(universe, left_entry)
    right = build_chain_data(universe, right_entry)

    if left["heavy_group"].n_atoms == 0 or right["heavy_group"].n_atoms == 0:
        return None
    if not left["residues"] or not right["residues"]:
        return None

    pair_id = f"{left_entry['display_label']}__{right_entry['display_label']}"
    residue_pair_stats = {}
    residue_pair_series = defaultdict(list)
    time_rows = []
    framewise_rows = []
    validated_event_rows = []
    analyzed_times_ns = []

    iterator = frame_indices
    if tqdm is not None:
        iterator = tqdm(frame_indices, desc=f"{left_entry['display_label']}:{right_entry['display_label']}", leave=False)

    for frame_counter, frame_index in enumerate(iterator):
        universe.trajectory[frame_index]
        analyzed_times_ns.append(float(universe.trajectory.time) / 1000.0)
        pos_left = left["heavy_group"].positions
        pos_right = right["heavy_group"].positions

        contact_pairs, distances = capped_distance(
            pos_left,
            pos_right,
            max_cutoff=cutoff_ang,
            return_distances=True,
        )

        frame_pair_min = {}
        min_distance_ang = None

        if len(contact_pairs) > 0:
            min_distance_ang = float(np.min(distances))
            for (idx_left, idx_right), dist in zip(contact_pairs, distances):
                res_left = left["atom_to_residue"][idx_left]
                res_right = right["atom_to_residue"][idx_right]
                if res_left < 0 or res_right < 0:
                    continue
                pair_key = (res_left, res_right)
                current = frame_pair_min.get(pair_key)
                if current is None or dist < current:
                    frame_pair_min[pair_key] = float(dist)
        else:
            full_distances = distance_array(pos_left, pos_right)
            min_distance_ang = float(np.min(full_distances)) if full_distances.size else math.nan

        for pair_key, dist in frame_pair_min.items():
            if pair_key not in residue_pair_series:
                residue_pair_series[pair_key] = [0] * frame_counter
            stats = residue_pair_stats.setdefault(
                pair_key,
                {"contact_frames": 0, "distance_sum": 0.0, "min_distance": float("inf")},
            )
            stats["contact_frames"] += 1
            stats["distance_sum"] += dist
            stats["min_distance"] = min(stats["min_distance"], dist)
            res_left = left["residues"][pair_key[0]]
            res_right = right["residues"][pair_key[1]]
            pair_id_value = make_pair_id(
                left_entry["display_label"],
                right_entry["display_label"],
                res_left["label"],
                res_right["label"],
            )
            interaction_class = classify_residue_pair(
                res_left["residue"], res_right["residue"], frame_index,
                float(universe.trajectory.time) / 1000.0,
            )
            for event in classify_protein_residue_pair_events(
                res_left["residue"], res_right["residue"], frame_index,
                float(universe.trajectory.time) / 1000.0,
            ):
                event = dict(event)
                event["Partner"] = res_left["label"]
                event["Component"] = pair_id
                event["InterfaceResidueA"] = res_left["label"]
                event["InterfaceResidueB"] = res_right["label"]
                validated_event_rows.append(event)
            framewise_rows.append(
                {
                    "Pair_ID": pair_id_value,
                    "Frame": int(frame_index),
                    "Time_ns": float(universe.trajectory.time) / 1000.0,
                    "ChainA_Label": left_entry["display_label"],
                    "ChainB_Label": right_entry["display_label"],
                    "ResidueA": res_left["label"],
                    "ResidueB": res_right["label"],
                    "ResidueA_Name": res_left["resname"],
                    "ResidueB_Name": res_right["resname"],
                    "ResidueA_Index": res_left["resnum"],
                    "ResidueB_Index": res_right["resnum"],
                    "MinDistance_A": float(dist),
                    "InteractionClass": interaction_class,
                }
            )

        for pair_key in residue_pair_stats:
            residue_pair_series[pair_key].append(1 if pair_key in frame_pair_min else 0)

        contact_count = len(frame_pair_min)
        time_rows.append(
            {
                "Time_ns": float(universe.trajectory.time) / 1000.0,
                "Frame": int(frame_index),
                "ChainA_Label": left_entry["display_label"],
                "ChainB_Label": right_entry["display_label"],
                "ContactCount": int(contact_count),
                "MinDistance_A": min_distance_ang,
                "InContact": int(contact_count > 0),
            }
        )

    frames_analyzed = len(frame_indices)
    residue_rows = []
    for pair_key, stats in residue_pair_stats.items():
        res_left = left["residues"][pair_key[0]]
        res_right = right["residues"][pair_key[1]]
        pair_id_value = make_pair_id(
            left_entry["display_label"],
            right_entry["display_label"],
            res_left["label"],
            res_right["label"],
        )
        residue_rows.append(
            {
                "Pair_ID": pair_id_value,
                "ChainA_Label": left_entry["display_label"],
                "ChainB_Label": right_entry["display_label"],
                "ResidueA": res_left["label"],
                "ResidueB": res_right["label"],
                "ResidueA_Name": res_left["resname"],
                "ResidueB_Name": res_right["resname"],
                "ResidueA_Index": res_left["resnum"],
                "ResidueB_Index": res_right["resnum"],
                "ContactFrames": int(stats["contact_frames"]),
                "Frames_Analyzed": int(frames_analyzed),
                "ContactFraction": float(stats["contact_frames"] / frames_analyzed) if frames_analyzed else 0.0,
                "MeanMinDistance_A": float(stats["distance_sum"] / stats["contact_frames"]),
                "MinObservedDistance_A": float(stats["min_distance"]),
                "InteractionClass": "Evidence in interface_residue_pair_framewise.csv",
                "_pair_key": pair_key,
                "_pair_id": pair_id,
                "_series_key": f"{pair_id}::{pair_key[0]}::{pair_key[1]}",
            }
        )

    in_contact_values = [row["InContact"] for row in time_rows]
    contact_counts = [row["ContactCount"] for row in time_rows]
    min_distances = [row["MinDistance_A"] for row in time_rows if not math.isnan(row["MinDistance_A"])]

    summary_row = {
        "ChainA_Label": left_entry["display_label"],
        "ChainB_Label": right_entry["display_label"],
        "ChainA_Type": left_entry["chain_type"],
        "ChainB_Type": right_entry["chain_type"],
        "Frames_Analyzed": int(frames_analyzed),
        "Frames_In_Contact": int(sum(in_contact_values)),
        "ContactFraction": float(np.mean(in_contact_values)) if in_contact_values else 0.0,
        "MeanContactsPerFrame": float(np.mean(contact_counts)) if contact_counts else 0.0,
        "MaxContactsInFrame": int(np.max(contact_counts)) if contact_counts else 0,
        "MeanMinDistance_A": float(np.mean(min_distances)) if min_distances else math.nan,
        "MinObservedDistance_A": float(np.min(min_distances)) if min_distances else math.nan,
        "ContactCutoff_A": float(cutoff_ang),
    }

    return {
        "summary_row": summary_row,
        "time_rows": time_rows,
        "residue_rows": residue_rows,
        "framewise_rows": framewise_rows,
        "validated_event_rows": validated_event_rows,
        "pair_binary_series": residue_pair_series,
        "times_ns": analyzed_times_ns,
    }


def build_node_summary(residue_df):
    if residue_df.empty:
        return pd.DataFrame(
            columns=[
                "Chain_Label",
                "Residue_Label",
                "Residue_Name",
                "Residue_Index",
                "Degree",
                "WeightedDegree",
                "MaxContactFraction",
                "MeanContactFraction",
            ]
        )

    node_stats = defaultdict(lambda: {"Chain_Label": None, "Residue_Name": None, "Residue_Index": None, "fractions": []})

    for _, row in residue_df.iterrows():
        left_key = row["ResidueA"]
        right_key = row["ResidueB"]

        node_stats[left_key]["Chain_Label"] = row["ChainA_Label"]
        node_stats[left_key]["Residue_Name"] = row["ResidueA_Name"]
        node_stats[left_key]["Residue_Index"] = row["ResidueA_Index"]
        node_stats[left_key]["fractions"].append(float(row["ContactFraction"]))

        node_stats[right_key]["Chain_Label"] = row["ChainB_Label"]
        node_stats[right_key]["Residue_Name"] = row["ResidueB_Name"]
        node_stats[right_key]["Residue_Index"] = row["ResidueB_Index"]
        node_stats[right_key]["fractions"].append(float(row["ContactFraction"]))

    rows = []
    for residue_key, stats in node_stats.items():
        fractions = stats["fractions"]
        rows.append(
            {
                "Chain_Label": stats["Chain_Label"],
                "Residue_Label": residue_key,
                "Residue_Name": stats["Residue_Name"],
                "Residue_Index": stats["Residue_Index"],
                "Degree": int(len(fractions)),
                "WeightedDegree": float(sum(fractions)),
                "MaxContactFraction": float(max(fractions)),
                "MeanContactFraction": float(sum(fractions) / len(fractions)),
            }
        )
    return pd.DataFrame(rows).sort_values(["WeightedDegree", "Degree"], ascending=[False, False])


def plot_placeholder(path, title, message):
    plt.figure(figsize=(8, 5))
    plt.axis("off")
    plt.title(title)
    plt.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def plot_chain_pair_heatmap(summary_df, outdir):
    path = outdir / "Interface_ChainPair_ContactFractions.png"
    if summary_df.empty:
        plot_placeholder(path, "Interface Chain-Pair Contact Fractions", "No chain pairs were analyzed.")
        return

    labels = sorted(set(summary_df["ChainA_Label"]).union(summary_df["ChainB_Label"]))
    matrix = pd.DataFrame(0.0, index=labels, columns=labels)
    for _, row in summary_df.iterrows():
        matrix.loc[row["ChainA_Label"], row["ChainB_Label"]] = row["ContactFraction"]
        matrix.loc[row["ChainB_Label"], row["ChainA_Label"]] = row["ContactFraction"]
    np.fill_diagonal(matrix.values, np.nan)

    plt.figure(figsize=(max(6, len(labels) * 1.2), max(5, len(labels) * 0.9)))
    sns.heatmap(matrix, annot=True, fmt=".2f", cmap="mako", cbar_kws={"label": "Contact fraction"})
    plt.title("Interface Chain-Pair Contact Fractions")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def plot_contact_timeseries(timeseries_df, outdir):
    path = outdir / "Interface_ContactCount_Timeseries.png"
    if timeseries_df.empty:
        plot_placeholder(path, "Interface Contact Count Over Time", "No time-series data were generated.")
        return

    plt.figure(figsize=(12, 6))
    for (chain_a, chain_b), group in timeseries_df.groupby(["ChainA_Label", "ChainB_Label"]):
        label = f"{chain_a} vs {chain_b}"
        plt.plot(group["Time_ns"], group["ContactCount"], lw=1.4, label=label)
    plt.xlabel("Time (ns)")
    plt.ylabel("Residue-residue contacts")
    plt.title("Interface Contact Count Over Time")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def plot_min_distance_timeseries(timeseries_df, outdir, cutoff_ang):
    path = outdir / "Interface_MinDistance_Timeseries.png"
    if timeseries_df.empty:
        plot_placeholder(path, "Interface Minimum Distance Over Time", "No minimum-distance data were generated.")
        return

    plt.figure(figsize=(12, 6))
    for (chain_a, chain_b), group in timeseries_df.groupby(["ChainA_Label", "ChainB_Label"]):
        label = f"{chain_a} vs {chain_b}"
        plt.plot(group["Time_ns"], group["MinDistance_A"], lw=1.4, label=label)
    plt.axhline(cutoff_ang, color="black", ls="--", lw=1.0, alpha=0.8)
    plt.xlabel("Time (ns)")
    plt.ylabel("Minimum inter-chain distance (Å)")
    plt.title("Interface Minimum Distance Over Time")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def plot_residue_pair_heatmap(residue_df, presence_payload, outdir):
    path = outdir / "Interface_ResiduePair_ContactHeatmap.png"
    if residue_df.empty:
        plot_placeholder(
            path,
            "Top Interface Residue Pairs Over Time",
            "No residue-residue contacts were observed."
        )
        return

    contact_matrix = np.asarray(presence_payload.get("contact_matrix", np.zeros((0, 0), dtype=np.uint8)), dtype=np.uint8)
    analyzed_times_ns = np.asarray(presence_payload.get("times_ns", np.array([], dtype=float)), dtype=float)
    pair_ids = [str(value) for value in np.asarray(presence_payload.get("pair_ids", np.array([], dtype=str)))]
    pair_id_to_index = {pair_id: index for index, pair_id in enumerate(pair_ids)}

    top_df = residue_df.sort_values(
        ["ContactFraction", "ContactFrames"],
        ascending=[False, False]
    ).head(80)
    rows = []
    row_labels = []
    for _, row in top_df.iterrows():
        pair_id = str(row.get("Pair_ID", ""))
        row_index = pair_id_to_index.get(pair_id)
        if row_index is None or row_index >= contact_matrix.shape[0]:
            continue
        row_labels.append(
            f"{row['ResidueA']} -- {row['ResidueB']}"
        )
        rows.append(contact_matrix[row_index])

    if not rows:
        plot_placeholder(
            path,
            "Top Interface Residue Pairs Over Time",
            "No residue-pair time series were available."
        )
        return

    heatmap_df = pd.DataFrame(
        rows,
        index=row_labels,
        columns=[f"{t:.3f}" for t in analyzed_times_ns]
    )

    n_rows, n_cols = heatmap_df.shape

    # Keep the full trajectory data, but do not let figure width scale with frame count.
    fig_width = 18.0
    fig_height = min(14.0, max(6.0, n_rows * 0.20))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        heatmap_df,
        cmap="rocket_r",
        cbar_kws={"label": "Contact present (0/1)"},
        vmin=0,
        vmax=1,
        xticklabels=False,
        ax=ax
    )

    # Show a small number of evenly spaced time labels while retaining all analyzed frames.
    if n_cols > 0:
        n_ticks = min(11, n_cols)
        tick_indices = np.linspace(0, n_cols - 1, n_ticks, dtype=int)
        ax.set_xticks(tick_indices + 0.5)

        if len(analyzed_times_ns) == n_cols:
            tick_labels = [f"{analyzed_times_ns[i]:.1f}" for i in tick_indices]
        else:
            tick_labels = [str(i) for i in tick_indices]

        ax.set_xticklabels(tick_labels, rotation=45, ha="right")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Residue pair")
    ax.set_title("Top Interface Residue-Pair Contacts Over Time")

    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_residue_network(residue_df, outdir, min_contact_frac, max_edges):
    path = outdir / "Interface_RIN_Network.png"
    if nx is None:
        plot_placeholder(path, "Interface Residue Interaction Network", "NetworkX is unavailable in this environment.")
        return
    filtered = residue_df[residue_df["ContactFraction"] >= float(min_contact_frac)].copy()
    filtered = filtered.sort_values(["ContactFraction", "ContactFrames"], ascending=[False, False]).head(int(max_edges))

    if filtered.empty:
        plot_placeholder(path, "Interface Residue Interaction Network", "No persistent edges passed the contact-fraction threshold.")
        return

    graph = nx.Graph()
    chain_labels = sorted(set(filtered["ChainA_Label"]).union(filtered["ChainB_Label"]))
    palette = sns.color_palette("tab10", n_colors=max(1, len(chain_labels)))
    color_map = {label: palette[i] for i, label in enumerate(chain_labels)}

    for _, row in filtered.iterrows():
        graph.add_node(row["ResidueA"], chain=row["ChainA_Label"], color=color_map[row["ChainA_Label"]])
        graph.add_node(row["ResidueB"], chain=row["ChainB_Label"], color=color_map[row["ChainB_Label"]])
        graph.add_edge(
            row["ResidueA"],
            row["ResidueB"],
            weight=float(row["ContactFraction"]),
            label=f"{int(round(row['ContactFraction'] * 100.0))}%",
        )

    pos = nx.spring_layout(graph, seed=42, weight="weight")
    widths = [1.5 + 6.0 * graph.edges[edge]["weight"] for edge in graph.edges]
    node_colors = [graph.nodes[node]["color"] for node in graph.nodes]

    plt.figure(figsize=(14, 10))
    nx.draw_networkx_nodes(graph, pos, node_color=node_colors, node_size=900, alpha=0.95)
    nx.draw_networkx_edges(graph, pos, width=widths, alpha=0.7, edge_color="#4a5568")
    nx.draw_networkx_labels(graph, pos, font_size=8)
    nx.draw_networkx_edge_labels(
        graph,
        pos,
        edge_labels={(u, v): data["label"] for u, v, data in graph.edges(data=True)},
        font_size=7,
    )
    plt.title("Residue-Level Interface Interaction Network")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=350)
    plt.close()


def plot_chain_network(summary_df, outdir):
    path = outdir / "Interface_ChainLevel_Network.png"
    if nx is None:
        plot_placeholder(path, "Interface Chain-Level Network", "NetworkX is unavailable in this environment.")
        return
    if summary_df.empty:
        plot_placeholder(path, "Interface Chain-Level Network", "No chain-pair summary data were generated.")
        return

    graph = nx.Graph()
    for _, row in summary_df.iterrows():
        graph.add_node(row["ChainA_Label"])
        graph.add_node(row["ChainB_Label"])
        graph.add_edge(
            row["ChainA_Label"],
            row["ChainB_Label"],
            weight=float(row["ContactFraction"]),
            label=f"{int(round(row['ContactFraction'] * 100.0))}%",
        )

    pos = nx.spring_layout(graph, seed=42, weight="weight")
    widths = [1.5 + 8.0 * graph.edges[edge]["weight"] for edge in graph.edges]
    plt.figure(figsize=(10, 8))
    nx.draw_networkx_nodes(graph, pos, node_color=sns.color_palette("Set2", n_colors=graph.number_of_nodes()), node_size=1800)
    nx.draw_networkx_edges(graph, pos, width=widths, alpha=0.75, edge_color="#2d3748")
    nx.draw_networkx_labels(graph, pos, font_size=10)
    nx.draw_networkx_edge_labels(
        graph,
        pos,
        edge_labels={(u, v): data["label"] for u, v, data in graph.edges(data=True)},
        font_size=9,
    )
    plt.title("Chain-Level Interface Network")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=350)
    plt.close()


def plot_top_residue_pairs_barplot(residue_df, outdir):
    path = outdir / "Interface_TopResiduePairs_Barplot.png"
    if residue_df.empty:
        plot_placeholder(path, "Top Interface Residue Pairs", "No residue-residue contacts were observed.")
        return

    top_df = residue_df.sort_values(["ContactFraction", "ContactFrames"], ascending=[False, False]).head(20).copy()
    top_df["PairLabel"] = top_df["ResidueA"] + " -- " + top_df["ResidueB"]
    plt.figure(figsize=(12, max(6, len(top_df) * 0.35)))
    sns.barplot(data=top_df, x="ContactFraction", y="PairLabel", palette="crest")
    plt.xlabel("Contact fraction")
    plt.ylabel("Residue pair")
    plt.title("Top Interface Residue Pairs by Persistence")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def interface_chain_sequences(universe, chain_pairs):
    """Return complete labelled residue sequences for the selected interface chains."""
    sequences = {}
    for left, right in chain_pairs:
        for entry in (left, right):
            label = entry["display_label"]
            if label in sequences:
                continue
            residues = atomgroup_for_entry(universe, entry).residues
            seen = set()
            sequence = []
            for residue in residues:
                key = (residue.ix, residue.resname, residue_number(residue))
                if key in seen:
                    continue
                seen.add(key)
                sequence.append({
                    "label": residue_label(label, residue),
                    "name": str(residue.resname).upper(),
                    "number": residue_number(residue),
                    "one": AA1.get(str(residue.resname).upper(), "X"),
                })
            sequences[label] = {
                "label": label,
                "type": entry.get("chain_type", "mixed_or_unknown"),
                "residues": sequence,
            }
    return sequences


def write_interactive_interface_rin_html(outdir, residue_df, events_df, frame_count, chain_pairs, universe):
    """Write an interactive, self-contained peptide/protein interface RIN explorer."""
    outdir = Path(outdir)
    sequences = interface_chain_sequences(universe, chain_pairs)
    records = []
    if events_df is not None and not events_df.empty:
        valid = events_df[events_df.get("ValidationStatus", "") == "validated"].copy()
        for col in ("InterfaceResidueA", "InterfaceResidueB", "InteractionClass", "Frame"):
            if col not in valid:
                valid[col] = pd.Series(dtype=object)
        typed = (valid.dropna(subset=["InterfaceResidueA", "InterfaceResidueB"])
                 .drop_duplicates(["InterfaceResidueA", "InterfaceResidueB", "InteractionClass", "Frame"])
                 .groupby(["InterfaceResidueA", "InterfaceResidueB", "InteractionClass"], as_index=False)["Frame"].nunique()
                 .rename(columns={"InterfaceResidueA": "a", "InterfaceResidueB": "b", "InteractionClass": "type", "Frame": "frames"}))
        for row in typed.itertuples(index=False):
            records.append({"a": row.a, "b": row.b, "type": row.type,
                            "fraction": round(float(row.frames) / max(1, frame_count), 6)})
    # Preserve every observed interface edge, including contacts whose exact
    # chemistry is deliberately reported only as proximity evidence.
    typed_pairs = {(edge["a"], edge["b"]) for edge in records}
    for row in residue_df.itertuples(index=False):
        pair = (str(row.ResidueA), str(row.ResidueB))
        if pair not in typed_pairs:
            records.append({"a": pair[0], "b": pair[1], "type": "Other contact",
                            "fraction": round(float(row.ContactFraction), 6)})

    pair_options = [
        {"id": f"{left['display_label']}__{right['display_label']}", "left": left["display_label"], "right": right["display_label"]}
        for left, right in chain_pairs
    ]
    # Peptide interfaces benefit from the full occupancy range; readers often
    # need to inspect a dispersed, lower-persistence patch along the complete
    # peptide sequence.  Protein-only interfaces begin at 50% to avoid a
    # visually unhelpful hairball of weak residue contacts.
    has_peptide = any(str(item.get("type", "")).casefold() == "peptide" for item in sequences.values())
    thresholds = tuple(step / 100 for step in range(10 if has_peptide else 50, 100, 10))
    default_threshold = 0.50 if has_peptide else 0.70
    threshold_buttons = "".join(
        f'<button data-t="{threshold:.2f}" class="{"active" if threshold == default_threshold else ""}">{int(threshold * 100)}% retained</button>'
        for threshold in thresholds
    )
    payload = json.dumps({
        "pairs": pair_options,
        "chains": sequences,
        "edges": records,
        "thresholds": thresholds,
        "default_threshold": default_threshold,
        "has_peptide": has_peptide,
    }, separators=(",", ":")).replace("</", "<\\/")
    title = "Interface RIN explorer"
    html = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title><style>
:root{{--bg:#05070d;--panel:#0e1625;--line:#2a405e;--text:#ecf5ff;--muted:#9bb0c9;--cyan:#4bd7ff;--pep:#ffcf5b;--hbond:#ff5578;--ionic:#c084fc;--hydrophobic:#7dff91;--other:#94a3b8}}*{{box-sizing:border-box}}body{{margin:0;background:radial-gradient(circle at 90% 0,rgba(75,215,255,.16),transparent 31%),var(--bg);color:var(--text);font-family:Inter,ui-sans-serif,system-ui,sans-serif}}main{{max-width:1500px;margin:auto;padding:28px}}.eyebrow{{color:var(--cyan);font-size:.73rem;font-weight:900;letter-spacing:.14em}}h1{{font-size:clamp(2rem,4vw,4.2rem);letter-spacing:-.06em;margin:.25rem 0}}p{{line-height:1.5;color:var(--muted);max-width:850px}}.controls{{display:flex;gap:9px;flex-wrap:wrap;margin:20px 0}}button{{border:1px solid var(--line);border-radius:10px;background:#101d30;color:var(--text);font:800 .88rem inherit;padding:10px 13px;cursor:pointer}}button:hover,button.active{{border-color:var(--cyan);background:#12425b}}.layout{{display:grid;grid-template-columns:minmax(0,1fr) 310px;gap:16px}}.panel{{background:linear-gradient(145deg,#111b2c,#080c14);border:1px solid var(--line);border-radius:20px}}.diagram{{padding:14px;overflow:auto}}svg{{display:block;width:100%;min-width:720px;height:650px}}.edge{{fill:none;stroke-width:3;stroke-dasharray:8 6;opacity:.85;cursor:pointer}}.edge:hover,.edge.selected{{stroke-width:6;opacity:1}}.hbond{{stroke:var(--hbond)}}.ionic{{stroke:var(--ionic)}}.hydrophobic{{stroke:var(--hydrophobic)}}.other{{stroke:var(--other)}}.node circle{{fill:#d8e6f4;stroke:#07101e;stroke-width:2.5}}.node.left circle{{fill:#ffdb76}}.node text{{font-size:11px;font-weight:900;text-anchor:middle;pointer-events:none;fill:#07101e}}.side{{padding:18px}}.metric{{font-size:2.5rem;font-weight:900;color:var(--cyan);letter-spacing:-.07em}}.detail{{border-top:1px solid var(--line);margin-top:14px;padding-top:14px;color:var(--muted);line-height:1.5}}.legend{{display:grid;gap:8px;margin-top:18px;color:var(--muted);font-size:.85rem}}.swatch{{display:inline-block;width:27px;border-top:4px dashed currentColor;vertical-align:middle;margin-right:7px}}.sequence-panel{{margin-top:16px;padding:18px}}.sequence-title{{font-weight:900}}.sequence{{display:flex;flex-wrap:wrap;gap:4px;margin-top:10px;max-height:170px;overflow:auto;padding:3px}}.aa{{min-width:25px;text-align:center;padding:5px 3px;border-radius:6px;background:#17243a;color:#dbeafe;font:800 .76rem ui-monospace,monospace;cursor:default}}.aa.contact{{background:var(--pep);color:#17110a;box-shadow:0 0 0 1px rgba(255,207,91,.55)}}.aa.selected{{outline:3px solid var(--cyan)}}.empty{{padding:30px;color:var(--muted)}}@media(max-width:900px){{main{{padding:18px}}.layout{{grid-template-columns:1fr}}.side{{display:grid;grid-template-columns:1fr 1fr;gap:12px}}.side h2,.legend{{grid-column:1/-1}}}}@media(max-width:520px){{.side{{display:block}}}}
</style></head><body><main><div class="eyebrow">PYMACS · INTERACTIVE INTERFACE RIN</div><h1>Persistent protein interface</h1><p>Choose an interface and persistence threshold. Edges retain their chemically validated class when available; gray edges are observed heavy-atom contacts without a more specific validated class. Peptide chains are shown below as their complete connected sequence, with contacting residues highlighted.</p><div class="controls" id="pairs"></div><div class="controls" id="thresholds"><button data-t=".95">95% retained</button><button data-t=".90" class="active">90% retained</button><button data-t=".80">80% retained</button></div><section class="layout"><div class="panel diagram"><svg id="network" viewBox="0 0 1000 650" role="img" aria-label="Interactive residue interface network"></svg></div><aside class="panel side"><h2>Visible interactions</h2><div><div class="metric" id="count">—</div><div id="label">≥ 90% occupancy</div></div><div class="detail" id="detail">Click an edge or residue to inspect its interaction type and occupancy.</div><div class="legend"><span style="color:var(--hbond)"><i class="swatch"></i>Hydrogen bond</span><span style="color:var(--ionic)"><i class="swatch"></i>Ionic</span><span style="color:var(--hydrophobic)"><i class="swatch"></i>Hydrophobic</span><span style="color:var(--other)"><i class="swatch"></i>Other contact</span></div></aside></section><section class="panel sequence-panel" id="sequences"></section></main><script>const DATA={payload};let pair=DATA.pairs[0]?.id,threshold=.90,selected=null;const svg=document.querySelector('#network'),count=document.querySelector('#count'),label=document.querySelector('#label'),detail=document.querySelector('#detail');const esc=s=>String(s).replace(/[&<>"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c]));const edgeClass=t=>/h-bond/i.test(t)?'hbond':/ionic/i.test(t)?'ionic':/hydrophobic/i.test(t)?'hydrophobic':'other';function current(){{const p=DATA.pairs.find(x=>x.id===pair);return {{p,edges:DATA.edges.filter(e=>e.a.startsWith(p.left+':')&&e.b.startsWith(p.right+':')&&e.fraction+1e-9>=threshold)}}}}function renderSequence(p,edges){{const contacts=new Set(edges.flatMap(e=>[e.a,e.b]));const chains=[p.left,p.right].map(id=>DATA.chains[id]).filter(Boolean);document.querySelector('#sequences').innerHTML=chains.map(c=>`<div class="sequence-title">${{esc(c.label)}} <span style="color:var(--muted);font-weight:600">(${{esc(c.type)}} · full sequence)</span></div><div class="sequence">${{c.residues.map(r=>`<span class="aa ${{contacts.has(r.label)?'contact':''}} ${{selected===r.label?'selected':''}}" title="${{esc(r.label)}}">${{esc(r.one)}}<sup>${{r.number}}</sup></span>`).join('')}}</div>`).join('')}}function render(){{const {{p,edges}}=current();if(!p){{svg.innerHTML='';return}}const left=[...new Set(edges.map(e=>e.a))],right=[...new Set(edges.map(e=>e.b))];const positions={{}};[left,right].forEach((nodes,side)=>nodes.forEach((node,i)=>positions[node]=[side?770:230,110+i*Math.min(480/Math.max(1,nodes.length-1),58)]));const lines=edges.map((e,i)=>{{const [x1,y1]=positions[e.a],[x2,y2]=positions[e.b],dy=(i%2?1:-1)*(18+i%3*10);return `<path class="edge ${{edgeClass(e.type)}} ${{selected===i?'selected':''}}" data-id="${{i}}" d="M ${{x1}} ${{y1}} Q 500 ${{(y1+y2)/2+dy}} ${{x2}} ${{y2}}"/>`}}).join('');const nodes=[...left,...right].map((n,i)=>{{const [x,y]=positions[n];return `<g class="node ${{left.includes(n)?'left':''}}" data-node="${{esc(n)}}"><circle cx="${{x}}" cy="${{y}}" r="40"/><text x="${{x}}" y="${{y-5}}">${{esc(n.split(':')[1]||n)}}</text><text x="${{x}}" y="${{y+11}}">${{esc(n.split(':')[0])}}</text></g>`}}).join('');svg.innerHTML=`<text x="230" y="45" text-anchor="middle" fill="#ffdb76" font-weight="900">${{esc(p.left)}}</text><text x="770" y="45" text-anchor="middle" fill="#dbeafe" font-weight="900">${{esc(p.right)}}</text>${{lines+nodes}}`;count.textContent=edges.length;label.textContent=`≥ ${{Math.round(threshold*100)}}% occupancy`;svg.querySelectorAll('.edge').forEach(el=>el.onclick=()=>{{const e=edges[+el.dataset.id];selected=+el.dataset.id;detail.innerHTML=`<b>${{esc(e.a)}}</b> ↔ <b>${{esc(e.b)}}</b><br>${{esc(e.type)}} · <b>${{(e.fraction*100).toFixed(1)}}%</b> occupancy`;render()}});svg.querySelectorAll('.node').forEach(el=>el.onclick=()=>{{selected=el.dataset.node;detail.innerHTML=`<b>${{esc(selected)}}</b><br>Highlighted in the complete chain sequence below.`;renderSequence(p,edges)}});if(!edges.length)detail.textContent='No interactions meet this threshold for the selected interface.';renderSequence(p,edges)}}document.querySelector('#pairs').innerHTML=DATA.pairs.map((p,i)=>`<button data-p="${{esc(p.id)}}" class="${{i?'':'active'}}">${{esc(p.left)}} ↔ ${{esc(p.right)}}</button>`).join('');document.querySelectorAll('[data-p]').forEach(b=>b.onclick=()=>{{pair=b.dataset.p;selected=null;document.querySelectorAll('[data-p]').forEach(x=>x.classList.toggle('active',x===b));render()}});document.querySelectorAll('[data-t]').forEach(b=>b.onclick=()=>{{threshold=+b.dataset.t;selected=null;document.querySelectorAll('[data-t]').forEach(x=>x.classList.toggle('active',x===b));render()}});render();</script></body></html>'''
    html = html.replace(
        '<div class="controls" id="thresholds"><button data-t=".95">95% retained</button><button data-t=".90" class="active">90% retained</button><button data-t=".80">80% retained</button></div>',
        f'<div class="controls" id="thresholds">{threshold_buttons}</div>'
        '<div class="controls" id="colour-controls" aria-label="Residue coloring mode"></div>',
    )
    html = html.replace(
        "</style>",
        ".colour-note{font-size:.8rem;color:var(--muted);margin:0 0 12px}.node circle{transition:fill .16s ease,stroke .16s ease}.node[data-side=left] circle{fill:#ffcf5b}.node[data-side=right] circle{fill:#4bd7ff}</style>",
        1,
    )
    # Exactly one chain may use residue-chemistry colours at a time.  The
    # opposite chain remains a fixed, high-contrast colour, preventing the
    # two interfacial surfaces from becoming visually indistinguishable.
    interactive_script = r'''<script>
const DATA=__PAYLOAD__;
let pair=DATA.pairs[0]?.id,threshold=__DEFAULT__,selected=null,colourMode='chain';
const svg=document.querySelector('#network'),count=document.querySelector('#count'),label=document.querySelector('#label'),detail=document.querySelector('#detail'),colourControls=document.querySelector('#colour-controls');
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const edgeClass=t=>/h-bond/i.test(t)?'hbond':/ionic/i.test(t)?'ionic':/hydrophobic/i.test(t)?'hydrophobic':'other';
const chemistry={hydrophobic:'#7dff91',positive:'#c084fc',negative:'#ff8b82',polar:'#70eaff',other:'#d8e6f4'};
function residueChemistry(name){const match=String(name).match(/:([A-Z]{3})/i),aa=(match?.[1]||'').toUpperCase();if(['ALA','VAL','LEU','ILE','MET','PHE','TRP','TYR','PRO'].includes(aa))return'hydrophobic';if(['LYS','ARG','HIS'].includes(aa))return'positive';if(['ASP','GLU'].includes(aa))return'negative';if(['SER','THR','ASN','GLN','CYS'].includes(aa))return'polar';return'other'}
function current(){const p=DATA.pairs.find(x=>x.id===pair);return{p,edges:p?DATA.edges.filter(e=>e.a.startsWith(p.left+':')&&e.b.startsWith(p.right+':')&&e.fraction+1e-9>=threshold):[]}}
function nodeColour(name,side){const typeMode=colourMode===`${side}-type`;return typeMode?chemistry[residueChemistry(name)]:(side==='left'?'#ffcf5b':'#4bd7ff')}
function renderSequence(p,edges){const contacts=new Set(edges.flatMap(e=>[e.a,e.b]));const chains=[['left',p.left],['right',p.right]].map(([side,id])=>[side,DATA.chains[id]]).filter(([,c])=>c);document.querySelector('#sequences').innerHTML='<p class="colour-note">Chain colors keep the two surfaces separate. “Color by chemistry” applies to only one chain at a time; the other remains uniform.</p>'+chains.map(([side,c])=>`<div class="sequence-title" style="color:${nodeColour(c.residues[0]?.label||'',side)}">${esc(c.label)} <span style="color:var(--muted);font-weight:600">(${esc(c.type)} · full sequence)</span></div><div class="sequence">${c.residues.map(r=>`<span class="aa ${contacts.has(r.label)?'contact':''} ${selected===r.label?'selected':''}" style="background:${nodeColour(r.label,side)};color:#07101e" title="${esc(r.label)} · ${residueChemistry(r.label)}">${esc(r.one)}<sup>${r.number}</sup></span>`).join('')}</div>`).join('')}
function render(){const {p,edges}=current();if(!p){svg.innerHTML='';return}const left=[...new Set(edges.map(e=>e.a))],right=[...new Set(edges.map(e=>e.b))],maxNodes=Math.max(1,left.length,right.length),height=Math.max(650,160+maxNodes*76),proteinRadius=Math.max(24,Math.min(40,44-(left.length+right.length)*.35)),isPeptide=side=>DATA.chains[side==='left'?p.left:p.right]?.type==='peptide',radiusFor=side=>isPeptide(side)?40:proteinRadius,fontFor=side=>Math.max(8,Math.round(11*radiusFor(side)/40)),top=100,bottom=height-56;svg.setAttribute('viewBox',`0 0 1000 ${height}`);svg.style.height=`${height}px`;const positions={};[left,right].forEach((nodes,side)=>nodes.forEach((node,i)=>positions[node]=[side?770:230,top+i*((bottom-top)/Math.max(1,nodes.length-1))]));const lines=edges.map((e,i)=>{const [x1,y1]=positions[e.a],[x2,y2]=positions[e.b],dy=(i%2?1:-1)*(18+i%3*10);return `<path class="edge ${edgeClass(e.type)} ${selected===i?'selected':''}" data-id="${i}" d="M ${x1} ${y1} Q 500 ${(y1+y2)/2+dy} ${x2} ${y2}"/>`}).join('');const nodes=[...left,...right].map(n=>{const side=left.includes(n)?'left':'right',[x,y]=positions[n],radius=radiusFor(side),font=fontFor(side);return `<g class="node" data-side="${side}" data-node="${esc(n)}"><circle cx="${x}" cy="${y}" r="${radius}"/><text x="${x}" y="${y-5}" style="font-size:${font}px">${esc(n.split(':')[1]||n)}</text><text x="${x}" y="${y+11}" style="font-size:${font}px">${esc(n.split(':')[0])}</text></g>`}).join('');svg.innerHTML=`<text x="230" y="45" text-anchor="middle" fill="#ffcf5b" font-weight="900">${esc(p.left)}</text><text x="770" y="45" text-anchor="middle" fill="#4bd7ff" font-weight="900">${esc(p.right)}</text>${lines+nodes}`;svg.querySelectorAll('.node').forEach(el=>{const side=el.dataset.side;el.querySelector('circle').style.fill=nodeColour(el.dataset.node,side)});count.textContent=edges.length;label.textContent=`≥ ${Math.round(threshold*100)}% occupancy`;svg.querySelectorAll('.edge').forEach(el=>el.onclick=()=>{const e=edges[+el.dataset.id];selected=+el.dataset.id;detail.innerHTML=`<b>${esc(e.a)}</b> ↔ <b>${esc(e.b)}</b><br>${esc(e.type)} · <b>${(e.fraction*100).toFixed(1)}%</b> occupancy`;render()});svg.querySelectorAll('.node').forEach(el=>el.onclick=()=>{selected=el.dataset.node;detail.innerHTML=`<b>${esc(selected)}</b><br>Highlighted in the complete chain sequence below.`;renderSequence(p,edges)});if(!edges.length)detail.textContent='No interactions meet this threshold for the selected interface.';renderSequence(p,edges)}
function renderColourControls(){const p=DATA.pairs.find(x=>x.id===pair);if(!p)return;colourControls.innerHTML=`<button data-colour="chain" class="${colourMode==='chain'?'active':''}">Chain colors</button><button data-colour="left-type" class="${colourMode==='left-type'?'active':''}">Color ${esc(p.left)} by chemistry</button><button data-colour="right-type" class="${colourMode==='right-type'?'active':''}">Color ${esc(p.right)} by chemistry</button>`;colourControls.querySelectorAll('[data-colour]').forEach(b=>b.onclick=()=>{colourMode=b.dataset.colour;renderColourControls();render()})}
document.querySelector('#pairs').innerHTML=DATA.pairs.map((p,i)=>`<button data-p="${esc(p.id)}" class="${i?'':'active'}">${esc(p.left)} ↔ ${esc(p.right)}</button>`).join('');document.querySelectorAll('[data-p]').forEach(b=>b.onclick=()=>{pair=b.dataset.p;selected=null;colourMode='chain';document.querySelectorAll('[data-p]').forEach(x=>x.classList.toggle('active',x===b));renderColourControls();render()});document.querySelectorAll('[data-t]').forEach(b=>b.onclick=()=>{threshold=+b.dataset.t;selected=null;document.querySelectorAll('[data-t]').forEach(x=>x.classList.toggle('active',x===b));render()});renderColourControls();render();
</script>'''.replace("__PAYLOAD__", payload).replace("__DEFAULT__", f"{default_threshold:.2f}")
    html = re.sub(r"<script>const DATA=.*?</script>", interactive_script, html, count=1, flags=re.DOTALL)
    (outdir / "Interface_RIN_Explorer.html").write_text(html, encoding="utf-8")


def write_interface_manifest(
    manifest_path: Path,
    args,
    input_files: Dict[str, str],
    chain_pairs,
    polymer_entries,
    total_frames: int,
    frame_indices: Sequence[int],
    numerical_outputs: Sequence[Path],
    plot_outputs: Sequence[Path],
    data_parameters: Dict[str, Any],
    plot_parameters: Dict[str, Any],
    *,
    legacy_adopted: bool = False,
) -> Dict[str, Any]:
    chain_metadata = []
    selected_labels = {
        left["display_label"]: left for left, _ in chain_pairs
    }
    selected_labels.update({right["display_label"]: right for _, right in chain_pairs})

    for entry in polymer_entries:
        if entry["display_label"] not in selected_labels:
            continue
        chain_metadata.append(
            {
                "display_label": entry["display_label"],
                "chain_type": entry.get("chain_type"),
                "label": entry.get("label"),
                "current_chain": entry.get("current_chain"),
                "source_chain": entry.get("source_chain"),
                "atom_range_1based": [int(entry["start"]) + 1, int(entry["end"]) + 1],
            }
        )

    manifest = record_stage_completion(
        manifest_path,
        stage_name="INTERFACE_RIN",
        analysis_script=Path(__file__).name,
        analysis_source_path=__file__,
        input_files=input_files,
        data_parameters=data_parameters,
        plot_parameters=plot_parameters,
        outputs=numerical_outputs,
        plot_outputs=plot_outputs,
        system="interface_rin",
        mode="protein_interface",
        legacy_adopted=legacy_adopted,
        extra_stage_fields={
            "total_trajectory_frames": int(total_frames),
            "analyzed_frame_count": int(len(frame_indices)),
            "selected_chain_pairs": [
                f"{left['display_label']}:{right['display_label']}" for left, right in chain_pairs
            ],
            "chain_metadata": chain_metadata,
        },
    )
    manifest["topology_fingerprint"] = manifest["input_fingerprints"].get("topology")
    manifest["trajectory_fingerprint"] = manifest["input_fingerprints"].get("trajectory")
    manifest["atomindex_fingerprint"] = manifest["input_fingerprints"].get("atomindex")
    manifest["total_trajectory_frames"] = int(total_frames)
    manifest["analyzed_frame_count"] = int(len(frame_indices))
    manifest["frame_step"] = int(args.frame_step)
    manifest["contact_cutoff_A"] = float(args.contact_cutoff)
    manifest["selected_chain_pairs"] = [
        f"{left['display_label']}:{right['display_label']}" for left, right in chain_pairs
    ]
    manifest["chain_metadata"] = chain_metadata
    manifest["numerical_output_files"] = describe_output_files(numerical_outputs, base_dir=manifest_path.parent)
    manifest["plot_output_files"] = describe_output_files(plot_outputs, base_dir=manifest_path.parent)
    write_json_atomic(manifest_path, manifest)
    return manifest


def main():
    start_time = time.perf_counter()
    args = build_parser().parse_args()
    if mda is None or capped_distance is None or distance_array is None:
        raise SystemExit(
            "MDAnalysis is required to run interface analysis. "
            "Install the MDAnalysis environment before executing 3C_Interface_RIN.py."
        )
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    progress = ProgressReporter(
        outdir, pipeline="interface_rin",
        stages=("load_inputs", "contact_scan", "write_tables", "render_figures", "write_manifest"),
    )
    progress.start("load_inputs", "Load trajectory and interface selections")
    manifest_path = outdir / "interface_analysis_manifest.json"
    numerical_outputs = interface_numerical_output_paths(outdir)
    plot_outputs = interface_plot_output_paths(outdir)

    stage_times = {}

    t0 = time.perf_counter()
    universe = load_universe(Path(args.topo), Path(args.traj))
    stage_times["load_universe_s"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    atomindex_entries = parse_atomindex_entries(args.atomindex)
    polymer_entries = detect_polymer_entries(apply_peptide_chain_hints(atomindex_entries, args.peptide_chains))
    stage_times["parse_atomindex_s"] = time.perf_counter() - t0
    progress.complete("load_inputs")

    print(f"🧬 Chains parsed from atomIndex.txt: {len(polymer_entries)}")
    for entry in polymer_entries:
        print(
            f"   • {entry['display_label']} ({entry['chain_type']}) "
            f"atoms {entry['start'] + 1}-{entry['end'] + 1}"
        )

    if len(polymer_entries) < 2:
        print("ℹ️ Fewer than two polymer chains were found in atomIndex.txt. Nothing to analyze.")
        return 0

    chain_pairs = parse_requested_chain_pairs(args.chain_pairs, polymer_entries)
    if not chain_pairs:
        print("ℹ️ No valid chain pairs were selected for analysis. Nothing to do.")
        return 0

    print(
        "🧬 Chain pairs analyzed: "
        + ", ".join(f"{left['display_label']}:{right['display_label']}" for left, right in chain_pairs)
    )
    print(f"🧬 Cutoff: {args.contact_cutoff:.2f} Å")
    print(f"🧬 Frame step: {args.frame_step}")

    frame_indices = iter_frame_indices(len(universe.trajectory), args.frame_step)
    print(f"🧬 Number of frames analyzed: {len(frame_indices)} / {len(universe.trajectory)} total")
    print(f"🧬 Output folder: {outdir}")

    if args.dry_run:
        print("🧪 Dry run requested; exiting before trajectory analysis.")
        return 0

    input_files = {
        "topology": str(Path(args.topo).resolve()),
        "trajectory": str(Path(args.traj).resolve()),
        "atomindex": str(Path(args.atomindex).resolve()),
    }
    data_parameters = {
        "contact_cutoff_A": float(args.contact_cutoff),
        "frame_step": int(args.frame_step),
        "chunk_size": int(args.chunk_size),
        "selected_chain_pairs": [
            f"{left['display_label']}:{right['display_label']}" for left, right in chain_pairs
        ],
        "total_trajectory_frames": int(len(universe.trajectory)),
        "analyzed_frame_count": int(len(frame_indices)),
        "chain_metadata": [
            {
                "display_label": entry["display_label"],
                "chain_type": entry.get("chain_type"),
                "atom_range_1based": [int(entry["start"]) + 1, int(entry["end"]) + 1],
            }
            for entry in polymer_entries
        ],
    }
    plot_parameters = interface_plot_parameters(args)
    data_signature = current_data_signature(input_files, data_parameters)

    cache_status = validate_stage_cache(
        manifest_path=manifest_path,
        stage_name="INTERFACE_RIN",
        required_outputs=numerical_outputs,
        current_data_signature=data_signature,
    )
    legacy_outputs_present = all(
        (outdir / name).exists()
        for name in [
            "interface_chain_pair_summary.csv",
            "interface_contact_timeseries.csv",
            "interface_residue_pair_contacts.csv",
            "interface_node_summary.csv",
        ]
    )

    if args.force_recompute:
        print("🔄 --force-recompute requested; ignoring any existing cache.")
    elif cache_status["reusable"]:
        cached = load_cached_interface_outputs(outdir)
        stage_record = cache_status.get("stage_record") or {}
        cached_plot_parameters = ((stage_record.get("parameters") or {}).get("plot_only") or {})
        plots_ready = all(path.exists() and path.stat().st_size > 0 for path in plot_outputs)

        if cache_status["legacy_adopt"]:
            print("♻️ Existing interface dataset found without a registered manifest.")
            print("📝 Adopting legacy interface outputs into interface_analysis_manifest.json.")
            write_interface_manifest(
                manifest_path,
                args,
                input_files,
                chain_pairs,
                polymer_entries,
                len(universe.trajectory),
                frame_indices,
                numerical_outputs,
                plot_outputs,
                data_parameters,
                plot_parameters,
                legacy_adopted=True,
            )
            cached_plot_parameters = plot_parameters

        if cached_plot_parameters == plot_parameters and plots_ready:
            print("✅ Existing interface dataset found")
            print("✅ Input trajectory/topology compatible")
            print("♻️ Reusing cached interface data")
            print("✅ Interface plots already present; no trajectory rescan needed.")
            progress.complete("contact_scan")
            progress.complete("write_tables")
            progress.complete("render_figures")
            progress.complete("write_manifest")
            progress.finish()
            return 0

        print("✅ Existing interface dataset found")
        print("✅ Input trajectory/topology compatible")
        print("♻️ Reusing cached interface data")
        if cached_plot_parameters != plot_parameters:
            print("🎨 Plot-only parameters changed; regenerating interface figures from cached data.")
        else:
            print("🎨 One or more interface figures are missing; regenerating figures from cached data.")

        t0 = time.perf_counter()
        progress.start("render_figures", "Regenerate interface figures from cached contact data")
        plot_chain_pair_heatmap(cached["summary_df"], outdir)
        plot_contact_timeseries(cached["timeseries_df"], outdir)
        plot_min_distance_timeseries(cached["timeseries_df"], outdir, args.contact_cutoff)
        plot_residue_pair_heatmap(cached["residue_df"], cached["presence_payload"], outdir)
        plot_residue_network(cached["residue_df"], outdir, args.min_contact_frac, args.max_edges)
        plot_chain_network(cached["summary_df"], outdir)
        plot_top_residue_pairs_barplot(cached["residue_df"], outdir)
        events_path = outdir / "ValidatedInterfaceInteractionEvents.csv"
        cached_events = pd.read_csv(events_path) if events_path.exists() else pd.DataFrame(columns=EVENT_COLUMNS)
        write_interactive_interface_rin_html(
            outdir, cached["residue_df"], cached_events, len(frame_indices), chain_pairs, universe
        )
        stage_times["plotting_s"] = time.perf_counter() - t0
        progress.complete("render_figures")

        runtime_lines = [
            f"input_topology={Path(args.topo).resolve()}",
            f"input_trajectory={Path(args.traj).resolve()}",
            f"atomindex_file={Path(args.atomindex).resolve()}",
            f"frames_total={len(universe.trajectory)}",
            f"frames_analyzed={len(frame_indices)}",
            f"frame_step={args.frame_step}",
            f"contact_cutoff_A={args.contact_cutoff}",
            f"min_contact_frac={args.min_contact_frac}",
            "chain_pairs_analyzed=" + ",".join(
                f"{left['display_label']}:{right['display_label']}" for left, right in chain_pairs
            ),
            f"plotting_s={stage_times['plotting_s']:.3f}",
            f"total_runtime_s={time.perf_counter() - start_time:.3f}",
        ]
        (outdir / "interface_runtime_summary.txt").write_text("\n".join(runtime_lines) + "\n")
        write_interface_manifest(
            manifest_path,
            args,
            input_files,
            chain_pairs,
            polymer_entries,
            len(universe.trajectory),
            frame_indices,
            numerical_outputs,
            plot_outputs,
            data_parameters,
            plot_parameters,
            legacy_adopted=bool(cache_status["legacy_adopt"]),
        )
        progress.complete("contact_scan")
        progress.complete("write_tables")
        progress.complete("write_manifest")
        progress.finish()
        return 0
    else:
        if legacy_outputs_present:
            print("⚠️ Legacy interface cache is incomplete for the current schema.")
            print("🔄 Backfilling interface temporal dataset.")
        elif cache_status["reason"] == "data_signature_mismatch":
            print("⚠️ Existing interface cache is incompatible with current data-affecting inputs.")
            print("🔄 Recomputing interface contact dataset.")
        else:
            print("🔄 No reusable interface cache found; running trajectory contact scan.")

    workdir = outdir / ".work"
    signature = _work_signature(data_signature, args)
    state_path = workdir / "progress.json"
    state = read_json_manifest(state_path) if state_path.exists() else None
    if args.force_recompute or not state or state.get("signature") != signature:
        if workdir.exists():
            shutil.rmtree(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        state = {
            "version": 1, "signature": signature, "chunk_size": int(args.chunk_size),
            "input_fingerprints": fingerprint_files(input_files),
            "contact_cutoff_A": float(args.contact_cutoff), "frame_step": int(args.frame_step),
            "selected_chain_pairs": data_parameters["selected_chain_pairs"],
            "peptide_chains": str(args.peptide_chains or ""), "completed_chunks": [],
        }
        write_json_atomic(state_path, state)
    t0 = time.perf_counter()
    progress.start("contact_scan", "Scan trajectory frames into bounded resumable shards")
    scan_interface_shards(universe, chain_pairs, frame_indices, args.contact_cutoff, workdir, state)
    stage_times["analyze_contacts_s"] = time.perf_counter() - t0
    progress.complete("contact_scan")
    analyzed_times_ns = [float(universe.trajectory[index].time) / 1000.0 for index in frame_indices]
    t0 = time.perf_counter()
    progress.start("write_tables", "Reduce interface shards into numerical outputs")
    summary_df, timeseries_df, residue_df, node_df, frame_index_df, presence_payload, validated_events_df = reduce_interface_shards(outdir, workdir, frame_indices, analyzed_times_ns, chain_pairs, args.contact_cutoff)
    stage_times["write_csv_s"] = time.perf_counter() - t0
    progress.complete("write_tables")

    runtime_lines = [
        f"input_topology={Path(args.topo).resolve()}",
        f"input_trajectory={Path(args.traj).resolve()}",
        f"atomindex_file={Path(args.atomindex).resolve()}",
        f"frames_total={len(universe.trajectory)}",
        f"frames_analyzed={len(frame_indices)}",
        f"frame_step={args.frame_step}",
        f"contact_cutoff_A={args.contact_cutoff}",
        f"min_contact_frac={args.min_contact_frac}",
        "chain_pairs_analyzed=" + ",".join(
            f"{left['display_label']}:{right['display_label']}" for left, right in chain_pairs
        ),
    ]
    runtime_lines.extend(f"{key}={value:.3f}" for key, value in stage_times.items())
    runtime_lines.append(f"total_runtime_s={time.perf_counter() - start_time:.3f}")

    (outdir / "interface_runtime_summary.txt").write_text("\n".join(runtime_lines) + "\n")

    t0 = time.perf_counter()
    progress.start("render_figures", "Render interface figures and explorer")
    plot_chain_pair_heatmap(summary_df, outdir)
    plot_contact_timeseries(timeseries_df, outdir)
    plot_min_distance_timeseries(timeseries_df, outdir, args.contact_cutoff)
    plot_residue_pair_heatmap(residue_df, presence_payload, outdir)
    plot_residue_network(residue_df, outdir, args.min_contact_frac, args.max_edges)
    plot_chain_network(summary_df, outdir)
    plot_top_residue_pairs_barplot(residue_df, outdir)
    write_interactive_interface_rin_html(
        outdir, residue_df, validated_events_df, len(frame_indices), chain_pairs, universe
    )
    stage_times["plotting_s"] = time.perf_counter() - t0
    progress.complete("render_figures")

    runtime_lines.extend(
        [
            f"write_csv_s={stage_times['write_csv_s']:.3f}",
            f"plotting_s={stage_times['plotting_s']:.3f}",
        ]
    )
    (outdir / "interface_runtime_summary.txt").write_text("\n".join(runtime_lines) + "\n")
    progress.start("write_manifest", "Write interface completion manifest")
    write_interface_manifest(
        manifest_path,
        args,
        input_files,
        chain_pairs,
        polymer_entries,
        len(universe.trajectory),
        frame_indices,
        numerical_outputs,
        plot_outputs,
        data_parameters,
        plot_parameters,
    )
    progress.complete("write_manifest")
    progress.finish()

    outputs = [
        "interface_chain_pair_summary.csv",
        "interface_contact_timeseries.csv",
        "interface_residue_pair_contacts.csv",
        "interface_node_summary.csv",
        "interface_frame_index.csv",
        "interface_residue_pair_framewise.csv.gz",
        "interface_residue_pair_presence.npz",
        "ValidatedInterfaceInteractionEvents.csv",
        "ValidatedInterfaceInteractionOccupancy.csv",
        "interface_analysis_manifest.json",
        "interface_runtime_summary.txt",
        "Interface_ChainPair_ContactFractions.png",
        "Interface_ContactCount_Timeseries.png",
        "Interface_MinDistance_Timeseries.png",
        "Interface_ResiduePair_ContactHeatmap.png",
        "Interface_RIN_Network.png",
        "Interface_ChainLevel_Network.png",
        "Interface_TopResiduePairs_Barplot.png",
        "Interface_RIN_Explorer.html",
    ]
    print("✅ Main outputs:")
    for output in outputs:
        print(f"   • {outdir / output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
