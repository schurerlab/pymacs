#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import math
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cbook as mpl_cbook
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, distance_array

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    tqdm = None


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
    parser.add_argument("--max-edges", type=int, default=100)
    parser.add_argument("--chain-pairs", default=None)
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
    return (
        entry.get("current_chain")
        or entry.get("source_chain")
        or entry.get("label")
        or f"Chain_{index}"
    )


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


def classify_residue_pair(resname_a, resname_b):
    a = (resname_a or "").upper()
    b = (resname_b or "").upper()
    if (a in POSITIVE and b in NEGATIVE) or (a in NEGATIVE and b in POSITIVE):
        return "ionic_or_salt_bridge_like"
    if a in HYDROPHOBIC and b in HYDROPHOBIC:
        return "hydrophobic"
    if a in AROMATIC and b in AROMATIC:
        return "aromatic"
    if a in POLAR and b in POLAR:
        return "polar"
    return "mixed_or_other"


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
        residue_rows.append(
            {
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
                "InteractionClass": classify_residue_pair(res_left["resname"], res_right["resname"]),
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


def write_csv_outputs(outdir, summary_df, timeseries_df, residue_df, node_df, runtime_lines):
    summary_df.to_csv(outdir / "interface_chain_pair_summary.csv", index=False)
    timeseries_df.to_csv(outdir / "interface_contact_timeseries.csv", index=False)
    residue_export_df = residue_df.drop(columns=["_pair_key", "_pair_id", "_series_key"], errors="ignore")
    residue_export_df.to_csv(outdir / "interface_residue_pair_contacts.csv", index=False)
    node_df.to_csv(outdir / "interface_node_summary.csv", index=False)
    (outdir / "interface_runtime_summary.txt").write_text("\n".join(runtime_lines) + "\n")


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


def plot_residue_pair_heatmap(residue_df, pair_binary_series, analyzed_times_ns, outdir):
    path = outdir / "Interface_ResiduePair_ContactHeatmap.png"
    if residue_df.empty:
        plot_placeholder(path, "Top Interface Residue Pairs Over Time", "No residue-residue contacts were observed.")
        return

    top_df = residue_df.sort_values(["ContactFraction", "ContactFrames"], ascending=[False, False]).head(80)
    rows = []
    row_labels = []
    for _, row in top_df.iterrows():
        series = pair_binary_series.get(row["_series_key"])
        if series is None:
            continue
        row_labels.append(f"{row['ResidueA']} -- {row['ResidueB']}")
        rows.append(series)

    if not rows:
        plot_placeholder(path, "Top Interface Residue Pairs Over Time", "No residue-pair time series were available.")
        return

    heatmap_df = pd.DataFrame(rows, index=row_labels, columns=[f"{t:.3f}" for t in analyzed_times_ns])
    plt.figure(figsize=(max(10, len(analyzed_times_ns) * 0.18), max(6, len(row_labels) * 0.25)))
    sns.heatmap(heatmap_df, cmap="rocket_r", cbar_kws={"label": "Contact present (0/1)"}, vmin=0, vmax=1)
    plt.xlabel("Time (ns)")
    plt.ylabel("Residue pair")
    plt.title("Top Interface Residue-Pair Contacts Over Time")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def plot_residue_network(residue_df, outdir, min_contact_frac, max_edges):
    path = outdir / "Interface_RIN_Network.png"
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


def main():
    start_time = time.perf_counter()
    args = build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    stage_times = {}

    t0 = time.perf_counter()
    universe = load_universe(Path(args.topo), Path(args.traj))
    stage_times["load_universe_s"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    atomindex_entries = parse_atomindex_entries(args.atomindex)
    polymer_entries = detect_polymer_entries(atomindex_entries)
    stage_times["parse_atomindex_s"] = time.perf_counter() - t0

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

    t0 = time.perf_counter()
    summary_rows = []
    timeseries_rows = []
    residue_rows = []
    pair_binary_series = {}
    analyzed_times_ns = None

    for left_entry, right_entry in chain_pairs:
        result = analyze_chain_pair_contacts(universe, left_entry, right_entry, frame_indices, args.contact_cutoff)
        if result is None:
            print(f"⚠️ Skipping pair {left_entry['display_label']}:{right_entry['display_label']} due to empty heavy-atom selections.")
            continue
        summary_rows.append(result["summary_row"])
        timeseries_rows.extend(result["time_rows"])
        residue_rows.extend(result["residue_rows"])
        for pair_key, binary_series in result["pair_binary_series"].items():
            pair_binary_series[f"{result['summary_row']['ChainA_Label']}__{result['summary_row']['ChainB_Label']}::{pair_key[0]}::{pair_key[1]}"] = binary_series
        if analyzed_times_ns is None:
            analyzed_times_ns = result["times_ns"]
    stage_times["analyze_contacts_s"] = time.perf_counter() - t0

    summary_columns = [
        "ChainA_Label",
        "ChainB_Label",
        "ChainA_Type",
        "ChainB_Type",
        "Frames_Analyzed",
        "Frames_In_Contact",
        "ContactFraction",
        "MeanContactsPerFrame",
        "MaxContactsInFrame",
        "MeanMinDistance_A",
        "MinObservedDistance_A",
        "ContactCutoff_A",
    ]
    summary_df = pd.DataFrame(summary_rows, columns=summary_columns)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(["ContactFraction", "MeanContactsPerFrame"], ascending=[False, False])

    timeseries_columns = [
        "Time_ns",
        "Frame",
        "ChainA_Label",
        "ChainB_Label",
        "ContactCount",
        "MinDistance_A",
        "InContact",
    ]
    timeseries_df = pd.DataFrame(timeseries_rows, columns=timeseries_columns)
    if not timeseries_df.empty:
        timeseries_df = timeseries_df.sort_values(["ChainA_Label", "ChainB_Label", "Frame"])

    residue_columns = [
        "ChainA_Label",
        "ChainB_Label",
        "ResidueA",
        "ResidueB",
        "ResidueA_Name",
        "ResidueB_Name",
        "ResidueA_Index",
        "ResidueB_Index",
        "ContactFrames",
        "Frames_Analyzed",
        "ContactFraction",
        "MeanMinDistance_A",
        "MinObservedDistance_A",
        "InteractionClass",
        "_pair_key",
        "_pair_id",
        "_series_key",
    ]
    residue_df = pd.DataFrame(residue_rows, columns=residue_columns)
    if not residue_df.empty:
        residue_df = residue_df.sort_values(["ContactFraction", "ContactFrames"], ascending=[False, False])

    node_df = build_node_summary(residue_df)

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

    t0 = time.perf_counter()
    write_csv_outputs(outdir, summary_df, timeseries_df, residue_df, node_df, runtime_lines)
    stage_times["write_csv_s"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    plot_chain_pair_heatmap(summary_df, outdir)
    plot_contact_timeseries(timeseries_df, outdir)
    plot_min_distance_timeseries(timeseries_df, outdir, args.contact_cutoff)
    plot_residue_pair_heatmap(residue_df, pair_binary_series, analyzed_times_ns or [], outdir)
    plot_residue_network(residue_df, outdir, args.min_contact_frac, args.max_edges)
    plot_chain_network(summary_df, outdir)
    plot_top_residue_pairs_barplot(residue_df, outdir)
    stage_times["plotting_s"] = time.perf_counter() - t0

    runtime_lines.extend(
        [
            f"write_csv_s={stage_times['write_csv_s']:.3f}",
            f"plotting_s={stage_times['plotting_s']:.3f}",
        ]
    )
    (outdir / "interface_runtime_summary.txt").write_text("\n".join(runtime_lines) + "\n")

    outputs = [
        "interface_chain_pair_summary.csv",
        "interface_contact_timeseries.csv",
        "interface_residue_pair_contacts.csv",
        "interface_node_summary.csv",
        "interface_runtime_summary.txt",
        "Interface_ChainPair_ContactFractions.png",
        "Interface_ContactCount_Timeseries.png",
        "Interface_MinDistance_Timeseries.png",
        "Interface_ResiduePair_ContactHeatmap.png",
        "Interface_RIN_Network.png",
        "Interface_ChainLevel_Network.png",
        "Interface_TopResiduePairs_Barplot.png",
    ]
    print("✅ Main outputs:")
    for output in outputs:
        print(f"   • {outdir / output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
