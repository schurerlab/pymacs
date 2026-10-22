#!/usr/bin/env python3
"""Non-destructive diagnostic audit of UMF607-Z historical Level-3 output.

This is intentionally an *endpoint-screen audit*, not a replacement for the
geometry rerun.  A trajectory rerun requires MDAnalysis in the selected
runtime; this script makes no chemical claims that the old CSV cannot support.
"""
from pathlib import Path
import shutil
import pandas as pd
import matplotlib.pyplot as plt

SOURCE = Path("/Users/jxs794/Library/CloudStorage/Box-Box/PROTAC_Degrader_DB_Project_2025/Kif20a/Manuscript Pieces/KIF20A_Final_Analysis/MolecularDynamics/UMF607-Z")
OUT = Path(__file__).resolve().parent / "umf607_z"
OUT.mkdir(parents=True, exist_ok=True)

old = SOURCE / "Analysis_Results" / "InteractionTypes_Framewise.csv"
if not old.exists():
    old = SOURCE / "Analysis_Results_before_UMF607-Z_20261009" / "InteractionTypes_Framewise.csv"
df = pd.read_csv(old)
atom = df.get("LigAtomName", pd.Series("", index=df.index)).fillna("").astype(str).str.upper()
typ = df.get("Type", pd.Series("", index=df.index)).fillna("").astype(str)
invalid_endpoint = atom.str.startswith(("H", "LP", "VS", "MWS"))
screened = df.loc[~invalid_endpoint].copy()

rows = []
for label, table in (("original", df), ("endpoint-screened", screened)):
    for interaction in ("Ionic", "Hydrophobic", "H-bond"):
        historical_name = "H-bond" if interaction == "H-bond" else interaction
        source_name = interaction
        subset = table[typ.loc[table.index].eq(source_name)]
        rows.append({"Dataset": label, "InteractionClass": historical_name, "Events": int(len(subset)),
                     "Residues": int(subset.get("Residue", pd.Series(dtype=str)).nunique())})
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "historical_endpoint_screen_summary.csv", index=False)

artifacts = pd.DataFrame({
    "Metric": ["Historical rows", "Hydrogen-anchored ionic rows", "Hydrogen-anchored hydrophobic rows",
               "Virtual-site ionic/hydrophobic rows", "Duplicate residue-atom-type-frame rows"],
    "Count": [
        len(df),
        int(((typ == "Ionic") & atom.str.startswith("H")).sum()),
        int(((typ == "Hydrophobic") & atom.str.startswith("H")).sum()),
        int(((typ.isin(["Ionic", "Hydrophobic"])) & atom.str.startswith(("LP", "VS", "MWS"))).sum()),
        int(df.duplicated([c for c in ["Time_ns", "Residue", "LigAtomName", "Type"] if c in df]).sum()),
    ],
})
artifacts.to_csv(OUT / "historical_artifact_counts.csv", index=False)

fig, ax = plt.subplots(figsize=(8, 4.5))
pivot = summary.pivot(index="InteractionClass", columns="Dataset", values="Events").reindex(["H-bond", "Hydrophobic", "Ionic"])
pivot.plot.bar(ax=ax, color=["#9AA0A6", "#31A354"])
ax.set_ylabel("Historical CSV rows")
ax.set_title("UMF607-Z: endpoint screening of historical labels")
ax.set_xlabel("")
ax.legend(title="Data state")
fig.tight_layout()
fig.savefig(OUT / "umf607_z_endpoint_screen_before_after.png", dpi=220)
plt.close(fig)

for candidate in (SOURCE / "Analysis_Results" / "NETWORX" / "UNK_expanded.png",
                  SOURCE / "Analysis_Results_before_UMF607-Z_20261009" / "NETWORX" / "UNK_expanded.png"):
    if candidate.exists():
        shutil.copy2(candidate, OUT / "historical_sid.png")
        break

(OUT / "README.md").write_text(
    "# UMF607-Z historical-output audit\n\n"
    "This directory is non-destructive. `historical_sid.png` is the historical diagram. "
    "The before/after plot removes only rows with hydrogen or virtual-site ligand endpoints; "
    "it is not a geometry-validated rerun and must not be used as a validated SID. "
    "The shared engine requires an MDAnalysis-capable runtime to evaluate the supplied XTC/PDB.\n"
)
