# Results dossier backfill playbook

This playbook refreshes an existing PyMACS results dossier after MD has
completed. It does **not** rerun molecular dynamics or modify trajectories.

## What the shared dashboard script can do

`5_PyMACS_Result_Dashboard.py` is a standalone post-processing publisher. Run
it from a completed simulation directory (or pass `--root`) after the normal
analysis stages have produced `Analysis_Results/`.

It reads already-published CSVs, final structures, `atomIndex.txt`, peptide
identity sidecars, and optional `Analysis_Results/Interface_RIN/` outputs. It
then writes or refreshes:

```text
Analysis_Results/PyMACS_Result_Dashboard.html
Analysis_Results/portal_results_manifest.json
Analysis_Results/PyMACS_Heatmap_Explorer.html   # when applicable
```

The dashboard refresh is inexpensive. It does **not** create missing science.
For a peptide interface, residue-node maps and interface persistence graphs
require the separate `3C_Interface_RIN.py` analysis to have access to the
completed `Final_Trajectory.pdb`, `Final_Trajectory.xtc`, and `atomIndex.txt`.

## Atom-index naming contract

The component label is the first field in `atomIndex.txt`; `current_chain` is
the chain identifier in the simulated coordinates. The report uses both.

```text
2wo3           1-2816 # current_chain=A chain_type=protein source_rescount=176
xEFN4_TOPdock 2817-3070 # current_chain=B chain_type=protein source_rescount=15
```

For a peptide run, the RIN launch passes chain B explicitly as the peptide
partner. The resulting dossier labels this pair as:

```text
2wo3 (receptor, chain A)
xEFN4_TOPdock (peptide, chain B)
```

This avoids both anonymous `Protein` labels and the incorrect treatment of a
standard-amino-acid peptide as a second generic protein chain.

## Kyle cache audit — 2026-10-28

Audit source: `/home/jxs794/lab-job-portal/jobs/md_results_index.json` on
Kyle. A dashboard is considered current when it contains the current shared
`Interface maps`/`interface-preview-grid` dossier marker.

| State | Count | Meaning |
|---|---:|---|
| Current dossier | 0 | No published cache dossier matched the current marker at audit time. |
| Legacy dossier | 49 | A report exists, but was generated before the current shared dossier. |
| Missing dossier | 56 | No `Analysis_Results/PyMACS_Result_Dashboard.html` exists. |

The inventory is a publication-cache audit, not a claim that every source run
still has its trajectory. Every run must pass the preflight below before an
LSF job is submitted.

## Generate a cross-reference inventory before a bulk backfill

Run the following on Kyle. It creates a tab-separated inventory that pairs
each cache result with the discovered Triton source directory by system name.
The script performs no writes on Triton and does not submit jobs. A blank
`triton_source` means that source data was not found under the searched roots;
`AMBIGUOUS` means more than one source directory had the same basename and
must be resolved manually.

```bash
python3 - <<'PY'
import csv, json, subprocess
from collections import defaultdict
from pathlib import Path

index = Path("/home/jxs794/lab-job-portal/jobs/md_results_index.json")
output = Path("/home/jxs794/dashboard_backfill_inventory.tsv")
cache_root = Path("/home/jxs794/lab-job-portal/jobs/published-md-results")
triton_roots = [
    "/projectnb/triton/home/jxs794/Molecular_Dynamics",
    "/projectnb/triton/home/jxs794/dashboard-refresh",
]

# One remote listing avoids launching one SSH connection per cache result.
find_cmd = "find " + " ".join(repr(p) for p in triton_roots) + \
    " -mindepth 1 -maxdepth 4 -type d -name Analysis_Results -printf '%h\\n' 2>/dev/null"
listed = subprocess.check_output(["ssh", "triton", find_cmd], text=True)
by_name = defaultdict(list)
for path in listed.splitlines():
    by_name[Path(path).name].append(path)

def cache_state(result):
    dashboard = cache_root / result["path"] / "Analysis_Results" / "PyMACS_Result_Dashboard.html"
    if not dashboard.is_file():
        return "missing"
    text = dashboard.read_text(errors="ignore")
    return "current" if "interface-preview-grid" in text and "Interface maps" in text else "legacy"

records = json.loads(index.read_text())
with output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["cache_path", "system", "dashboard_state", "triton_source", "match_state"])
    for result in sorted(records, key=lambda row: row["path"]):
        system = Path(result["path"]).name
        matches = by_name[system]
        source = matches[0] if len(matches) == 1 else ("AMBIGUOUS" if matches else "")
        match_state = "matched" if len(matches) == 1 else ("ambiguous" if matches else "not_found")
        writer.writerow([result["path"], system, cache_state(result), source, match_state])
print(output)
PY

column -ts $'\t' ~/dashboard_backfill_inventory.tsv | less -S
```

## Ephrin peptide priority batch

All 24 Ephrin systems map directly between the same system name on Triton and
Kyle:

```text
Triton source: /projectnb/triton/home/jxs794/Molecular_Dynamics/<system>
Kyle cache:    /home/jxs794/lab-job-portal/jobs/published-md-results/triton/Ephrin-analysis-recovery-20261027/<system>
```

| Systems | Cache status at audit | Backfill action |
|---|---|---|
| `1kgy-xEFN[1-4]_{2dock,3dock,TOPdock}` | 12 missing dossiers | Preflight, then dashboard refresh; run Interface RIN only when source trajectory files are present. |
| `2wo3-xEFN[1-4]_{2dock,3dock}` | 11 missing dossiers | Same. |
| `2wo3-xEFN4_TOPdock` | Legacy dossier | Regenerate immediately; it has peptide RMSD/RMSF but its published cache lacked Interface RIN outputs at audit time. |

## Preflight: no mutation

Run this before creating a job for any system:

```bash
RUN=/projectnb/triton/home/jxs794/Molecular_Dynamics/2wo3-xEFN4_TOPdock
test -f "$RUN/Final_Trajectory.pdb" && echo "final PDB: yes" || echo "final PDB: no"
test -f "$RUN/Final_Trajectory.xtc" && echo "final XTC: yes" || echo "final XTC: no"
test -f "$RUN/atomIndex.txt" && echo "atom index: yes" || echo "atom index: no"
test -d "$RUN/Analysis_Results" && echo "analysis directory: yes" || echo "analysis directory: no"
sed -n '1,12p' "$RUN/atomIndex.txt"
```

Eligibility rules:

1. **Dashboard-only refresh:** requires `Analysis_Results/`; optionally uses
   the final PDB for the 3D viewer.
2. **Interface RIN refresh:** also requires `Final_Trajectory.pdb`,
   `Final_Trajectory.xtc`, and at least two entries in `atomIndex.txt`.
3. Never run either job against a directory with a running MD/analysis job.

## LSF job tiers

Use two tiers. This keeps inexpensive dashboard publication separate from
trajectory scanning.

### Tier 1 — dashboard-only refresh

This is appropriate for every eligible completed result. It does not change
the trajectory or replace CSV science outputs.

```bash
#!/usr/bin/env bash
#BSUB -J pymacs_dashboard_<system>
#BSUB -q gpu
#BSUB -n 2
#BSUB -W 00:30
#BSUB -oo dashboard_refresh.%J.out
#BSUB -eo dashboard_refresh.%J.err
set -euo pipefail

RUN=/projectnb/triton/home/jxs794/Molecular_Dynamics/<system>
PYMACS=/projectnb/triton/home/jxs794/Pymacs
cd "$RUN"
test -d Analysis_Results
source /projectnb/cheminfo/share/pymacs-lab/lab-tools/triton/envs/mdanalysis/bin/activate
python "$PYMACS/5_PyMACS_Result_Dashboard.py" --root "$RUN"
```

### Tier 2 — peptide interface RIN plus dashboard refresh

Submit this only after Tier 1 preflight confirms the final PDB/XTC and atom
index exist. The current shared `3A` runner passes the peptide-chain role into
Interface RIN automatically for new analyses. For a historical backfill,
invoke the RIN directly with the explicit peptide chain from `atomIndex.txt`.

```bash
#!/usr/bin/env bash
#BSUB -J pymacs_interface_<system>
#BSUB -q gpu
#BSUB -n 8
#BSUB -W 12:00
#BSUB -oo interface_refresh.%J.out
#BSUB -eo interface_refresh.%J.err
set -euo pipefail

RUN=/projectnb/triton/home/jxs794/Molecular_Dynamics/<system>
PYMACS=/projectnb/triton/home/jxs794/Pymacs
PEPTIDE_CHAIN=B                    # read this from atomIndex.txt; do not guess
cd "$RUN"
test -f Final_Trajectory.pdb
test -f Final_Trajectory.xtc
test -f atomIndex.txt
source /projectnb/cheminfo/share/pymacs-lab/lab-tools/triton/envs/mdanalysis/bin/activate
python "$PYMACS/3C_Interface_RIN.py" \
  --topo Final_Trajectory.pdb --traj Final_Trajectory.xtc \
  --atomindex atomIndex.txt --outdir Analysis_Results/Interface_RIN \
  --peptide-chains "$PEPTIDE_CHAIN" --contact-cutoff 4.0 \
  --min-contact-frac 0.10 --frame-step 1 --headless
python "$PYMACS/5_PyMACS_Result_Dashboard.py" --root "$RUN"
```

## Cache publication checklist

After a job is complete, publish only derived outputs to the matching Kyle
cache directory. Do not copy or overwrite trajectories, topology inputs, or
historical raw output folders.

1. Confirm the output dossier contains the expected peptide/receptor labels.
2. Confirm `Interface_RIN/Interface_RIN_Explorer.html` and its CSVs exist for
   Tier 2 jobs.
3. Copy the refreshed HTML, manifest, generated ligand/heatmap assets, and
   `Interface_RIN/` directory where present.
4. Atomically update the Kyle result index with the new asset list.
5. Visually inspect at least one 1kgy and one 2wo3 dossier before bulk
publication.

## Broader stale-cache groups

The same audit identified stale/missing dossier candidates beyond Ephrin:

- Joey/Angela ranked-ligand recovery jobs (`joey/*` and `angela/*`).
- Historical Angela ranked docks.
- Historical HDGF/LEDGF results.
- Historical MacroH2A peptide models.
- Historical Nahuel apo/ligand systems.

These should be handled only after the Ephrin pilot validates the copied cache
workflow. Their source directory roots vary, so do not infer a Triton path
from a Kyle cache path without the preflight check.
