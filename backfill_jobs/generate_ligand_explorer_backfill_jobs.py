#!/usr/bin/env python3
"""Generate isolated, analysis-only 3B + dashboard LSF jobs for ligands.

The supplied manifest is validated before this generator is used.  Each job
works in a new staging directory, leaves its completed source run untouched,
and verifies that the interactive artifact was actually created.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


PYMACS = "/projectnb/triton/home/jxs794/Pymacs"
DEFAULT_STAGE_ROOT = f"{PYMACS}/backfill_staging/ligand_explorer_20261029"


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")[-72:]


def render(entry: dict[str, str], stage_root: str) -> str:
    name = slug(entry["cache"])
    stage = f"{stage_root}/{name}"
    return f'''#!/usr/bin/env bash
#BSUB -J networx_{name[:42]}
#BSUB -L /bin/bash
#BSUB -q normal
#BSUB -P brd
#BSUB -W 02:00
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -oo {stage_root}/logs/{name}.%J.out
#BSUB -eo {stage_root}/logs/{name}.%J.err

# Analysis-only recovery: no MD is started and SOURCE is never modified.
set -eo pipefail
SOURCE={entry["source"]!r}
STAGE={stage!r}
test -d "$SOURCE"
test -s "$SOURCE/Analysis_Results/InteractionTypes_Framewise.csv"
test ! -e "$STAGE"
mkdir -p "$STAGE"
cp -a "$SOURCE/Analysis_Results" "$STAGE/"
find "$SOURCE" -maxdepth 1 -type f \\
  \\( -name '*.cgenff.mol2' -o -name '*.mol2' -o -name '*.pdb' -o -name '*.gro' -o -name 'atomIndex*.txt' -o -name '.chainmap.json' \\) \\
  -exec cp -a {{}} "$STAGE/" \\
  \\;

source /etc/bashrc
module purge
module load gcc-runtime/13.4.0-none-none-gh6rpjx
module load gromacs/2025.1-gcc-13.4.0-6vq7xfo
source /projectnb/triton/home/jxs794/miniforge3/etc/profile.d/conda.sh
conda activate mdanalysis
export OMP_NUM_THREADS=2
cd "$STAGE"

ligand="$(find . -maxdepth 1 -type f -name '*.cgenff.mol2' -printf '%f\\n' | head -1 | sed 's/\\.cgenff\\.mol2$//')"
test -n "$ligand"
python {PYMACS}/3B_NETWORX.py --ligand "$ligand" --filter 0.10 --output "$ligand" --mode ligand
test -s "Analysis_Results/NETWORX/${{ligand}}_interactive.html"
python {PYMACS}/5_PyMACS_Result_Dashboard.py --root .
grep -q 'Interactive ligand interaction explorer' Analysis_Results/PyMACS_Result_Dashboard.html
'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--stage-root", default=DEFAULT_STAGE_ROOT)
    args = parser.parse_args()
    entries = json.loads(args.manifest.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for entry in entries:
        path = args.output_dir / f"{slug(entry['cache'])}.lsf"
        path.write_text(render(entry, args.stage_root), encoding="utf-8")
        path.chmod(0o755)
    print(f"wrote {len(entries)} jobs to {args.output_dir}")


if __name__ == "__main__":
    main()
