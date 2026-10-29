#!/usr/bin/env python3
"""Generate isolated, analysis-only LSF jobs for final ligand dashboard backfill."""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "final_ligand_scorecard_backfill.json"
OUT = ROOT / "final_ligand_scorecard_jobs"
STAGE_ROOT = "/projectnb/triton/home/jxs794/Pymacs/backfill_staging/final_ligand_scorecard_20261029"


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")[-72:]


def render(entry: dict[str, str]) -> str:
    name = slug(entry["cache"])
    stage = f"{STAGE_ROOT}/{name}"
    return f'''#!/usr/bin/env bash
#BSUB -J scorecard_{name[:45]}
#BSUB -L /bin/bash
#BSUB -q normal
#BSUB -P brd
#BSUB -W 02:00
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -oo {stage}/backfill.%J.out
#BSUB -eo {stage}/backfill.%J.err

# Analysis-only staging: no MD is started and the source directory is never written.
set -eo pipefail
SOURCE={entry["source"]!r}
STAGE={stage!r}
rm -rf "$STAGE"
mkdir -p "$STAGE"
test -d "$SOURCE"
test -d "$SOURCE/Analysis_Results"
cp -a "$SOURCE/Analysis_Results" "$STAGE/"
find "$SOURCE" -maxdepth 1 -type f \
  \( -name '*.cgenff.mol2' -o -name '*.mol2' -o -name '*.pdb' -o -name '*.gro' -o -name 'atomIndex*.txt' -o -name '.chainmap.json' \) \
  -exec cp -a {{}} "$STAGE/" \;

source /etc/bashrc
set -u
module purge
module load gcc-runtime/13.4.0-none-none-gh6rpjx
module load gromacs/2025.1-gcc-13.4.0-6vq7xfo
source /projectnb/triton/home/jxs794/miniforge3/etc/profile.d/conda.sh
conda activate mdanalysis
export OMP_NUM_THREADS=2
cd "$STAGE"

ligand="$(find . -maxdepth 1 -type f -name '*.cgenff.mol2' -printf '%f\\n' | head -1 | sed 's/\\.cgenff\\.mol2$//')"
test -n "$ligand"
if test -s Analysis_Results/InteractionTypes_Framewise.csv; then
  python /projectnb/triton/home/jxs794/Pymacs/3B_NETWORX.py --ligand "$ligand" --filter 0.10 --output "$ligand" --mode ligand
else
  echo "No InteractionTypes_Framewise.csv: preserving available analyses; explorer cannot be regenerated." >&2
fi
python /projectnb/triton/home/jxs794/Pymacs/5_PyMACS_Result_Dashboard.py --root .
grep -q 'Binding scorecard' Analysis_Results/PyMACS_Result_Dashboard.html
'''


def main() -> None:
    entries = json.loads(MANIFEST.read_text())
    OUT.mkdir(exist_ok=True)
    for old in OUT.glob("*.lsf"):
        old.unlink()
    for entry in entries:
        path = OUT / f"{slug(entry['cache'])}.lsf"
        path.write_text(render(entry))
        path.chmod(0o755)
    print(f"wrote {len(entries)} jobs to {OUT}")


if __name__ == "__main__":
    main()
