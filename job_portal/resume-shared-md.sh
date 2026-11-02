#!/usr/bin/env bash
# Submit a GROMACS continuation from a PyMACS campaign stored in Cheminfo's
# shared lab area. Run this script while logged in as the intended submitter.
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /projectnb/cheminfo/share/pymacs-lab/campaigns/.../systems/<system>" >&2
  exit 2
fi

system_dir="$1"
[[ "$system_dir" == /projectnb/cheminfo/share/pymacs-lab/* ]] || {
  echo "Refusing: system directory must be inside /projectnb/cheminfo/share/pymacs-lab" >&2
  exit 2
}
[[ -s "$system_dir/md_0_1.tpr" && -s "$system_dir/md_0_1.cpt" ]] || {
  echo "Refusing: md_0_1.tpr and a nonempty md_0_1.cpt are required." >&2
  exit 2
}

job_name="$(basename "$system_dir")-continue"
job_name="${job_name,,}"

bsub -J "$job_name" -q normal -P brd -W 144:00 -n 16 \
  -gpu 'num=1' -R 'span[hosts=1]' \
  -oo "$system_dir/${job_name}.%J.out" \
  -eo "$system_dir/${job_name}.%J.err" \
  bash -lc "
    source /etc/bashrc
    module purge
    module load gcc-runtime/13.4.0-none-none-gh6rpjx gromacs/2025.1-gcc-13.4.0-6vq7xfo
    cd '$system_dir'
    export OMP_NUM_THREADS=16
    gmx_mpi mdrun -v -deffnm md_0_1 -cpi md_0_1.cpt -append -pin on -pinoffset 0 -ntomp 16 -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update cpu
  "
