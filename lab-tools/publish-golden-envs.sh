#!/usr/bin/env bash
# Lab-maintainer utility: clone the validated source environments into the
# Cheminformatics shared release. Run it on each relevant cluster only.
set -euo pipefail

ROOT="${PYMACS_LAB_TOOLS_ROOT:-/projectnb/cheminfo/share/pymacs-lab/lab-tools}"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
clone_environment() {
    local source_prefix="$1" destination="$2"
    [[ -d "$source_prefix/conda-meta" ]] || die "Source environment is missing: $source_prefix"
    [[ ! -e "$destination" ]] || die "Destination already exists: $destination (publish a new release instead of overwriting it)"
    "$CONDA_EXE" create -y --prefix "$destination" --clone "$source_prefix"
}

case "$(uname -m)" in
  ppc64le)
    CONDA_EXE="${CONDA_EXE:-$HOME/miniforge3/bin/conda}"
    [[ -x "$CONDA_EXE" ]] || die "Set CONDA_EXE to the Triton Conda executable."
    install -d -m 2770 "$ROOT/triton/envs"
    clone_environment "$HOME/miniforge3/envs/cgenff" "$ROOT/triton/envs/cgenff"
    clone_environment "$HOME/miniforge3/envs/mdanalysis" "$ROOT/triton/envs/mdanalysis"
    clone_environment "$HOME/miniforge3/envs/vina_env" "$ROOT/triton/envs/vina_env"
    chmod -R g-w "$ROOT/triton/envs"
    ;;
  x86_64)
    CONDA_EXE="${CONDA_EXE:-$HOME/miniconda3/bin/conda}"
    [[ -x "$CONDA_EXE" ]] || die "Set CONDA_EXE to the Pegasus Conda executable."
    install -d -m 2770 "$ROOT/pegasus/envs"
    clone_environment "$HOME/miniconda3/envs/vina_env" "$ROOT/pegasus/envs/vina_env"
    chmod -R g-w "$ROOT/pegasus/envs"
    ;;
  *) die "Unsupported platform: $(uname -m)" ;;
esac
