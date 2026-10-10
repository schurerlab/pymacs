#!/usr/bin/env bash
# Create PyMACS's two Leela environments with micromamba or mamba.
# GROMACS is intentionally supplied by Leela's CUDA-11.5 installation, not Conda.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SETUP_ENV="pymacs-leela-setup"
ANALYSIS_ENV="pymacs-leela-analysis"

if [ -x "$HOME/bin/micromamba" ]; then
    # Leela's default installation is configured by ~/.bashrc, but this
    # installer must also work from a non-interactive shell.
    MAMBA_BIN="$HOME/bin/micromamba"
elif command -v micromamba >/dev/null 2>&1; then
    MAMBA_BIN="$(command -v micromamba)"
elif command -v mamba >/dev/null 2>&1; then
    MAMBA_BIN="$(command -v mamba)"
else
    echo "ERROR: micromamba or mamba must be available before running this installer." >&2
    exit 1
fi

env_exists() {
    "$MAMBA_BIN" run -n "$1" python -c "import sys; assert sys.version_info >= (3, 9)" >/dev/null 2>&1
}

create_if_missing() {
    local environment_name="$1"
    local environment_file="$2"
    if env_exists "$environment_name"; then
        echo "Environment '$environment_name' already exists; leaving it unchanged."
    else
        echo "Creating '$environment_name' from $(basename "$environment_file") ..."
        "$MAMBA_BIN" env create -f "$environment_file"
    fi
}

create_if_missing "$SETUP_ENV" "$SCRIPT_DIR/environment_leela_setup.yml"
create_if_missing "$ANALYSIS_ENV" "$SCRIPT_DIR/environment_leela_analysis.yml"

CONFIG_DIR="$HOME/.config/pymacs"
HELPER_SOURCE="$SCRIPT_DIR/hpc_profiles/leela/pymacs-leela-env.sh"
HELPER_DESTINATION="$CONFIG_DIR/leela-env.sh"
mkdir -p "$CONFIG_DIR"
chmod 700 "$CONFIG_DIR"
install -m 600 "$HELPER_SOURCE" "$HELPER_DESTINATION"

if ! grep -q "PYMACS LEELA ENVIRONMENT" "$HOME/.bashrc" 2>/dev/null; then
    {
        printf '\n# BEGIN PYMACS LEELA ENVIRONMENT\n'
        printf '[ -r "$HOME/.config/pymacs/leela-env.sh" ] && . "$HOME/.config/pymacs/leela-env.sh"\n'
        printf '# END PYMACS LEELA ENVIRONMENT\n'
    } >> "$HOME/.bashrc"
fi

cat <<'EOF'

Leela environments are ready.

Open a new shell, or run:
  source "$HOME/.bashrc"
  pymacs-leela

For CPU trajectory analysis:
  pymacs-leela-analysis

This installer deliberately does not install GROMACS or private Funnel credentials.
EOF
