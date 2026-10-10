#!/usr/bin/env bash
# Source this file from ~/.bashrc on Leela. It does not contain any secret.
# A site administrator may set PYMACS_LEELA_ROOT (for example /data/pymacs)
# to use the shared PyMACS installation instead of a per-user installation.

PYMACS_LEELA_ROOT="${PYMACS_LEELA_ROOT:-}"

_pymacs_leela_mamba() {
    if [ -n "$PYMACS_LEELA_ROOT" ] && [ -x "$PYMACS_LEELA_ROOT/bin/micromamba" ]; then
        printf '%s\n' "$PYMACS_LEELA_ROOT/bin/micromamba"
        return 0
    fi
    if command -v micromamba >/dev/null 2>&1; then
        command -v micromamba
    elif command -v mamba >/dev/null 2>&1; then
        command -v mamba
    else
        echo "PyMACS Leela setup: micromamba or mamba is not available." >&2
        return 1
    fi
}

_pymacs_leela_activate() {
    local environment_name="$1"
    local mamba_bin
    local mamba_command
    mamba_bin="$(_pymacs_leela_mamba)" || return 1
    eval "$("$mamba_bin" shell hook --shell bash)"
    # The shell hook defines a function named after the executable.  Calling
    # the binary path here would launch a subprocess and cannot change PATH.
    mamba_command="$(basename "$mamba_bin")"
    "$mamba_command" activate "$environment_name"
}

_pymacs_leela_gromacs() {
    export CUDA_HOME="/usr/local/cuda-11.5"
    local gromacs_root="${PYMACS_LEELA_GROMACS_ROOT:-$HOME/gromacs-2021.5-install}"
    if [ -x "$gromacs_root/bin/gmx" ]; then
        export PATH="$gromacs_root/bin:$CUDA_HOME/bin:$PATH"
        export LD_LIBRARY_PATH="$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
    else
        echo "PyMACS Leela setup: expected CUDA-11.5 GROMACS was not found at $gromacs_root/bin/gmx." >&2
        return 1
    fi
}

pymacs-leela() {
    local setup_environment="pymacs-leela-setup"
    if [ -n "$PYMACS_LEELA_ROOT" ]; then
        setup_environment="$PYMACS_LEELA_ROOT/envs/pymacs-leela-setup"
    fi
    _pymacs_leela_activate "$setup_environment" || return 1
    _pymacs_leela_gromacs || return 1
    if [ -r "$HOME/.config/pymacs/cgenff-funnel.env" ]; then
        . "$HOME/.config/pymacs/cgenff-funnel.env"
    fi
    export PYMACS_CGENFF_TOKEN_FILE="$HOME/.config/pymacs/cgenff-funnel.token"
}

pymacs-leela-analysis() {
    local analysis_environment="pymacs-leela-analysis"
    if [ -n "$PYMACS_LEELA_ROOT" ]; then
        analysis_environment="$PYMACS_LEELA_ROOT/envs/pymacs-leela-analysis"
    fi
    _pymacs_leela_activate "$analysis_environment" || return 1
}
