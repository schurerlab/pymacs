# System-wide Bash profile fragment for a shared Leela PyMACS installation.
# Installed by a Leela administrator at /etc/profile.d/pymacs-leela.sh.
# No credentials belong in this file.
export PYMACS_LEELA_ROOT="/data/pymacs"
export PYMACS_LEELA_GROMACS_ROOT="$PYMACS_LEELA_ROOT/gromacs-2021.5"

if [ -r "$PYMACS_LEELA_ROOT/share/pymacs-leela-env.sh" ]; then
    . "$PYMACS_LEELA_ROOT/share/pymacs-leela-env.sh"
fi
