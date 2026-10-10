#!/usr/bin/env bash
# Install an already-built PyMACS Leela stack for all Leela users.
# Run this as root (for example: sudo bash install_leela_shared.sh).
# Private Kyle Funnel URLs and tokens are deliberately excluded.
set -euo pipefail

if [ "${EUID}" -ne 0 ]; then
    echo "ERROR: run this administrator installer with sudo." >&2
    exit 1
fi

SOURCE_USER="${SUDO_USER:-${USER}}"
SOURCE_HOME="$(getent passwd "$SOURCE_USER" | awk -F: '{print $6}')"
SOURCE_ROOT="${PYMACS_LEELA_SOURCE_ROOT:-$SOURCE_HOME}"
DESTINATION_ROOT="${PYMACS_LEELA_ROOT:-/data/pymacs}"
GROUP_NAME="${PYMACS_LEELA_GROUP:-cheminfo}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

require_directory() {
    if [ ! -d "$1" ]; then
        echo "ERROR: missing required source directory: $1" >&2
        exit 1
    fi
}

require_directory "$SOURCE_ROOT/micromamba/envs/pymacs-leela-setup"
require_directory "$SOURCE_ROOT/micromamba/envs/pymacs-leela-analysis"
require_directory "$SOURCE_ROOT/gromacs-2021.5-install"
if [ ! -x "$SOURCE_ROOT/bin/micromamba" ]; then
    echo "ERROR: missing micromamba at $SOURCE_ROOT/bin/micromamba" >&2
    exit 1
fi

install -d -o root -g "$GROUP_NAME" -m 0755 \
    "$DESTINATION_ROOT" "$DESTINATION_ROOT/bin" "$DESTINATION_ROOT/envs" "$DESTINATION_ROOT/share"
cp -a "$SOURCE_ROOT/bin/micromamba" "$DESTINATION_ROOT/bin/micromamba"
cp -a "$SOURCE_ROOT/gromacs-2021.5-install" "$DESTINATION_ROOT/gromacs-2021.5"
cp -a "$SOURCE_ROOT/micromamba/envs/pymacs-leela-setup" "$DESTINATION_ROOT/envs/"
cp -a "$SOURCE_ROOT/micromamba/envs/pymacs-leela-analysis" "$DESTINATION_ROOT/envs/"
install -o root -g "$GROUP_NAME" -m 0644 \
    "$SCRIPT_DIR/hpc_profiles/leela/pymacs-leela-env.sh" "$DESTINATION_ROOT/share/pymacs-leela-env.sh"
install -o root -g root -m 0644 \
    "$SCRIPT_DIR/hpc_profiles/leela/pymacs-leela-profile.sh" /etc/profile.d/pymacs-leela.sh
chown -R root:"$GROUP_NAME" "$DESTINATION_ROOT"
chmod -R go-w "$DESTINATION_ROOT"
find "$DESTINATION_ROOT" -type d -exec chmod 0755 {} +
chmod 0755 "$DESTINATION_ROOT/bin/micromamba" "$DESTINATION_ROOT/gromacs-2021.5/bin/gmx"

echo "Shared Leela PyMACS installation is ready at $DESTINATION_ROOT."
echo "No Funnel token or private Funnel URL was installed."
