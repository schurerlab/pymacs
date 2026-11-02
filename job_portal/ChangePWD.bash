#!/usr/bin/env bash
# Change one approved PyMACS portal user's password on Kyle.
set -euo pipefail

APP_ROOT="/home/jxs794/lab-job-portal"
cd "$APP_ROOT"

name="${1:-}"
if [[ -z "$name" ]]; then
  read -r -p "Portal account first name: " name
fi
[[ "$name" =~ ^[A-Za-z][A-Za-z-]{0,39}$ ]] || {
  printf 'Use a first name containing letters or hyphens only.\n' >&2
  exit 2
}

set -a
. ./.env
set +a

exec env PYTHONPATH="$HOME/lab-job-portal/vendor" python3 manage_accounts.py set-password "$name"
