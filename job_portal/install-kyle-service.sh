#!/usr/bin/env bash
# Run interactively on Kyle once. It prompts for sudo only for OS packages and
# installing the service definition; the application runs as jxs794.
set -euo pipefail
APP_ROOT="/home/jxs794/lab-job-portal"

sudo apt-get update
sudo apt-get install -y python3-venv
python3 -m venv "$APP_ROOT/.venv"
"$APP_ROOT/.venv/bin/pip" install --upgrade pip
"$APP_ROOT/.venv/bin/pip" install -r "$APP_ROOT/requirements.txt"
sudo install -m 0644 "$APP_ROOT/lab-job-portal.service" /etc/systemd/system/lab-job-portal.service
sudo systemctl daemon-reload
sudo systemctl enable --now lab-job-portal
sudo systemctl status --no-pager lab-job-portal
echo "Portal: https://kyle.rove-vernier.ts.net:8788/"
