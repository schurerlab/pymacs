#!/usr/bin/env bash
# Run from this repository once Randy is reachable through the VPN.
set -euo pipefail
REMOTE_HOST="${1:-randy}"
REMOTE_DIR="/opt/lab-job-portal"
scp -r job_portal "${REMOTE_HOST}:/tmp/lab-job-portal"
ssh "$REMOTE_HOST" "sudo install -d -o jxs794 -g staff ${REMOTE_DIR} /var/lib/lab-job-portal/jobs && sudo rsync -a --delete /tmp/lab-job-portal/ ${REMOTE_DIR}/ && sudo chown -R jxs794:staff ${REMOTE_DIR} /var/lib/lab-job-portal && python3 -m venv ${REMOTE_DIR}/.venv && ${REMOTE_DIR}/.venv/bin/pip install -r ${REMOTE_DIR}/requirements.txt && sudo install -m 0644 ${REMOTE_DIR}/lab-job-portal.service /etc/systemd/system/lab-job-portal.service && sudo systemctl daemon-reload && sudo systemctl enable --now lab-job-portal"
echo "Create /etc/lab-job-portal.env from lab-job-portal.env.example, then run:"
echo "  sudo systemctl restart lab-job-portal"
echo "Expose only after local verification, e.g.: tailscale serve --https=443 http://127.0.0.1:8788"
