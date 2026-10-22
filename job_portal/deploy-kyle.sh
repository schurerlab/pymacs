#!/usr/bin/env bash
# Deploy the portal code without overwriting Kyle's credentials, runtime, or
# persistent job/user records.  The remote .env is intentionally managed only
# on Kyle.
set -euo pipefail

REMOTE_HOST="${1:-100.81.66.47}"
REMOTE_DIR="/home/jxs794/lab-job-portal"

rsync -az --delete \
  --exclude '.env' \
  --exclude 'vendor/' \
  --exclude 'jobs/' \
  --exclude 'portal.log' \
  --exclude 'portal.pid' \
  ./ "$REMOTE_HOST:$REMOTE_DIR/"

ssh "$REMOTE_HOST" '
  set -e
  cd /home/jxs794/lab-job-portal
  test -f .env
  pkill -f "^python3 portal.py$" || true
  set -a; . ./.env; set +a
  nohup env PYTHONPATH="$HOME/lab-job-portal/vendor" python3 portal.py > portal.log 2>&1 &
  echo $! > portal.pid
  sleep 1
  curl -fsS http://127.0.0.1:8788/api/submitters >/dev/null
'

echo "Portal deployed and responding on Kyle."
