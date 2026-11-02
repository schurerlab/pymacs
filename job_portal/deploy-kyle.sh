#!/usr/bin/env bash
# Deploy the portal code without overwriting Kyle's credentials, runtime, or
# persistent job/user records.  The remote .env is intentionally managed only
# on Kyle.
set -euo pipefail

REMOTE_HOST="${1:-100.81.66.47}"
REMOTE_DIR="/home/jxs794/lab-job-portal"

rsync -az --delete \
  --exclude '.env' \
  --exclude '.venv/' \
  --exclude 'vendor/' \
  --exclude 'jobs/' \
  --exclude 'portal.log' \
  --exclude 'portal.pid' \
  ./ "$REMOTE_HOST:$REMOTE_DIR/"

ssh "$REMOTE_HOST" '
  set -e
  cd /home/jxs794/lab-job-portal
  test -f .env
  set -a; . ./.env; set +a
  if test -s portal.pid; then
    old_pid=$(cat portal.pid)
    if kill -0 "$old_pid" 2>/dev/null; then kill "$old_pid"; fi
  fi
  # Cover a stale PID file without matching the deploy shell itself.
  pkill -f "[p]ython.*portal.py" || true
  runtime=python3
  if test -x .venv/bin/python && .venv/bin/python -c "import flask" >/dev/null 2>&1; then
    runtime=.venv/bin/python
  fi
  nohup env PYTHONPATH="$HOME/lab-job-portal/vendor" "$runtime" portal.py > portal.log 2>&1 < /dev/null &
  echo $! > portal.pid
  sleep 1
  curl -fsS "http://127.0.0.1:${PORT:-8788}/api/submitters" >/dev/null
'

echo "Portal deployed and responding on Kyle."
