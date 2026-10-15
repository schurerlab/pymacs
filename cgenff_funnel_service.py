#!/usr/bin/env python3
"""Kyle-side, token-gated CGenFF service for PyMACS.

Run behind Tailscale Funnel. The application deliberately listens only on
127.0.0.1, accepts a single MOL2 parameterization operation, and never accepts
shell commands, file paths, or executable names from a caller.
"""

from __future__ import annotations

import argparse
import base64
import hmac
import os
import re
import secrets
import subprocess
import tempfile
import threading
import time
import json
from pathlib import Path
from urllib.request import Request, urlopen

from flask import Flask, jsonify, request


MAX_MOL2_BYTES = 2 * 1024 * 1024
LIGAND_CODE_RE = re.compile(r"^[A-Z0-9]{1,4}$")


def _analytics_enabled() -> bool:
    return os.environ.get("PYMACS_ANALYTICS_ENABLED", "").lower() in {"1", "true", "yes", "on"}


def _duration_bucket(seconds: float) -> str:
    return "under_10s" if seconds < 10 else "10_to_30s" if seconds < 30 else "30_to_120s" if seconds < 120 else "over_120s"


def _emit_analytics(event_type: str, feature: str = "", status: str = "") -> None:
    """Send only a coarse operational counter; never include request or ligand data."""
    url = os.environ.get("PYMACS_ANALYTICS_RANDY_URL", "").strip().rstrip("/")
    token = os.environ.get("PYMACS_ANALYTICS_RANDY_TOKEN", "").strip()
    if not _analytics_enabled() or not url or not token:
        return
    payload = {"event_id": secrets.token_urlsafe(18), "event_type": event_type,
               "visitor_id": secrets.token_urlsafe(18), "session_id": secrets.token_urlsafe(18),
               "path": "/v1/parameterize", "referrer": "service", "device": "desktop",
               "browser": "other", "feature": feature, "status": status}
    def send() -> None:
        try:
            req = Request(url + "/event", data=json.dumps(payload).encode(), headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json", "User-Agent": "pymacs-cgenff-analytics/1.0"}, method="POST")
            urlopen(req, timeout=2.0).close()
        except Exception:
            pass
    threading.Thread(target=send, daemon=True).start()


def read_secret(path: Path) -> str:
    secret = path.expanduser().read_text(encoding="utf-8").strip()
    if len(secret) < 24:
        raise ValueError("The Funnel token must be at least 24 characters.")
    return secret


def create_app(*, token_file: Path, cgenff_exec: Path, timeout_seconds: int) -> Flask:
    token = read_secret(token_file)
    if not cgenff_exec.is_file() or not os.access(cgenff_exec, os.X_OK):
        raise ValueError(f"CGenFF executable is unavailable: {cgenff_exec}")
    app = Flask(__name__)
    app.config["MAX_CONTENT_LENGTH"] = 3 * 1024 * 1024

    @app.post("/v1/parameterize")
    def parameterize():
        started = time.monotonic()
        supplied = request.headers.get("Authorization", "").removeprefix("Bearer ").strip()
        if not supplied or not hmac.compare_digest(supplied, token):
            _emit_analytics("cgenff_request_failed", "unauthorized", "failure")
            return jsonify(error="unauthorized"), 401
        _emit_analytics("cgenff_request_started")
        data = request.get_json(silent=True)
        if not isinstance(data, dict):
            _emit_analytics("cgenff_request_failed", "invalid_input", "failure")
            return jsonify(error="expected a JSON request"), 400
        ligand_code = str(data.get("ligand_code", "")).upper().strip()
        if not LIGAND_CODE_RE.fullmatch(ligand_code):
            _emit_analytics("cgenff_request_failed", "invalid_input", "failure")
            return jsonify(error="ligand_code must contain 1-4 uppercase letters or digits"), 400
        try:
            mol2_bytes = base64.b64decode(str(data.get("mol2_b64", "")), validate=True)
        except ValueError:
            _emit_analytics("cgenff_request_failed", "invalid_input", "failure")
            return jsonify(error="mol2_b64 is not valid base64"), 400
        if not mol2_bytes or len(mol2_bytes) > MAX_MOL2_BYTES:
            _emit_analytics("cgenff_request_failed", "invalid_input", "failure")
            return jsonify(error="MOL2 exceeds the allowed size"), 413
        if b"@<TRIPOS>MOLECULE" not in mol2_bytes.upper() or b"\x00" in mol2_bytes:
            _emit_analytics("cgenff_request_failed", "invalid_input", "failure")
            return jsonify(error="input is not a valid text MOL2"), 400

        with tempfile.TemporaryDirectory(prefix="pymacs-cgenff-") as temp_dir:
            workdir = Path(temp_dir)
            mol2_path = workdir / f"{ligand_code}.cgenff.mol2"
            mol2_path.write_bytes(mol2_bytes)
            try:
                completed = subprocess.run(
                    [str(cgenff_exec), mol2_path.name],
                    cwd=workdir,
                    stdin=subprocess.DEVNULL,
                    capture_output=True,
                    timeout=timeout_seconds,
                    check=False,
                )
            except subprocess.TimeoutExpired:
                _emit_analytics("cgenff_request_failed", "timeout", "failure")
                return jsonify(error="CGenFF timed out"), 504
            stdout, stderr = completed.stdout, completed.stderr
        if completed.returncode != 0 or b"skipped molecule" in stderr.lower() or b"unfulfilled valence" in stderr.lower():
            _emit_analytics("cgenff_request_failed", "execution", "failure")
            return jsonify(
                error="CGenFF could not parameterize this molecule",
                log_b64=base64.b64encode(stderr[-65536:]).decode("ascii"),
            ), 422
        if not stdout.strip():
            _emit_analytics("cgenff_request_failed", "execution", "failure")
            return jsonify(error="CGenFF returned an empty stream file"), 502
        _emit_analytics("cgenff_request_completed", _duration_bucket(time.monotonic() - started), "success")
        log = b"CGenFF completed successfully.\n" + stderr[-65536:]
        return jsonify(
            request_id=secrets.token_hex(8),
            str_b64=base64.b64encode(stdout).decode("ascii"),
            log_b64=base64.b64encode(log).decode("ascii"),
        )

    return app


def main() -> None:
    parser = argparse.ArgumentParser(description="Token-gated local CGenFF API for Tailscale Funnel.")
    parser.add_argument("--token-file", required=True, type=Path)
    parser.add_argument("--cgenff-exec", required=True, type=Path)
    parser.add_argument("--port", type=int, default=8787)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()
    app = create_app(token_file=args.token_file, cgenff_exec=args.cgenff_exec, timeout_seconds=args.timeout)
    try:
        from waitress import serve
    except ImportError as exc:
        raise SystemExit("Install requirements_cgenff_funnel_service.txt before starting the service.") from exc
    serve(app, host="127.0.0.1", port=args.port, threads=4)


if __name__ == "__main__":
    main()
