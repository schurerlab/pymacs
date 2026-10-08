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
from pathlib import Path

from flask import Flask, jsonify, request


MAX_MOL2_BYTES = 2 * 1024 * 1024
LIGAND_CODE_RE = re.compile(r"^[A-Z0-9]{1,4}$")


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
        supplied = request.headers.get("Authorization", "").removeprefix("Bearer ").strip()
        if not supplied or not hmac.compare_digest(supplied, token):
            return jsonify(error="unauthorized"), 401
        data = request.get_json(silent=True)
        if not isinstance(data, dict):
            return jsonify(error="expected a JSON request"), 400
        ligand_code = str(data.get("ligand_code", "")).upper().strip()
        if not LIGAND_CODE_RE.fullmatch(ligand_code):
            return jsonify(error="ligand_code must contain 1-4 uppercase letters or digits"), 400
        try:
            mol2_bytes = base64.b64decode(str(data.get("mol2_b64", "")), validate=True)
        except ValueError:
            return jsonify(error="mol2_b64 is not valid base64"), 400
        if not mol2_bytes or len(mol2_bytes) > MAX_MOL2_BYTES:
            return jsonify(error="MOL2 exceeds the allowed size"), 413
        if b"@<TRIPOS>MOLECULE" not in mol2_bytes.upper() or b"\x00" in mol2_bytes:
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
                return jsonify(error="CGenFF timed out"), 504
            stdout, stderr = completed.stdout, completed.stderr
        if completed.returncode != 0 or b"skipped molecule" in stderr.lower() or b"unfulfilled valence" in stderr.lower():
            return jsonify(
                error="CGenFF could not parameterize this molecule",
                log_b64=base64.b64encode(stderr[-65536:]).decode("ascii"),
            ), 422
        if not stdout.strip():
            return jsonify(error="CGenFF returned an empty stream file"), 502
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
