#!/usr/bin/env python3
"""Small, dependency-free client for the PyMACS CGenFF Funnel service.

The token authenticates one narrowly scoped operation: parameterizing one
already-normalized MOL2 file.  It never grants shell or SSH access to Kyle.
"""

from __future__ import annotations

import base64
import json
import re
from pathlib import Path
from typing import Final
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


MAX_MOL2_BYTES: Final = 2 * 1024 * 1024
LIGAND_CODE_RE: Final = re.compile(r"^[A-Z0-9]{1,4}$")


class CGenFFFunnelError(RuntimeError):
    """Raised when a remote CGenFF request cannot safely be completed."""


def read_token(token_file: str | Path) -> str:
    path = Path(token_file).expanduser()
    try:
        token = path.read_text(encoding="utf-8").strip()
    except OSError as exc:
        raise CGenFFFunnelError(f"Cannot read CGenFF Funnel token file: {path}") from exc
    if len(token) < 24:
        raise CGenFFFunnelError(f"CGenFF Funnel token in {path} is missing or too short.")
    return token


def parameterize_mol2(
    *,
    endpoint: str,
    token_file: str | Path,
    ligand_code: str,
    mol2_path: str | Path,
    timeout_seconds: int = 180,
) -> tuple[str, str]:
    """Submit one MOL2 and return ``(stream_file, service_log)`` text.

    The API is intentionally synchronous. CGenFF jobs are small and this keeps
    Triton's setup process deterministic: a successful response means the
    corresponding ``.str`` can immediately be passed to charmm2gmx.
    """
    ligand_code = ligand_code.upper().strip()
    if not LIGAND_CODE_RE.fullmatch(ligand_code):
        raise CGenFFFunnelError("Ligand code must contain 1–4 uppercase letters or digits.")
    mol2 = Path(mol2_path)
    try:
        mol2_bytes = mol2.read_bytes()
    except OSError as exc:
        raise CGenFFFunnelError(f"Cannot read prepared MOL2: {mol2}") from exc
    if not mol2_bytes or len(mol2_bytes) > MAX_MOL2_BYTES:
        raise CGenFFFunnelError(f"Prepared MOL2 must be between 1 byte and {MAX_MOL2_BYTES} bytes.")
    if b"@<TRIPOS>MOLECULE" not in mol2_bytes.upper():
        raise CGenFFFunnelError("Prepared ligand file is not a valid MOL2 (missing MOLECULE section).")

    url = endpoint.rstrip("/") + "/v1/parameterize"
    token = read_token(token_file)
    payload = json.dumps(
        {
            "ligand_code": ligand_code,
            "mol2_b64": base64.b64encode(mol2_bytes).decode("ascii"),
        }
    ).encode("utf-8")
    request = Request(
        url,
        data=payload,
        method="POST",
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
    )
    try:
        with urlopen(request, timeout=timeout_seconds) as response:  # nosec B310: configured HTTPS Funnel URL
            raw = response.read()
    except HTTPError as exc:
        try:
            detail = json.loads(exc.read().decode("utf-8")).get("error", exc.reason)
        except (UnicodeDecodeError, json.JSONDecodeError):
            detail = exc.reason
        raise CGenFFFunnelError(f"CGenFF Funnel rejected the request ({exc.code}): {detail}") from exc
    except URLError as exc:
        raise CGenFFFunnelError(f"Cannot reach the CGenFF Funnel: {exc.reason}") from exc

    try:
        response_data = json.loads(raw.decode("utf-8"))
        stream = base64.b64decode(response_data["str_b64"], validate=True).decode("utf-8")
        service_log = base64.b64decode(response_data.get("log_b64", ""), validate=True).decode("utf-8")
    except (KeyError, TypeError, ValueError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise CGenFFFunnelError("CGenFF Funnel returned an invalid response.") from exc
    if not stream.strip():
        raise CGenFFFunnelError("CGenFF Funnel returned an empty stream file.")
    return stream, service_log
