from __future__ import annotations

import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence

import numpy as np


SCHEMA_VERSION = 1
SMALL_FILE_SHA256_LIMIT = 5 * 1024 * 1024


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _normalize_signature_value(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _normalize_signature_value(value[key]) for key in sorted(value)}
    if isinstance(value, (list, tuple)):
        return [_normalize_signature_value(item) for item in value]
    if isinstance(value, set):
        return [_normalize_signature_value(item) for item in sorted(value, key=lambda item: str(item))]
    return value


def compute_signature(payload: Any) -> str:
    normalized = _normalize_signature_value(payload)
    serialized = json.dumps(normalized, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def sha256_file(path: os.PathLike[str] | str, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint_file(
    path: os.PathLike[str] | str,
    *,
    hash_small_file: bool = True,
    small_file_limit: int = SMALL_FILE_SHA256_LIMIT,
) -> Dict[str, Any]:
    resolved = Path(path).expanduser()
    record: Dict[str, Any] = {
        "path": str(resolved.resolve(strict=False)),
        "name": resolved.name,
        "exists": resolved.exists(),
    }
    if not resolved.exists():
        return record

    stat = resolved.stat()
    record.update(
        {
            "size_bytes": int(stat.st_size),
            "mtime_ns": int(stat.st_mtime_ns),
        }
    )
    if hash_small_file and resolved.is_file() and stat.st_size <= int(small_file_limit):
        record["sha256"] = sha256_file(resolved)
    return record


def fingerprint_files(
    files: Mapping[str, os.PathLike[str] | str],
    *,
    hash_small_files: bool = True,
    small_file_limit: int = SMALL_FILE_SHA256_LIMIT,
) -> Dict[str, Dict[str, Any]]:
    return {
        str(name): fingerprint_file(
            path,
            hash_small_file=hash_small_files,
            small_file_limit=small_file_limit,
        )
        for name, path in files.items()
    }


def file_is_nonempty(path: os.PathLike[str] | str) -> bool:
    candidate = Path(path)
    return candidate.exists() and candidate.is_file() and candidate.stat().st_size > 0


def inspect_required_outputs(paths: Sequence[os.PathLike[str] | str]) -> Dict[str, Any]:
    normalized = [Path(path) for path in paths]
    missing = [str(path) for path in normalized if not path.exists()]
    empty = [str(path) for path in normalized if path.exists() and path.is_file() and path.stat().st_size <= 0]
    present = [str(path) for path in normalized if file_is_nonempty(path)]
    return {
        "present": present,
        "missing": missing,
        "empty": empty,
        "valid": not missing and not empty,
    }


def require_nonempty_files(paths: Sequence[os.PathLike[str] | str]) -> None:
    status = inspect_required_outputs(paths)
    if status["valid"]:
        return
    problems = []
    if status["missing"]:
        problems.append("missing: " + ", ".join(status["missing"]))
    if status["empty"]:
        problems.append("empty: " + ", ".join(status["empty"]))
    raise FileNotFoundError("; ".join(problems))


def read_json_manifest(path: os.PathLike[str] | str) -> Dict[str, Any]:
    candidate = Path(path)
    if not candidate.exists():
        return {}
    with candidate.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    return data if isinstance(data, dict) else {}


def write_json_atomic(path: os.PathLike[str] | str, payload: Mapping[str, Any]) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=str(target.parent),
        prefix=f".{target.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temp_name = handle.name
    os.replace(temp_name, target)
    return target


def save_npz_atomic(path: os.PathLike[str] | str, **arrays: Any) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "wb",
        dir=str(target.parent),
        prefix=f".{target.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        np.savez_compressed(handle, **arrays)
        temp_name = handle.name
    os.replace(temp_name, target)
    return target


def load_npz_dict(path: os.PathLike[str] | str, *, allow_legacy_object_arrays: bool = False) -> Dict[str, np.ndarray]:
    """Load a PyMACS NPZ cache without accepting pickles by default.

    Older, locally generated interface-presence caches stored pandas string
    columns as object arrays.  They are accepted only by the explicit legacy
    path after the caller has validated the cache manifest.
    """
    with np.load(Path(path), allow_pickle=allow_legacy_object_arrays) as handle:
        return {key: handle[key] for key in handle.files}


def _maybe_relpath(path: Path, base_dir: Path) -> str:
    try:
        return str(path.resolve().relative_to(base_dir.resolve()))
    except Exception:
        return str(path.resolve())


def describe_output_files(
    paths: Sequence[os.PathLike[str] | str],
    *,
    base_dir: Optional[os.PathLike[str] | str] = None,
) -> list[Dict[str, Any]]:
    base = Path(base_dir) if base_dir is not None else None
    descriptions = []
    for raw_path in paths:
        path = Path(raw_path)
        record = {
            "path": _maybe_relpath(path, base) if base is not None else str(path.resolve(strict=False)),
            "exists": path.exists(),
        }
        if path.exists():
            stat = path.stat()
            record["size_bytes"] = int(stat.st_size)
            record["mtime_ns"] = int(stat.st_mtime_ns)
        descriptions.append(record)
    return descriptions


def _default_manifest(
    *,
    analysis_script: str,
    analysis_source_path: os.PathLike[str] | str,
    input_files: Mapping[str, os.PathLike[str] | str],
    system: Optional[str] = None,
    mode: Optional[str] = None,
) -> Dict[str, Any]:
    resolved_inputs = {name: str(Path(path).resolve(strict=False)) for name, path in input_files.items()}
    manifest: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": utc_now_iso(),
        "analysis_script": analysis_script,
        "analysis_script_version": fingerprint_file(analysis_source_path, hash_small_file=True),
        "input_files": resolved_inputs,
        "input_fingerprints": fingerprint_files(input_files),
        "system": system,
        "mode": mode,
        "stages": {},
    }
    return manifest


def record_stage_completion(
    manifest_path: os.PathLike[str] | str,
    *,
    stage_name: str,
    analysis_script: str,
    analysis_source_path: os.PathLike[str] | str,
    input_files: Mapping[str, os.PathLike[str] | str],
    data_parameters: Optional[Mapping[str, Any]] = None,
    plot_parameters: Optional[Mapping[str, Any]] = None,
    outputs: Optional[Sequence[os.PathLike[str] | str]] = None,
    plot_outputs: Optional[Sequence[os.PathLike[str] | str]] = None,
    system: Optional[str] = None,
    mode: Optional[str] = None,
    legacy_adopted: bool = False,
    status: str = "complete",
    extra_stage_fields: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    manifest_file = Path(manifest_path)
    manifest = read_json_manifest(manifest_file)
    if not manifest:
        manifest = _default_manifest(
            analysis_script=analysis_script,
            analysis_source_path=analysis_source_path,
            input_files=input_files,
            system=system,
            mode=mode,
        )

    manifest["schema_version"] = SCHEMA_VERSION
    manifest["generated_utc"] = utc_now_iso()
    manifest["analysis_script"] = analysis_script
    manifest["analysis_script_version"] = fingerprint_file(analysis_source_path, hash_small_file=True)
    manifest["input_files"] = {name: str(Path(path).resolve(strict=False)) for name, path in input_files.items()}
    manifest["input_fingerprints"] = fingerprint_files(input_files)
    if system is not None:
        manifest["system"] = system
    if mode is not None:
        manifest["mode"] = mode

    data_parameters = dict(data_parameters or {})
    plot_parameters = dict(plot_parameters or {})
    outputs = list(outputs or [])
    plot_outputs = list(plot_outputs or [])

    stage_record: Dict[str, Any] = {
        "stage": stage_name,
        "status": status,
        "legacy_adopted": bool(legacy_adopted),
        "data_signature": compute_signature(
            {
                "inputs": manifest["input_fingerprints"],
                "parameters": data_parameters,
            }
        ),
        "inputs": manifest["input_fingerprints"],
        "parameters": {
            "data_affecting": data_parameters,
            "plot_only": plot_parameters,
        },
        "outputs": describe_output_files(outputs, base_dir=manifest_file.parent),
        "plot_outputs": describe_output_files(plot_outputs, base_dir=manifest_file.parent),
        "completed_utc": utc_now_iso(),
    }
    if extra_stage_fields:
        stage_record.update(dict(extra_stage_fields))

    manifest.setdefault("stages", {})
    manifest["stages"][stage_name] = stage_record
    manifest["parameters"] = {
        "data_affecting": data_parameters,
        "plot_only": plot_parameters,
    }
    manifest["outputs"] = {
        "numerical": describe_output_files(outputs, base_dir=manifest_file.parent),
        "plots": describe_output_files(plot_outputs, base_dir=manifest_file.parent),
    }
    write_json_atomic(manifest_file, manifest)
    return manifest


def validate_stage_cache(
    *,
    manifest_path: os.PathLike[str] | str,
    stage_name: str,
    required_outputs: Sequence[os.PathLike[str] | str],
    current_data_signature: Optional[str] = None,
    checkpoint_path: Optional[os.PathLike[str] | str] = None,
) -> Dict[str, Any]:
    checkpoint_exists = True if checkpoint_path is None else Path(checkpoint_path).exists()
    output_status = inspect_required_outputs(required_outputs)
    manifest = read_json_manifest(manifest_path)
    stage_record = manifest.get("stages", {}).get(stage_name) if manifest else None

    if not checkpoint_exists and checkpoint_path is not None:
        return {
            "reusable": False,
            "legacy_adopt": False,
            "reason": "checkpoint_missing",
            "outputs": output_status,
            "manifest": manifest,
            "stage_record": stage_record,
        }

    if not output_status["valid"]:
        return {
            "reusable": False,
            "legacy_adopt": False,
            "reason": "required_outputs_missing",
            "outputs": output_status,
            "manifest": manifest,
            "stage_record": stage_record,
        }

    if stage_record and stage_record.get("status") == "complete":
        cached_signature = stage_record.get("data_signature")
        if current_data_signature and cached_signature != current_data_signature:
            return {
                "reusable": False,
                "legacy_adopt": False,
                "reason": "data_signature_mismatch",
                "outputs": output_status,
                "manifest": manifest,
                "stage_record": stage_record,
            }
        return {
            "reusable": True,
            "legacy_adopt": False,
            "reason": "manifest_match",
            "outputs": output_status,
            "manifest": manifest,
            "stage_record": stage_record,
        }

    return {
        "reusable": True,
        "legacy_adopt": True,
        "reason": "legacy_outputs_complete",
        "outputs": output_status,
        "manifest": manifest,
        "stage_record": stage_record,
    }
