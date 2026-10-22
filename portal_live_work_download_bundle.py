#!/usr/bin/env python3
"""Internal lab portal for staging and submitting docking and PyMACS jobs.

The service deliberately accepts archives only, safely extracts them, copies
them to an HPC login host, then invokes a fixed scheduler command.  It never
passes browser input through a shell.
"""
from __future__ import annotations

import datetime as dt
import base64
import csv
import fcntl
import json
import os
import re
import shutil
import shlex
import sqlite3
import subprocess
import sys
import tempfile
import tarfile
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.error import URLError, HTTPError
from urllib.request import Request, urlopen
import uuid
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from flask import Flask, abort, flash, jsonify, redirect, render_template, request, send_file, send_from_directory, url_for
from werkzeug.utils import secure_filename

MAX_ARCHIVE_BYTES = 2 * 1024 * 1024 * 1024
JOB_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$")
SUBMITTER_NAME = re.compile(r"^[A-Za-z][A-Za-z -]{0,39}$")
EMAIL = re.compile(r"^[^\s@]+@[^\s@]+\.[^\s@]+$")
JOB_ID = re.compile(r"Job <(\d+)>")
STRUCTURE_SUFFIXES = {".pdb", ".cif", ".mmcif"}
RCSB_ID = re.compile(r"^[A-Za-z0-9]{4}$")
RESULT_FILE_SUFFIXES = {".csv", ".html", ".pdf", ".png"}
RESULT_TRAJECTORY_NAMES = {
    "Final_Trajectory.pdb", "Final_Trajectory.xtc",
    "binding_pocket_only.pdb", "binding_pocket_only.xtc",
}
RESULT_FILE_LIMIT = 2 * 1024 * 1024 * 1024
DEFAULT_LEGACY_MD_COLLECTIONS = (
    {
        "root": "/projectnb/triton/home/jxs794/Molecular_Dynamics/Angela-OCT26",
        "campaign": "Angela · October 2026",
        "submitter": "Angela",
        "protein_name": "BTK",
        "ligand_name": "LIG",
    },
)

# Migrated from the former SouthPark Servers inventory.  These are names only:
# connection details remain in Kyle's SSH configuration, never in the portal.
LAB_MACHINES = (
    {"name": "Amy", "host": "amy", "role": "Lab server"},
    {"name": "Bender", "host": "bender", "role": "Lab server"},
    {"name": "Butters", "host": "butters", "role": "Workstation"},
    {"name": "Cartman", "host": "cartman", "role": "Workstation"},
    {"name": "Jimmy", "host": "jimmy", "role": "Workstation"},
    {"name": "Kenny", "host": "kenny", "role": "Workstation"},
    {"name": "Kyle", "host": "kyle", "role": "Portal host"},
    {"name": "Leela", "host": "leela", "role": "Lab server"},
    {"name": "MrHankey", "host": "mrhankey", "role": "Workstation"},
    {"name": "Pegasus", "host": "pegasus", "role": "HPC cluster"},
    {"name": "Randy", "host": "randy", "role": "Workstation"},
    {"name": "Stan", "host": "stan", "role": "Workstation"},
    {"name": "Timmy", "host": "timmy", "role": "Workstation"},
    {"name": "Token", "host": "token", "role": "Workstation"},
    {"name": "Towelie", "host": "towelie", "role": "Workstation"},
    {"name": "Triton", "host": "triton", "role": "HPC cluster"},
)


class PortalError(RuntimeError):
    pass


@dataclass(frozen=True)
class Settings:
    stage_root: Path
    pegasus_remote_root: str
    triton_remote_root: str
    pegasus_host: str
    triton_host: str
    pymacs_root: str
    cgenff_url: str
    cgenff_token_file: Path
    ssh_bin: str = "ssh"
    scp_bin: str = "scp"
    # Reviewed lab copy used for every submitted docking campaign. Uploaded
    # ZIPs provide scientific data only; they never provide executable code.
    autodock_root: str = "/home/jxs794/INST11-Docking"
    pegasus_docking_root: str = "/nethome/jxs794/Docking_Experiments/IE1-Docking"
    pegasus_results_root: str = "/projectnb/cheminfo/home/jxs794/AutoDockVina-Results"
    md_results_root: str = "/projectnb/cheminfo/home/jxs794/PyMACS-Results"
    legacy_md_collections: tuple[dict[str, str], ...] = ()
    legacy_md_auto_restart: bool = False


def legacy_md_collections_from_env() -> tuple[dict[str, str], ...]:
    """Read explicitly configured, pre-portal MD collections.

    A collection is read-only.  It exists to surface completed legacy runs,
    never to infer a scheduler job or alter the original project.
    """
    raw = os.environ.get("LAB_PORTAL_LEGACY_MD_COLLECTIONS")
    configured: object = DEFAULT_LEGACY_MD_COLLECTIONS
    if raw:
        try:
            configured = json.loads(raw)
        except json.JSONDecodeError:
            configured = ()
    if not isinstance(configured, (list, tuple)):
        return ()
    items: list[dict[str, str]] = []
    for value in configured:
        if not isinstance(value, dict):
            continue
        root = str(value.get("root", "")).rstrip("/")
        if not root.startswith("/projectnb/"):
            continue
        items.append({
            "root": root,
            "campaign": str(value.get("campaign") or Path(root).name),
            "submitter": str(value.get("submitter") or "External"),
            "protein_name": str(value.get("protein_name") or ""),
            "ligand_name": str(value.get("ligand_name") or "LIG"),
        })
    return tuple(items)


def settings_from_env() -> Settings:
    missing = [key for key in ("LAB_PORTAL_PEGASUS_REMOTE_ROOT", "LAB_PORTAL_TRITON_REMOTE_ROOT", "LAB_PORTAL_PYMACS_ROOT") if not os.environ.get(key)]
    if missing:
        raise PortalError("Missing required configuration: " + ", ".join(missing))
    return Settings(
        stage_root=Path(os.environ.get("LAB_PORTAL_STAGE_ROOT", "/var/lib/lab-job-portal/jobs")),
        pegasus_remote_root=os.environ["LAB_PORTAL_PEGASUS_REMOTE_ROOT"].rstrip("/"),
        triton_remote_root=os.environ["LAB_PORTAL_TRITON_REMOTE_ROOT"].rstrip("/"),
        pegasus_host=os.environ.get("LAB_PORTAL_PEGASUS_HOST", "pegasus"),
        triton_host=os.environ.get("LAB_PORTAL_TRITON_HOST", "triton"),
        pymacs_root=os.environ["LAB_PORTAL_PYMACS_ROOT"].rstrip("/"),
        cgenff_url=os.environ.get("LAB_PORTAL_CGENFF_URL", "http://127.0.0.1:8787/v1/parameterize"),
        cgenff_token_file=Path(os.environ.get("LAB_PORTAL_CGENFF_TOKEN_FILE", "/home/jxs794/.config/pymacs/cgenff-funnel.token")),
        ssh_bin=os.environ.get("LAB_PORTAL_SSH_BIN", "ssh"),
        scp_bin=os.environ.get("LAB_PORTAL_SCP_BIN", "scp"),
        autodock_root=os.environ.get("LAB_PORTAL_AUTODOCK_ROOT", "/home/jxs794/INST11-Docking"),
        pegasus_docking_root=os.environ.get("LAB_PORTAL_PEGASUS_DOCKING_ROOT", "/nethome/jxs794/Docking_Experiments/IE1-Docking"),
        pegasus_results_root=os.environ.get("LAB_PORTAL_PEGASUS_RESULTS_ROOT", "/projectnb/cheminfo/home/jxs794/AutoDockVina-Results").rstrip("/"),
        md_results_root=os.environ.get("LAB_PORTAL_MD_RESULTS_ROOT", "/projectnb/cheminfo/home/jxs794/PyMACS-Results").rstrip("/"),
        legacy_md_collections=legacy_md_collections_from_env(),
        legacy_md_auto_restart=os.environ.get("LAB_PORTAL_LEGACY_MD_AUTORESTART", "").strip().lower() in {"1", "true", "yes"},
    )


def safe_extract(archive: Path, destination: Path) -> None:
    """Extract a normal ZIP without allowing path traversal, links, or zip bombs."""
    with zipfile.ZipFile(archive) as source:
        infos = source.infolist()
        if len(infos) > 20_000:
            raise PortalError("Archive contains too many files (limit: 20,000).")
        total = sum(info.file_size for info in infos)
        if total > 10 * 1024 * 1024 * 1024:
            raise PortalError("Archive expands beyond the 10 GB safety limit.")
        for info in infos:
            path = PurePosixPath(info.filename)
            if path.is_absolute() or ".." in path.parts or not path.parts:
                raise PortalError("Archive contains an unsafe path.")
            # UNIX symlink bit in external attributes.
            if (info.external_attr >> 16) & 0o170000 == 0o120000:
                raise PortalError("Archive may not contain symlinks.")
        source.extractall(destination)


def project_root(extracted: Path) -> Path:
    children = [item for item in extracted.iterdir() if item.name != "__MACOSX"]
    if len(children) == 1 and children[0].is_dir():
        return children[0]
    return extracted


def inspect_docking_package(root: Path) -> dict:
    """Describe the actual ligand inputs shipped by an AutoDock package."""
    ignored = {"__MACOSX", ".git", "logs", "results"}
    inputs = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or any(part in ignored or part.startswith(".") for part in path.relative_to(root).parts):
            continue
        suffix = path.suffix.lower()
        count, headers = 0, []
        if suffix == ".sdf":
            count = sum(1 for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip() == "$$$$")
        elif suffix in {".smi", ".smiles"}:
            count = sum(1 for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip())
        elif suffix == ".csv":
            if "center" in path.name.lower():
                continue
            with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
                rows = list(csv.reader(handle)); headers = rows[0] if rows else []; count = sum(1 for row in rows[1:] if any(cell.strip() for cell in row))
        else:
            continue
        inputs.append({"path": path.relative_to(root).as_posix(), "type": suffix[1:], "ligands": count, "headers": headers})
    return {"inputs": inputs, "total_ligands": sum(item["ligands"] for item in inputs),
            "has_splitter": bool(list(root.rglob("0_LIGSPLIT.py"))), "has_confgen": bool(list(root.rglob("1B_confgen_batch.py")))}


def _safe_relative(value: str) -> Path:
    path = PurePosixPath(value or "")
    if not path.parts or path.is_absolute() or ".." in path.parts:
        raise PortalError("Choose a ligand input from the uploaded package.")
    return Path(*path.parts)


def _positive_int(value: str, label: str, *, default: int, minimum: int = 1) -> int:
    try:
        number = int(value or default)
    except ValueError as exc:
        raise PortalError(f"{label} must be a whole number.") from exc
    if number < minimum:
        raise PortalError(f"{label} must be at least {minimum}.")
    return number


def _number(value: str, label: str, *, default: float, minimum: float = 0) -> float:
    try:
        number = float(value or default)
    except ValueError as exc:
        raise PortalError(f"{label} must be numeric.") from exc
    if number < minimum:
        raise PortalError(f"{label} must be at least {minimum}.")
    return number


def _safe_lsf_value(value: str, label: str, default: str) -> str:
    cleaned = (value or default).strip()
    if not re.fullmatch(r"[A-Za-z0-9_.:-]+", cleaned):
        raise PortalError(f"{label} contains unsupported characters.")
    return cleaned


def _shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def build_docking_campaign(settings: Settings, root: Path, form) -> dict[str, object]:
    """Turn an AutoDock data ZIP into trusted per-batch LSF jobs on Kyle."""
    report = inspect_docking_package(root)
    selected = _safe_relative(form.get("ligand_input", ""))
    input_path = (root / selected).resolve()
    if root.resolve() not in input_path.parents or not input_path.is_file() or input_path.suffix.lower() not in {".csv", ".sdf", ".smi", ".smiles"}:
        raise PortalError("Choose one raw CSV, SDF, SMI, or SMILES input from the ZIP.")
    centers = sorted(root.rglob("vina_centers*.csv"))
    receptors = [path for path in root.rglob("*.pdbqt") if "Ligands_PDBQT" not in path.as_posix()]
    if not centers or not receptors:
        raise PortalError("The docking ZIP needs prepared receptor .pdbqt file(s) and vina_centers.csv.")
    for script in ("0_LIGSPLIT.py", "1_ConformerGeneration.py", "3_Complete_batch_docking.py", "ligand_manifest.py", "ligand_naming.py"):
        source = Path(settings.autodock_root) / script
        if not source.is_file():
            raise PortalError(f"Kyle is missing the reviewed docking helper: {source}")
        shutil.copy2(source, root / script)

    batch_size = _positive_int(form.get("batch_size", ""), "Ligands per batch", default=500)
    batch_prefix = _safe_lsf_value(form.get("batch_prefix", "batch"), "Batch prefix", "batch")
    subprocess.run([sys.executable, str(root / "0_LIGSPLIT.py"), "--input", str(input_path), "--output", str(root / "batches"), "--group-size", str(batch_size), "--prefix", batch_prefix, "--yes"], check=True, text=True, capture_output=True, timeout=3600)
    batches = sorted(path for path in (root / "batches").glob(f"{batch_prefix}_*") if path.is_dir())
    if not batches:
        raise PortalError("The ligand splitter did not create any batches.")
    if input_path.suffix.lower() == ".csv" and (not form.get("csv_smiles_col", "").strip() or not form.get("csv_id_col", "").strip()):
        raise PortalError("Choose both the CSV SMILES and compound-ID columns.")

    opts = {
        "confgen_poses": _positive_int(form.get("confgen_poses", ""), "Conformers per ligand", default=64),
        "confgen_workers": _positive_int(form.get("confgen_workers", ""), "ConfGen workers", default=16),
        "vina_poses": _positive_int(form.get("vina_poses", ""), "Vina poses", default=20),
        "vina_workers": _positive_int(form.get("vina_workers", ""), "Vina workers", default=8),
        "queue": _safe_lsf_value(form.get("queue", "normal"), "Queue", "normal"),
        "project": _safe_lsf_value(form.get("project", "brd"), "Project", "brd"),
        "confgen_walltime": _safe_lsf_value(form.get("confgen_walltime", "48:00"), "ConfGen walltime", "48:00"),
        "vina_walltime": _safe_lsf_value(form.get("vina_walltime", "240:00"), "Vina walltime", "240:00"),
        "csv_smiles_col": form.get("csv_smiles_col", "").strip(), "csv_id_col": form.get("csv_id_col", "").strip(),
        "enumerate_protomers": form.get("enumerate_protomers") == "on",
        "ph_min": _number(form.get("ph_min", ""), "Minimum pH", default=6.8),
        "ph_max": _number(form.get("ph_max", ""), "Maximum pH", default=7.4),
        "ph_precision": _number(form.get("ph_precision", ""), "pH precision", default=0.5, minimum=0.01),
        "max_protomers": _positive_int(form.get("max_protomers", ""), "Maximum protomers", default=4),
        "max_tautomers": _positive_int(form.get("max_tautomers", ""), "Maximum tautomers", default=4),
        "max_transforms": _positive_int(form.get("max_transforms", ""), "Maximum tautomer transforms", default=200),
    }
    if opts["ph_max"] < opts["ph_min"]:
        raise PortalError("Maximum pH must be greater than or equal to minimum pH.")
    state_flags = f" --ph-min {opts['ph_min']} --ph-max {opts['ph_max']} --ph-precision {opts['ph_precision']} --max-protomers {opts['max_protomers']} --max-tautomers {opts['max_tautomers']} --max-transforms {opts['max_transforms']}"
    if opts["enumerate_protomers"]:
        state_flags += " --enumerate-protomers"
    receptor_dir = receptors[0].parent.relative_to(root).as_posix()
    centers_csv = centers[0].relative_to(root).as_posix()
    for position, batch in enumerate(batches, start=1):
        job_tag = f"{batch_prefix}-{position:03d}"
        # LSF opens -oo/-eo before it starts the shell script, so this must
        # exist while the portal stages the campaign—not inside the job body.
        (batch / "logs").mkdir(exist_ok=True)
        files = sorted(path for path in batch.iterdir() if path.is_file() and path.suffix.lower() in {".csv", ".sdf", ".smi", ".smiles"})
        if len(files) != 1:
            raise PortalError(f"Batch {batch.name} has no unambiguous ligand input.")
        source = files[0]
        if source.suffix.lower() == ".csv":
            command = f'"$PYBIN" "$ROOT/1_ConformerGeneration.py" --mode 1 --csv {_shell_quote(source.name)} --smiles-col {_shell_quote(opts["csv_smiles_col"])} --id-col {_shell_quote(opts["csv_id_col"])}'
        elif source.suffix.lower() == ".sdf":
            command = '"$PYBIN" "$ROOT/1_ConformerGeneration.py" --mode 2 --folder . --filetype sdf'
        else:
            command = '"$PYBIN" "$ROOT/1_ConformerGeneration.py" --mode 2 --folder . --filetype smiles'
        header = lambda stage, workers, walltime: f"#!/bin/bash\n#BSUB -J portal-{stage}-{job_tag}\n#BSUB -q {opts['queue']}\n#BSUB -P {opts['project']}\n#BSUB -n {workers}\n#BSUB -R \"rusage[mem=8000]\"\n#BSUB -W {walltime}\n#BSUB -oo logs/{stage}-%J.out\n#BSUB -eo logs/{stage}-%J.err\nset -euo pipefail\n"
        # LSF starts the job in the submitted batch directory.  ``$0`` is
        # merely the shell when a job is fed through ``bsub < file``, so use
        # that execution CWD rather than trying to infer the script path.
        prelude = 'BATCH_DIR="$PWD"\nROOT="$(cd "$BATCH_DIR/../.." && pwd)"\ncd "$BATCH_DIR"\nmkdir -p logs\nsource /nethome/jxs794/miniconda3/etc/profile.d/conda.sh\nconda activate vina_env\nPYBIN="$(command -v python3 || command -v python)"\n'
        confgen_handoff = (
            command
            + f" --num-confs {opts['confgen_poses']} --workers {opts['confgen_workers']}{state_flags}\n"
            # Each batch is a clean directory, so ConfGen must create exactly
            # one PDBQT directory.  Record that concrete path for its paired
            # Vina job instead of having Vina guess from similarly named dirs.
            + 'mapfile -t PDBQT_DIRS < <(find "$BATCH_DIR" -maxdepth 1 -type d -name "*_Ligands_PDBQT_*" -print | sort)\n'
            + '[ "${#PDBQT_DIRS[@]}" -eq 1 ] || { echo "Expected exactly one ConfGen PDBQT output directory; found ${#PDBQT_DIRS[@]}" >&2; exit 2; }\n'
            + 'printf "%s\\n" "${PDBQT_DIRS[0]}" > confgen_pdbqt_dir.txt\n'
        )
        (batch / "run_confgen.lsf").write_text(header("conf", opts["confgen_workers"], opts["confgen_walltime"]) + prelude + confgen_handoff)
        vina = header("vina", opts["vina_workers"], opts["vina_walltime"]) + prelude
        vina += 'HANDOFF_FILE="$BATCH_DIR/confgen_pdbqt_dir.txt"\n[ -s "$HANDOFF_FILE" ] || { echo "ConfGen handoff file is missing or empty" >&2; exit 2; }\nLIGAND_DIR="$(< "$HANDOFF_FILE")"\n[ -d "$LIGAND_DIR" ] || { echo "ConfGen PDBQT directory does not exist: $LIGAND_DIR" >&2; exit 2; }\nfind "$LIGAND_DIR" -maxdepth 1 -type f -name "*.pdbqt" -print -quit | grep -q . || { echo "ConfGen PDBQT directory is empty: $LIGAND_DIR" >&2; exit 2; }\n'
        # This trusted runtime file is staged after the final Cheminfo campaign
        # path is known.  Do not rely on bsub inheriting dispatcher variables.
        vina += f'RUNTIME_CONFIG="$ROOT/portal_runtime.env"\n[ -s "$RUNTIME_CONFIG" ] || {{ echo "Portal runtime configuration is missing" >&2; exit 2; }}\nsource "$RUNTIME_CONFIG"\n: "${{CAMPAIGN_RESULTS_ROOT:?Portal results root is missing}}"\nSHARD_DIR="$CAMPAIGN_RESULTS_ROOT/shards/{batch.name}"\nmkdir -p "$SHARD_DIR"\nexport VINA_SCRATCH_ROOT="/tmp"\n"$PYBIN" "$ROOT/3_Complete_batch_docking.py" --receptors "$ROOT/{receptor_dir}" --ligands "$LIGAND_DIR" --centers_csv "$ROOT/{centers_csv}" --poses {opts["vina_poses"]} --reserve_cores 0 --compact-results --results-dir "$SHARD_DIR" --compact-scratch-root "$VINA_SCRATCH_ROOT"\n'
        (batch / "run_vina.lsf").write_text(vina)
    manifest = {"schema_version": 1, "input": selected.as_posix(), "batch_size": batch_size, "batches": [batch.name for batch in batches], "options": opts}
    (root / "campaign.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def require_md_project(root: Path) -> Path:
    configs = list(root.rglob("pymacs_run.json"))
    if len(configs) != 1:
        raise PortalError("MD ZIP must contain exactly one PyMACS project with pymacs_run.json.")
    project = configs[0].parent
    if not any(project.glob("*.pdb")) and not any(project.glob("*.cif")) and not any(project.glob("*.mmcif")):
        raise PortalError("The PyMACS project needs its input PDB or mmCIF structure.")
    # The portal supplies its reviewed PyMACS scripts on Triton. Do not execute
    # Python uploaded through a browser under the shared submission account.
    prohibited = [path.name for path in project.rglob("*.py")]
    if prohibited:
        raise PortalError("MD ZIPs may contain configuration and input data only; do not include Python scripts.")
    return project


def fetch_rcsb_pdb(pdb_id: str, destination: Path) -> str:
    """Download coordinates from RCSB, preferring PDB then falling back to mmCIF."""
    identifier = pdb_id.strip().upper()
    if not RCSB_ID.fullmatch(identifier):
        raise PortalError("RCSB PDB ID must be exactly four letters/numbers.")
    failures = []
    for suffix, valid in (
        (".pdb", lambda data: data.startswith((b"HEADER", b"ATOM", b"HETATM"))),
        (".cif", lambda data: data.lstrip().startswith(b"data_")),
    ):
        try:
            with urlopen(f"https://files.rcsb.org/download/{identifier}{suffix}", timeout=30) as response:
                content = response.read(50 * 1024 * 1024 + 1)
        except URLError as exc:
            failures.append(f"{suffix}: {exc.reason}")
            continue
        if content and len(content) <= 50 * 1024 * 1024 and valid(content):
            filename = f"{identifier}{suffix}"
            destination.joinpath(filename).write_bytes(content)
            return filename
        failures.append(f"{suffix}: unusable coordinate content")
    raise PortalError(
        f"Could not download usable PDB or mmCIF coordinates for {identifier} from RCSB "
        f"({' ; '.join(failures)})."
    )


def default_chain_map(structure: Path, system_name: str = "") -> str | None:
    """Give headless PyMACS a deterministic name for every detected PDB chain."""
    if structure.suffix.lower() != ".pdb":
        return None
    chains: list[str] = []
    for line in structure.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM"):
            chain = line[21:22].strip()
            if chain and chain not in chains:
                chains.append(chain)
    # PyMACS carries these labels into chain-specific analysis figures and
    # interface graphs.  Prefixing them avoids the generic A/B labels while
    # retaining the unambiguous source chain identifier.
    prefix = re.sub(r"[^A-Za-z0-9_.-]+", "_", system_name).strip("_.-")
    return ",".join(f"{chain}:{prefix}_{chain}" if prefix else f"{chain}:{chain}" for chain in chains) or None


def inspect_structure_file(structure: Path) -> dict[str, object]:
    """Return only the choices a scientist needs before starting PyMACS."""
    chains: dict[str, set[tuple[str, str]]] = {}
    hetero: dict[str, int] = {}
    hetero_instances: dict[tuple[str, str, str], int] = {}
    waters = {"HOH", "WAT", "SOL", "H2O"}
    if structure.suffix.lower() in {".cif", ".mmcif"}:
        lines = structure.read_text(encoding="utf-8", errors="replace").splitlines()
        try:
            start = next(index for index, line in enumerate(lines) if line.strip() == "_atom_site.group_PDB")
        except StopIteration:
            return {"format": "MMCIF", "chains": [], "hetero": [],
                    "note": "This mmCIF does not contain an atom-site coordinate table."}
        headers = []
        index = start
        while index < len(lines) and lines[index].startswith("_atom_site."):
            headers.append(lines[index].strip().split(".", 1)[1])
            index += 1
        for line in lines[index:]:
            if not line.strip() or line.startswith("#"):
                if line.startswith("#"):
                    break
                continue
            try:
                values = shlex.split(line)
            except ValueError:
                continue
            if len(values) < len(headers):
                continue
            atom = dict(zip(headers, values))
            record = atom.get("group_PDB", "")
            residue = atom.get("auth_comp_id") or atom.get("label_comp_id", "")
            chain = atom.get("auth_asym_id") or atom.get("label_asym_id") or "(blank)"
            sequence = atom.get("auth_seq_id") or atom.get("label_seq_id") or ""
            if record == "ATOM" and residue:
                chains.setdefault(chain, set()).add((sequence, residue))
            elif record == "HETATM" and residue and residue not in waters:
                hetero[residue] = hetero.get(residue, 0) + 1
                key = (residue, chain, sequence)
                hetero_instances[key] = hetero_instances.get(key, 0) + 1
        return {"format": "MMCIF", "chains": [{"id": key, "residues": len(value)} for key, value in chains.items()],
                "hetero": [{"code": key, "atoms": value} for key, value in sorted(hetero.items())],
                "hetero_instances": [{"code": code, "chain": chain, "resid": resid, "atoms": atoms,
                                      "token": f"{code}@{chain}:{resid}"} for (code, chain, resid), atoms in sorted(hetero_instances.items())],
                "note": ""}
    for line in structure.read_text(encoding="utf-8", errors="replace").splitlines():
        record, residue = line[:6].strip(), line[17:20].strip().upper()
        if record == "ATOM":
            chain = line[21:22].strip() or "(blank)"
            chains.setdefault(chain, set()).add((line[22:26].strip(), residue))
        elif record == "HETATM" and residue and residue not in waters:
            hetero[residue] = hetero.get(residue, 0) + 1
            chain = line[21:22].strip() or "(blank)"
            resid = line[22:26].strip()
            key = (residue, chain, resid)
            hetero_instances[key] = hetero_instances.get(key, 0) + 1
    return {"format": "PDB", "chains": [{"id": key, "residues": len(value)} for key, value in chains.items()],
            "hetero": [{"code": key, "atoms": value} for key, value in sorted(hetero.items())],
            "hetero_instances": [{"code": code, "chain": chain, "resid": resid, "atoms": atoms,
                                  "token": f"{code}@{chain}:{resid}"} for (code, chain, resid), atoms in sorted(hetero_instances.items())],
            "note": ""}


def structure_preview(structure: Path, maximum_atoms: int = 900) -> list[list[float | int]]:
    """Small coordinate payload for the portal's local, rotatable PDB preview."""
    if structure.suffix.lower() != ".pdb":
        return []
    atoms: list[list[float | int]] = []
    for line in structure.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        # Hetero atoms are shown in green, protein atoms in blue.
        atoms.append([x, y, z, 1 if line.startswith("HETATM") else 0])
        if len(atoms) >= maximum_atoms:
            break
    return atoms


def discover_md_structures(root: Path) -> list[dict[str, object]]:
    """Find each independently runnable structure in an uploaded MD bundle."""
    ignored = {"__MACOSX", ".git", "logs", "results"}
    found = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in STRUCTURE_SUFFIXES:
            continue
        relative = path.relative_to(root)
        if any(part in ignored or part.startswith(".") for part in relative.parts):
            continue
        report = inspect_structure_file(path)
        found.append({"path": relative.as_posix(), "filename": path.name, **report,
                      "preview": structure_preview(path),
                      # 3Dmol supports both PDB and mmCIF for reasonably
                      # sized structures; larger files remain reviewable
                      # without a browser preview.
                      "pdb_text": path.read_text(encoding="utf-8", errors="replace")
                      if path.stat().st_size <= 2 * 1024 * 1024 else "",
                      "viewer_format": "pdb" if path.suffix.lower() == ".pdb" else "cif"})
    if not found:
        raise PortalError("No PDB, CIF, or mmCIF structures were found in this upload.")
    if len(found) > 100:
        raise PortalError("An MD campaign may contain at most 100 structures.")
    return found


def write_md_config(root: Path, filename: str, form) -> None:
    """Write the same per-folder PyMACS configuration created by its Q&A."""
    if Path(filename).suffix.lower() not in STRUCTURE_SUFFIXES:
        raise PortalError("MD input must be a .pdb, .cif, or .mmcif file.")
    system_type = form.get("system_type", "protein")
    if system_type not in {"protein", "peptide", "ligand", "protac", "biological"}:
        raise PortalError("Choose a supported system type.")
    system_name = form.get("system_name", "").strip()
    if not system_name:
        system_name = Path(filename).stem
    if len(system_name) > 80 or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9 _.()'/-]{0,79}", system_name):
        raise PortalError("System name must be 1–80 letters, numbers, spaces, or . _ ( ) ' / -.")
    ligand = form.get("ligand", "").strip().upper()
    ligand_name = form.get("ligand_name", "").strip()
    component_instances = form.get("component_instances", "").strip().upper()
    if system_type in {"ligand", "protac"} and not re.fullmatch(r"[A-Z0-9]{1,4}", ligand):
        raise PortalError("Ligand and PROTAC MD need a 1–4 character residue code.")
    try:
        length_ns, box_distance = float(form.get("length_ns", "100")), float(form.get("box_distance_nm", "1.0"))
        threads = int(form.get("threads", "16"))
        pocket_cutoff, contact_cutoff = float(form.get("pocket_cutoff", "5.0")), float(form.get("contact_cutoff", "4.0"))
        min_contact_fraction = float(form.get("min_contact_fraction", "0.10"))
        interface_contact_cutoff = float(form.get("interface_contact_cutoff", "4.0"))
        interface_min_contact_fraction = float(form.get("interface_min_contact_fraction", "0.10"))
        interface_frame_step, interface_max_edges = int(form.get("interface_frame_step", "1")), int(form.get("interface_max_edges", "100"))
    except ValueError as exc:
        raise PortalError("Length, box distance, and threads must be numeric.") from exc
    if not 0 < length_ns <= 10_000 or not 0 < box_distance <= 10 or not 1 <= threads <= 64:
        raise PortalError("Use 0–10,000 ns, 0–10 nm box distance, and 1–64 threads.")
    if not 0 < pocket_cutoff <= 20 or not 0 < contact_cutoff <= 20 or not 0 <= min_contact_fraction <= 1 or not 0 < interface_contact_cutoff <= 20 or not 0 <= interface_min_contact_fraction <= 1 or not 1 <= interface_frame_step <= 1000 or not 1 <= interface_max_edges <= 10000:
        raise PortalError("Analysis cutoffs, fractions, frame step, or maximum edges are out of range.")
    box = form.get("box_type", "dodecahedron")
    if box not in {"dodecahedron", "cubic", "octahedron", "triclinic"}:
        raise PortalError("Unsupported box type.")
    setup_mode = "ligand" if system_type in {"ligand", "protac"} else "protein"
    chain_map = form.get("chain_map", "").strip() or default_chain_map(root / filename, system_name)
    if len(ligand_name) > 80:
        raise PortalError("Ligand name must be at most 80 characters.")
    config = {
        "schema_version": 1, "input": filename,
        "metadata": {"system_name": system_name, "ligand_name": ligand_name or None,
                     "replicate_id": form.get("replicate_id", "").strip() or None,
                     "purpose_notes": form.get("purpose_notes", "").strip() or None,
                     "portal_assumptions": {"temperature_k": 300, "pressure_bar": 1.0,
                         "force_field": "CHARMM36", "water_model": "PyMACS/GROMACS default",
                         "ion_settings": "PyMACS setup default; no user-selectable salt concentration",
                         "seed": "GROMACS default/random; PyMACS does not expose a seed"}},
        "setup": {"mode": setup_mode, "chain_map": chain_map,
                  "keep_chains": form.get("keep_chains", "").strip() or None, "ligand": ligand or None,
                  "cofactors": form.get("cofactors", "").strip().upper() or None,
                  "component_instances": component_instances or None,
                  "remove_input_waters": form.get("remove_waters") == "on",
                  "remove_input_ions": form.get("remove_ions") == "on",
                  "remove_input_solvents": form.get("remove_solvents") == "on",
                  "box_type": box, "box_distance_nm": box_distance, "cgenff_backend": "local"},
        "simulation": {"mode": system_type, "length_ns": length_ns, "threads": threads, "compute": "GPU"},
        "analysis": {"enabled": form.get("analysis") == "on", "figurebook": form.get("figurebook") == "on",
                     "threads": 8, "chunk_size": 1000,
                     "pocket_cutoff": pocket_cutoff, "contact_cutoff": contact_cutoff,
                     "min_contact_fraction": min_contact_fraction,
                     "interface_contact_cutoff": interface_contact_cutoff,
                     "interface_min_contact_fraction": interface_min_contact_fraction,
                     "interface_frame_step": interface_frame_step, "interface_max_edges": interface_max_edges},
    }
    if system_type in {"ligand", "protac"}:
        config["analysis"]["compound_name"] = ligand_name or ligand
    if system_type == "protac":
        config["analysis"]["protac_frame_step"] = 10
        target_chains = [chain.strip() for chain in form.get("protac_target_chains", "").split(",") if chain.strip()]
        try:
            target_chain_names = json.loads(form.get("protac_target_chain_names", "{}"))
        except (TypeError, ValueError):
            target_chain_names = {}
        if not isinstance(target_chain_names, dict):
            target_chain_names = {}
        target_chain_names = {str(chain): str(name).strip() for chain, name in target_chain_names.items()
                              if str(chain) in target_chains and str(name).strip()}
        config["metadata"]["protac"] = {
            "ligase_chain": form.get("protac_ligase_chain", "").strip() or None,
            "ligase_name": form.get("protac_ligase_name", "").strip() or None,
            "target_chains": target_chains,
            "target_assembly_name": form.get("protac_target_name", "").strip() or None,
            "target_chain_names": target_chain_names,
        }
        config["analysis"]["target_assembly_name"] = config["metadata"]["protac"]["target_assembly_name"]
    root.joinpath("pymacs_run.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")


def build_md_campaign(settings: Settings, root: Path, selected: list[dict[str, object]]) -> dict[str, object]:
    """Create one isolated PyMACS project per selected uploaded structure."""
    if not selected:
        raise PortalError("Choose at least one structure for the MD campaign.")
    projects = []
    systems_root = root / "systems"
    for position, item in enumerate(selected, start=1):
        source_value = str(item.get("path", ""))
        source = (root / _safe_relative(source_value)).resolve()
        if root.resolve() not in source.parents or not source.is_file() or source.suffix.lower() not in STRUCTURE_SUFFIXES:
            raise PortalError("One selected MD structure is no longer present in the uploaded ZIP.")
        requested_name = str(item.get("system_name", "")).strip() or source.stem
        stem = re.sub(r"[^A-Za-z0-9_.-]+", "-", requested_name).strip(".-") or "structure"
        project_id = f"system-{position:03d}-{stem[:48]}"
        project = systems_root / project_id
        project.mkdir(parents=True, exist_ok=False)
        structure_name = secure_filename(source.name)
        shutil.copy2(source, project / structure_name)
        write_md_config(project, structure_name, item)
        if item.get("system_type") in {"ligand", "protac"}:
            prepare_ligand_cgenff(settings, project, str(item.get("ligand", "")))
        try:
            target_chain_names = json.loads(str(item.get("protac_target_chain_names", "{}")))
        except (TypeError, ValueError):
            target_chain_names = {}
        if not isinstance(target_chain_names, dict):
            target_chain_names = {}
        projects.append({"id": project_id, "input": structure_name,
                         "system_name": requested_name, "ligand_name": str(item.get("ligand_name", "")).strip() or None,
                         "replicate_id": str(item.get("replicate_id", "")).strip() or None,
                         "analysis": item.get("analysis") == "on",
                         "protac": {"ligase_chain": str(item.get("protac_ligase_chain", "")).strip() or None,
                                    "target_chains": [chain.strip() for chain in str(item.get("protac_target_chains", "")).split(",") if chain.strip()],
                                    "target_assembly_name": str(item.get("protac_target_name", "")).strip() or None,
                                    "target_chain_names": target_chain_names}
                         if item.get("system_type") == "protac" else None})
    manifest = {"schema_version": 1, "projects": projects}
    (root / "md_campaign.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def prepare_ligand_cgenff(settings: Settings, project: Path, ligand_code: str) -> None:
    """Create the CGenFF pair PyMACS reuses when the project reaches Triton."""
    code = ligand_code.upper()
    structure = next((item for item in project.iterdir() if item.suffix.lower() == ".pdb"), None)
    if structure is None:
        raise PortalError("Kyle-side ligand preparation currently requires a PDB structure.")
    ligand_pdb = project / f"{code}_extracted.pdb"
    atoms = [line for line in structure.read_text(encoding="utf-8", errors="replace").splitlines()
             if line.startswith("HETATM") and line[17:20].strip().upper() == code]
    if not atoms:
        raise PortalError(f"Could not find HETATM residue {code} in {structure.name} for CGenFF preparation.")
    ligand_pdb.write_text("\n".join(atoms) + "\nEND\n", encoding="utf-8")
    mol2 = project / f"{code}.cgenff.mol2"
    completed = subprocess.run(["obabel", str(ligand_pdb), "-O", str(mol2), "-h", "--title", code],
                               text=True, capture_output=True, timeout=120)
    if completed.returncode or not mol2.is_file() or not mol2.stat().st_size:
        raise PortalError("Open Babel could not make the ligand MOL2: " + (completed.stderr or completed.stdout)[-500:])
    try:
        token = settings.cgenff_token_file.read_text(encoding="utf-8").strip()
        payload = json.dumps({"ligand_code": code, "mol2_b64": base64.b64encode(mol2.read_bytes()).decode("ascii")}).encode()
        request_obj = Request(settings.cgenff_url, data=payload, method="POST",
                              headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json"})
        with urlopen(request_obj, timeout=180) as response:
            reply = json.loads(response.read().decode("utf-8"))
        stream = base64.b64decode(reply["str_b64"], validate=True)
    except (OSError, KeyError, ValueError, URLError, HTTPError) as exc:
        raise PortalError(f"Kyle CGenFF preparation failed: {exc}") from exc
    if not stream.strip() or b"RESI" not in stream.upper():
        raise PortalError("Kyle CGenFF returned an invalid stream file.")
    project.joinpath(f"{code}.str").write_bytes(stream)


def run(command: list[str], *, cwd: Path | None = None, timeout: int = 120) -> str:
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True, timeout=timeout)
    if completed.returncode:
        detail = (completed.stderr or completed.stdout or "command failed").strip()
        raise PortalError(detail[-2000:])
    return completed.stdout.strip()


def shell_quote(value: str) -> str:
    # Remote commands are fixed except generated server-side path components.
    return "'" + value.replace("'", "'\\\"'\\\"'") + "'"


def submitter_slug(name: str) -> str:
    cleaned = "-".join(name.strip().lower().split())
    if not SUBMITTER_NAME.fullmatch(name.strip()):
        raise PortalError("Enter your first name using letters, spaces, or hyphens only.")
    return re.sub(r"[^a-z-]", "", cleaned)


def record_submission(settings: Settings, submitter: str, kind: str, name: str, result: dict[str, str]) -> None:
    record = {"submitted_at": dt.datetime.now(dt.timezone.utc).isoformat(), "submitter": submitter,
              "kind": kind, "project": name, **result}
    with (settings.stage_root / "submissions.jsonl").open("a", encoding="utf-8") as ledger:
        ledger.write(json.dumps(record) + "\n")


def submission_records(settings: Settings) -> list[dict[str, object]]:
    """Return the portal's durable submission ledger, newest campaign first."""
    path = settings.stage_root / "submissions.jsonl"
    records: list[dict[str, object]] = []
    try:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    item = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if isinstance(item, dict) and item.get("remote") and item.get("kind") in {"docking", "md"}:
                    records.append(item)
    except OSError:
        pass
    return sorted(records, key=lambda item: str(item.get("submitted_at", "")), reverse=True)[:100]


def published_md_results_root(settings: Settings) -> Path:
    """Kyle-local, browser-friendly copies of finished MD result bundles."""
    return settings.stage_root / "published-md-results"


def md_results_index_path(settings: Settings) -> Path:
    return settings.stage_root / "md_results_index.json"


def read_md_results_index(settings: Settings) -> dict[str, object]:
    try:
        data = json.loads(md_results_index_path(settings).read_text(encoding="utf-8"))
        return data if isinstance(data, dict) and isinstance(data.get("systems"), list) else {"systems": []}
    except (OSError, json.JSONDecodeError):
        return {"systems": []}


def safe_result_relative(value: str) -> PurePosixPath | None:
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or any(part in {"", ".", ".."} for part in path.parts):
        return None
    return path


def remote_md_result_inventory(settings: Settings, campaigns: list[dict[str, str]]) -> list[dict[str, object]]:
    """List final, shareable PyMACS artifacts without scanning from Kyle.

    Analysis completion is represented by a figurebook.  Only then are the
    result files considered stable enough to publish or archive.
    """
    if not campaigns:
        return []
    payload = base64.b64encode(json.dumps(campaigns).encode()).decode()
    script = r'''import base64, json
from pathlib import Path

campaigns = json.loads(base64.b64decode(PAYLOAD))
suffixes = {".csv", ".html", ".pdf", ".png", ".svg", ".mol2", ".sdf", ".mol"}
trajectory_names = {"Final_Trajectory.pdb", "Final_Trajectory.xtc", "binding_pocket_only.pdb", "binding_pocket_only.xtc"}
limit = 2 * 1024 * 1024 * 1024

def projects(root):
    manifest = root / "md_campaign.json"
    try:
        entries = json.loads(manifest.read_text()).get("projects", [])
        found = [root / "systems" / str(entry.get("id", "")) for entry in entries]
        return [item for item in found if item.is_dir()]
    except Exception:
        return [root] if root.is_dir() else []

items = []
for campaign in campaigns:
    root = Path(campaign["remote"])
    for project in projects(root):
        config = {}
        try:
            config = json.loads((project / "pymacs_run.json").read_text())
        except Exception:
            pass
        metadata = config.get("metadata", {}) if isinstance(config, dict) else {}
        files = []
        seen = set()
        for candidate in project.iterdir() if project.is_dir() else []:
            if candidate.is_file() and (candidate.name in trajectory_names or candidate.suffix.lower() in suffixes) and candidate.stat().st_size <= limit:
                seen.add(candidate)
        for output in (project / "Analysis_Results",):
            if output.is_dir():
                for candidate in output.rglob("*"):
                    if candidate.is_file() and (candidate.suffix.lower() in suffixes or candidate.name in trajectory_names) and candidate.stat().st_size <= limit:
                        seen.add(candidate)
        for candidate in sorted(seen):
            try:
                rel = candidate.relative_to(project).as_posix()
                files.append({"relative": rel, "size": candidate.stat().st_size})
            except (OSError, ValueError):
                pass
        ready = any(entry["relative"].endswith(("MD_ANALYSIS_FIGUREBOOK.pdf", "PROTAC_MD_ANALYSIS_FIGUREBOOK.pdf")) for entry in files)
        if ready:
            items.append({"campaign_key": campaign["key"], "project": project.name,
                          "display_name": str(metadata.get("system_name") or project.name),
                          "ligand_name": str(metadata.get("ligand_name") or ""), "files": files})
print(json.dumps(items))
'''.replace("PAYLOAD", repr(payload), 1)
    encoded = base64.b64encode(script.encode()).decode()
    try:
        output = run([settings.ssh_bin, "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                      settings.triton_host, f"printf %s {shell_quote(encoded)} | base64 -d | python3"], timeout=90)
        parsed = json.loads(output)
        return parsed if isinstance(parsed, list) else []
    except (PortalError, ValueError, json.JSONDecodeError, subprocess.TimeoutExpired):
        return []


def remote_legacy_md_inventory(settings: Settings) -> list[dict[str, object]]:
    """Discover completed legacy MD runs using their GROMACS completion marker.

    These folders predate portal submission.  They are intentionally marked as
    external and are not auto-archived: a large raw trajectory must never be
    copied merely because someone opened a dashboard page.
    """
    collections = list(settings.legacy_md_collections)
    if not collections:
        return []
    payload = base64.b64encode(json.dumps(collections).encode()).decode()
    script = r'''import base64, json
from pathlib import Path

collections = json.loads(base64.b64decode(PAYLOAD))
portable_names = {"selected_input_curated.pdb", "selected_input_source.pdb", "md_0_1.gro", "pymacs_run.json"}
limit = 2 * 1024 * 1024 * 1024
items = []
for collection_index, collection in enumerate(collections):
    root = Path(collection["root"])
    try:
        projects = sorted(item for item in root.iterdir() if item.is_dir() and not item.name.startswith("."))
    except OSError:
        continue
    for project in projects:
        run_log = project / "md_0_1.log"
        try:
            completed = "Finished mdrun" in run_log.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            completed = False
        if not completed:
            continue
        files = []
        for candidate in project.iterdir():
            try:
                if candidate.name == "selected_input_source.pdb" and (project / "selected_input_curated.pdb").is_file():
                    continue
                if candidate.is_file() and (candidate.name in portable_names or candidate.suffix.lower() in {".mol2", ".sdf", ".mol"}) and candidate.stat().st_size <= limit:
                    files.append({"relative": candidate.name, "size": candidate.stat().st_size})
            except OSError:
                pass
        # If a legacy run is analysed later, expose its portable outputs too.
        output = project / "Analysis_Results"
        if output.is_dir():
            for candidate in output.rglob("*"):
                try:
                    if candidate.is_file() and (candidate.suffix.lower() in {".csv", ".html", ".pdf", ".png", ".svg"} or candidate.name in {"Final_Trajectory.pdb", "Final_Trajectory.xtc", "binding_pocket_only.pdb", "binding_pocket_only.xtc"}) and candidate.stat().st_size <= limit:
                        files.append({"relative": candidate.relative_to(project).as_posix(), "size": candidate.stat().st_size})
                except OSError:
                    pass
        metadata = {}
        try:
            metadata = json.loads((project / "pymacs_run.json").read_text()).get("metadata", {})
        except Exception:
            pass
        ready = any(file["relative"].endswith(("MD_ANALYSIS_FIGUREBOOK.pdf", "PROTAC_MD_ANALYSIS_FIGUREBOOK.pdf")) for file in files)
        items.append({
            "collection_key": str(collection_index), "project": project.name, "project_remote": str(project),
            "display_name": str(metadata.get("system_name") or f"{collection.get('protein_name', '')} · {project.name}"),
            "ligand_name": str(metadata.get("ligand_name") or collection.get("ligand_name") or "LIG"),
            "files": files, "analysis_ready": ready,
        })
print(json.dumps(items))
'''.replace("PAYLOAD", repr(payload), 1)
    encoded = base64.b64encode(script.encode()).decode()
    try:
        output = run([settings.ssh_bin, "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                      settings.triton_host, f"printf %s {shell_quote(encoded)} | base64 -d | python3"], timeout=90)
        parsed = json.loads(output)
        return parsed if isinstance(parsed, list) else []
    except (PortalError, ValueError, json.JSONDecodeError, subprocess.TimeoutExpired):
        return []


def start_full_md_archive(settings: Settings, campaign: dict[str, str]) -> str:
    """Start a complete, non-destructive Cheminfo mirror after analysis is final."""
    submitter = submitter_slug(campaign["submitter"])
    remote = campaign["remote"]
    archive = f"{settings.md_results_root}/{submitter}/{Path(remote).name}"
    lock = archive + "/.portal-full-archive.lock"
    complete = archive + "/.portal-full-archive.complete"
    archive_script = (
        f"if rsync -aH --partial --append-verify {shell_quote(remote + '/')} {shell_quote(archive + '/')}; then "
        f"touch {shell_quote(complete)}; fi; rmdir {shell_quote(lock)}"
    )
    command = (
        f"mkdir -p {shell_quote(archive)}; "
        f"if [ ! -f {shell_quote(complete)} ] && mkdir {shell_quote(lock)} 2>/dev/null; then "
        f"nohup bash -lc {shell_quote(archive_script)} >/dev/null 2>&1 </dev/null & fi"
    )
    try:
        run([settings.ssh_bin, "-o", "BatchMode=yes", settings.triton_host, command], timeout=30)
    except (PortalError, subprocess.TimeoutExpired):
        pass
    return archive


def summarize_published_md_results(destination: Path, files: list[dict[str, object]]) -> dict[str, object]:
    """Extract small, presentation-safe facts from Kyle-local analysis outputs.

    This intentionally reads only the published copy.  The results page must
    remain independent from Triton and never open a trajectory or scan a
    remote project while a user is waiting for the page.
    """
    summary: dict[str, object] = {"qx": None}
    q_csv = next((entry for entry in files if str(entry.get("relative", "")).endswith("native_contact_q.csv")), None)
    if not isinstance(q_csv, dict):
        return summary
    relative = safe_result_relative(str(q_csv.get("relative", "")))
    if relative is None:
        return summary
    path = destination.joinpath(*relative.parts)
    try:
        values: list[float] = []
        cutoff = None
        pairs = None
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                try:
                    values.append(float(str(row.get("Q_native_contact_fraction", ""))))
                except (TypeError, ValueError):
                    continue
                if cutoff is None:
                    try:
                        cutoff = float(str(row.get("ContactCutoff_A", "")))
                    except (TypeError, ValueError):
                        pass
                if pairs is None:
                    try:
                        pairs = int(float(str(row.get("NativeHeavyAtomContactPairs", ""))))
                    except (TypeError, ValueError):
                        pass
        if values:
            summary["qx"] = {
                "mean": round(sum(values) / len(values), 4),
                "final": round(values[-1], 4),
                "frames": len(values),
                "cutoff_angstrom": cutoff,
                "reference_pairs": pairs,
            }
    except (OSError, UnicodeError, csv.Error):
        pass
    return summary


def refresh_md_results(settings: Settings) -> None:
    """Publish final analysis bundles on Kyle and trigger the full Cheminfo archive.

    The browser reads only the local JSON index and Kyle-local copies.  It
    never waits on Triton or opens a remote trajectory during a page load.
    """
    records = [
        {"key": str(index), "remote": str(record["remote"]), "submitter": str(record.get("submitter", "")),
         "campaign": str(record.get("project", "")), "external": False}
        for index, record in enumerate(submission_records(settings))
        if record.get("kind") == "md" and record.get("remote") and record.get("submitter")
    ]
    portal_inventory = remote_md_result_inventory(settings, records)
    legacy_records = [
        {"key": f"legacy-{index}", "remote": collection["root"], "submitter": collection["submitter"],
         "campaign": collection["campaign"], "external": True, "protein_name": collection["protein_name"]}
        for index, collection in enumerate(settings.legacy_md_collections)
    ]
    legacy_inventory = remote_legacy_md_inventory(settings)
    # A remote listing error is represented as an empty list.  Never turn a
    # transient Triton/SSH outage into an empty Kyle results library: keep the
    # last successfully published index until the next healthy refresh.
    if (records and not portal_inventory) or (settings.legacy_md_collections and not legacy_inventory):
        return
    inventory = list(portal_inventory)
    for item in legacy_inventory:
        collection_key = str(item.pop("collection_key", ""))
        item["campaign_key"] = f"legacy-{collection_key}"
        item["external"] = True
        inventory.append(item)
    records.extend(legacy_records)
    record_by_key = {item["key"]: item for item in records}
    root = published_md_results_root(settings)
    root.mkdir(parents=True, exist_ok=True)
    systems: list[dict[str, object]] = []
    archived: set[str] = set()
    for item in inventory:
        campaign = record_by_key.get(str(item.get("campaign_key", "")))
        if campaign is None:
            continue
        external = bool(campaign.get("external") or item.get("external"))
        archive = "Not yet requested for this external legacy run."
        if not external and campaign["key"] not in archived:
            archive = start_full_md_archive(settings, campaign)
            archived.add(campaign["key"])
        source_root = str(item.get("project_remote") or campaign["remote"])
        destination_base = root / submitter_slug(campaign["submitter"]) / Path(campaign["remote"]).name
        if external:
            destination_base = root / "external" / submitter_slug(campaign["submitter"]) / Path(campaign["remote"]).name
        destination = destination_base / str(item["project"])
        copied = []
        for source in item.get("files", []):
            if not isinstance(source, dict):
                continue
            relative = safe_result_relative(str(source.get("relative", "")))
            if relative is None or int(source.get("size", 0)) > RESULT_FILE_LIMIT:
                continue
            target = destination.joinpath(*relative.parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists() or target.stat().st_size != int(source.get("size", 0)):
                partial = target.with_name(target.name + ".partial")
                try:
                    run([settings.scp_bin, "-p", f"{settings.triton_host}:{source_root}/{relative.as_posix()}", str(partial)], timeout=900)
                    partial.replace(target)
                except (PortalError, OSError, subprocess.TimeoutExpired):
                    partial.unlink(missing_ok=True)
                    continue
            copied.append({"relative": relative.as_posix(), "size": target.stat().st_size})
        if copied:
            presentation = summarize_published_md_results(destination, copied)
            systems.append({"submitter": campaign["submitter"], "campaign": campaign["campaign"],
                            "system": item.get("display_name") or item["project"], "project": item["project"],
                            "ligand_name": item.get("ligand_name", ""), "archive": archive,
                            "published_at": dt.datetime.now(dt.timezone.utc).isoformat(), "files": copied,
                            "directory": destination.relative_to(root).as_posix(), "external": external,
                            "analysis_ready": bool(item.get("analysis_ready", not external)),
                            "source_path": source_root, **presentation})
    payload = {"generated_at": dt.datetime.now(dt.timezone.utc).isoformat(), "systems": systems}
    target = md_results_index_path(settings)
    temporary = target.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload), encoding="utf-8")
    temporary.replace(target)


def restart_explicit_legacy_timeouts(settings: Settings) -> None:
    """Resume only legacy projects that were explicitly opted in on Triton.

    The remote watcher requires a valid GROMACS checkpoint and LSF's timeout
    signal.  It has no authority over ordinary portal or legacy folders unless
    they contain a reviewed ``.portal-auto-resume.json`` marker.
    """
    if not settings.legacy_md_auto_restart:
        return
    watcher = f"{settings.pymacs_root}/portal_legacy_resume_watch.py"
    for collection in settings.legacy_md_collections:
        root = str(collection.get("root", "")).rstrip("/")
        if not root.startswith("/projectnb/"):
            continue
        command = f"python3 {shell_quote(watcher)} --root {shell_quote(root)} --execute"
        try:
            run([settings.ssh_bin, "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                 settings.triton_host, command], timeout=120)
        except (PortalError, subprocess.TimeoutExpired):
            # A monitor outage must never interrupt result publication or a
            # browser request.  The next five-minute cycle will retry safely.
            pass


def start_md_results_worker(settings: Settings) -> None:
    """Keep publication out of the request path; completed results refresh every 5 min."""
    def worker() -> None:
        while True:
            try:
                refresh_md_results(settings)
            except Exception:
                pass
            threading.Event().wait(300)
    threading.Thread(target=worker, name="md-results-publisher", daemon=True).start()


def start_legacy_resume_worker(settings: Settings) -> None:
    """Poll explicit legacy timeout recovery independently every five minutes."""
    if not settings.legacy_md_auto_restart:
        return
    def worker() -> None:
        while True:
            try:
                restart_explicit_legacy_timeouts(settings)
            except Exception:
                pass
            threading.Event().wait(300)
    threading.Thread(target=worker, name="legacy-md-timeout-recovery", daemon=True).start()


def remote_monitor(settings: Settings, host: str, campaigns: list[dict[str, object]]) -> dict[str, object]:
    """Ask one scheduler host for progress on its already-staged campaigns.

    The remote program receives only server-recorded paths/job IDs, emits JSON,
    and never invokes bkill, bsub, or modifies campaign data.
    """
    payload = base64.b64encode(json.dumps(campaigns).encode()).decode()
    script = '''import base64, json, re, subprocess
from pathlib import Path

campaigns = json.loads(base64.b64decode(PAYLOAD))
states = {{}}
try:
    result = subprocess.run(["bjobs", "-a", "-noheader", "-o", "jobid stat job_name"], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=15)
    for line in result.stdout.splitlines():
        fields = line.split(maxsplit=2)
        if len(fields) >= 2:
            states[fields[0]] = {{"state": fields[1], "name": fields[2] if len(fields) == 3 else ""}}
except Exception:
    pass

def job_state(job_id):
    return states.get(str(job_id), {{}}).get("state", "UNKNOWN")

def log_progress(project, target_ns):
    logs = sorted(project.rglob("md_0_1.log"), key=lambda p: p.stat().st_mtime if p.exists() else 0)
    if not logs:
        return None
    try:
        text = logs[-1].read_text(errors="replace")[-50000:]
    except OSError:
        return None
    values = re.findall(r"^\\s*([0-9]+)\\s+([0-9]+(?:\\.[0-9]+)?)\\s*$", text, re.M)
    if not values:
        return None
    step, time_ps = values[-1]
    percent = min(99.9, max(0.0, float(time_ps) / (float(target_ns) * 1000) * 100)) if target_ns else 0.0
    return {{"step": int(step), "time_ps": float(time_ps), "percent": round(percent, 1), "log": str(logs[-1])}}

def md_campaign(campaign):
    root = Path(campaign["remote"])
    ids = [value.strip() for value in str(campaign.get("md_jobs", "")).split(",") if value.strip()]
    try:
        manifest = json.loads((root / "md_campaign.json").read_text())
        projects = [root / "systems" / str(item["id"]) for item in manifest.get("projects", [])]
    except Exception:
        projects = [root]
    systems = []
    for index, project in enumerate(projects):
        try:
            config = json.loads((project / "pymacs_run.json").read_text())
            target_ns = float(config.get("simulation", {{}}).get("length_ns", 0))
        except Exception:
            target_ns = 0.0
        job_id = ids[index] if index < len(ids) else ""
        state = job_state(job_id)
        progress = log_progress(project, target_ns)
        if state == "DONE":
            progress = {{"percent": 100.0}}
        systems.append({{"name": project.name, "job_id": job_id, "state": state, "target_ns": target_ns, "progress": progress}})
    done = sum(1 for item in systems if item["state"] == "DONE")
    return {{"systems": systems, "done": done, "total": len(systems)}}

def docking_campaign(campaign):
    root = Path(campaign["remote"])
    results = Path(str(campaign.get("results", "")))
    try:
        manifest = json.loads((root / "campaign.json").read_text())
        batches = manifest.get("batches", [])
    except Exception:
        manifest = {{}}
        batches = []
    if not batches:
        return {"batches": [], "done": 0, "total": 0, "progress": 0,
                "message": "Campaign files are unavailable (likely an archived or deleted test)."}
    options = manifest.get("options", {})

    def input_count(folder):
        inputs = list(folder.glob("*.csv")) + list(folder.glob("*.sdf")) + list(folder.glob("*.smi")) + list(folder.glob("*.smiles"))
        if not inputs:
            return 0
        try:
            source = inputs[0]
            if source.suffix.lower() == ".csv":
                return max(0, len(source.read_text(errors="replace").splitlines()) - 1)
            if source.suffix.lower() == ".sdf":
                return source.read_text(errors="replace").count("$$$$")
            return sum(1 for line in source.read_text(errors="replace").splitlines() if line.strip())
        except OSError:
            return 0

    def center_count():
        # The launch workflow writes the target centers at campaign root.  Do
        # not recursively walk every generated ligand/result file just to
        # locate this small, stable input file.
        centers = sorted(root.glob("vina_centers*.csv"))
        if not centers:
            return 0
        try:
            return max(0, len(centers[0].read_text(errors="replace").splitlines()) - 1)
        except OSError:
            return 0

    state_multiplier = (int(options.get("max_protomers", 4)) * int(options.get("max_tautomers", 4))) if options.get("enumerate_protomers") else 1
    pdbqt_per_input = int(options.get("confgen_poses", 64)) * state_multiplier
    docking_targets = center_count()
    job_ids = {}
    try:
        for line in (root / "job_ids.jsonl").read_text().splitlines():
            entry = json.loads(line)
            if isinstance(entry, dict) and entry.get("batch"):
                job_ids[str(entry["batch"])] = entry
    except Exception:
        pass

    def vina_log_info(folder):
        logs = sorted((folder / "logs").glob("vina-*.out"), key=lambda p: p.stat().st_mtime if p.exists() else 0)
        if not logs:
            return "", None
        try:
            text = logs[-1].read_text(errors="replace")[-100000:]
        except OSError:
            return "", None
        matches = re.findall(r"\\]\s+([0-9]+)%\s+\\(([0-9]+)/([0-9]+)\\)", text)
        progress = None
        if matches:
            percent, done, total = matches[-1]
            progress = {{"percent": int(percent), "done": int(done), "total": int(total)}}
        if "Exited with exit code 0" in text:
            return "DONE", progress or {{"percent": 100, "done": 0, "total": 0}}
        if "Exited with exit code" in text:
            return "EXIT", progress
        return "", progress

    def confgen_pdbqt_info(folder):
        """Read the runner's own generated-PDBQT tally without walking NFS."""
        logs = sorted((folder / "logs").glob("conf-*.out"), key=lambda p: p.stat().st_mtime if p.exists() else 0)
        if not logs:
            return 0
        try:
            text = logs[-1].read_text(errors="replace")
        except OSError:
            return 0
        completed = re.findall(r"Done!\\s+([0-9]+)\\s+poses written", text)
        if completed:
            return int(completed[-1])
        # An in-progress runner emits one line for every successfully written
        # state; summing those reported conformers is an inexpensive live count.
        return sum(int(value) for value in re.findall(r"\\[STATE\\].*?:\\s*([0-9]+)/[0-9]+\\s+conformers generated", text))

    items = []
    for batch in batches:
        batch = str(batch)
        folder = root / "batches" / batch
        conf_done = (folder / "confgen_pdbqt_dir.txt").is_file()
        ligand_inputs = input_count(folder)
        # Expensive directory counts are refreshed by Kyle's background cache,
        # never while somebody is loading the board.
        pdbqt_count = 0
        pdbqt_maximum = ligand_inputs * pdbqt_per_input
        vina_done = False
        ids = job_ids.get(batch, {})
        job_batch = batch.replace("_", "-")
        # The handoff file is written only after ConfGen has successfully
        # produced its PDBQT directory.  It remains trustworthy after LSF has
        # retired the original job record, so prefer it over an aged-out ID.
        conf_state = "DONE" if conf_done else (job_state(ids.get("confgen")) if ids.get("confgen") else next((row["state"] for row in states.values() if row.get("name") == f"portal-conf-{job_batch}"), "UNKNOWN"))
        vina_state = job_state(ids.get("vina")) if ids.get("vina") else next((row["state"] for row in states.values() if row.get("name") == f"portal-vina-{job_batch}"), "DONE" if vina_done else "UNKNOWN")
        log_state, vina_progress = vina_log_info(folder)
        if vina_progress and vina_progress.get("total"):
            pass
        elif pdbqt_count and docking_targets:
            vina_progress = {"percent": 0, "done": 0, "total": pdbqt_count * docking_targets}
        if vina_state == "UNKNOWN" and log_state:
            vina_state = log_state
        if vina_state == "DONE":
            vina_done = True
            vina_progress = vina_progress or {{"percent": 100, "done": 0, "total": 0}}
        items.append({{"name": batch, "confgen": conf_state, "vina": vina_state, "confgen_done": conf_done or conf_state == "DONE", "vina_done": vina_done, "input_count": ligand_inputs, "pdbqt": {{"produced": pdbqt_count, "maximum": pdbqt_maximum, "percent": round(pdbqt_count / pdbqt_maximum * 100, 1) if pdbqt_maximum else 0}}, "vina_progress": vina_progress}})
    pdbqt_produced = sum(item["pdbqt"]["produced"] for item in items)
    pdbqt_maximum = sum(item["pdbqt"]["maximum"] for item in items)
    return {{"batches": items, "done": sum(1 for item in items if item["vina_done"]), "total": len(items), "pdbqt": {{"produced": pdbqt_produced, "maximum": pdbqt_maximum, "percent": round(pdbqt_produced / pdbqt_maximum * 100, 1) if pdbqt_maximum else 0}}}}

output = {{}}
for campaign in campaigns:
    key = str(campaign["key"])
    try:
        output[key] = docking_campaign(campaign) if campaign["kind"] == "docking" else md_campaign(campaign)
    except Exception as error:
        output[key] = {{"error": str(error)}}
print(json.dumps(output))
'''
    script = script.replace("PAYLOAD", repr(payload), 1).replace("{{", "{").replace("}}", "}")
    encoded = base64.b64encode(script.encode()).decode()
    try:
        output = run([settings.ssh_bin, "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                      "-o", "ConnectionAttempts=1", host,
                      f"printf %s {shell_quote(encoded)} | base64 -d | python3"], timeout=90)
        parsed = json.loads(output)
        return parsed if isinstance(parsed, dict) else {}
    except (PortalError, ValueError, json.JSONDecodeError, subprocess.TimeoutExpired):
        return {}


def pdbqt_cache_path(settings: Settings) -> Path:
    return settings.stage_root / "campaign_progress.sqlite3"


def prepare_pdbqt_cache(settings: Settings) -> None:
    settings.stage_root.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(pdbqt_cache_path(settings)) as connection:
        connection.execute("""
            CREATE TABLE IF NOT EXISTS pdbqt_counts (
                remote TEXT NOT NULL,
                batch TEXT NOT NULL,
                produced INTEGER NOT NULL,
                vina_completed INTEGER NOT NULL DEFAULT 0,
                vina_total INTEGER NOT NULL DEFAULT 0,
                vina_successes INTEGER NOT NULL DEFAULT 0,
                vina_failures INTEGER NOT NULL DEFAULT 0,
                checked_at TEXT NOT NULL,
                PRIMARY KEY (remote, batch)
            )
        """)
        columns = {row[1] for row in connection.execute("PRAGMA table_info(pdbqt_counts)")}
        if "vina_completed" not in columns:
            connection.execute("ALTER TABLE pdbqt_counts ADD COLUMN vina_completed INTEGER NOT NULL DEFAULT 0")
        if "vina_total" not in columns:
            connection.execute("ALTER TABLE pdbqt_counts ADD COLUMN vina_total INTEGER NOT NULL DEFAULT 0")
        if "vina_successes" not in columns:
            connection.execute("ALTER TABLE pdbqt_counts ADD COLUMN vina_successes INTEGER NOT NULL DEFAULT 0")
        if "vina_failures" not in columns:
            connection.execute("ALTER TABLE pdbqt_counts ADD COLUMN vina_failures INTEGER NOT NULL DEFAULT 0")


def cached_pdbqt_counts(settings: Settings) -> dict[tuple[str, str], tuple[int, int, int, int, int, str]]:
    prepare_pdbqt_cache(settings)
    with sqlite3.connect(pdbqt_cache_path(settings)) as connection:
        rows = connection.execute("SELECT remote, batch, produced, vina_completed, vina_total, vina_successes, vina_failures, checked_at FROM pdbqt_counts").fetchall()
    return {(str(remote), str(batch)): (int(produced), int(vina_completed), int(vina_total), int(vina_successes), int(vina_failures), str(checked_at))
            for remote, batch, produced, vina_completed, vina_total, vina_successes, vina_failures, checked_at in rows}


def refresh_pdbqt_cache(settings: Settings) -> None:
    """Run the expensive Pegasus directory count outside request handling."""
    campaigns = [
        {"key": str(index), "remote": str(record.get("remote", "")), "results": str(record.get("results", ""))}
        for index, record in enumerate(submission_records(settings))
        if record.get("kind") == "docking" and record.get("remote")
    ]
    if not campaigns:
        return
    payload = base64.b64encode(json.dumps(campaigns).encode()).decode()
    script = '''import base64, json
from pathlib import Path

campaigns = json.loads(base64.b64decode(PAYLOAD))
output = {}
for campaign in campaigns:
    root = Path(campaign["remote"])
    counts = {}
    try:
        batches = json.loads((root / "campaign.json").read_text()).get("batches", [])
        with (root / "vina_centers.csv").open(newline="") as handle:
            center_rows = list(__import__("csv").DictReader(handle))
        receptor_dir = root / "Receptors"
        valid_receptors = sum(
            1 for row in center_rows
            if (receptor_dir / row["PDB_ID"]).exists() or any((receptor_dir / (row["PDB_ID"] + ext)).exists() for ext in (".pdbqt", ".pdb", ".mol2"))
        )
        for batch in batches:
            folder = root / "batches" / str(batch)
            produced = sum(
                1 for pdbqt_dir in folder.glob("*_Ligands_PDBQT_*")
                for path in pdbqt_dir.glob("*.pdbqt") if path.is_file()
            )
            shard = Path(campaign.get("results", "")) / "shards" / str(batch)
            def completed_tasks(path):
                try:
                    import csv
                    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
                        return {(row.get("ResolvedReceptor", ""), row.get("Ligand", ""), row.get("Seed", ""))
                                for row in csv.DictReader(handle)
                                if row.get("ResolvedReceptor") and row.get("Ligand")}
                except OSError:
                    return set()
            successes = completed_tasks(shard / "compact_results.csv")
            failures = completed_tasks(shard / "compact_failures.csv")
            failures -= successes
            completed = successes | failures
            counts[str(batch)] = {"produced": produced, "vina_completed": len(completed),
                                  "vina_total": produced * valid_receptors,
                                  "vina_successes": len(successes), "vina_failures": len(failures)}
    except Exception:
        pass
    output[str(campaign["key"])] = counts
print(json.dumps(output))
'''.replace("PAYLOAD", repr(payload), 1)
    encoded = base64.b64encode(script.encode()).decode()
    try:
        output = run([settings.ssh_bin, "-o", "BatchMode=yes", settings.pegasus_host,
                      f"printf %s {shell_quote(encoded)} | base64 -d | python3"], timeout=900)
        parsed = json.loads(output)
    except (PortalError, ValueError, json.JSONDecodeError, subprocess.TimeoutExpired):
        return
    if not isinstance(parsed, dict):
        return
    checked_at = dt.datetime.now(dt.timezone.utc).isoformat()
    remotes = {item["key"]: item["remote"] for item in campaigns}
    rows = [
        (remote, str(batch), int(values.get("produced", 0)), int(values.get("vina_completed", 0)),
         int(values.get("vina_total", 0)), int(values.get("vina_successes", 0)), int(values.get("vina_failures", 0)), checked_at)
        for key, batches in parsed.items() if isinstance(batches, dict)
        for remote in [remotes.get(str(key), "")]
        for batch, values in batches.items() if remote and isinstance(values, dict)
    ]
    if rows:
        prepare_pdbqt_cache(settings)
        with sqlite3.connect(pdbqt_cache_path(settings)) as connection:
            connection.executemany("""
                INSERT INTO pdbqt_counts (remote, batch, produced, vina_completed, vina_total, vina_successes, vina_failures, checked_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(remote, batch) DO UPDATE SET
                    produced = excluded.produced, checked_at = excluded.checked_at
                    , vina_completed = excluded.vina_completed, vina_total = excluded.vina_total,
                    vina_successes = excluded.vina_successes, vina_failures = excluded.vina_failures
            """, rows)


def start_pdbqt_cache_worker(settings: Settings) -> None:
    """Refresh persistent PDBQT counts every five minutes in one daemon thread."""
    def worker() -> None:
        while True:
            try:
                refresh_pdbqt_cache(settings)
            except Exception:
                pass
            threading.Event().wait(300)
    threading.Thread(target=worker, name="pdbqt-cache-refresh", daemon=True).start()


def monitored_submissions(settings: Settings) -> list[dict[str, object]]:
    """Merge the submission ledger with current read-only scheduler evidence."""
    records = submission_records(settings)
    requests: dict[str, list[dict[str, object]]] = {settings.pegasus_host: [], settings.triton_host: []}
    prepared = []
    for index, record in enumerate(records):
        result = {key: value for key, value in record.items() if key not in {"submitted_at", "submitter", "kind", "project"}}
        kind = str(record["kind"])
        key = str(index)
        remote = str(result.get("remote", ""))
        request_item = {"key": key, "kind": kind, "remote": remote, "results": str(result.get("results", "")),
                        "md_jobs": str(result.get("md jobs", result.get("md", "")))}
        requests[settings.pegasus_host if kind == "docking" else settings.triton_host].append(request_item)
        prepared.append({"key": key, "submitted_at": record.get("submitted_at"), "submitter": record.get("submitter"),
                         "kind": kind, "project": record.get("project"), "remote": remote})
    evidence = {}
    for host, entries in requests.items():
        if entries:
            evidence.update(remote_monitor(settings, host, entries))
    pdbqt_counts = cached_pdbqt_counts(settings)
    for item in prepared:
        item["monitor"] = evidence.get(item.pop("key"), {"error": "Status host did not respond."})
        monitor = item["monitor"]
        if item["kind"] == "docking" and isinstance(monitor, dict):
            produced_total = 0
            maximum_total = 0
            for batch in monitor.get("batches", []):
                pdbqt = batch.get("pdbqt", {})
                cached = pdbqt_counts.get((str(item["remote"]), str(batch.get("name", ""))))
                if cached:
                    pdbqt["produced"], vina_completed, vina_total, vina_successes, vina_failures, pdbqt["checked_at"] = cached
                    maximum = int(pdbqt.get("maximum", 0))
                    pdbqt["percent"] = round(int(pdbqt["produced"]) / maximum * 100, 1) if maximum else 0
                    if vina_total:
                        batch["vina_progress"] = {"done": vina_completed, "total": vina_total,
                                                  "successes": vina_successes, "failures": vina_failures,
                                                  "percent": round(vina_completed / vina_total * 100, 1)}
                produced_total += int(pdbqt.get("produced", 0))
                maximum_total += int(pdbqt.get("maximum", 0))
            if monitor.get("batches"):
                monitor["pdbqt"] = {"produced": produced_total, "maximum": maximum_total,
                                    "percent": round(produced_total / maximum_total * 100, 1) if maximum_total else 0}
    return prepared


def job_monitor_cache_path(settings: Settings) -> Path:
    """Local, atomically-written live-status snapshot for the campaign board."""
    return settings.stage_root / "job_monitor_snapshot.json"


def read_job_monitor_cache(settings: Settings) -> dict[str, object] | None:
    try:
        cached = json.loads(job_monitor_cache_path(settings).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(cached, dict) or not isinstance(cached.get("campaigns"), list):
        return None
    return cached


def refresh_job_monitor_cache(settings: Settings) -> None:
    """Poll remote schedulers out of band, never on the browser request path."""
    snapshot = {
        "campaigns": monitored_submissions(settings),
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    settings.stage_root.mkdir(parents=True, exist_ok=True)
    target = job_monitor_cache_path(settings)
    temporary = target.with_suffix(".tmp")
    temporary.write_text(json.dumps(snapshot), encoding="utf-8")
    temporary.replace(target)


def start_job_monitor_cache_worker(settings: Settings) -> None:
    """Keep the campaign board responsive with a 30-second status snapshot."""
    def worker() -> None:
        while True:
            try:
                refresh_job_monitor_cache(settings)
            except Exception:
                pass
            threading.Event().wait(30)
    threading.Thread(target=worker, name="job-monitor-cache-refresh", daemon=True).start()


def remember_submitter(settings: Settings, submitter: str, email: str = "") -> None:
    """Persist friendly names on Kyle for portal autocomplete and folder ownership."""
    path = settings.stage_root / "submitters.json"
    settings.stage_root.mkdir(parents=True, exist_ok=True)
    with path.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        handle.seek(0)
        try:
            names = json.load(handle)
        except json.JSONDecodeError:
            names = []
        users = [{"name": item, "email": ""} if isinstance(item, str) else item for item in names]
        canonical = submitter.strip()
        match = next((item for item in users if str(item.get("name", "")).lower() == canonical.lower()), None)
        if match is None:
            users.append({"name": canonical, "email": email.strip(), "created_at": dt.datetime.now(dt.timezone.utc).isoformat()})
            users.sort(key=lambda item: str(item["name"]).casefold())
        elif email.strip():
            match["email"] = email.strip()
        if users != names:
            handle.seek(0); handle.truncate()
            json.dump(users, handle)
        fcntl.flock(handle, fcntl.LOCK_UN)


def saved_submitters(settings: Settings) -> list[dict[str, str]]:
    try:
        names = json.loads((settings.stage_root / "submitters.json").read_text(encoding="utf-8"))
        users = [{"name": item, "email": ""} if isinstance(item, str) else item for item in names]
        return [{"name": str(item.get("name", "")), "email": str(item.get("email", ""))}
                for item in users if SUBMITTER_NAME.fullmatch(str(item.get("name", "")).strip())]
    except (OSError, json.JSONDecodeError):
        return []


def reserve_project_number(settings: Settings, remote_root: str, submitter: str, project: str, label: str) -> int:
    """Reserve a stable, human-readable run number across portal restarts."""
    counters = settings.stage_root / "campaign_counters.json"
    key = f"{remote_root}:{submitter_slug(submitter)}:{project}:{label}"
    settings.stage_root.mkdir(parents=True, exist_ok=True)
    with counters.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        handle.seek(0)
        try:
            data = json.load(handle)
        except json.JSONDecodeError:
            data = {}
        number = int(data.get(key, 0)) + 1
        data[key] = number
        handle.seek(0); handle.truncate(); json.dump(data, handle, indent=2)
        fcntl.flock(handle, fcntl.LOCK_UN)
    return number


def stage_project(settings: Settings, root: Path, host: str, remote_root: str, job_name: str, submitter: str, label: str) -> str:
    number = reserve_project_number(settings, remote_root, submitter, job_name, label)
    remote = f"{remote_root}/{submitter_slug(submitter)}/{job_name}/{label}-{number:03d}"
    run([settings.ssh_bin, "-o", "BatchMode=yes", host, "mkdir", "-p", remote])
    # Create a server-side archive instead of relying on scp's ambiguous
    # handling of a source directory's `/.` suffix.
    with tempfile.NamedTemporaryFile(prefix="lab-portal-stage-", suffix=".tar", delete=False) as handle:
        archive = Path(handle.name)
    try:
        with tarfile.open(archive, "w") as bundle:
            for item in root.iterdir():
                bundle.add(item, arcname=item.name, recursive=True)
        run([settings.scp_bin, str(archive), f"{host}:{remote}/project.tar"])
        run([settings.ssh_bin, "-o", "BatchMode=yes", host,
             f"cd {shell_quote(remote)} && tar -xf project.tar && rm project.tar"])
    finally:
        archive.unlink(missing_ok=True)
    return remote


def submit_docking(settings: Settings, root: Path, name: str, submitter: str) -> dict[str, str]:
    manifest_path = root / "campaign.json"
    if not manifest_path.is_file():
        raise PortalError("The docking campaign was not prepared.")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        batches = [str(item) for item in manifest["batches"]]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise PortalError("The docking campaign manifest is invalid.") from exc
    remote = stage_project(settings, root, settings.pegasus_host, settings.pegasus_remote_root, name, submitter, "campaign")
    results_root = f"{settings.pegasus_results_root}/{submitter_slug(submitter)}/{name}/{Path(remote).name}"
    runtime = f"CAMPAIGN_RESULTS_ROOT={shell_quote(results_root)}\n"
    runtime_encoded = base64.b64encode(runtime.encode()).decode()
    helpers = " ".join(shell_quote(item) for item in ("3_Complete_batch_docking.py", "ligand_manifest.py", "ligand_naming.py"))
    runtime_command = (
        f"set -e; mkdir -p {shell_quote(results_root)}; "
        f"for helper in {helpers}; do test -f {shell_quote(settings.pegasus_docking_root)}/$helper || {{ echo \"Missing reviewed compact helper: $helper\" >&2; exit 2; }}; cp -f {shell_quote(settings.pegasus_docking_root)}/$helper {shell_quote(remote)}/$helper; done; "
        f"printf %s {shell_quote(runtime_encoded)} | base64 -d > {shell_quote(remote)}/portal_runtime.env; chmod 600 {shell_quote(remote)}/portal_runtime.env"
    )
    run([settings.ssh_bin, "-o", "BatchMode=yes", settings.pegasus_host, runtime_command])
    # Submit every batch immediately. LSF performs the actual scheduling;
    # each Vina job remains dependent on its matching ConfGen completion.
    lines = ["#!/bin/bash", "set -euo pipefail", f"ROOT={shell_quote(remote)}", "for batch in " + " ".join(shell_quote(batch) for batch in batches) + "; do", "  vina_script=\"$ROOT/batches/$batch/run_vina.lsf\"", "  grep -Fq -- \"--compact-results\" \"$vina_script\" && grep -Fq -- \"--results-dir\" \"$vina_script\" && grep -Fq -- \"portal_runtime.env\" \"$vina_script\" || { echo \"Refusing non-compact Vina script: $vina_script\" >&2; exit 2; }", "  out=$(cd \"$ROOT/batches/$batch\" && bsub < run_confgen.lsf)", "  conf_id=$(printf '%s' \"$out\" | sed -n 's/.*<\\([0-9][0-9]*\\)>.*/\\1/p')", "  [ -n \"$conf_id\" ] || { echo \"Could not parse ConfGen submission: $out\" >&2; exit 2; }", "  vina_out=$(cd \"$ROOT/batches/$batch\" && bsub -w \"done($conf_id)\" < run_vina.lsf)", "  vina_id=$(printf '%s' \"$vina_out\" | sed -n 's/.*<\\([0-9][0-9]*\\)>.*/\\1/p')", "  [ -n \"$vina_id\" ] || { echo \"Could not parse Vina submission: $vina_out\" >&2; exit 2; }", "  printf '{\"batch\":\"%s\",\"confgen\":\"%s\",\"vina\":\"%s\"}\\n' \"$batch\" \"$conf_id\" \"$vina_id\" >> \"$ROOT/job_ids.jsonl\"", "done"]
    dispatcher = "\n".join(lines) + "\n"
    encoded = base64.b64encode(dispatcher.encode()).decode()
    command = f"cd {shell_quote(remote)} && printf %s {shell_quote(encoded)} | base64 -d > submit_campaign.sh && chmod 700 submit_campaign.sh && nohup ./submit_campaign.sh > dispatcher.log 2>&1 & echo dispatcher-started"
    output = run([settings.ssh_bin, "-o", "BatchMode=yes", settings.pegasus_host, command])
    if "dispatcher-started" not in output:
        raise PortalError("Pegasus did not start the campaign dispatcher: " + output)
    return {"remote": remote, "results": results_root, "campaign": "dispatcher started", "batches": str(len(batches))}


def write_docking_answers(root: Path, form) -> None:
    """Persist the portal equivalents of 0_LIGSPLIT/1B prompts with the package."""
    answers = {key: form.get(key, "") for key in ("ligand_input", "batch_size", "batch_prefix", "csv_smiles_col", "csv_id_col", "confgen_poses", "confgen_workers", "vina_poses", "vina_workers", "queue", "project", "confgen_walltime", "vina_walltime")}
    answers["enumerate_protomers"] = bool(form.get("enumerate_protomers"))
    root.joinpath("portal_docking_answers.json").write_text(json.dumps(answers, indent=2) + "\n", encoding="utf-8")


def submit_md(settings: Settings, root: Path, name: str, submitter: str) -> dict[str, str]:
    project = require_md_project(root)
    relative = project.relative_to(root).as_posix()
    remote = stage_project(settings, root, settings.triton_host, settings.triton_remote_root, name, submitter, "run")
    workdir = remote if relative == "." else f"{remote}/{relative}"
    # Setup is intentionally a separate explicit phase: PyMACS only renders a
    # production LSF job after required setup outputs are verified.
    command = (
        f"source /etc/bashrc && module load gcc-runtime/13.4.0-none-none-gh6rpjx gromacs/2025.1-gcc-13.4.0-6vq7xfo && "
        f"source /projectnb/triton/home/jxs794/miniforge3/etc/profile.d/conda.sh && conda activate cgenff && "
        f"cp -a {shell_quote(settings.pymacs_root + '/.')} {shell_quote(workdir)}/ && "
        f"cd {shell_quote(workdir)} && python pymacs_run.py validate && "
        f"python pymacs_run.py setup && python pymacs_run.py submit --with-analysis"
    )
    output = run([settings.ssh_bin, "-o", "BatchMode=yes", settings.triton_host, command])
    ids = JOB_ID.findall(output)
    if not ids:
        raise PortalError("Triton completed setup but did not return an LSF job ID: " + output)
    return {"remote": workdir, "md": ids[0], "analysis": ids[1] if len(ids) > 1 else "not requested"}


def submit_md_campaign(settings: Settings, root: Path, name: str, submitter: str) -> dict[str, str]:
    """Stage and submit every independently configured MD system in a campaign."""
    try:
        projects = json.loads((root / "md_campaign.json").read_text(encoding="utf-8"))["projects"]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise PortalError("The MD campaign manifest is invalid.") from exc
    if not isinstance(projects, list) or not projects:
        raise PortalError("The MD campaign contains no projects.")
    remote = stage_project(settings, root, settings.triton_host, settings.triton_remote_root, name, submitter, "campaign")
    lines = ["#!/bin/bash", "set -euo pipefail",
             "source /etc/bashrc", "module load gcc-runtime/13.4.0-none-none-gh6rpjx gromacs/2025.1-gcc-13.4.0-6vq7xfo",
             "source /projectnb/triton/home/jxs794/miniforge3/etc/profile.d/conda.sh", "conda activate cgenff",
             f"ROOT={shell_quote(remote)}"]
    for project in projects:
        project_id = str(project.get("id", ""))
        if not re.fullmatch(r"system-[0-9]{3}-[A-Za-z0-9_.-]{1,48}", project_id):
            raise PortalError("The MD campaign manifest contains an invalid project ID.")
        submit_flag = " --with-analysis" if project.get("analysis") else ""
        lines += [f"WORKDIR=\"$ROOT/systems/{project_id}\"",
                  f"cp -a {shell_quote(settings.pymacs_root + '/.')} \"$WORKDIR\"/",
                  "cd \"$WORKDIR\"", "python pymacs_run.py validate", "python pymacs_run.py setup",
                  f"python pymacs_run.py submit{submit_flag}"]
    command = "\n".join(lines) + "\n"
    output = run([settings.ssh_bin, "-o", "BatchMode=yes", settings.triton_host, command])
    ids = JOB_ID.findall(output)
    if len(ids) < len(projects):
        raise PortalError("Triton did not return an MD job ID for every system: " + output[-1000:])
    return {"remote": remote, "systems": str(len(projects)), "md jobs": ", ".join(ids)}


def machine_status_path(settings: Settings) -> Path:
    return settings.stage_root / "machine_status.json"


def read_machine_status(settings: Settings) -> dict[str, object]:
    try:
        payload = json.loads(machine_status_path(settings).read_text(encoding="utf-8"))
        return payload if isinstance(payload, dict) else {}
    except (OSError, ValueError):
        return {}


def probe_machine(settings: Settings, machine: dict[str, str]) -> dict[str, str]:
    """Perform a harmless, non-interactive SSH availability check from Kyle."""
    try:
        completed = subprocess.run(
            [settings.ssh_bin, "-o", "BatchMode=yes", "-o", "ConnectTimeout=20",
             "-o", "StrictHostKeyChecking=yes", machine["host"], "printf portal-ssh-ok"],
            text=True, capture_output=True, timeout=25,
        )
        detail = (completed.stderr or completed.stdout or "").strip().replace("\n", " ")
        if completed.returncode == 0 and "portal-ssh-ok" in completed.stdout:
            state, detail = "online", "SSH check succeeded"
        elif "Host key verification failed" in detail:
            state, detail = "host_key", "Host key needs administrator review"
        elif "Permission denied" in detail:
            state, detail = "auth", "Host reached, but Kyle's SSH key is not authorized"
        elif "Could not resolve hostname" in detail:
            state, detail = "unconfigured", "No resolvable SSH host configuration"
        else:
            state, detail = "offline", "Host is unreachable or SSH is unavailable"
    except subprocess.TimeoutExpired:
        state, detail = "offline", "SSH connection timed out"
    return {**machine, "state": state, "detail": detail[:220]}


def refresh_machine_status(settings: Settings) -> None:
    checked_at = dt.datetime.now(dt.timezone.utc).isoformat()
    with ThreadPoolExecutor(max_workers=8, thread_name_prefix="machine-probe") as pool:
        futures = [pool.submit(probe_machine, settings, machine) for machine in LAB_MACHINES]
        machines = [future.result() for future in as_completed(futures)]
    machines.sort(key=lambda item: item["name"].casefold())
    settings.stage_root.mkdir(parents=True, exist_ok=True)
    temporary = machine_status_path(settings).with_suffix(".tmp")
    temporary.write_text(json.dumps({"checked_at": checked_at, "machines": machines}), encoding="utf-8")
    temporary.replace(machine_status_path(settings))


def start_machine_status_worker(settings: Settings) -> None:
    """Probe the inventory out of band; page loads only read the cached result."""
    def worker() -> None:
        while True:
            try:
                refresh_machine_status(settings)
            except Exception:
                pass
            threading.Event().wait(300)
    threading.Thread(target=worker, name="machine-status-refresh", daemon=True).start()


def create_app(settings: Settings | None = None) -> Flask:
    app = Flask(__name__)
    app.config["MAX_CONTENT_LENGTH"] = MAX_ARCHIVE_BYTES
    app.config["SECRET_KEY"] = os.environ.get("LAB_PORTAL_SECRET", os.urandom(32))
    app.config["SETTINGS"] = settings or settings_from_env()
    prepare_pdbqt_cache(app.config["SETTINGS"])
    start_pdbqt_cache_worker(app.config["SETTINGS"])
    start_job_monitor_cache_worker(app.config["SETTINGS"])
    start_machine_status_worker(app.config["SETTINGS"])
    start_md_results_worker(app.config["SETTINGS"])
    start_legacy_resume_worker(app.config["SETTINGS"])

    @app.get("/")
    def index():
        return render_template("index.html")

    @app.get("/machines")
    def machines():
        return render_template("machines.html", status=read_machine_status(app.config["SETTINGS"]))

    @app.get("/results")
    def results():
        index = read_md_results_index(app.config["SETTINGS"])
        for system in index.get("systems", []):
            if isinstance(system, dict):
                files = system.get("files", [])
                system["images"] = [item for item in files if isinstance(item, dict) and str(item.get("relative", "")).lower().endswith(".png")]
                # Every PNG exported by a finished analysis is potentially
                # meaningful.  Do not hide a valid scientific figure simply
                # because its filename does not match a preset category.
                system["graphs"] = system["images"]
        return render_template("results.html", results=index)

    @app.get("/results/system/<path:directory>")
    def result_detail(directory: str):
        """Serve one published result dossier without reopening Triton."""
        requested = safe_result_relative(directory)
        if requested is None:
            abort(404)
        index = read_md_results_index(app.config["SETTINGS"])
        for system in index.get("systems", []):
            if not isinstance(system, dict) or str(system.get("directory", "")) != requested.as_posix():
                continue
            files = system.get("files", [])
            # A generated single-system PyMACS report is the primary results
            # experience: it contains the interactive, explained gallery and
            # opens entirely from Kyle-local files.  Preserve this dossier
            # fallback for running/legacy systems that have no report yet.
            report_relative = "Analysis_Results/PyMACS_Result_Dashboard.html"
            # The publisher index is intentionally refreshed out of band.  A
            # freshly generated report may already be safely present on Kyle
            # in the short interval before that cache records it.  Check this
            # exact trusted local path as well, so a scientist never has to
            # wait through a dossier just because the five-minute inventory
            # cycle has not caught up.
            local_report = published_md_results_root(app.config["SETTINGS"]).joinpath(
                *requested.parts, *PurePosixPath(report_relative).parts
            )
            report = next(
                (item for item in files
                 if isinstance(item, dict)
                 and str(item.get("relative", "")) == report_relative),
                None,
            )
            if report or local_report.is_file():
                return redirect(url_for(
                    "result_file",
                    relative=requested.as_posix() + "/" + report_relative,
                    inline="1",
                ))
            system["images"] = [item for item in files if isinstance(item, dict) and str(item.get("relative", "")).lower().endswith(".png")]
            system["graphs"] = system["images"]
            # CSVs stay Kyle-local too.  The dossier loads them only when a
            # scientist asks for an interactive plot, keeping initial page
            # loads as fast as the static results library.
            system["csv_files"] = [item for item in files if isinstance(item, dict) and str(item.get("relative", "")).lower().endswith(".csv")]
            return render_template("result_detail.html", system=system)
        abort(404)

    @app.get("/results/files/<path:relative>")
    def result_file(relative: str):
        safe_relative = safe_result_relative(relative)
        if safe_relative is None:
            abort(404)
        root = published_md_results_root(app.config["SETTINGS"]).resolve()
        candidate = root.joinpath(*safe_relative.parts).resolve()
        if root not in candidate.parents or not candidate.is_file():
            abort(404)
        # Interactive PyMACS result dashboards are generated by our own
        # trusted analysis step and are intentionally opened in-browser.
        inline = request.args.get("inline") == "1" and candidate.suffix.lower() in {".png", ".pdf", ".html"}
        return send_from_directory(root, safe_relative.as_posix(), as_attachment=not inline,
                                   download_name=candidate.name)

    @app.get("/results/download/<path:directory>")
    def download_result_data(directory: str):
        """Deliver a compact, Kyle-local result package for one published system.

        Full trajectories and generated HTML reports are intentionally omitted:
        this route is for a scientist-friendly handoff of structures, figures,
        tables, chemistry, and run metadata rather than a raw project archive.
        """
        requested = safe_result_relative(directory)
        if requested is None:
            abort(404)
        index = read_md_results_index(app.config["SETTINGS"])
        system = next((item for item in index.get("systems", [])
                       if isinstance(item, dict)
                       and str(item.get("directory", "")) == requested.as_posix()), None)
        if system is None:
            abort(404)

        root = published_md_results_root(app.config["SETTINGS"]).resolve()
        system_root = root.joinpath(*requested.parts).resolve()
        if root not in system_root.parents or not system_root.is_dir():
            abort(404)

        files = [item for item in system.get("files", []) if isinstance(item, dict)]
        has_curated_pdb = any(str(item.get("relative", "")) == "selected_input_curated.pdb" for item in files)
        allowed_suffixes = {".csv", ".png", ".svg", ".mol", ".mol2", ".sdf"}
        allowed_names = {
            "selected_input_curated.pdb", "selected_input_source.pdb",
            "Final_Trajectory.pdb", "binding_pocket_only.pdb", "pymacs_run.json",
        }
        selected: list[tuple[Path, PurePosixPath]] = []
        for item in files:
            relative_text = str(item.get("relative", ""))
            relative = safe_result_relative(relative_text)
            if relative is None:
                continue
            name = relative.name
            if name == "selected_input_source.pdb" and has_curated_pdb:
                continue
            if relative.suffix.lower() not in allowed_suffixes and name not in allowed_names:
                continue
            candidate = system_root.joinpath(*relative.parts).resolve()
            if system_root not in candidate.parents or not candidate.is_file():
                continue
            selected.append((candidate, relative))

        safe_stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(system.get("system") or requested.name)).strip("._") or "pymacs_result"
        bundles = root / "_download_bundles"
        bundles.mkdir(mode=0o755, exist_ok=True)
        bundle = bundles / f"{requested.as_posix().replace('/', '__')}.zip"
        newest_source = max((source.stat().st_mtime for source, _ in selected), default=0)
        if not bundle.is_file() or bundle.stat().st_mtime < newest_source:
            temporary = bundle.with_suffix(".zip.new")
            with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
                prefix = f"{safe_stem}_lightweight_results/"
                archive.writestr(
                    prefix + "README.txt",
                    "PyMACS lightweight results package\n\n"
                    "Includes published plots/tables, shareable PDB structures, ligand chemistry, and run metadata.\n"
                    "Large raw trajectories, checkpoints, and the interactive HTML report remain available on the lab portal/archive.\n",
                )
                for source, relative in selected:
                    archive.write(source, prefix + relative.as_posix())
            temporary.replace(bundle)
        return send_file(bundle, as_attachment=True, conditional=True,
                         download_name=f"{safe_stem}_lightweight_results.zip", max_age=0)

    @app.get("/api/submitters")
    def submitters():
        return jsonify(names=saved_submitters(app.config["SETTINGS"]))

    @app.get("/api/jobs")
    def jobs():
        cached = read_job_monitor_cache(app.config["SETTINGS"])
        if cached is not None:
            return jsonify(**cached, refreshing=False)
        return jsonify(campaigns=[], generated_at=None, refreshing=True)

    @app.post("/api/inspect-docking")
    def inspect_docking():
        upload = request.files.get("project")
        if not upload or not upload.filename or not upload.filename.lower().endswith(".zip"):
            return jsonify(error="Upload the AutoDock Vina package ZIP first."), 400
        settings = app.config["SETTINGS"]
        with tempfile.TemporaryDirectory(prefix="lab-dock-inspect-", dir=settings.stage_root) as work:
            archive = Path(work) / secure_filename(upload.filename); upload.save(archive)
            try:
                safe_extract(archive, Path(work) / "project")
                report = inspect_docking_package(project_root(Path(work) / "project"))
                if not report["inputs"]:
                    raise PortalError("No SDF, SMILES, or CSV ligand input was found in this package.")
                return jsonify(report)
            except (PortalError, zipfile.BadZipFile) as exc:
                return jsonify(error=str(exc)), 400

    @app.post("/api/inspect-structure")
    def inspect_structure():
        upload, rcsb_id = request.files.get("structure"), request.form.get("rcsb_id", "").strip()
        if (not upload or not upload.filename) and not rcsb_id:
            return jsonify(error="Upload one structure/ZIP or enter an RCSB PDB ID."), 400
        settings = app.config["SETTINGS"]
        settings.stage_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="lab-inspect-", dir=settings.stage_root) as work:
            root = Path(work); structure: Path
            try:
                if rcsb_id:
                    structure = root / fetch_rcsb_pdb(rcsb_id, root)
                else:
                    incoming = root / secure_filename(upload.filename); upload.save(incoming)
                    if incoming.suffix.lower() == ".zip":
                        extracted = root / "extracted"; extracted.mkdir(); safe_extract(incoming, extracted)
                        candidates = [p for p in extracted.rglob("*") if p.is_file() and p.suffix.lower() in STRUCTURE_SUFFIXES]
                        if len(candidates) != 1: raise PortalError("ZIP must contain exactly one PDB, CIF, or mmCIF structure.")
                        structure = candidates[0]
                    else:
                        structure = incoming
                    if structure.suffix.lower() not in STRUCTURE_SUFFIXES: raise PortalError("Use PDB, CIF, or mmCIF input.")
                return jsonify({"filename": structure.name, **inspect_structure_file(structure)})
            except (PortalError, zipfile.BadZipFile) as exc:
                return jsonify(error=str(exc)), 400

    @app.post("/api/inspect-md-campaign")
    def inspect_md_campaign():
        upload, rcsb_id = request.files.get("structures"), request.form.get("rcsb_id", "").strip()
        if (not upload or not upload.filename) and not rcsb_id:
            return jsonify(error="Upload a PDB/mmCIF file or ZIP of structures first."), 400
        settings = app.config["SETTINGS"]
        settings.stage_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="lab-md-inspect-", dir=settings.stage_root) as work:
            root = Path(work) / "project"; root.mkdir()
            try:
                if rcsb_id:
                    fetch_rcsb_pdb(rcsb_id, root)
                else:
                    incoming = Path(work) / secure_filename(upload.filename); upload.save(incoming)
                    if incoming.suffix.lower() == ".zip":
                        safe_extract(incoming, root)
                    elif incoming.suffix.lower() in STRUCTURE_SUFFIXES:
                        shutil.move(str(incoming), root / incoming.name)
                    else:
                        raise PortalError("Use PDB, CIF, mmCIF, or a ZIP containing those files.")
                return jsonify(systems=discover_md_structures(root))
            except (PortalError, zipfile.BadZipFile) as exc:
                return jsonify(error=str(exc)), 400

    @app.post("/submit/<kind>")
    def submit(kind: str):
        if kind not in {"docking", "md"}:
            return "Unknown job type", 404
        name, submitter, email = request.form.get("job_name", "").strip(), request.form.get("submitter", "").strip(), request.form.get("submitter_email", "").strip()
        upload = request.files.get("project")
        if not JOB_NAME.fullmatch(name):
            flash("Job name: letters, numbers, dots, dashes, and underscores only (max 80).", "error")
            return redirect(url_for("index"))
        try:
            submitter_slug(submitter)
        except PortalError as exc:
            flash(str(exc), "error"); return redirect(url_for("index"))
        if email and not EMAIL.fullmatch(email):
            flash("Enter a valid email address or leave it blank.", "error"); return redirect(url_for("index"))
        remember_submitter(app.config["SETTINGS"], submitter, email)
        if not upload or not upload.filename or not upload.filename.lower().endswith(".zip"):
            flash("Upload one .zip project folder.", "error")
            return redirect(url_for("index"))
        settings = app.config["SETTINGS"]
        settings.stage_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="lab-portal-", dir=settings.stage_root) as work:
            workdir = Path(work)
            archive = workdir / secure_filename(upload.filename)
            upload.save(archive)
            if archive.stat().st_size > MAX_ARCHIVE_BYTES:
                flash("Archive exceeds the 2 GB upload limit.", "error")
                return redirect(url_for("index"))
            extracted = workdir / "project"
            extracted.mkdir()
            try:
                safe_extract(archive, extracted)
                root = project_root(extracted)
                if kind == "docking":
                    write_docking_answers(root, request.form)
                    build_docking_campaign(settings, root, request.form)
                result = submit_docking(settings, root, name, submitter) if kind == "docking" else submit_md(settings, root, name, submitter)
                record_submission(settings, submitter, kind, name, result)
            except (PortalError, zipfile.BadZipFile) as exc:
                flash(str(exc), "error")
                return redirect(url_for("index"))
        return render_template("submitted.html", kind=kind, result=result)

    @app.post("/submit/md-wizard")
    def submit_md_wizard():
        name, submitter, email, upload = request.form.get("job_name", "").strip(), request.form.get("submitter", "").strip(), request.form.get("submitter_email", "").strip(), request.files.get("structure")
        if not JOB_NAME.fullmatch(name):
            flash("Job name: letters, numbers, dots, dashes, and underscores only (max 80).", "error")
            return redirect(url_for("index"))
        try:
            submitter_slug(submitter)
        except PortalError as exc:
            flash(str(exc), "error"); return redirect(url_for("index"))
        if email and not EMAIL.fullmatch(email):
            flash("Enter a valid email address or leave it blank.", "error"); return redirect(url_for("index"))
        remember_submitter(app.config["SETTINGS"], submitter, email)
        rcsb_id = request.form.get("rcsb_id", "").strip()
        if (not upload or not upload.filename) and not rcsb_id:
            flash("Choose a structure file/ZIP or enter an RCSB PDB ID.", "error")
            return redirect(url_for("index"))
        settings = app.config["SETTINGS"]
        settings.stage_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="lab-portal-", dir=settings.stage_root) as work:
            root = Path(work) / "project"; root.mkdir()
            try:
                if rcsb_id:
                    structure_name = fetch_rcsb_pdb(rcsb_id, root)
                else:
                    incoming = Path(work) / secure_filename(upload.filename)
                    upload.save(incoming)
                    if incoming.suffix.lower() == ".zip":
                        safe_extract(incoming, root)
                        candidates = [p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in STRUCTURE_SUFFIXES]
                        if len(candidates) != 1:
                            raise PortalError("MD ZIP must contain exactly one PDB, CIF, or mmCIF structure.")
                        structure = candidates[0]
                        if structure.parent != root:
                            target = root / secure_filename(structure.name)
                            shutil.move(str(structure), target); structure = target
                    else:
                        structure = root / incoming.name
                        shutil.move(str(incoming), structure)
                    structure_name = structure.name
                write_md_config(root, structure_name, request.form)
                if request.form.get("system_type") in {"ligand", "protac"}:
                    prepare_ligand_cgenff(settings, root, request.form.get("ligand", ""))
                result = submit_md(settings, root, name, submitter)
                record_submission(settings, submitter, "md", name, result)
            except (PortalError, zipfile.BadZipFile) as exc:
                flash(str(exc), "error"); return redirect(url_for("index"))
        return render_template("submitted.html", kind="MD", result=result)

    @app.post("/submit/md-campaign")
    def submit_md_campaign_wizard():
        name = request.form.get("job_name", "").strip()
        submitter, email = request.form.get("submitter", "").strip(), request.form.get("submitter_email", "").strip()
        upload, rcsb_id = request.files.get("structures"), request.form.get("rcsb_id", "").strip()
        if not JOB_NAME.fullmatch(name):
            flash("Job name: letters, numbers, dots, dashes, and underscores only (max 80).", "error")
            return redirect(url_for("index"))
        try:
            submitter_slug(submitter)
        except PortalError as exc:
            flash(str(exc), "error"); return redirect(url_for("index"))
        if email and not EMAIL.fullmatch(email):
            flash("Enter a valid email address or leave it blank.", "error"); return redirect(url_for("index"))
        if (not upload or not upload.filename) and not rcsb_id:
            flash("Upload a structure ZIP first.", "error"); return redirect(url_for("index"))
        try:
            requested = json.loads(request.form.get("systems_json", "[]"))
            if not isinstance(requested, list):
                raise ValueError
        except ValueError:
            flash("The per-system MD settings were invalid. Inspect the upload again.", "error"); return redirect(url_for("index"))
        settings = app.config["SETTINGS"]
        settings.stage_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="lab-md-campaign-", dir=settings.stage_root) as work:
            root = Path(work) / "project"; root.mkdir()
            try:
                if rcsb_id:
                    fetch_rcsb_pdb(rcsb_id, root)
                else:
                    incoming = Path(work) / secure_filename(upload.filename); upload.save(incoming)
                    if incoming.suffix.lower() == ".zip":
                        safe_extract(incoming, root)
                    elif incoming.suffix.lower() in STRUCTURE_SUFFIXES:
                        shutil.move(str(incoming), root / incoming.name)
                    else:
                        raise PortalError("Use PDB, CIF, mmCIF, or a ZIP containing those files.")
                available = {str(item["path"]) for item in discover_md_structures(root)}
                chosen = []
                for item in requested:
                    if not isinstance(item, dict) or str(item.get("path", "")) not in available:
                        raise PortalError("One selected structure does not match the uploaded ZIP.")
                    normalized = {key: item.get(key, "") for key in ("path", "system_type", "system_name", "ligand", "ligand_name", "keep_chains", "cofactors", "box_type", "box_distance_nm", "length_ns", "threads", "chain_map", "replicate_id", "purpose_notes", "pocket_cutoff", "contact_cutoff", "min_contact_fraction", "interface_contact_cutoff", "interface_min_contact_fraction", "interface_frame_step", "interface_max_edges", "protac_ligase_chain", "protac_ligase_name", "protac_target_chains", "protac_target_name")}
                    normalized["protac_target_chain_names"] = json.dumps(item.get("protac_target_chain_names", {}))
                    for flag in ("remove_waters", "remove_ions", "remove_solvents", "analysis", "figurebook"):
                        normalized[flag] = "on" if item.get(flag) else ""
                    chosen.append(normalized)
                if len(chosen) != len(available):
                    raise PortalError("Configure every detected structure before submitting the MD campaign.")
                manifest = build_md_campaign(settings, root, chosen)
                result = submit_md_campaign(settings, root, name, submitter)
                record_submission(settings, submitter, "md", name, result)
            except (PortalError, zipfile.BadZipFile) as exc:
                flash(str(exc), "error"); return redirect(url_for("index"))
        return render_template("submitted.html", kind="MD campaign", result=result)

    return app


if __name__ == "__main__":
    app = create_app()
    try:
        from waitress import serve
    except ImportError:
        # Kyle's minimal OS image may not have pip yet. This fallback keeps the
        # VPN-only pilot usable; production service deployment uses Waitress.
        app.run(host="127.0.0.1", port=int(os.environ.get("PORT", "8788")))
    else:
        serve(app, host="127.0.0.1", port=int(os.environ.get("PORT", "8788")))
