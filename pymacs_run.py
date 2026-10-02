#!/usr/bin/env python3
"""Per-folder PyMACS run configuration and execution helper.

Run this script from one simulation folder.  ``configure`` asks questions,
stores the answers in ``pymacs_run.json``, and later commands replay those
answers without asking again.  The optional ``triton`` profile renders and
submits an LSF GPU job; it is never selected implicitly.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


CONFIG_NAME = "pymacs_run.json"
ROOT = Path(__file__).resolve().parent
STRUCTURE_SUFFIXES = {".pdb", ".cif", ".mmcif"}
GENERATED_NAMES = {"protein.pdb", "lig_h.pdb", "lig_raw.pdb"}


class ConfigError(RuntimeError):
    """Raised when a run folder cannot safely be acted on."""


def config_path(folder: Path) -> Path:
    return folder / CONFIG_NAME


def write_config(folder: Path, config: dict[str, Any]) -> None:
    config["schema_version"] = 1
    config["configured_at"] = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    config_path(folder).write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")


def load_config(folder: Path) -> dict[str, Any]:
    path = config_path(folder)
    if not path.is_file():
        raise ConfigError(f"No {CONFIG_NAME} found. Run 'python pymacs_run.py configure' first.")
    try:
        config = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ConfigError(f"{path.name} is not valid JSON: {exc}") from exc
    validate_config(config, folder)
    return config


def candidate_structures(folder: Path) -> list[Path]:
    candidates = []
    for path in sorted(folder.iterdir()):
        if not path.is_file() or path.suffix.lower() not in STRUCTURE_SUFFIXES:
            continue
        if path.name in GENERATED_NAMES or path.name.startswith(("selected_input_", "protein_processed")):
            continue
        candidates.append(path)
    return candidates


def inspect_pdb(path: Path) -> tuple[list[str], list[str]]:
    """Return polymer chain IDs and non-water HETATM residue names for PDB input."""
    if path.suffix.lower() != ".pdb":
        return [], []
    chains, hetero = set(), set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "ATOM":
                chains.add(line[21:22].strip() or "(blank)")
            elif record == "HETATM":
                residue = line[17:20].strip().upper()
                if residue and residue not in {"HOH", "WAT", "SOL"}:
                    hetero.add(residue)
    return sorted(chains), sorted(hetero)


def inspect_cif(path: Path) -> tuple[list[str], list[str]]:
    """Read chain IDs and non-water components from a standard mmCIF atom loop.

    This deliberately uses only the Python standard library so configuration
    works before either PyMACS conda environment has been activated.
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for loop_index, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        headers, cursor = [], loop_index + 1
        while cursor < len(lines) and lines[cursor].lstrip().startswith("_"):
            headers.append(lines[cursor].split()[0])
            cursor += 1
        if "_atom_site.group_PDB" not in headers:
            continue
        indices = {name: headers.index(name) for name in headers}
        chain_key = "_atom_site.auth_asym_id" if "_atom_site.auth_asym_id" in indices else "_atom_site.label_asym_id"
        residue_key = "_atom_site.label_comp_id"
        chains, hetero, tokens = set(), set(), []
        while cursor < len(lines):
            current = lines[cursor].strip()
            if not current or current.startswith("#") or current == "loop_" or current.startswith("_"):
                break
            tokens.extend(shlex.split(current))
            while len(tokens) >= len(headers):
                row, tokens = tokens[:len(headers)], tokens[len(headers):]
                record = row[indices["_atom_site.group_PDB"]]
                chain = row[indices[chain_key]]
                residue = row[indices[residue_key]].upper()
                if record == "ATOM" and chain not in {".", "?"}:
                    chains.add(chain)
                elif record == "HETATM" and residue not in {"HOH", "WAT", "SOL"}:
                    hetero.add(residue)
            cursor += 1
        return sorted(chains), sorted(hetero)
    return [], []


def inspect_structure(path: Path) -> tuple[list[str], list[str]]:
    return inspect_pdb(path) if path.suffix.lower() == ".pdb" else inspect_cif(path)


def ask(prompt: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{prompt}{suffix}: ").strip()
    return value or (default or "")


def ask_yes_no(prompt: str, default: bool = True) -> bool:
    answer = ask(f"{prompt} [{'Y/n' if default else 'y/N'}]").lower()
    return default if not answer else answer in {"y", "yes"}


def ask_choice(prompt: str, choices: list[tuple[str, str]], default: int = 1) -> str:
    print(f"\n{prompt}")
    for index, (_, description) in enumerate(choices, 1):
        print(f"  {index}) {description}")
    while True:
        answer = ask("Choose", str(default))
        if answer.isdigit() and 1 <= int(answer) <= len(choices):
            return choices[int(answer) - 1][0]
        print(f"Please enter a number from 1 to {len(choices)}.")


def configure(folder: Path) -> None:
    structures = candidate_structures(folder)
    if not structures:
        raise ConfigError("No .pdb, .cif, or .mmcif structure was found in this folder.")
    print("\nPyMACS per-folder configuration")
    print("Each folder describes exactly one structure and one simulation.\n")
    for index, path in enumerate(structures, 1):
        print(f"  {index}) {path.name}")
    while True:
        choice = ask("Choose input structure", "1")
        if choice.isdigit() and 1 <= int(choice) <= len(structures):
            structure = structures[int(choice) - 1]
            break
        print("Choose one of the numbered structures.")

    chains, hetero = inspect_structure(structure)
    if chains:
        print(f"\nDetected polymer chains: {', '.join(chains)}")
    if hetero:
        print(f"Detected non-water HETATM residues: {', '.join(hetero)}")

    system_type = ask_choice(
        "What system are you preparing?",
        [
            ("protein", "Protein–protein complex"),
            ("peptide", "Protein–peptide complex"),
            ("ligand", "Protein–small-molecule complex"),
            ("biological", "RNA/DNA/protein biological assembly (advanced)"),
        ],
    )
    selected_chains = chains[:]
    if len(chains) > 1 and not ask_yes_no("Keep all detected polymer chains", True):
        while True:
            raw_chains = ask("Chains to retain (comma-separated IDs)", ",".join(chains))
            selected_chains = [item.strip() for item in raw_chains.split(",") if item.strip()]
            unknown = sorted(set(selected_chains) - set(chains))
            if selected_chains and not unknown:
                break
            print(f"Use one or more detected chain IDs only: {', '.join(chains)}")
    chain_map = ""
    if selected_chains:
        names = []
        for chain in selected_chains:
            name = ask(f"Name for chain {chain}", chain)
            names.append(f"{chain}:{name}")
        chain_map = ",".join(names)

    ligand = cofactors = ""
    if system_type == "ligand":
        if hetero:
            print(f"Detected possible ligand residues: {', '.join(hetero)}")
        ligand = ask("Ligand residue code (for example LIG)").upper()
        if not ligand:
            raise ConfigError("A ligand residue code is required for a ligand workflow.")
        cofactors = ask("Optional retained cofactors (comma-separated; press Enter for none)").upper()

    box_type = ask_choice(
        "Choose simulation box shape:",
        [
            ("dodecahedron", "Dodecahedron — recommended for most proteins/complexes; compact and water-efficient."),
            ("cubic", "Cubic — simple conventional box; usually contains more water."),
            ("octahedron", "Octahedron — compact alternative for globular systems."),
            ("triclinic", "Triclinic — advanced general-purpose geometry."),
        ],
    )
    distance = float(ask("Solute-to-box-edge distance in nm", "1.0"))
    if distance <= 0:
        raise ConfigError("Box distance must be greater than zero.")
    length = float(ask("Production length in ns", "100"))
    if length <= 0:
        raise ConfigError("Production length must be greater than zero.")
    analysis_enabled = ask_yes_no("Run chunked Step 3 analysis automatically after MD", True)
    analysis = {"enabled": analysis_enabled, "figurebook": False, "threads": 8, "chunk_size": 1000}
    if analysis_enabled:
        analysis["figurebook"] = ask_yes_no("Generate the Step 4 figurebook PDF", True)
        analysis["threads"] = int(ask("Analysis CPU threads", "8"))
        analysis["chunk_size"] = int(ask("Frames per analysis chunk", "1000"))
        if system_type in {"protein", "peptide"}:
            if not ask_yes_no("Analyze all chain interfaces", True):
                analysis["interface_chain_pairs"] = ask("Chain pairs to analyze (for example A:B)")
            analysis["interface_contact_cutoff"] = float(ask("Interface contact cutoff in Å", "4.0"))
            analysis["interface_min_contact_fraction"] = float(ask("Minimum interface contact fraction", "0.10"))
            analysis["interface_frame_step"] = int(ask("Analyze every Nth interface frame", "1"))
            analysis["interface_max_edges"] = int(ask("Maximum interface-network edges", "100"))
        elif system_type == "ligand":
            analysis["pocket_cutoff"] = float(ask("Ligand pocket cutoff in Å", "5.0"))
            analysis["contact_cutoff"] = float(ask("Ligand contact cutoff in Å", "4.0"))
            analysis["min_contact_fraction"] = float(ask("Minimum persistent-contact fraction", "0.10"))
            analysis["compound_name"] = ask("Compound display name", ligand)

    config = {
        "input": structure.name,
        "setup": {
            "mode": "ligand" if system_type == "ligand" else "protein",
            "chain_map": chain_map or None,
            "keep_chains": ",".join(selected_chains) if set(selected_chains) != set(chains) else None,
            "ligand": ligand or None,
            "cofactors": cofactors or None,
            "remove_input_waters": ask_yes_no("Remove crystallographic waters", True),
            "remove_input_ions": ask_yes_no("Remove deposited ions", True),
            "box_type": box_type,
            "box_distance_nm": distance,
        },
        "simulation": {
            "mode": system_type,
            "length_ns": length,
            "threads": 16,
            "compute": "GPU",
        },
        "analysis": analysis,
    }
    print("\nReady to write:")
    print(f"  Input: {structure.name}")
    print(f"  Mode: {system_type}")
    print(f"  Box: {box_type}, {distance:g} nm")
    print(f"  Production: {length:g} ns")
    if not ask_yes_no(f"Write {CONFIG_NAME}", True):
        print("Configuration cancelled.")
        return
    write_config(folder, config)
    print(f"\nWrote {config_path(folder)}")
    print("Next: review with 'python pymacs_run.py validate', then run setup.")


def validate_config(config: dict[str, Any], folder: Path) -> None:
    for key in ("input", "setup", "simulation"):
        if key not in config:
            raise ConfigError(f"Configuration is missing required field '{key}'.")
    input_path = folder / str(config["input"])
    if not input_path.is_file() or input_path.suffix.lower() not in STRUCTURE_SUFFIXES:
        raise ConfigError(f"Configured input is missing or unsupported: {input_path.name}")
    setup, simulation = config["setup"], config["simulation"]
    if setup.get("mode") not in {"protein", "ligand"}:
        raise ConfigError("setup.mode must be protein or ligand.")
    if simulation.get("mode") not in {"protein", "peptide", "ligand", "biological"}:
        raise ConfigError("simulation.mode is not supported.")
    if setup.get("mode") == "ligand" and not setup.get("ligand"):
        raise ConfigError("A ligand workflow requires setup.ligand.")
    if setup.get("box_type") not in {"cubic", "dodecahedron", "octahedron", "triclinic"}:
        raise ConfigError("setup.box_type is not a supported GROMACS box type.")
    if float(setup.get("box_distance_nm", 0)) <= 0 or float(simulation.get("length_ns", 0)) <= 0:
        raise ConfigError("Box distance and production length must be positive.")


def setup_command(folder: Path, config: dict[str, Any]) -> list[str]:
    script = folder / "1_AutomateGromacs_MPI.py"
    if not script.is_file():
        raise ConfigError("1_AutomateGromacs_MPI.py is missing. Install/copy PyMACS into this run folder first.")
    setup = config["setup"]
    command = [sys.executable, str(script), "--gmx-bin", "gmx_mpi", "--pdb", str(config["input"]), "--mode", setup["mode"],
               "--box-type", setup["box_type"], "--box-distance", str(setup["box_distance_nm"])]
    if setup.get("chain_map"):
        command += ["--chain-map", setup["chain_map"]]
    if setup.get("keep_chains"):
        command += ["--keep-chains", setup["keep_chains"]]
    if setup.get("ligand"):
        command += ["--ligand", setup["ligand"]]
    if setup.get("cofactors"):
        command += ["--cofactors", setup["cofactors"]]
    if setup.get("remove_input_waters"):
        command.append("--remove-input-waters")
    if setup.get("remove_input_ions"):
        command.append("--remove-input-ions")
    return command


def simulation_command(config: dict[str, Any], resume: bool = False) -> list[str]:
    sim, setup = config["simulation"], config["setup"]
    command = ["python", "2_AutomateGromacs_MPI.py", "--gmx-bin", "gmx_mpi", "--mode", sim["mode"],
               "--ns", str(sim["length_ns"]), "--compute", sim.get("compute", "GPU"),
               "--ntomp", str(sim.get("threads", 16)), "--external-mpi", "--headless"]
    if sim["mode"] == "ligand":
        command += ["--ligand", setup["ligand"]]
        if setup.get("cofactors"):
            command += ["--cofactors", setup["cofactors"]]
    if resume:
        command += ["--resume", "--production_only"]
    return command


def analysis_command(config: dict[str, Any]) -> list[str]:
    sim, setup, analysis = config["simulation"], config["setup"], config.get("analysis", {})
    command = ["python", "3A_AutomateGromacs_chunks_MPI.py", "--gmx-bin", "gmx_mpi", "--mode", sim["mode"], "--headless",
               "--threads", str(analysis.get("threads", 8)), "--chunk-size", str(analysis.get("chunk_size", 1000))]
    if sim["mode"] == "ligand":
        command += ["--ligand", setup["ligand"], "--pocket-cutoff", str(analysis.get("pocket_cutoff", 5.0)),
                    "--contact_cutoff", str(analysis.get("contact_cutoff", 4.0)), "--min_contact_frac", str(analysis.get("min_contact_fraction", 0.10))]
        if analysis.get("compound_name"):
            command += ["--compound-name", analysis["compound_name"]]
    elif sim["mode"] in {"protein", "peptide"}:
        command += ["--interface-contact-cutoff", str(analysis.get("interface_contact_cutoff", 4.0)),
                    "--interface-min-contact-frac", str(analysis.get("interface_min_contact_fraction", 0.10)),
                    "--interface-frame-step", str(analysis.get("interface_frame_step", 1)),
                    "--interface-max-edges", str(analysis.get("interface_max_edges", 100))]
        if analysis.get("interface_chain_pairs"):
            command += ["--interface-chain-pairs", analysis["interface_chain_pairs"]]
    return command


def append_log(folder: Path, message: str) -> None:
    with (folder / "pymacs_commands.log").open("a", encoding="utf-8") as handle:
        handle.write(f"{dt.datetime.now().astimezone().isoformat(timespec='seconds')} {message}\n")


def execute_setup(folder: Path, config: dict[str, Any], dry_run: bool) -> None:
    command = setup_command(folder, config)
    display = shlex.join(command)
    print(display)
    if dry_run:
        return
    append_log(folder, f"setup: {display}")
    completed = subprocess.run(command, cwd=folder)
    if completed.returncode:
        raise SystemExit(completed.returncode)
    required = ["topol.top", "solv_ions.gro", "em.tpr"]
    missing = [name for name in required if not (folder / name).is_file()]
    if missing:
        raise ConfigError(f"Step 1 exited successfully but required outputs are missing: {', '.join(missing)}")
    print("Setup completed and required outputs are present.")


def load_profile(name: str) -> dict[str, Any]:
    path = ROOT / "hpc_profiles" / f"{name}.json"
    if not path.is_file():
        raise ConfigError(f"Unknown execution profile '{name}'. Available profiles are in {ROOT / 'hpc_profiles'}.")
    return json.loads(path.read_text(encoding="utf-8"))


def resolved_conda_root(profile: dict[str, Any]) -> str:
    """Resolve the submitting user's Conda path before LSF changes $HOME."""
    return os.path.expanduser(os.path.expandvars(str(profile["conda_root"])))


def render_lsf(folder: Path, config: dict[str, Any], profile_name: str, resume: bool) -> Path:
    profile = load_profile(profile_name)
    if profile.get("scheduler") != "lsf":
        raise ConfigError(f"Profile '{profile_name}' is not an LSF profile.")
    for name in ("topol.top", "solv_ions.gro", "em.tpr"):
        if not (folder / name).is_file():
            raise ConfigError(f"Cannot submit MD: required setup output is missing: {name}")
    template = (ROOT / "hpc_profiles" / profile["template"]).read_text(encoding="utf-8")
    job_stem = re.sub(r"[^a-z0-9._-]+", "_", folder.name.lower()).strip("_") or "pymacs_md"
    values = {
        "job_name": f"{job_stem}_{'resume' if resume else 'md'}",
        "workdir": str(folder.resolve()),
        "walltime": profile["walltime"],
        "threads": int(config["simulation"].get("threads", profile["threads"])),
        "gpus": profile["gpus"],
        "queue": profile["queue"],
        "project": profile.get("project", "brd"),
        "gromacs_module": profile["gromacs_module"],
        "gcc_module": profile["gcc_module"],
        "conda_root": resolved_conda_root(profile),
        "conda_env": profile["conda_env"],
        "command": shlex.join(simulation_command(config, resume=resume)),
    }
    destination = folder / f"run_{job_stem}_{'resume' if resume else 'md'}.lsf"
    destination.write_text(template.format(**values), encoding="utf-8")
    return destination


def submit_lsf(folder: Path, script: Path, dependency: str | None = None, dry_run: bool = False) -> str | None:
    if dry_run:
        print(script.read_text(encoding="utf-8"))
        return None
    bsub = shutil.which("bsub")
    if not bsub:
        raise ConfigError("bsub was not found. Submit from Triton or use --dry-run to inspect the generated script.")
    command = [bsub] + (["-w", dependency] if dependency else [])
    completed = subprocess.run(command, cwd=folder, input=script.read_text(encoding="utf-8"), text=True, capture_output=True)
    print(completed.stdout, end="")
    if completed.returncode:
        print(completed.stderr, file=sys.stderr, end="")
        raise SystemExit(completed.returncode)
    match = re.search(r"Job <(\d+)>", completed.stdout)
    return match.group(1) if match else None


def render_analysis_lsf(folder: Path, config: dict[str, Any], profile_name: str) -> Path:
    profile = load_profile(profile_name)
    template = (ROOT / "hpc_profiles" / profile["analysis_template"]).read_text(encoding="utf-8")
    job_stem = re.sub(r"[^a-z0-9._-]+", "_", folder.name.lower()).strip("_") or "pymacs_analysis"
    command = shlex.join(analysis_command(config))
    if config.get("analysis", {}).get("figurebook", True):
        command += " && python 4PDF4MD.py"
    values = {"job_name": f"{job_stem}_analysis", "workdir": str(folder.resolve()), "walltime": profile["analysis_walltime"],
              "threads": int(config.get("analysis", {}).get("threads", profile["analysis_threads"])), "queue": profile["queue"], "project": profile.get("project", "brd"),
              "gromacs_module": profile["gromacs_module"], "gcc_module": profile["gcc_module"],
              "conda_root": resolved_conda_root(profile), "conda_env": profile["conda_env"], "command": command}
    destination = folder / f"run_{job_stem}_analysis.lsf"
    destination.write_text(template.format(**values), encoding="utf-8")
    return destination


def submit_analysis(folder: Path, config: dict[str, Any], profile: str, after_job: str | None, dry_run: bool) -> None:
    if not config.get("analysis", {}).get("enabled", True):
        raise ConfigError("Analysis is disabled in pymacs_run.json. Reconfigure or enable analysis before submitting.")
    script = render_analysis_lsf(folder, config, profile)
    print(f"Generated {script.name}")
    jobid = submit_lsf(folder, script, f"done({after_job})" if after_job else None, dry_run)
    if jobid:
        append_log(folder, f"submitted analysis profile={profile} script={script.name} jobid={jobid} after={after_job or 'none'}")


def submit(folder: Path, config: dict[str, Any], profile: str, resume: bool, dry_run: bool, with_analysis: bool = False) -> None:
    script = render_lsf(folder, config, profile, resume)
    print(f"Generated {script.name}")
    jobid = submit_lsf(folder, script, dry_run=dry_run)
    if jobid:
        append_log(folder, f"submitted profile={profile} resume={resume} script={script.name} jobid={jobid}")
        if with_analysis:
            submit_analysis(folder, config, profile, jobid, dry_run=False)
    elif with_analysis and dry_run:
        print("\n--- Dependent analysis script ---")
        submit_analysis(folder, config, profile, "MD_JOBID", dry_run=True)


def status(folder: Path) -> None:
    log = folder / "pymacs_commands.log"
    print(f"Run folder: {folder}")
    print(f"Configuration: {'present' if config_path(folder).is_file() else 'missing'}")
    for name in ("topol.top", "solv_ions.gro", "em.tpr", "md_0_1.tpr", "md_0_1.cpt"):
        print(f"  {name}: {'present' if (folder / name).is_file() else 'missing'}")
    if log.is_file():
        print("\nRecent activity:")
        print("\n".join(log.read_text(encoding="utf-8").splitlines()[-8:]))


def main() -> None:
    parser = argparse.ArgumentParser(description="Configure and run one PyMACS simulation folder.")
    subparsers = parser.add_subparsers(dest="action", required=True)
    subparsers.add_parser("configure", help="Interactive Q&A that writes pymacs_run.json.")
    subparsers.add_parser("validate", help="Validate the saved configuration and input structure.")
    setup_parser = subparsers.add_parser("setup", help="Run PyMACS Step 1 from saved answers.")
    setup_parser.add_argument("--dry-run", action="store_true")
    for name, help_text in (("submit", "Render and submit a new Triton MD job."), ("resume", "Render and submit a Triton checkpoint-resume job.")):
        command = subparsers.add_parser(name, help=help_text)
        command.add_argument("--profile", default="triton")
        command.add_argument("--dry-run", action="store_true")
        if name == "submit":
            command.add_argument("--with-analysis", action="store_true", help="Submit chunked Step 3/4 after successful MD.")
    analysis_parser = subparsers.add_parser("submit-analysis", help="Render and submit the chunked Step 3/4 CPU job.")
    analysis_parser.add_argument("--profile", default="triton")
    analysis_parser.add_argument("--after", default=None, help="Run only after this LSF job succeeds.")
    analysis_parser.add_argument("--dry-run", action="store_true")
    subparsers.add_parser("status", help="Show local setup/checkpoint status.")
    args = parser.parse_args()
    folder = Path.cwd()
    try:
        if args.action == "configure":
            configure(folder)
            return
        config = load_config(folder)
        if args.action == "validate":
            print(f"Configuration is valid: {config_path(folder)}")
            print("\nNext: run PyMACS Step 1 setup:")
            print("  python pymacs_run.py setup")
            print("\nOn Triton, activate the cgenff environment and load gmx_mpi first.")
        elif args.action == "setup":
            execute_setup(folder, config, args.dry_run)
        elif args.action in {"submit", "resume"}:
            submit(folder, config, args.profile, args.action == "resume", args.dry_run, getattr(args, "with_analysis", False))
        elif args.action == "submit-analysis":
            submit_analysis(folder, config, args.profile, args.after, args.dry_run)
        else:
            status(folder)
    except ConfigError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc


if __name__ == "__main__":
    main()
