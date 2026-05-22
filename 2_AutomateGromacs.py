#!/usr/bin/env python3
"""
====================================================================
🏗️  Automated GROMACS MD Pipeline — Step 2: Equilibration + Production
====================================================================

Author:  Joseph-Michael Schulz (University of Miami)
Created: 2025-10-08
Version: 2.3 (GPU-Adaptive Build)
Dependencies:
    • GROMACS 2022+ (compiled with CUDA support)
    • Python 3.9+
    • MDAnalysis (optional, for binding pocket extraction)
    • NumPy
    • NVIDIA CUDA drivers installed and visible via `nvidia-smi`

--------------------------------------------------------------------
📖  OVERVIEW
--------------------------------------------------------------------
This script automates **equilibration and production molecular dynamics (MD)**
runs for protein, ligand, and PROTAC systems prepared in Step 1 of the workflow.

It dynamically configures:
• GPU and CPU threading
• MDP files (`tc-grps`, `nsteps`)
• Topology inclusion for ligands
• Ligand restraint generation
• Energy minimization → NVT → NPT → Production stages
• Optional post-run pocket extraction via MDAnalysis

The script supports three main simulation modes:
1️⃣  `protein` — single-protein systems (no ligand)
2️⃣  `ligand`  — small-molecule or peptide bound to a protein
3️⃣  `protac`  — ternary complex systems (E3 ligase + PROTAC + target)

--------------------------------------------------------------------
⚙️  KEY FEATURES
--------------------------------------------------------------------
• **Automatic GPU detection** (via `nvidia-smi`)
• **Thread scaling** via detected CPU cores
• **Dynamic tc-grps** updates to include ligand groups
• **Automatic nsteps scaling** (500,000 steps per ns)
• **Automatic topology patching** for ligand restraints
• **Energy minimization, NVT, NPT, and production MD** in series
• **Pocket builder** — extracts binding-site trajectories (optional)

--------------------------------------------------------------------
📦  REQUIRED INPUT FILES
--------------------------------------------------------------------
Expected to exist in the working directory (output of Step 1):
    topol.top
    em.mdp, nvt.mdp, npt.mdp, md.mdp
    protein_processed.gro / complex.gro
    ions.mdp
    forcefield directory (e.g., charmm36-mar2023.ff)
Optional (for ligand/protac):
    LIG.cgenff.mol2, LIG.str, LIG.gro, posre_LIG.itp, index_LIG.ndx

--------------------------------------------------------------------
🧠  USAGE
--------------------------------------------------------------------
Run directly from your system or cluster terminal after Step 1 setup:

    python 2_AutomateMD.py [options]

Examples:

▶ Protein-only simulation (default 50 ns):
    $ python 2_AutomateMD.py --mode protein

▶ Ligand-bound system (custom ligand and runtime):
    $ python 2_AutomateMD.py --mode ligand --ligand PTC --ns 25

▶ PROTAC ternary system (automatically includes atomIndex.txt):
    $ python 2_AutomateMD.py --mode protac --ligand PTX --ns 100

▶ Specify GPU manually (default = auto-detect):
    $ python 2_AutomateMD.py --mode ligand --ligand HDG --gpu 1

--------------------------------------------------------------------
🧩  OPTIONAL: BUILD POCKET TRAJECTORY
--------------------------------------------------------------------
After MD completion, extract binding-site residues interacting
with the ligand using:

    >>> from 2_AutomateMD import build_binding_pocket
    >>> build_binding_pocket("./", ligand_code="PTC", cutoff_ang=5.0)

Outputs:
    binding_pocket_only.pdb
    binding_pocket_only.xtc

--------------------------------------------------------------------
🧰  OUTPUTS
--------------------------------------------------------------------
Energy minimization:
    em.gro, em.log, em.edr
NVT/NPT equilibration:
    nvt.gro, npt.gro, nvt.cpt, npt.cpt, nvt.log, npt.log
Production MD:
    md_0_1.tpr, md_0_1.trr, md_0_1.xtc, md_0_1.log, md_0_1.edr
Optional (Ligand systems):
    posre_LIG.itp, index_LIG.ndx
Optional (Pocket extraction):
    binding_pocket_only.pdb, binding_pocket_only.xtc

--------------------------------------------------------------------
🔬  NOTES
--------------------------------------------------------------------
• Uses single precision and GPU-accelerated nonbonded + PME computations.
• Adjusts tc-grps automatically for two-temperature coupling groups:
    `Protein_LIGAND` and `Water_and_ions`
• Designed for reproducible multi-GPU use on HPC systems (e.g., Pegasus).
• All parameters editable via the MDP templates in the same directory.
• Will exit safely on any command failure with logged STDERR.

--------------------------------------------------------------------
🧭  FUTURE EXTENSIONS
--------------------------------------------------------------------
• Add adaptive sampling (REMD, REST2)
• Integrate with slurm/LSF submission scheduler
• Implement real-time energy drift monitoring
• Export JSON summaries for DockQ/PLIP analysis pipelines

====================================================================
"""


import os
import subprocess
import re
import multiprocessing
import argparse
import time
import multiprocessing, os, subprocess, time
import shutil
import os
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def parse_args():
    p = argparse.ArgumentParser(
        description="Automated GROMACS Step 2: MD equilibration + production run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--mode", choices=["ligand", "protein", "peptide", "protac"], required=False)
    p.add_argument("--ligand", type=str, help="3-letter ligand code (required for ligand or protac modes)")
    p.add_argument("--ns", type=float, default=None, help="Production length in nanoseconds (default 50)")
    p.add_argument("--gpu", "--gpu-id", dest="gpu", type=int, help="GPU ID to expose to GROMACS when running in GPU mode.")
    p.add_argument("--gpu-ids", type=str, default=None, help="Expert override for the exact gmx mdrun -gpu_id value (for example 0 or 01).")
    p.add_argument("--compute", choices=["CPU", "GPU"], default="CPU", help="Computation mode")
    p.add_argument("--no-gpu", action="store_true", help="Force CPU-only execution even if GPUs are available.")
    p.add_argument("--threads", type=int, default=None, help="Legacy alias for OpenMP thread count (same role as --ntomp).")
    p.add_argument("--ntomp", type=int, default=None, help="Number of OpenMP threads for gmx mdrun.")
    p.add_argument("--ntmpi", type=int, default=1, help="Number of MPI ranks passed to gmx mdrun.")
    p.add_argument("--pinoffset", type=int, default=None, help="Starting CPU core index for thread pinning")
    p.add_argument("--index-file", type=str, default=None, help="Use a user-supplied GROMACS index file instead of auto-generating index.ndx.")
    p.add_argument("--mdp-dir", type=str, default="MDPs", help="Directory searched for MDP templates when local copies are not present.")
    p.add_argument("--em-mdp", type=str, default="em.mdp", help="Energy-minimization MDP template.")
    p.add_argument("--nvt-mdp", type=str, default="nvt.mdp", help="NVT equilibration MDP template.")
    p.add_argument("--npt-mdp", type=str, default="npt.mdp", help="NPT equilibration MDP template.")
    p.add_argument("--md-mdp", type=str, default="md.mdp", help="Production MD MDP template.")
    p.add_argument(
        "--equilibration-plan",
        type=str,
        default=None,
        help="Optional JSON file describing a custom ordered equilibration plan. If omitted, PyMACS runs the default NVT then NPT sequence.",
    )
    p.add_argument("--headless", action="store_true",
                   help="Non-interactive mode; skip all prompts and use provided flags")
    p.add_argument("--production_only", action="store_true",
               help="Skip EM/NVT/NPT and run production only (resume if checkpoint exists).")
    p.add_argument("--resume", action="store_true",
                help="If md checkpoint exists, resume production from it (non-interactive).")
    p.add_argument("--force_restart", action="store_true",
                help="Do NOT resume even if checkpoint exists (start pipeline normally).")


    args = p.parse_args()
    if args.no_gpu:
        args.compute = "CPU"
        args.gpu = -1
    if args.ntomp is not None and args.ntomp <= 0:
        p.error("--ntomp must be a positive integer.")
    if args.ntmpi <= 0:
        p.error("--ntmpi must be a positive integer.")

    # Interactive fallback only when not headless
    if not args.headless and not args.mode:
        print("\nPlease select the simulation mode:")
        print("  1. ligand")
        print("  2. protein")
        print("  3. peptide")
        print("  4. protac")
        choice = input("Enter choice (1–4): ").strip()
        args.mode = {"1": "ligand", "2": "protein", "3": "peptide", "4": "protac"}.get(choice, "protein")
        print(f"🧠 Selected mode: {args.mode}")

    return args


def append_mdrun_log(directory, message):
    with open(os.path.join(directory, "mdrun.log"), "a", encoding="utf-8") as handle:
        handle.write(f"{message}\n")


def resolve_existing_path(base_directory, candidate):
    if not candidate:
        return None
    paths = [candidate]
    if not os.path.isabs(candidate):
        paths.extend([
            os.path.join(base_directory, candidate),
            os.path.join(SCRIPT_DIR, candidate),
        ])
    seen = set()
    for path in paths:
        abspath = os.path.abspath(path)
        if abspath in seen:
            continue
        seen.add(abspath)
        if os.path.exists(abspath):
            return abspath
    return None


def resolve_mdp_template(base_directory, mdp_name_or_path, mdp_dir):
    basename = os.path.basename(mdp_name_or_path) if mdp_name_or_path else None
    candidates = []
    if mdp_name_or_path:
        candidates.extend([
            mdp_name_or_path,
            os.path.join(base_directory, mdp_name_or_path),
            os.path.join(SCRIPT_DIR, mdp_name_or_path),
        ])
    if basename:
        candidates.extend([
            os.path.join(base_directory, basename),
            os.path.join(base_directory, mdp_dir, basename),
            os.path.join(SCRIPT_DIR, basename),
            os.path.join(SCRIPT_DIR, mdp_dir, basename),
        ])

    seen = set()
    for path in candidates:
        if not path:
            continue
        abspath = os.path.abspath(path)
        if abspath in seen:
            continue
        seen.add(abspath)
        if os.path.exists(abspath):
            return abspath
    return None


def prepare_local_mdp_template(directory, source_path, local_name):
    dest_path = os.path.join(directory, local_name)
    if os.path.abspath(source_path) != os.path.abspath(dest_path):
        shutil.copyfile(source_path, dest_path)
    return dest_path


def prepare_standard_mdp_files(directory, args):
    requested = {
        "em.mdp": args.em_mdp,
        "nvt.mdp": args.nvt_mdp,
        "npt.mdp": args.npt_mdp,
        "md.mdp": args.md_mdp,
    }
    local_paths = {}
    for local_name, requested_path in requested.items():
        resolved = resolve_mdp_template(directory, requested_path, args.mdp_dir)
        if not resolved:
            raise FileNotFoundError(
                f"Required MDP template '{requested_path}' for {local_name} was not found. "
                f"Check --mdp-dir or provide an explicit path."
            )
        local_paths[local_name] = prepare_local_mdp_template(directory, resolved, local_name)
        append_mdrun_log(directory, f"MDP {local_name}: source={resolved} local={local_paths[local_name]}")
    return local_paths


def load_equilibration_plan(plan_path, directory):
    if not plan_path:
        return None
    resolved = resolve_existing_path(directory, plan_path)
    if not resolved:
        raise FileNotFoundError(f"Equilibration plan JSON not found: {plan_path}")
    with open(resolved, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, list) or not payload:
        raise ValueError("Equilibration plan must be a non-empty JSON list of stage objects.")
    normalized = []
    for idx, item in enumerate(payload, start=1):
        if not isinstance(item, dict):
            raise ValueError(f"Equilibration plan entry {idx} must be a JSON object.")
        name = str(item.get("name", "")).strip()
        mdp = str(item.get("mdp", "")).strip()
        if not name or not mdp:
            raise ValueError(f"Equilibration plan entry {idx} must define both 'name' and 'mdp'.")
        normalized.append({"name": name, "mdp": mdp})
    return resolved, normalized


def system_has_virtual_sites(tpr_path: str) -> bool:
    """
    Detect whether the system contains virtual sites by inspecting `gmx dump`.
    """
    try:
        out = subprocess.check_output(
            f"gmx dump -s {tpr_path} | grep -i virtual",
            shell=True,
            text=True
        )
        return bool(out.strip())
    except subprocess.CalledProcessError:
        return False

GPU_PROFILES = [
    # 🚀 SAFE for ALL systems (virtual sites OK)
    {
        "name": "gpu_safe",
        "flags": "-nb gpu -pme gpu -bonded gpu -update cpu"
    },

    # 🚀 FULL GPU (ONLY if no virtual sites)
    {
        "name": "gpu_full",
        "flags": "-nb gpu -pme gpu -bonded gpu -update gpu"
    },
]



def run_mdrun_with_fallback(
    base_cmd,
    cwd,
    pin_offset,
    threads,
    gpu_id,
    use_gpu=True,
    tpr_check=None,
    ntmpi=1,
    gpu_id_arg=None,
):
    """
    Try best GPU profile first; if virtual sites exist, force CPU update.
    Fall back to CPU-only if all GPU attempts fail.
    """

    has_vs = False
    if tpr_check:
        tpr_path = os.path.join(cwd, tpr_check)
        if os.path.exists(tpr_path):
            has_vs = system_has_virtual_sites(tpr_path)
            if has_vs:
                print("⚠️ Virtual sites detected → forcing -update cpu (no GPU update).")

    # Profile order:
    #   - virtual sites: only safe
    #   - no virtual sites: full first, then safe
    if has_vs:
        profiles = [p for p in GPU_PROFILES if "-update cpu" in p["flags"]]
    else:
        profiles = []
        # prefer full if available
        for p in GPU_PROFILES:
            if "-update gpu" in p["flags"]:
                profiles.append(p)
        for p in GPU_PROFILES:
            if "-update cpu" in p["flags"]:
                profiles.append(p)

    # --------------------
    # GPU attempts
    # --------------------
    if use_gpu and gpu_id is not None and gpu_id >= 0:
        for profile in profiles:
            gpu_selector = gpu_id_arg if gpu_id_arg else "0"
            cmd = (
                f"{base_cmd} "
                f"-pin on -pinoffset {pin_offset} "
                f"-ntmpi {ntmpi} -ntomp {threads} "
                f"-gpu_id {gpu_selector} {profile['flags']}"
            )

            print(f"\n🚀 Trying {profile['name']} → {cmd}")
            append_mdrun_log(cwd, f"mdrun attempt ({profile['name']}): {cmd}")
            rc = run_command_check_rc(cmd, cwd=cwd)

            if rc == 0:
                print(f"✅ Success using {profile['name']}")
                append_mdrun_log(cwd, f"mdrun success ({profile['name']}): {cmd}")
                return

            print(f"⚠️ {profile['name']} failed, trying next fallback...")
            append_mdrun_log(cwd, f"mdrun failed ({profile['name']}): {cmd}")

    # --------------------
    # CPU fallback
    # --------------------
    print("\n🧯 Falling back to CPU-only execution")
    cpu_cmd = (
        f"{base_cmd} "
        f"-pin on -pinoffset {pin_offset} "
        f"-ntmpi {ntmpi} -ntomp {threads}"
    )

    append_mdrun_log(cwd, f"mdrun CPU fallback: {cpu_cmd}")
    rc = run_command_check_rc(cpu_cmd, cwd=cwd)
    if rc != 0:
        print("❌ CPU fallback also failed. Aborting.")
        exit(rc)

    print("✅ CPU-only execution successful")
    append_mdrun_log(cwd, f"mdrun CPU success: {cpu_cmd}")



def run_command(command, cwd=None, input_text=None):
    """Runs a shell command in the given directory and exits on failure."""
    print(f"\n🔹 Running command: {command}")
    result = subprocess.run(command, shell=True, cwd=cwd, text=True, capture_output=True, input=input_text)
    if result.stdout:
        print(f"📜 STDOUT:\n{result.stdout}")
    if result.returncode != 0:
        print(f"❌ Error running command: {command}")
        print(f"🔺 STDERR: {result.stderr}")
        exit(1)


import subprocess
import os


import shutil
import re
import subprocess
import os

def read_tpr_nsteps_dt(tpr_path: str):
    """
    Pull dt (ps) and nsteps from a TPR by streaming `gmx dump`,
    stopping as soon as we find both.
    """
    cmd = f"gmx dump -s {tpr_path}"
    proc = subprocess.Popen(cmd, shell=True, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

    dt_ps = None
    nsteps = None

    for line in proc.stdout:
        if dt_ps is None and "delta_t" in line:
            m = re.search(r"delta_t\s*=\s*([0-9.]+)", line)
            if m:
                dt_ps = float(m.group(1))

        if nsteps is None and "nsteps" in line:
            m = re.search(r"nsteps\s*=\s*(\d+)", line)
            if m:
                nsteps = int(m.group(1))

        if dt_ps is not None and nsteps is not None:
            proc.kill()
            break

    proc.wait()
    return nsteps, dt_ps


def tpr_total_ns(tpr_path: str):
    nsteps, dt_ps = read_tpr_nsteps_dt(tpr_path)
    if nsteps is None or dt_ps is None:
        return None
    # time = nsteps * dt (ps) -> convert ps to ns
    return (nsteps * dt_ps) / 1000.0


def extend_tpr_to_target_ns(directory: str, tpr_name: str, target_ns: float):
 
    tpr_path = os.path.join(directory, tpr_name)
    current_ns = tpr_total_ns(tpr_path)

    if current_ns is None:
        print("⚠️ Could not read current TPR length; skipping extension.")
        return

    if target_ns <= current_ns + 1e-6:
        print(f"ℹ️ TPR already set to ~{current_ns:.3f} ns (>= {target_ns} ns). No extension needed.")
        return

    extend_ns = target_ns - current_ns
    extend_ps = extend_ns * 1000.0

    tmp_tpr = tpr_path + ".tmp"
    print(f"🧩 Extending TPR: {current_ns:.3f} ns → {target_ns:.3f} ns (extend {extend_ns:.3f} ns)")

    run_command_cpu(
        f"gmx convert-tpr -s {tpr_name} -extend {extend_ps:.3f} -o {os.path.basename(tmp_tpr)}",
        cwd=directory
    )
    shutil.move(tmp_tpr, tpr_path)
    print("✅ TPR extension complete.")



def run_command_check_rc(command, cwd=None, input_text=None):
    """
    Special version of run_command:
    - DOES NOT exit on failure
    - Returns the command's exit code
    - Captures stdout/stderr for debugging
    """
    print(f"\n🔸 Running (check_rc): {command}")
    result = subprocess.run(
        command, shell=True, cwd=cwd, text=True,
        capture_output=True, input=input_text
    )
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr)
    return result.returncode  # <-- IMPORTANT



def load_chainmap_json(directory="."):
    """
    Load .chainmap.json if present.
    Returns dict like {"A": "H2B", "B": "H2A", "C": "IE1 Peptide"} or None.
    """
    path = os.path.join(directory, ".chainmap.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r") as f:
            data = json.load(f)
        if isinstance(data, dict) and data:
            return data
    except Exception as e:
        print(f"⚠️ Could not read {path}: {e}")
    return None

def parse_atomindex_ranges(directory="."):
    """
    Parse atomIndex.txt lines like:
      H2B 1-1485 # ...
      H2A 1486-3197 # ...
      IE1 Peptide 3198-3411 # ...

    Allows names with spaces.
    Returns:
      [{"name":"H2B","start":1,"end":1485,"length":1485}, ...]
    """
    path = os.path.join(directory, "atomIndex.txt")
    if not os.path.exists(path):
        return []

    entries = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            # Grab the first numeric range anywhere in the line
            m = re.search(r"(\d+)-(\d+)", s)
            if not m:
                continue

            start = int(m.group(1))
            end = int(m.group(2))

            # Everything before the range is the name
            name = s[:m.start()].strip()
            if not name:
                continue

            entries.append({
                "name": name,
                "start": start,
                "end": end,
                "length": end - start + 1
            })

    return entries

def detect_peptide_from_chainmap(directory="."):
    """
    Standalone peptide detector for Step 2 only.
    Uses:
      1) .chainmap.json
      2) atomIndex.txt
      3) shortest chain heuristic
      4) protein.pdb fallback if needed
    """
    chainmap = load_chainmap_json(directory)
    atom_entries = parse_atomindex_ranges(directory)

    if not chainmap:
        print("⚠️ No .chainmap.json found — peptide detection unavailable.")
        return None

    if atom_entries and len(chainmap) == len(atom_entries):
        chain_ids = list(chainmap.keys())

        combined = []
        for i, cid in enumerate(chain_ids):
            entry = atom_entries[i]
            combined.append({
                "chain_id": cid,
                "name": chainmap[cid],
                "atom_range": f"{entry['start']}-{entry['end']}",
                "length": entry["length"]
            })

        peptide_entry = min(combined, key=lambda x: x["length"])
        peptide_name = peptide_entry["name"].strip()
        peptide_group_name = re.sub(r"[^A-Za-z0-9_]", "", peptide_name.replace(" ", "_")) or "Peptide"

        return {
            "chain_id": peptide_entry["chain_id"],
            "name": peptide_name,
            "group_name": peptide_group_name,
            "atom_range": peptide_entry["atom_range"]
        }

    print("⚠️ atomIndex.txt parsing mismatch — falling back to protein.pdb residue counts...")

    protein_pdb = os.path.join(directory, "protein.pdb")
    if not os.path.exists(protein_pdb):
        print("⚠️ protein.pdb not found — cannot fall back.")
        return None

    from collections import OrderedDict
    chain_resids = OrderedDict()

    with open(protein_pdb, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain_id = line[21].strip() or "A"
            try:
                resid = int(line[22:26].strip())
            except ValueError:
                continue
            chain_resids.setdefault(chain_id, set()).add(resid)

    if not chain_resids:
        return None

    shortest_chain = min(chain_resids.items(), key=lambda kv: len(kv[1]))[0]
    peptide_name = chainmap.get(shortest_chain, "Peptide").strip()
    peptide_group_name = re.sub(r"[^A-Za-z0-9_]", "", peptide_name.replace(" ", "_")) or "Peptide"

    return {
        "chain_id": shortest_chain,
        "name": peptide_name,
        "group_name": peptide_group_name,
        "atom_range": None
    }

    

def save_peptide_chainmap(directory=".", peptide_info=None):
    """
    Save a simple peptide cache discovered entirely in Step 2.
    """
    if not peptide_info:
        return

    path = os.path.join(directory, ".chainmap.peptide.json")
    try:
        with open(path, "w") as f:
            json.dump(peptide_info, f, indent=2)
        print(f"💾 Saved peptide metadata → {path}")
    except Exception as e:
        print(f"⚠️ Could not save {path}: {e}")


def maybe_resume_production_only(directory, args, gpu_id, pin_offset, threads, use_gpu, target_ns,
                                 md_deffnm="md_0_1", ntmpi=1, gpu_id_arg=None):
    md_tpr = os.path.join(directory, f"{md_deffnm}.tpr")
    md_cpt = os.path.join(directory, f"{md_deffnm}.cpt")
    md_prev_cpt = os.path.join(directory, f"{md_deffnm}_prev.cpt")

    # Prefer md_0_1.cpt, fallback to md_0_1_prev.cpt if needed
    cpt_to_use = md_cpt if os.path.exists(md_cpt) else (md_prev_cpt if os.path.exists(md_prev_cpt) else None)

    if args.force_restart:
        return False

    if not (os.path.exists(md_tpr) and cpt_to_use):
        return False

    # Decide whether to resume
    if args.production_only or args.resume:
        do_resume = True
    elif args.headless:
        # sensible default for headless runs: resume if possible
        do_resume = True
    else:
        ans = input(
            f"\n♻️ Found production checkpoint ({os.path.basename(cpt_to_use)}).\n"
            f"Resume production ONLY (skip EM/NVT/NPT)? [Y/n]: "
        ).strip().lower()
        do_resume = (ans in ("", "y", "yes"))

    if not do_resume:
        return False

    print("\n🚀 Resuming PRODUCTION only...")
    # If user requested a longer run, extend the TPR so it doesn’t stop early
    if target_ns is not None:
        extend_tpr_to_target_ns(directory, f"{md_deffnm}.tpr", float(target_ns))

    base_cmd = f"gmx mdrun -v -deffnm {md_deffnm} -cpi {os.path.basename(cpt_to_use)} -append"

    # Reuse your existing GPU→GPU-lite→CPU fallback runner
    run_mdrun_with_fallback(
        base_cmd=base_cmd,
        cwd=directory,
        pin_offset=pin_offset,
        threads=threads,
        gpu_id=gpu_id,
        use_gpu=use_gpu,
        tpr_check=f"{md_deffnm}.tpr",
        ntmpi=ntmpi,
        gpu_id_arg=gpu_id_arg,
    )


    print("✅ Production resume complete.\n")
    return True



def run_command_cpu(command, cwd=None, input_text=None):
    """Run a shell command (CPU version). Exits on failure."""
    print(f"\n🔹 [CPU] Running command: {command}")
    result = subprocess.run(command, shell=True, cwd=cwd, text=True, capture_output=True, input=input_text)
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        print(f"❌ Error running command: {command}")
        print(f"🔺 STDERR: {result.stderr}")
        exit(result.returncode)


def run_command_gpu(command, cwd=None, input_text=None, gpu_id=0, threads=None):
    """Run a GPU-accelerated GROMACS command without overwriting the existing CUDA isolation."""
    threads = threads or os.cpu_count()

    # Do NOT reset CUDA_VISIBLE_DEVICES here — it is already set globally in environment setup
    os.environ["OMP_NUM_THREADS"] = str(threads)

    print(f"\n🚀 [GPU] Running on visible GPU 0 (physical GPU mapped via CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES')})")
    print(f"🔹 Command: {command}")

    process = subprocess.Popen(
        command,
        shell=True,
        cwd=cwd,
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        universal_newlines=True
    )
    if input_text:
        process.stdin.write(input_text)
        process.stdin.close()
    for line in process.stdout:
        print(line, end="")
    process.wait()

    if process.returncode != 0:
        print(f"❌ GPU command failed: {command}")
        exit(process.returncode)



def run_command_auto(command, cwd=None, input_text=None, use_gpu=True, gpu_id=0):
    if use_gpu:
        run_command_gpu(command, cwd=cwd, input_text=input_text, gpu_id=gpu_id)
    else:
        run_command_cpu(command, cwd=cwd, input_text=input_text)



def modify_topology_for_ligand(topol_path, ligand_code):
    """Ensures ligand position restraints are correctly added to topol.top (if ligand present)."""
    with open(topol_path, "r") as f:
        lines = f.readlines()
    if ligand_code:
        for i, line in enumerate(lines):
            if f'#include "{ligand_code.lower()}.itp"' in line:
                insert_index = i + 1
                break
        else:
            print(f"⚠️ `#include \"{ligand_code.lower()}.itp\"` not found in topol.top. Skipping modification.")
            return
        posres_block = [
            "; Ligand position restraints\n",
            "#ifdef POSRES\n",
            f'#include "posre_{ligand_code.lower()}.itp"\n',
            "#endif\n"
        ]
        if not any(f'#include "posre_{ligand_code.lower()}.itp"' in line for line in lines):
            lines[insert_index:insert_index] = posres_block
        with open(topol_path, "w") as f:
            f.writelines(lines)
        print(f"✅ Updated {topol_path} with ligand position restraints.")
    else:
        print("ℹ️ No ligand present; no topology modifications needed.")


def build_binding_pocket(directory: str, ligand_code: str, traj="Final_Trajectory.xtc",
                        topo_tpr="md_0_1.tpr", cutoff_ang=5.0,
                        out_pdb="binding_pocket_only.pdb", out_xtc="binding_pocket_only.xtc"):
    """
    Create pocket-only topology+trajectory containing protein residues within `cutoff_ang` Å of the ligand
    across the whole trajectory.
    """
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import distances
        import numpy as np

        traj_path = os.path.join(directory, traj)
        tpr_path  = os.path.join(directory, topo_tpr)

        if not os.path.exists(traj_path):
            print(f"❌ Pocket build: trajectory not found: {traj_path}")
            return
        if not ligand_code:
            print("ℹ️ Pocket build skipped (no ligand code provided).")
            return

        print(f"🔍 Pocket build: loading {topo_tpr} + {traj} ...")
        try:
            u = mda.Universe(tpr_path, traj_path)
        except Exception:
            print("⚠️ TPR not supported, using GRO fallback...")
            u = mda.Universe(os.path.join(directory, "md_0_1.gro"), traj_path)

        # ====== main logic ======True
        protein = u.select_atoms("protein")
        ligand  = u.select_atoms(f"resname {ligand_code}")

        if ligand.n_atoms == 0:
            print(f"❌ Pocket build: no atoms for ligand resname '{ligand_code}'. Skipping.")
            return

        print(f"📏 Finding protein residues within {cutoff_ang:.1f} Å of ligand {ligand_code} over all frames...")
        close_resids = set()
        for ts in u.trajectory:
            D = distances.distance_array(ligand.positions, protein.positions)  # Å
            for res in protein.residues:
                if np.any(D[:, res.atoms.indices] < cutoff_ang):
                    close_resids.add(res.resid)

        if not close_resids:
            print("⚠️ Pocket build: no contacting residues found; writing ligand-only pocket.")
            selection = f"resname {ligand_code}"
        else:
            res_sel = " or ".join([f"resid {r}" for r in sorted(close_resids)])
            selection = f"resname {ligand_code} or ({res_sel})"

        pocket_atoms = u.select_atoms(selection)
        print(f"🧬 Pocket atoms: {pocket_atoms.n_atoms}")

        # Write reduced trajectory
        out_xtc_path = os.path.join(directory, out_xtc)
        out_pdb_path = os.path.join(directory, out_pdb)

        print(f"💾 Writing pocket trajectory → {out_xtc_path}")
        with mda.Writer(out_xtc_path, pocket_atoms.n_atoms) as w:
            for ts in u.trajectory:
                w.write(pocket_atoms)

        # Write a single-frame PDB topology (first frame)
        u.trajectory[0]
        print(f"💾 Writing pocket topology → {out_pdb_path}")
        pocket_atoms.write(out_pdb_path)

        print("✅ Pocket build complete.")

    except ImportError as e:
        print(f"❌ Pocket build: missing dependency → {e}. Install MDAnalysis in this environment.")


def print_numa_summary():
    print("\n🧩 NUMA/GPU core mapping summary:")
    try:
        gpu_list = subprocess.check_output(
            "nvidia-smi --query-gpu=name --format=csv,noheader",
            shell=True, text=True
        ).strip().splitlines()
    except Exception:
        gpu_list = []
    for i, name in enumerate(gpu_list):
        cores = get_numa_cores_for_gpu(i)
        print(f"  GPU {i} ({name.strip()}): cores {cores[0]}–{cores[-1]} ({len(cores)} cores)")



def set_gpu_env(gpu_id, threads):
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    os.environ["OMP_NUM_THREADS"] = str(threads)

def clear_gpu_env():
    for k in [
        "CUDA_VISIBLE_DEVICES",
        "GMX_FORCE_GPU_ID",
        "GMX_GPU_DD_COMMS",
        "GMX_GPU_PME_DECOMPOSITION"
    ]:
        os.environ.pop(k, None)




def get_numa_cores_for_gpu(gpu_id: int):
    """
    Return a list of CPU cores local to a GPU (NUMA-aware if possible).

    Behavisor:
    • Uses `lscpu -e=CPU,NODE` to map physical cores to NUMA nodes.
    • If NUMA nodes ≥ number of GPUs → assign GPU N → Node N.
    • If fewer NUMA nodes than GPUs → split all cores evenly.
    • If NUMA info unavailable → fallback to even split across GPUs.
    """

    import subprocess, multiprocessing

    total_cores = multiprocessing.cpu_count()

    # --- Detect GPU count ---
    try:
        gpu_list = subprocess.check_output(
            "nvidia-smi --query-gpu=name --format=csv,noheader",
            shell=True, text=True
        ).strip().splitlines()
        gpu_count = len([g for g in gpu_list if g])
    except Exception:
        gpu_count = 1

    # --- Attempt NUMA-aware mapping ---
    try:
        output = subprocess.check_output("lscpu -e=CPU,NODE", shell=True, text=True)
        node_map = {}
        for line in output.splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 2 and parts[1] != "-":
                cpu, node = int(parts[0]), int(parts[1])
                node_map.setdefault(node, []).append(cpu)

        # Sort nodes for deterministic order
        node_ids = sorted(node_map.keys())

        # If NUMA nodes ≥ GPUs, map 1:1
        if gpu_id < len(node_ids):
            mapped_cores = node_map[node_ids[gpu_id]]
            if mapped_cores:
                return mapped_cores

        # Otherwise fall through to even-split fallback
    except Exception:
        pass  # will fallback below

    # --- Fallback: even split across GPUs ---
    per_gpu = max(1, total_cores // max(1, gpu_count))
    start = gpu_id * per_gpu
    end = min(start + per_gpu, total_cores)
    return list(range(start, end))


def set_env(gpu_id:int, threads:int, offset:int):
    """Prepare environment vars for proper GPU/CPU pinning."""
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    env["OMP_NUM_THREADS"] = str(threads)
    core_range = f"{offset}-{offset + threads - 1}"
    env["OMP_PLACES"] = f"{{{core_range}}}"
    env["OMP_PROC_BIND"] = "close"
    env["GOMP_CPU_AFFINITY"] = core_range
    env["GMX_GPU_DD_COMMS"] = "CUDA"
    env["GMX_GPU_PME_DECOMPOSITION"] = "CUDA"
    return env




def detect_ligand_from_cgenff():
    """
    Detect ligand name from *.cgenff.mol2 file.
    Returns ligand code (uppercase) or None.
    """
    mol2_files = [
        f for f in os.listdir(".")
        if f.lower().endswith(".cgenff.mol2")
    ]

    if not mol2_files:
        return None

    if len(mol2_files) > 1:
        print("⚠️ Multiple CGenFF MOL2 files detected:")
        for i, f in enumerate(mol2_files):
            print(f"  [{i}] {f}")
        try:
            idx = int(input("Select ligand MOL2 index [default 0]: ").strip() or "0")
            fname = mol2_files[idx]
        except Exception:
            fname = mol2_files[0]
    else:
        fname = mol2_files[0]

    ligand = fname.split(".cgenff.mol2")[0]
    return ligand.upper()


def index_has_group(index_path, group_name):
    """Return True if index.ndx contains a group with the exact given name."""
    if not os.path.exists(index_path):
        return False

    with open(index_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                name = line[1:-1].strip()
                if name == group_name:
                    return True
    return False


def standardize_mdp_groups(mdp_path, coupling_group_str):
    """
    Force tc-grps, comm-mode, and comm-grps to match the expected group layout.
    This prevents stale values from previous ligand/protein runs.
    """
    if not os.path.exists(mdp_path):
        print(f"❌ ERROR: Missing MDP file: {mdp_path}")
        exit(1)

    with open(mdp_path, "r") as f:
        lines = f.readlines()

    filtered = []
    found_tc = False
    found_comm_mode = False
    found_comm_grps = False

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("tc-grps"):
            filtered.append(f"tc-grps                 = {coupling_group_str}\n")
            found_tc = True
        elif stripped.startswith("comm-mode"):
            filtered.append("comm-mode               = Linear\n")
            found_comm_mode = True
        elif stripped.startswith("comm-grps"):
            filtered.append(f"comm-grps               = {coupling_group_str}\n")
            found_comm_grps = True
        else:
            filtered.append(line)

    if not found_tc:
        filtered.append(f"tc-grps                 = {coupling_group_str}\n")
    if not found_comm_mode:
        filtered.append("comm-mode               = Linear\n")
    if not found_comm_grps:
        filtered.append(f"comm-grps               = {coupling_group_str}\n")

    with open(mdp_path, "w") as f:
        f.writelines(filtered)


def set_md_nsteps(md_path, simulation_time_ns):
    if not os.path.exists(md_path):
        print(f"❌ ERROR: md.mdp not found: {md_path}")
        exit(1)

    nsteps = int(round(simulation_time_ns * 500000))
    with open(md_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    updated = []
    found_nsteps = False
    for line in lines:
        if line.strip().startswith("nsteps"):
            updated.append(f"nsteps                  = {nsteps}\n")
            found_nsteps = True
        else:
            updated.append(line)
    if not found_nsteps:
        updated.append(f"nsteps                  = {nsteps}\n")

    with open(md_path, "w", encoding="utf-8") as handle:
        handle.writelines(updated)
    return nsteps


def maybe_use_custom_index(directory, args):
    if not args.index_file:
        return False
    resolved = resolve_existing_path(directory, args.index_file)
    if not resolved:
        print(f"❌ User-supplied index file not found: {args.index_file}")
        exit(1)
    dest_path = os.path.join(directory, "index.ndx")
    backup_path = os.path.join(directory, "index.user_supplied.ndx")
    shutil.copyfile(resolved, dest_path)
    if os.path.abspath(resolved) != os.path.abspath(backup_path):
        shutil.copyfile(resolved, backup_path)
    print(f"📋 Using user-supplied index file: {resolved}")
    append_mdrun_log(directory, f"Index mode: user-supplied ({resolved})")
    return True


def run_equilibration_stage(stage_name, mdp_path, start_gro, ref_gro, prev_cpt, directory,
                            pin_offset, available_threads, gpu_id, use_gpu, ntmpi, gpu_id_arg,
                            coupling_group_str):
    local_name = f"{stage_name}.mdp"
    local_mdp = prepare_local_mdp_template(directory, mdp_path, local_name)
    standardize_mdp_groups(local_mdp, coupling_group_str)
    append_mdrun_log(directory, f"Equilibration stage {stage_name}: mdp={mdp_path} local={local_mdp}")

    grompp_cmd = (
        f"gmx grompp -f {local_name} -c {os.path.basename(start_gro)} "
        f"-r {os.path.basename(ref_gro)} "
    )
    if prev_cpt and os.path.exists(prev_cpt):
        grompp_cmd += f"-t {os.path.basename(prev_cpt)} "
    grompp_cmd += f"-p topol.top -n index.ndx -o {stage_name}.tpr -maxwarn 2"

    print(f"🧪 Running custom equilibration stage '{stage_name}'")
    append_mdrun_log(directory, f"grompp ({stage_name}): {grompp_cmd}")
    run_command_cpu(grompp_cmd, cwd=directory)
    run_mdrun_with_fallback(
        base_cmd=f"gmx mdrun -v -deffnm {stage_name}",
        cwd=directory,
        pin_offset=pin_offset,
        threads=available_threads,
        gpu_id=gpu_id,
        use_gpu=use_gpu,
        ntmpi=ntmpi,
        gpu_id_arg=gpu_id_arg,
    )
    return os.path.join(directory, f"{stage_name}.gro"), os.path.join(directory, f"{stage_name}.cpt")


def setup_md(directory, gpu_id, ligand_code, simulation_type, simulation_time_ns, use_gpu, args, peptide_info=None):
    """
    Runs full MD simulation setup on a specified GPU, with dynamic MDP configuration.
    Includes automatic tc-grps update, nsteps scaling, ligand topology handling,
    and indexed group management.
    """

    # ============================================================
    # ##0##  INITIAL ENVIRONMENT + GPU SELECTION + NUMA-AWARE THREADING
    # ============================================================
    print("\n🧩 [STEP 0] Configuring environment and CPU/GPU threading ...")

    available_cores = multiprocessing.cpu_count()
    print(f"🧠 Detected {available_cores} total CPU cores on system.")
    gpu_id_arg = args.gpu_ids if getattr(args, "gpu_ids", None) else ("0" if gpu_id is not None and gpu_id >= 0 else None)

    # --- Detect GPUs ---
    try:
        gpu_list = subprocess.check_output(
            "nvidia-smi --query-gpu=name --format=csv,noheader",
            shell=True, text=True
        ).strip().splitlines()
        gpu_list = [g.strip() for g in gpu_list if g.strip()]
        gpu_count = len(gpu_list)
    except Exception:
        gpu_list, gpu_count = [], 0

    # --- Interactive GPU selection ---
    headless_mode = getattr(args, "headless", False)
    if not headless_mode and gpu_count > 0:
        print(f"\n🧠 Detected {gpu_count} GPU(s):")
        for i, name in enumerate(gpu_list):
            print(f"  [{i}] {name}")
        gpu_choice = input(f"\nSelect GPU number (0–{gpu_count-1}) [default 0]: ").strip()
        gpu_id = int(gpu_choice) if gpu_choice.isdigit() and int(gpu_choice) < gpu_count else 0
        print(f"🎯 Selected GPU {gpu_id}: {gpu_list[gpu_id]}")
    elif gpu_count > 0:
        gpu_id = args.gpu if getattr(args, "gpu", None) is not None else 0
        print(f"🤖 Headless mode: using GPU {gpu_id} ({gpu_list[gpu_id]})")
    else:
        gpu_id = -1
        print("⚠️ No GPUs detected — running in CPU mode.")

    # --- NUMA-aware core detection ---
    try:
        cores = get_numa_cores_for_gpu(gpu_id) if gpu_id >= 0 else list(range(available_cores))
    except Exception:
        cores = list(range(available_cores))

    if not cores:
        per_gpu = max(1, available_cores // max(1, gpu_count or 1))
        start = gpu_id * per_gpu
        end = min(start + per_gpu, available_cores)
        cores = list(range(start, end))

    suggested_threads = len(cores)
    suggested_offset = cores[0]
    print(f"📊 Suggested NUMA mapping for GPU {gpu_id}: cores {suggested_offset}–{cores[-1]} ({suggested_threads} cores)")

    # --- Interactive override ---
    if not headless_mode:
        accept = input(f"Use this mapping for GPU {gpu_id}? [Y/n]: ").strip().lower()
        if accept in ["n", "no"]:
            try:
                threads_in = input(f"Enter number of threads [default {suggested_threads}]: ").strip()
                threads = int(threads_in) if threads_in else suggested_threads
                offset_in = input(f"Enter pin offset [default {suggested_offset}]: ").strip()
                pin_offset = int(offset_in) if offset_in else suggested_offset
            except Exception:
                print("⚠️ Invalid input, reverting to suggested values.")
                threads, pin_offset = suggested_threads, suggested_offset
        else:
            threads, pin_offset = suggested_threads, suggested_offset
    else:
        threads, pin_offset = suggested_threads, suggested_offset
        print(f"🤖 Headless mode: auto-accepted NUMA mapping → threads={threads}, pin offset={pin_offset}")

    # --- Apply CLI overrides ---
    if getattr(args, "ntomp", None):
        threads = int(args.ntomp)
    elif getattr(args, "threads", None):
        threads = int(args.threads)
    if getattr(args, "pinoffset", None):
        pin_offset = int(args.pinoffset)

    available_threads = max(1, threads)

    # --- Environment setup (GPU isolation like 2A_AutoGMXrestart.py) ---
    if gpu_id is not None and gpu_id >= 0:
        # Map the chosen physical GPU (e.g., 1 or 2) into a single visible device (index 0)
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        if args.gpu_ids:
            os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_ids)
            os.environ.pop("GMX_FORCE_GPU_ID", None)
            print(f"🎯 Using user-specified GPU visibility string: {args.gpu_ids}")
        else:
            os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
            os.environ["GMX_FORCE_GPU_ID"] = "0"  # GROMACS now always sees it as device 0
            print(f"🎯 Forcing isolation: physical GPU {gpu_id} → visible GPU 0")
    else:
        # No GPU explicitly chosen: expose all GPUs
        try:
            visible = subprocess.check_output(
                "nvidia-smi --query-gpu=index --format=csv,noheader",
                shell=True, text=True
            ).strip().replace("\n", ",")
            os.environ["CUDA_VISIBLE_DEVICES"] = visible
            print(f"💻 No GPU ID specified; exposing all GPUs: {visible}")
        except Exception:
            os.environ["CUDA_VISIBLE_DEVICES"] = "0"
            print("⚠️ Could not query nvidia-smi; defaulting to GPU 0")

    # Set thread binding
    os.environ["OMP_NUM_THREADS"] = str(available_threads)



    print(f"\n🧠 Environment configured:")
    print(f"   • Detected total cores: {available_cores}")
    if gpu_id >= 0:
        print(f"   • GPU {gpu_id}: {gpu_list[gpu_id]} (cores {suggested_offset}–{cores[-1]})")
    else:
        print("   • CPU-only mode active.")
    print(f"   • Using {available_threads} OpenMP thread(s), pin offset {pin_offset}")
    print(f"   • Using {args.ntmpi} MPI rank(s)")
    print(f"   • CUDA_VISIBLE_DEVICES={gpu_id if gpu_id >= 0 else 'None'}")
    print(f"   • OMP_NUM_THREADS={available_threads}")
    if gpu_id_arg:
        print(f"   • gmx mdrun -gpu_id {gpu_id_arg}")

    append_mdrun_log(directory, f"Compute mode: {'GPU' if use_gpu else 'CPU'}")
    append_mdrun_log(directory, f"Resources: ntmpi={args.ntmpi} ntomp={available_threads} pinoffset={pin_offset} gpu_id={gpu_id} gpu_id_arg={gpu_id_arg}")




    # after your Step 0 environment config (once pin_offset and available_threads are known)
    if maybe_resume_production_only(
        directory=directory,
        args=args,
        gpu_id=gpu_id,
        pin_offset=pin_offset,
        threads=available_threads,
        use_gpu=use_gpu,
        target_ns=simulation_time_ns,
        md_deffnm="md_0_1",
        ntmpi=args.ntmpi,
        gpu_id_arg=gpu_id_arg,
    ):
        return  # <- IMPORTANT: skip steps 1–7 entirely


    # ============================================================
    # 🧩 1. Update MDP Files (tc-grps + comm-grps + nsteps)
    # ============================================================
    print("\n🧩 [STEP 1] Updating MDP files (tc-grps, comm-grps, nsteps) ...")

    try:
        mdp_files = prepare_standard_mdp_files(directory, args)
    except FileNotFoundError as exc:
        print(f"❌ {exc}")
        exit(1)

    # Decide the exact thermostat / COM groups based on system type
    if ligand_code:
        coupling_group_str = f"Protein_{ligand_code} Water_and_ions"
    else:
        coupling_group_str = "Protein Water_and_ions"

    for mdp in ["md.mdp", "nvt.mdp", "npt.mdp"]:
        standardize_mdp_groups(mdp_files[mdp], coupling_group_str)

    print(f"🛠️ Standardized tc-grps/comm-grps to: {coupling_group_str}")

    nsteps = set_md_nsteps(mdp_files["md.mdp"], simulation_time_ns)

    print(f"🧮 Set nsteps = {nsteps} (~{simulation_time_ns} ns).")
    print("✅ Finished MDP standardization.\n")

    # ============================================================
    # 🧱 2. Topology Preparation
    # ============================================================
    print("🧱 [STEP 2] Preparing topology ...")
    topol_path = os.path.join(directory, "topol.top")
    if ligand_code:
        print(f"🔧 Modifying topology to include ligand {ligand_code} ...")
        modify_topology_for_ligand(topol_path, ligand_code)
        print(f"✅ Ligand '{ligand_code}' integrated into topology.")
    else:
        print("ℹ️ No ligand present; skipping topology modification.")

    # ============================================================
    # ⚡ 3. Energy Minimization
    # ============================================================
    print("\n⚡ [STEP 3] Starting energy minimization ...")

    run_mdrun_with_fallback(
        base_cmd="gmx mdrun -v -deffnm em",
        cwd=directory,
        pin_offset=pin_offset,
        threads=available_threads,
        gpu_id=gpu_id,
        use_gpu=use_gpu,
        ntmpi=args.ntmpi,
        gpu_id_arg=gpu_id_arg,
    )

    print("✅ Energy minimization complete.\n")


    


    # ============================================================
    # 🧬 4. Ligand Restraints (generate if missing)
    # ============================================================
    if ligand_code:
        posre_file = f"posre_{ligand_code.lower()}.itp"
        ligand_gro = f"{ligand_code.lower()}.gro"
        index_file = f"index_{ligand_code.lower()}.ndx"

        if not os.path.exists(os.path.join(directory, posre_file)):
            print(f"🧪 Generating position restraints for ligand: {ligand_code}")

            if os.path.exists(os.path.join(directory, ligand_gro)):

                # ------------------------------
                # TRY 1 →  H* wildcard``
                # TRY 2 →  explicit H
                # TRY 3 →  keep full ligand
                # ------------------------------

                cmd1 = (
                    f"printf \"0 & ! a H*\\nq\\n\" | "
                    f"gmx make_ndx -f {ligand_gro} -o {index_file}"
                )

                cmd2 = (
                    f"printf \"0 & ! a H\\nq\\n\" | "
                    f"gmx make_ndx -f {ligand_gro} -o {index_file}"
                )

                cmd3 = (
                    f"printf \"0\\nq\\n\" | "
                    f"gmx make_ndx -f {ligand_gro} -o {index_file}"
                )

                print("🔹 Attempt 1: Removing hydrogens using wildcard H* ...")
                rc1 = run_command_check_rc(cmd1, cwd=directory)

                if rc1 != 0:
                    print("⚠️ Attempt 1 failed — trying explicit hydrogen selector...")

                    rc2 = run_command_check_rc(cmd2, cwd=directory)

                    if rc2 != 0:
                        print("⚠️ Attempt 2 failed — using full ligand group...")

                        rc3 = run_command_check_rc(cmd3, cwd=directory)

                        if rc3 != 0:
                            print("❌ All index generation attempts failed.")
                            index_file = None


                # If we successfully generated an index → create restraints
                if index_file and os.path.exists(os.path.join(directory, index_file)):
                    print("🔧 Creating ligand position restraints (genrester)...")
                    run_command(
                        f"echo '3' | gmx genrestr -f {ligand_gro} -n {index_file} "
                        f"-o {posre_file} -fc 1000 1000 1000",
                        cwd=directory
                    )
                    print(f"✅ Created {posre_file}")

            else:
                print(f"⚠️ Ligand GRO file not found ({ligand_gro}). Skipping posre generation.")
        else:
            print(f"✅ Using existing ligand position restraints: {posre_file}")






    # ============================================================
    # 🧩 5. Build Master Index File (Protein / Ligand / PROTAC)
    # ============================================================
    print("\n🧩 [STEP 5] Building clean index.ndx ...")

    custom_index_in_use = maybe_use_custom_index(directory, args)

    # ------------------------------------------------------------
    # 5A — Clean old index files
    # ------------------------------------------------------------
    cleanup_files = ["index.ndx", "index_default.ndx"]
    if ligand_code:
        cleanup_files.append(f"index_{ligand_code.lower()}.ndx")

    if not custom_index_in_use:
        for f in cleanup_files:
            fp = os.path.join(directory, f)
            if os.path.exists(fp):
                os.remove(fp)

    # ------------------------------------------------------------
    # 5B — Ensure ligand restraint file exists if ligand is present
    # ------------------------------------------------------------
    if ligand_code:
        lig = ligand_code.upper()
        posre = os.path.join(directory, f"posre_{lig.lower()}.itp")
        lig_gro = os.path.join(directory, f"{lig.lower()}.gro")

        if not os.path.exists(posre):
            print(f"⚠️ posre_{lig.lower()}.itp missing — generating fallback restraints...")

            if not os.path.exists(lig_gro):
                print(f"❌ FATAL: {lig_gro} not found — cannot create restraints.")
                exit(1)

            print(f"🔍 Scanning {lig_gro} for heavy atoms...")
            heavy_atoms = []

            with open(lig_gro, "r") as g:
                for line in g:
                    if len(line) < 20:
                        continue
                    resname = line[5:10].strip()
                    if resname != lig:
                        continue
                    atom_name = line[10:15].strip()
                    atom_id = int(line[15:20])
                    if not atom_name.startswith("H"):
                        heavy_atoms.append(atom_id)

            print(f"🔎 Found {len(heavy_atoms)} heavy atoms")

            if heavy_atoms:
                with open(posre, "w") as out:
                    out.write("[ position_restraints ]\n")
                    out.write("; atom  type    fx    fy    fz\n")
                    for a in heavy_atoms:
                        out.write(f"{a:5d}   1   1000   1000   1000\n")
                print(f"✅ Wrote fallback restraint file: {posre}")
            else:
                print("⚠️ No heavy atoms found — posre will be empty.")
        else:
            print(f"✔ Using existing ligand restraints: {posre}")

    # ------------------------------------------------------------
    # 5C — MODE-SPECIFIC INDEX BUILD
    # ------------------------------------------------------------
    simtype = simulation_type.upper().strip()

    if custom_index_in_use:
        if (not ligand_code) and (simtype != "PROTAC"):
            expected_group = "Protein"
            expected_water_group = "Water_and_ions"
        else:
            expected_group = f"Protein_{ligand_code.upper()}"
            expected_water_group = "Water_and_ions"
        print("📋 Skipping index auto-generation because --index-file was provided.")
    elif (not ligand_code) and (simtype != "PROTAC"):
        # ========================================================
        # PROTEIN-ONLY MODE
        # ========================================================
        print("🧬 Protein-only system → building default index.ndx")

        run_command(
            "gmx make_ndx -f em.gro -o index.ndx",
            cwd=directory,
            input_text="q\n"
        )

        expected_group = "Protein"
        expected_water_group = "Water_and_ions"
        print("✅ index.ndx built (protein only)")

    elif simtype == "PROTAC":
        # ========================================================
        # PROTAC MODE
        # ========================================================
        print("🧬 PROTAC system → building custom multi-group index.ndx")

        if not ligand_code:
            print("❌ FATAL: PROTAC mode requires a ligand_code.")
            exit(1)

        atom_file = os.path.join(directory, "atomIndex.txt")
        if not os.path.exists(atom_file):
            print("❌ FATAL: atomIndex.txt missing for PROTAC mode.")
            exit(1)

        run_command(
            "gmx make_ndx -f em.gro -o index.ndx",
            cwd=directory,
            input_text=(
                "1 | 13\n"
                f"name 21 Protein_{ligand_code.upper()}\n"
                + "".join(
                    f"{line.strip()}\nname {22+i} "
                    f"{['Ligase','Target_A','Target_B','Target_C','Target_D'][i] if i < 5 else f'Target_{i}'}\n"
                    for i, line in enumerate(open(atom_file))
                    if line.strip()
                )
                + "q\n"
            )
        )

        expected_group = f"Protein_{ligand_code.upper()}"
        expected_water_group = "Water_and_ions"
        print(f"✅ PROTAC index.ndx built with merged group: {expected_group}")

    elif ligand_code:
        # ========================================================
        # LIGAND MODE
        # KEEP OLD LOGIC THAT PREVIOUSLY WORKED
        # ========================================================
        expected_group = f"Protein_{ligand_code.upper()}"
        expected_water_group = "Water_and_ions"

        print(f"💊 Ligand-bound system → building merged index group {expected_group}")

        ndx_input = (
            "1 | 13\n"
            f"name 21 {expected_group}\n"
            "q\n"
        )

        run_command(
            "gmx make_ndx -f em.gro -o index.ndx",
            cwd=directory,
            input_text=ndx_input
        )

        print(f"✅ index.ndx built cleanly: {expected_group}")
        print("   ✔ Old ligand merge preserved")
        print("   ✔ Protein-only logic kept separate")
        print("   ✔ tc-grps can safely use this group")

    else:
        print("❌ FATAL: Could not determine simulation mode for Step 5.")
        print(f"   ligand_code={ligand_code}")
        print(f"   simulation_type={simulation_type}")
        exit(1)

    # ------------------------------------------------------------
    # 5D — Final validation
    # ------------------------------------------------------------
    index_path = os.path.join(directory, "index.ndx")
    if not os.path.exists(index_path):
        print("❌ FATAL: index.ndx was not created.")
        exit(1)

    missing_groups = []

    if not index_has_group(index_path, expected_group):
        missing_groups.append(expected_group)

    if not index_has_group(index_path, expected_water_group):
        missing_groups.append(expected_water_group)

    if missing_groups:
        print(f"❌ FATAL: Missing expected index group(s): {', '.join(missing_groups)}")
        print("📋 Groups present in index.ndx:")
        with open(index_path, "r") as f:
            for line in f:
                s = line.strip()
                if s.startswith("[") and s.endswith("]"):
                    print(f"   {s}")
        exit(1)

    print(f"✅ Verified Step 5 index groups exist: {expected_group}, {expected_water_group}\n")

    # ============================================================
    # 🌡️ 6–7. Equilibration
    # ============================================================
    if ligand_code:
        expected_group = f"Protein_{ligand_code.upper()}"
        coupling_group_str = f"{expected_group} Water_and_ions"
    else:
        expected_group = "Protein"
        coupling_group_str = "Protein Water_and_ions"

    print(f"🧭 Expected equilibration coupling/index group: {expected_group}")

    index_path = os.path.join(directory, "index.ndx")
    if not os.path.exists(index_path):
        print("❌ FATAL: index.ndx not found before equilibration grompp.")
        exit(1)

    if not index_has_group(index_path, expected_group):
        print(f"❌ FATAL: Expected index group '{expected_group}' not found in index.ndx")
        print("   The equilibration MDP file and index file are out of sync.")
        print("\n📋 Groups found in index.ndx:")
        with open(index_path, "r") as f:
            for line in f:
                s = line.strip()
                if s.startswith("[") and s.endswith("]"):
                    print(f"   {s}")
        exit(1)

    if not index_has_group(index_path, "Water_and_ions"):
        print("❌ FATAL: Expected index group 'Water_and_ions' not found in index.ndx")
        print("   The default solvent/ion group is missing from the generated or supplied index.")
        print("\n📋 Groups found in index.ndx:")
        with open(index_path, "r") as f:
            for line in f:
                s = line.strip()
                if s.startswith("[") and s.endswith("]"):
                    print(f"   {s}")
        exit(1)

    print(f"✅ Verified index groups exist: {expected_group}, Water_and_ions")

    equilibration_start_gro = os.path.join(directory, "em.gro")
    equilibration_ref_gro = os.path.join(directory, "em.gro")
    equilibration_cpt = None

    plan_info = None
    if args.equilibration_plan:
        try:
            plan_info = load_equilibration_plan(args.equilibration_plan, directory)
        except Exception as exc:
            print(f"❌ Invalid equilibration plan: {exc}")
            exit(1)

    if plan_info:
        plan_path, plan_rows = plan_info
        print(f"🧭 [STEP 6] Running custom equilibration plan from {plan_path}")
        append_mdrun_log(directory, f"Equilibration plan: {plan_path}")
        for item in plan_rows:
            resolved_stage_mdp = resolve_mdp_template(directory, item["mdp"], args.mdp_dir)
            if not resolved_stage_mdp:
                print(f"❌ Could not locate equilibration stage MDP: {item['mdp']}")
                exit(1)
            equilibration_start_gro, equilibration_cpt = run_equilibration_stage(
                stage_name=item["name"],
                mdp_path=resolved_stage_mdp,
                start_gro=equilibration_start_gro,
                ref_gro=equilibration_ref_gro,
                prev_cpt=equilibration_cpt,
                directory=directory,
                pin_offset=pin_offset,
                available_threads=available_threads,
                gpu_id=gpu_id,
                use_gpu=use_gpu,
                ntmpi=args.ntmpi,
                gpu_id_arg=gpu_id_arg,
                coupling_group_str=coupling_group_str,
            )
        print("✅ Custom equilibration plan complete.\n")
    else:
        print("🌡️ [STEP 6] Starting NVT equilibration ...")
        append_mdrun_log(directory, "grompp (nvt): gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2")
        run_command_cpu(
            "gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2",
            cwd=directory
        )
        run_mdrun_with_fallback(
            base_cmd="gmx mdrun -v -deffnm nvt",
            cwd=directory,
            pin_offset=pin_offset,
            threads=available_threads,
            gpu_id=gpu_id,
            use_gpu=use_gpu,
            ntmpi=args.ntmpi,
            gpu_id_arg=gpu_id_arg,
        )
        print("✅ NVT equilibration complete.\n")

        print("💧 [STEP 7] Starting NPT equilibration ...")
        append_mdrun_log(directory, "grompp (npt): gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2")
        run_command_cpu(
            "gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2",
            cwd=directory
        )
        run_mdrun_with_fallback(
            base_cmd="gmx mdrun -v -deffnm npt",
            cwd=directory,
            pin_offset=pin_offset,
            threads=available_threads,
            gpu_id=gpu_id,
            use_gpu=use_gpu,
            ntmpi=args.ntmpi,
            gpu_id_arg=gpu_id_arg,
        )
        equilibration_start_gro = os.path.join(directory, "npt.gro")
        equilibration_cpt = os.path.join(directory, "npt.cpt")
        print("✅ NPT equilibration complete.\n")


    # ============================================================
    # 🧠 8. Production MD (resume from checkpoint if present)
    # ============================================================
    print("🧠 [STEP 8] Starting production MD ...")

    md_deffnm = "md_0_1"
    md_tpr = os.path.join(directory, f"{md_deffnm}.tpr")
    md_cpt = os.path.join(directory, f"{md_deffnm}.cpt")

    # build / rebuild TPR if needed
    if not (os.path.exists(md_tpr) and os.path.exists(md_cpt)):
        append_mdrun_log(
            directory,
            f"grompp (production): gmx grompp -f md.mdp -c {os.path.basename(equilibration_start_gro)} "
            + (f"-t {os.path.basename(equilibration_cpt)} " if equilibration_cpt and os.path.exists(equilibration_cpt) else "")
            + f"-p topol.top -n index.ndx -o {md_deffnm}.tpr -maxwarn 2"
        )
        run_command_cpu(
            f"gmx grompp -f md.mdp -c {os.path.basename(equilibration_start_gro)} "
            + (f"-t {os.path.basename(equilibration_cpt)} " if equilibration_cpt and os.path.exists(equilibration_cpt) else "")
            + "-p topol.top -n index.ndx "
            f"-o {md_deffnm}.tpr -maxwarn 2",
            cwd=directory
        )

    # base mdrun (resume if checkpoint exists)
    base_cmd = f"gmx mdrun -v -deffnm {md_deffnm}"
    if os.path.exists(md_cpt) and os.path.exists(md_tpr):
        print(f"♻️ Found checkpoint: {md_deffnm}.cpt → resuming with -cpi/-append")
        base_cmd += f" -cpi {md_deffnm}.cpt -append"

    # GPU→fallback runner (virtual-site aware)
    run_mdrun_with_fallback(
        base_cmd=base_cmd,
        cwd=directory,
        pin_offset=pin_offset,
        threads=available_threads,
        gpu_id=gpu_id,
        use_gpu=use_gpu,
        tpr_check=f"{md_deffnm}.tpr",
        ntmpi=args.ntmpi,
        gpu_id_arg=gpu_id_arg,
    )

    print("✅ Production MD complete.\n")







if __name__ == "__main__":
    args = parse_args()
    print_numa_summary()

    # GPU selection
    try:
        gpu_list = subprocess.check_output(
            "nvidia-smi --query-gpu=name --format=csv,noheader", shell=True
        ).decode().strip().split('\n')
        gpu_list = [g for g in gpu_list if g]
        gpu_count = len(gpu_list)
    except Exception:
        gpu_list = []
        gpu_count = 0

    if args.compute.upper() == "GPU":
        if args.gpu is not None:
            gpu_id = args.gpu
        elif args.headless:
            # Headless GPU but no gpu id given -> default 0 if available, else fail to CPU
            gpu_id = 0 if gpu_count > 0 else -1
        else:
            if gpu_count == 0:
                print("❌ No GPU detected. Ensure NVIDIA drivers and CUDA are properly installed.")
                exit(1)
            elif gpu_count == 1:
                print("\n🎯 Detected 1 GPU. Automatically selecting GPU 0.")
                gpu_id = 0
            else:
                print(f"\n🧠 Detected {gpu_count} GPUs:")
                for idx, name in enumerate(gpu_list):
                    print(f"  [{idx}] {name}")
                gpu_choice = input(f"\nEnter GPU number (0-{gpu_count - 1}): ").strip()
                if not gpu_choice.isdigit() or int(gpu_choice) not in range(gpu_count):
                    print("❌ Invalid GPU selection.")
                    exit(1)
                gpu_id = int(gpu_choice)
    else:
        gpu_id = -1  # CPU mode


    simulation_directory = os.getcwd()
    print(f"\n📂 Running MD setup in: {simulation_directory} on GPU {gpu_id}")

    # Simulation type
    sim_map = {
        "1": "Ligand:Protein",
        "2": "Protein:Protein",
        "3": "Peptide:Protein",
        "4": "PROTAC",
        "5": "Just Protein"
    }

    if args.mode:
        mode_map = {
            "ligand": "Ligand:Protein",
            "protein": "Just Protein",
            "peptide": "Peptide:Protein",
            "protac": "PROTAC"
        }
        simulation_type = mode_map.get(args.mode.lower(), args.mode)

    else:
        print("\nPlease select the simulation type:")
        print("  1. Ligand:Protein")
        print("  2. Protein:Protein")
        print("  3. Peptide:Protein")
        print("  4. PROTAC")
        print("  5. Just Protein")
        sim_choice = input("Enter choice (1-5): ").strip()
        if sim_choice not in ["1", "2", "3", "4", "5"]:
            print("❌ Invalid selection. Exiting.")
            exit(1)
        simulation_type = sim_map[sim_choice]
    
    peptide_info = None

    if simulation_type.lower() in ["peptide:protein", "peptide"]:
        peptide_info = detect_peptide_from_chainmap(simulation_directory)

        if peptide_info:
            print("\n🧬 Detected peptide from existing chain map:")
            print(f"   • Chain ID:   {peptide_info['chain_id']}")
            print(f"   • Name:       {peptide_info['name']}")
            print(f"   • Atom range: {peptide_info['atom_range']}")
            print(f"   • Group name: {peptide_info['group_name']}")

            if not args.headless:
                ans = input("Use this peptide assignment? [Y/n]: ").strip().lower()
                if ans in ["n", "no"]:
                    print("❌ Peptide mode requires a valid peptide assignment from .chainmap.json + atomIndex.txt")
                    exit(1)

            save_peptide_chainmap(simulation_directory, peptide_info)

        else:
            print("❌ Could not detect peptide automatically.")
            print("   Peptide mode requires:")
            print("   • .chainmap.json")
            print("   • atomIndex.txt")
            print("   • shortest chain logic to identify peptide")
            exit(1)

    # ─── Ligand code (auto-detect + confirm) ──────────────────
    ligand_code = None

    if simulation_type.lower() in ["ligand:protein", "ligand", "protac"]:

        # CLI always wins
        if args.ligand:
            ligand_code = args.ligand.upper()
            print(f"🧬 Using ligand from CLI: {ligand_code}")

        else:
            detected = detect_ligand_from_cgenff()

            if detected:
                ans = input(
                    f"🔍 Detected ligand '{detected}' from CGenFF files. "
                    "Use this ligand? [Y/n]: "
                ).strip().lower()

                if ans in ("", "y", "yes"):
                    ligand_code = detected
                else:
                    ligand_code = input(
                        "Enter the 3-letter ligand code: "
                    ).strip().upper()

            else:
                ligand_code = input(
                    "Enter the 3-letter ligand code: "
                ).strip().upper()

        if not ligand_code:
            print("❌ Ligand code is required for this simulation type.")
            exit(1)

    print(f"🧬 Final ligand selection: {ligand_code}")


    # Production length
    if args.ns is not None:
        simulation_time_ns = float(args.ns)
    else:
        length_input = input("How long should the production MD run (in nanoseconds)? [Press Enter for 50 ns]: ").strip()
        if length_input == "":
            simulation_time_ns = 50.0
        else:
            try:
                simulation_time_ns = float(length_input)
            except ValueError:
                print("ℹ️ Invalid input. Defaulting to 50 ns.")
                simulation_time_ns = 50.0
    if simulation_time_ns <= 0:
        print("ℹ️ Non-positive simulation time provided. Defaulting to 50 ns.")
        simulation_time_ns = 50.0
    # Determine compute mode
    use_gpu = (args.compute.upper() == "GPU")
    print(f"⚙️ Compute mode selected: {'GPU-accelerated' if use_gpu else 'CPU-only'}")


    setup_md(
        simulation_directory,
        gpu_id,
        ligand_code,
        simulation_type,
        simulation_time_ns,
        use_gpu,
        args,
        peptide_info=peptide_info
    )






