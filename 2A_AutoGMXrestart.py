#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=======================================================================
💠 2B_RestartAndFinalizeMD — Resume Production MD + Postprocessing
=======================================================================

Resumes production MD (md_0_1) from checkpoint or starts fresh if needed,
then automatically:
✅ wraps/recenters trajectory (Protein + Ligand only)
✅ extracts binding pocket (pocket PDB/XTC)
✅ writes Final_Trajectory.pdb for visualization
Skips downstream Python analyses.

Author: Joseph-Michael Schulz (PhD Candidate, UM BMB)
Version: 1.3 — Oct 2025
"""

import os
import subprocess
import argparse
from datetime import datetime
from pymacs_component_utils import (
    build_polymer_group_name,
    get_all_retained_components,
    index_has_group,
    load_component_registry,
    load_system_registry,
    merged_component_group_name,
    parse_component_list,
)

# ============================================================
# 🧾 Logging & Utilities
# ============================================================
import os, subprocess, psutil

def set_env(gpu_id=None, threads=32, cpu_range=None):
    """
    Configure environment for isolated GROMACS GPU execution.
    Ensures only the selected GPU is visible to this process.
    """
    env = os.environ.copy()

    # --- Hard reset all CUDA context variables ---
    for k in ["CUDA_VISIBLE_DEVICES", "CUDA_DEVICE_ORDER", "GMX_FORCE_GPU_ID"]:
        env.pop(k, None)

    # --- Map the selected GPU into its own sandbox ---
    if gpu_id is not None:
        env["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)  # expose only that GPU
        env["GMX_FORCE_GPU_ID"] = "0"              # GROMACS sees it as GPU 0
        log(f"🎯 Forcing isolation: physical GPU {gpu_id} → visible GPU 0")
    else:
        log("💻 Running in CPU-only mode.")

    # --- Threading and affinity settings ---
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PLACES"] = "cores"
    env["OMP_PROC_BIND"] = "close"
    if cpu_range:
        env["GOMP_CPU_AFFINITY"] = cpu_range
    env["GMX_GPU_DD_COMMS"] = "CUDA"
    env["GMX_GPU_PME_DECOMPOSITION"] = "CUDA"

    return env



def log(msg):
    ts = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    print(f"{ts} {msg}", flush=True)

def run(cmd, cwd=None, env=None, shell=True, check=True):
    log(f"🔹 Running: {cmd}")
    res = subprocess.run(cmd, cwd=cwd, env=env, shell=shell, text=True,
                        capture_output=True)
    if res.stdout:
        print(res.stdout.strip())
    if res.stderr and res.returncode != 0:
        print(res.stderr.strip())
    if check and res.returncode != 0:
        raise SystemExit(f"❌ Command failed: {cmd}")
    return res


def is_update_gpu_failure(stderr: str) -> bool:
    """
    Detect GROMACS failures related to GPU update incompatibility
    (e.g., virtual sites).
    """
    if not stderr:
        return False
    patterns = [
        "Update task can not run on the GPU",
        "Virtual sites are not supported",
        "decideWhetherToUseGpuForUpdate"
    ]
    return any(p in stderr for p in patterns)




def prompt_if_none(value, message, cast=str, default=None):
    if value is not None:
        return value
    s = input(f"{message}" + (f" [default: {default}] " if default is not None else " ")).strip()
    if s == "" and default is not None:
        return default
    try:
        return cast(s)
    except Exception:
        log(f"⚠️ Invalid input; using default {default}")
        return default

# ============================================================
# ⚙️ 1. Resume MD
# ============================================================


def resume_md(env, use_gpu=False, gpu_id=0, threads=32):
    if not os.path.exists("md_0_1.tpr"):
        log("❌ md_0_1.tpr not found — cannot start production MD.")
        return

    base_cmd = (
        f"gmx mdrun -v -deffnm md_0_1 "
        f"-pin off -ntmpi 1 -ntomp {threads}"
    )

    if os.path.exists("md_0_1.cpt"):
        base_cmd += " -cpi md_0_1.cpt"
        log("⏯️  Resuming MD from checkpoint (md_0_1.cpt)")
    else:
        log("▶️  Starting fresh MD run (no checkpoint found)")

    if not use_gpu:
        log("🧠 Launching CPU-only MD run.")
        run(base_cmd, env=env, check=True)
        return

    # --- GPU attempt #1: update on GPU ---
    gpu_cmd = (
        base_cmd +
        " -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu"
    )

    log("🚀 Attempting GPU MD run with GPU update")
    run("nvidia-smi -L", env=env, shell=True, check=False)

    res = run(gpu_cmd, env=env, check=False)

    # --- If it failed, inspect why ---
    if res.returncode != 0 and is_update_gpu_failure(res.stderr):
        log("⚠️ GPU update not supported (virtual sites detected)")
        log("🔁 Retrying with CPU-based updates (hybrid GPU/CPU mode)")

        cpu_update_cmd = (
            base_cmd +
            " -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update cpu"
        )

        run(cpu_update_cmd, env=env, check=True)
        log("✅ MD resumed successfully with CPU updates.")
        return

    # --- Other failure → abort ---
    if res.returncode != 0:
        raise SystemExit("❌ MD run failed for an unexpected reason.")

    log("✅ Production MD completed successfully.")




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
        log("⚠️ Multiple CGenFF MOL2 files detected:")
        for i, f in enumerate(mol2_files):
            log(f"  [{i}] {f}")
        idx = prompt_if_none(
            None,
            "Select ligand MOL2 index:",
            int,
            default=0
        )
        fname = mol2_files[idx]
    else:
        fname = mol2_files[0]

    ligand = fname.split(".cgenff.mol2")[0]
    return ligand.upper()



# ============================================================
# 📦 2. Wrap / Center Final Trajectory (Protein + Ligand only)
# ============================================================

def resolve_component_context(cli_ligand=None, cli_cofactors=None):
    registry = load_component_registry(".") or {}
    system_registry = load_system_registry(".") or {}
    ligand_code = (cli_ligand or registry.get("primary_ligand") or None)
    ligand_code = ligand_code.upper() if ligand_code else None
    cofactors = parse_component_list(cli_cofactors) if cli_cofactors else registry.get("cofactors", [])
    cofactors = [code for code in cofactors if code != ligand_code]
    all_components = get_all_retained_components(ligand_code, cofactors)
    merged_group = merged_component_group_name(all_components, fallback_ligand=ligand_code, protein_only=not all_components)
    polymer_group = build_polymer_group_name(system_registry, all_components, override="auto")
    return ligand_code, cofactors, all_components, (polymer_group if system_registry.get("has_nucleic_acid") else merged_group)


def wrap_trajectory(env, ligand_code=None, component_group_name=None):
    """Recenter & wrap the trajectory (protein + ligand only)."""
    if os.path.exists("Final_Trajectory.xtc"):
        log("ℹ️ Final_Trajectory.xtc already exists; skipping trjconv.")
        return
    if not os.path.exists("md_0_1.xtc") or not os.path.exists("md_0_1.tpr"):
        log("⚠️ Missing trajectory or TPR; cannot wrap.")
        return

    group_name = component_group_name or (f"Protein_{ligand_code}" if ligand_code else "Protein")
    if os.path.exists("index.ndx") and not index_has_group("index.ndx", group_name):
        fallback = f"Protein_{ligand_code}" if ligand_code else "Protein"
        group_name = fallback if index_has_group("index.ndx", fallback) else "Protein"
    log(f"📦 Re-centering & wrapping trajectory → Final_Trajectory.xtc (group: {group_name})")

    # Use index group from index.ndx (Protein + Ligand)
    if os.path.exists("index.ndx") and index_has_group("index.ndx", group_name):
        cmd = (
            f"bash -lc \"printf '{group_name}\\n{group_name}\\n' | "
            "gmx trjconv -s md_0_1.tpr -f md_0_1.xtc "
            "-o Final_Trajectory.xtc -n index.ndx -pbc mol -center\""
        )
    else:
        # fallback: just center on protein
        cmd = (
            "bash -lc \"printf 'Protein\\nProtein\\n' | gmx trjconv ..."
            "gmx trjconv -s md_0_1.tpr -f md_0_1.xtc "
            "-o Final_Trajectory.xtc -pbc mol -center\""
        )

    run(cmd, env=env, shell=True, check=True)

    # Create a PDB of the first frame for visualization
    if os.path.exists("Final_Trajectory.xtc"):
        log("🧩 Extracting first frame → Final_Trajectory.pdb")
        pdb_cmd = (
            "bash -lc \"printf '0\\n' | "
            "gmx trjconv -s md_0_1.tpr -f Final_Trajectory.xtc "
            "-o Final_Trajectory.pdb -dump 0\""
        )
        run(pdb_cmd, env=env, shell=True, check=True)
        if os.path.exists("Final_Trajectory.pdb"):
            log("✅ Final_Trajectory.pdb written successfully (first frame).")
        else:
            log("⚠️ PDB file not generated — check gmx trjconv output.")
    else:
        log("⚠️ Trajectory wrapping step did not produce output.")

# ============================================================
# 🧩 3. Build Binding Pocket
# ============================================================

def build_binding_pocket(ligand_code, cutoff_ang=5.0):
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import distances
        import numpy as np
    except Exception as e:
        log(f"❌ MDAnalysis not available: {e}")
        return

    if not ligand_code:
        log("ℹ️ No ligand code specified; skipping pocket extraction.")
        return
    if not os.path.exists("Final_Trajectory.xtc") or not os.path.exists("md_0_1.tpr"):
        log("⚠️ Missing trajectory/TPR; cannot extract pocket.")
        return

    log(f"🧬 Building binding-pocket subsystem for {ligand_code} (cutoff = {cutoff_ang:.1f} Å)")
    try:
        u = mda.Universe("md_0_1.tpr", "Final_Trajectory.xtc")
    except Exception:
        u = mda.Universe("md_0_1.gro", "Final_Trajectory.xtc")

    protein = u.select_atoms("protein and not name LP*")
    ligand  = u.select_atoms(f"resname {ligand_code} and not name LP*")

    if ligand.n_atoms == 0:
        log(f"❌ No atoms found for ligand {ligand_code}.")
        return

    close_resids = set()
    for ts in u.trajectory[::5]:  # sample every 5th frame for speed
        D = distances.distance_array(ligand.positions, protein.positions)
        for res in protein.residues:
            if np.any(D[:, res.atoms.indices] < cutoff_ang):
                close_resids.add(res.resid)

    if not close_resids:
        sel = f"resname {ligand_code}"
        log("⚠️ No residues within cutoff; saving ligand-only selection.")
    else:
        sel = f"resname {ligand_code} or (" + " or ".join(f"resid {r}" for r in sorted(close_resids)) + ")"

    pocket = u.select_atoms(sel)
    log(f"🧩 Pocket atoms: {pocket.n_atoms}")
    with mda.Writer("binding_pocket_only.xtc", pocket.n_atoms) as W:
        for ts in u.trajectory[::5]:
            W.write(pocket)
    pocket.write("binding_pocket_only.pdb")
    log("✅ Binding-pocket files written (binding_pocket_only.pdb / .xtc).")

# ============================================================
# 🏁 MAIN
# ============================================================

def main():
    import psutil  # for CPU core detection

    ap = argparse.ArgumentParser(description="Resume production MD and rebuild wrapped trajectory + pocket.")
    ap.add_argument("--ligand", default=None, help="Ligand 3-letter code (e.g., GDP)")
    ap.add_argument("--cofactors", default=None, help="Comma-separated retained cofactors/context components (e.g. CLR,HEM).")
    ap.add_argument("--gpu", type=int, default=None, help="GPU ID to use (e.g., 0)")
    ap.add_argument("--threads", type=int, default=None, help="Number of CPU threads to use")
    ap.add_argument("--offset", type=int, default=None, help="CPU core offset (e.g., 64 means use cores 65–128)")
    ap.add_argument("--pocket-cutoff", type=float, default=None, help="Å cutoff for pocket extraction")
    ap.add_argument("--cpu-only", action="store_true", help="Force CPU-only run even if GPU available")
    args = ap.parse_args()

    # ─── 1️⃣ Ligand Detection / Confirmation ─────────────────
    ligand, cofactors, all_components, merged_group = resolve_component_context(args.ligand, args.cofactors)

    if ligand is None:
        detected = detect_ligand_from_cgenff()

        if detected:
            ans = input(
                f"🔍 Detected ligand '{detected}' from CGenFF files. "
                "Use this ligand? [Y/n]: "
            ).strip().lower()
            if ans in ("", "y", "yes"):
                ligand = detected
            else:
                ligand = prompt_if_none(
                    None,
                    "Enter 3-letter ligand code:",
                    str
                )
        else:
            log("⚠️ No CGenFF ligand detected.")
            ligand = prompt_if_none(
                None,
                "Enter 3-letter ligand code:",
                str,
                default=None
            )

    log(f"🧬 Using ligand: {ligand}")
    log(f"🧬 Retained cofactors: {', '.join(cofactors) if cofactors else 'None'}")



    # ─── 2️⃣ CPU Threads + Offset ─────────────────────────
    total_cores = psutil.cpu_count(logical=True)
    default_threads = args.threads or total_cores
    threads = prompt_if_none(args.threads,
                            f"Enter number of CPU threads [max {total_cores}]:",
                            int, default=default_threads)
    if threads > total_cores:
        log(f"⚠️ Requested {threads} > available {total_cores}; using {total_cores}.")
        threads = total_cores

    offset = prompt_if_none(args.offset, f"Enter CPU offset (0–{total_cores - 1}):", int, default=0)
    if offset + threads > total_cores:
        log(f"⚠️ Offset+threads ({offset}+{threads}) exceeds total {total_cores}; adjusting offset.")
        offset = max(0, total_cores - threads)

    cpu_range = f"{offset}-{offset + threads - 1}"
    log(f"🧠 Using CPU cores {cpu_range} ({threads} threads).")

    # ─── 3️⃣ GPU Selection ───────────────────────────────
    use_gpu = not args.cpu_only
    gpu_id = args.gpu

    if use_gpu:
        try:
            gpu_list = subprocess.check_output(
                "nvidia-smi --query-gpu=name --format=csv,noheader", shell=True
            ).decode().strip().split('\n')
            if not gpu_list or gpu_list == ['']:
                log("⚠️ No GPU detected — switching to CPU-only mode.")
                use_gpu = False
            else:
                log("🧠 Available GPUs:")
                for i, name in enumerate(gpu_list):
                    log(f"  [{i}] {name}")
                gpu_id = prompt_if_none(gpu_id, "Select GPU ID:", int, default=0)
        except Exception as e:
            log(f"⚠️ GPU detection failed ({e}); running on CPU.")
            use_gpu = False

    env = set_env(gpu_id if use_gpu else None, threads)

    # Explicitly pin to CPU range if desired
    env["OMP_PLACES"] = "cores"
    env["OMP_PROC_BIND"] = "close"
    env["GOMP_CPU_AFFINITY"] = cpu_range  # works for OpenMP / GROMACS

    # ─── 4️⃣ Main Steps ──────────────────────────────────
    # Map selected GPU to CUDA_VISIBLE_DEVICES so only that one is seen
    if use_gpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        os.environ["GMX_FORCE_GPU_ID"] = "0"
        gpu_id = 0  # make GROMACS think it’s device 0

    resume_md(env, use_gpu=use_gpu, gpu_id=gpu_id, threads=threads)
    wrap_trajectory(env, ligand_code=ligand, component_group_name=merged_group)
    build_binding_pocket(ligand, cutoff_ang=args.pocket_cutoff or 5.0)
    log("✅ Restart + Finalization complete.")




if __name__ == "__main__":
    main()
