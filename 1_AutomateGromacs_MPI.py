"""
====================================================================
⚙️  1_AutomateGromacs.py — Automated GROMACS MD Setup & Run Pipeline
====================================================================

This script automates the setup and execution of molecular dynamics (MD)
simulations involving one or more proteins and an optional small-molecule
ligand. It supports both **interactive** (guided prompts) and **headless**
(non-interactive, fully automated) operation.

--------------------------------------------------------------
🚀 QUICKSTART SUMMARY
--------------------------------------------------------------
# Protein–ligand system (headless):
python 1_AutomateGromacs.py \
    --pdb Model_2.pdb \
    --chain-map "A:Rap1B,B:Rap1GAP" \
    --ligand GDP

# Protein–protein system (no ligand):
python 1_AutomateGromacs.py \
    --pdb Model_2.pdb \
    --chain-names "Rap1B,Rap1GAP" \
    --mode protein

--------------------------------------------------------------
🧭 INTERACTIVE MODE
--------------------------------------------------------------
Run with no arguments to launch guided mode:
    python 1_AutomateGromacs.py

You will be prompted to:
  • Select a PDB file
  • Indicate whether the system includes a ligand
  • Enter the ligand’s 3-letter code (if applicable)
  • Assign custom names to protein chains

--------------------------------------------------------------
🤖 HEADLESS MODE (Recommended for pipelines)
--------------------------------------------------------------
Skip all prompts by providing arguments:

  --pdb <FILENAME>             Select input PDB file directly
  --chain-map "A:Name,B:Name"  Explicit chain naming (recommended)
  --chain-names "Name1,Name2"  Shorthand naming (ordered by chain order)
  --ligand GDP                 Ligand 3-letter code (optional)
  --mode [ligand|protein]      Skip system-type prompt entirely

🧠 NOTE:
If `--ligand` is provided, the script will automatically
switch to *ligand mode* — specifying `--mode ligand` is optional.

Only one of `--chain-map` or `--chain-names` is required.

--------------------------------------------------------------
🧩 EXAMPLES
--------------------------------------------------------------
# 1️⃣ Full headless run (Rap1B : Rap1GAP + GDP)
python 1_AutomateGromacs.py \
    --pdb Model_2.pdb \
    --chain-map "A:Rap1B,B:Rap1GAP" \
    --ligand GDP

# 2️⃣ Protein-only run (Rap1B : Rap1GAP)
python 1_AutomateGromacs.py \
    --pdb Model_2.pdb \
    --chain-names "Rap1B,Rap1GAP" \
    --mode protein

# 3️⃣ Minimal ligand-only quick test
python 1_AutomateGromacs.py \
    --pdb test.pdb \
    --ligand ATP

--------------------------------------------------------------
🧠 NOTES
--------------------------------------------------------------
- Detected chain IDs are cached in `.chainmap.json` for reuse.
- Ligand topologies are automatically generated via SILCSBio CGenFF.
- Ion concentration defaults to 0.15 M (override with):
      export GMX_SALT_M=0.10
- Compatible with CHARMM36 forcefields and GROMACS ≥2023.
- All setup steps are logged to `mdrun.log` in the working directory.
- A runtime banner automatically summarizes the chosen mode,
  e.g. “🚀 Running headless: Rap1B–Rap1GAP + GDP [ligand mode]”.

====================================================================


--------------------------------------------------------------
📦 OUTPUT SUMMARY
--------------------------------------------------------------

After a successful run, the following files and folders will be generated
in your working directory (`Model_X/` or equivalent):

🧬 STRUCTURE & TOPOLOGY FILES
─────────────────────────────
protein.pdb                 → Ligand-stripped protein (from original PDB)
protein_processed.gro       → Processed protein structure (via pdb2gmx)
topol.top                   → Master topology including forcefield and ligand entries
atomIndex.txt               → Chain index mapping (e.g., Rap1B 1–2773, Rap1GAP 2774–9550)
.chainmap.json              → Cached chain-to-name map for future reuse

💊 LIGAND FILES (only if --ligand provided)
──────────────────────────────────────────
<lig>.cgenff.mol2           → Ligand structure for CGenFF parameterization
<lig>.str                   → CGenFF output parameter stream
<lig>.itp / <lig>.prm       → Ligand topology and parameters (auto-generated)
<lig>_ini.pdb               → Ligand reference coordinates (for pose grafting)
<lig>_pose_match.pdb        → Pose-aligned ligand in correct coordinates
<lig>.gro                   → Final GROMACS-format ligand coordinates

⚙️ SYSTEM ASSEMBLY FILES
────────────────────────
complex.gro                 → Combined protein + ligand structure
newbox.gro                  → Simulation box (1.0 nm buffer, cubic)
solv.gro                    → Solvated box with water molecules
ions.tpr                    → Pre-ionization input structure
solv_ions.gro               → Final solvated, neutralized system with 0.15 M NaCl

🧠 DIAGNOSTIC & LOG FILES
──────────────────────────
mdrun.log                   → Master run log (appends key runtime info)
grompp.mdp / ions.mdp       → Input parameter files used for setup
*.itp / *.prm               → Topology includes generated by pdb2gmx and CGenFF
_water.ndx (if used)        → Backup water-only index for genion
index.ndx                   → System index groups used by GROMACS
index_default.ndx           → Default group definitions before merging
chainmap.json               → Cached chain-to-name mapping for reruns

✅ FINAL PREPARED INPUT
────────────────────────
em.tpr                      → Ready-to-run energy minimization input file

This marks the completion of the **Setup Phase**.
The next script (e.g., `2_AutomateMD.py` or `2_AutomateGromacsRun.py`)
continues with equilibration (NVT, NPT) and production MD.

====================================================================
"""
# === END OF 1_AutomateGromacs.py HEADER ===


import os
import subprocess
# Define number of threads for GROMACS execution
import multiprocessing
import shutil  # For file operations
import shlex

import re


import argparse


VALID_BOX_TYPES = {"cubic", "dodecahedron", "octahedron"}
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GMX_BIN = None


def resolve_gmx_binary(requested="auto"):
    explicit_requested = requested not in (None, "", "auto")
    env_requested = os.environ.get("PYMACS_GMX_BIN")

    if explicit_requested:
        candidates = [requested]
    elif env_requested:
        candidates = [env_requested]
    else:
        candidates = ["gmx_mpi", "gmx-mpi", "gmx"]

    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved

    if explicit_requested:
        raise FileNotFoundError(f"Requested GROMACS executable was not found: {requested}")
    if env_requested:
        raise FileNotFoundError(f"PYMACS_GMX_BIN is set but not callable: {env_requested}")
    raise FileNotFoundError(
        "No GROMACS executable found. Tried: gmx_mpi, gmx-mpi, gmx. "
        "Use --gmx-bin or set PYMACS_GMX_BIN."
    )


def gmx_cmd(subcommand, args=""):
    base = f"{shlex.quote(GMX_BIN)} {subcommand}"
    return f"{base} {args}".strip()


def print_gmx_preflight(directory=None):
    print(f"🧬 Using GROMACS executable: {GMX_BIN}")
    try:
        result = subprocess.run(
            f"{shlex.quote(GMX_BIN)} --version",
            shell=True,
            text=True,
            capture_output=True,
            cwd=directory,
        )
    except Exception as exc:
        raise RuntimeError(f"Failed to execute '{GMX_BIN} --version': {exc}") from exc

    if result.returncode != 0:
        stderr = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(
            f"GROMACS executable is not callable: {GMX_BIN}\n{stderr}"
        )

    version_line = next(
        (line.strip() for line in result.stdout.splitlines() if line.strip()),
        None,
    )
    if version_line:
        print(f"🧪 GROMACS version: {version_line}")
    if directory:
        append_setup_log(directory, f"GROMACS executable: {GMX_BIN}")
        if version_line:
            append_setup_log(directory, f"GROMACS version: {version_line}")


parser = argparse.ArgumentParser(
    description="Automate MD setup with optional chain naming, PDB preflight checks, and ligand support.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument("--mode", choices=["ligand", "protein"], help="Skip interactive mode selection (headless mode).")
parser.add_argument("--chain-names", type=str, help="Comma-separated custom chain names (e.g. Rap1B,Rap1GAP)")
parser.add_argument("--chain-map", type=str, help="Explicit chainID:name pairs (e.g. A:Rap1B,B:Rap1GAP)")
parser.add_argument("--ligand", "-l", type=str, help="Ligand residue code used for non-covalent ligand workflows.")
parser.add_argument("--pdb", type=str, help="Input PDB filename (e.g., Model_2.pdb)")
parser.add_argument(
    "--gmx-bin",
    default="auto",
    help="GROMACS executable to use: auto, gmx, gmx_mpi, gmx-mpi, or an explicit path.",
)
parser.add_argument(
    "--box-type",
    choices=sorted(VALID_BOX_TYPES),
    default="cubic",
    help="GROMACS box shape passed to gmx editconf -bt.",
)
parser.add_argument(
    "--box-distance",
    type=float,
    default=1.0,
    help="Distance in nm between the solute and the periodic box edge.",
)
parser.add_argument(
    "--strict-pdb-validation",
    action="store_true",
    help="Escalate selected PDB preflight warnings into errors before setup starts.",
)
parser.add_argument(
    "--allow-covalent-ligand",
    action="store_true",
    help="Allow ligand workflows to continue after LINK/CONECT records suggest a covalent attachment.",
)
parser.add_argument(
    "--mdp-dir",
    default="MDPs",
    help="Directory searched for setup-stage MDP templates when local copies are not present.",
)
parser.add_argument(
    "--ions-mdp",
    default="ions.mdp",
    help="Path or filename for the ion-preparation MDP template.",
)
parser.add_argument(
    "--em-mdp",
    default="em.mdp",
    help="Path or filename for the energy-minimization MDP template used at the end of Step 1.",
)

args = parser.parse_args()
if args.box_distance <= 0:
    parser.error("--box-distance must be greater than 0.0 nm.")
try:
    GMX_BIN = resolve_gmx_binary(args.gmx_bin)
    print_gmx_preflight()
except (FileNotFoundError, RuntimeError) as exc:
    raise SystemExit(f"❌ {exc}")

# ============================================================
# 🚀 Runtime Summary Banner
# ============================================================
if args.mode or args.ligand or args.chain_names or args.chain_map:
    chain_display = ""
    if args.chain_map:
        # Convert "A:Rap1B,B:Rap1GAP" → "Rap1B–Rap1GAP"
        chain_display = "–".join([p.split(":")[1] for p in args.chain_map.split(",")])
    elif args.chain_names:
        chain_display = "–".join(args.chain_names.split(","))
    else:
        chain_display = "UnnamedChains"

    if args.mode == "ligand" or args.ligand:
        mode_label = "ligand mode"
        ligand_label = args.ligand or "UnknownLigand"
        banner = f"🚀 Running headless: {chain_display} + {ligand_label} [{mode_label}]"
    elif args.mode == "protein":
        mode_label = "protein-only mode"
        banner = f"🚀 Running headless: {chain_display} [{mode_label}]"
    else:
        banner = f"🚀 Running headless: {chain_display} [auto mode]"

    print(f"\n{banner}\n")
    # Optional: also log it to mdrun.log for traceability
    with open(os.path.join(os.getcwd(), "mdrun.log"), "a") as log:
        log.write(f"{banner}\n")


def append_setup_log(directory, message):
    log_path = os.path.join(directory, "mdrun.log")
    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write(f"{message}\n")


def resolve_existing_path(base_directory, candidate):
    if not candidate:
        return None
    options = []
    if os.path.isabs(candidate):
        options.append(candidate)
    else:
        options.extend([
            os.path.join(base_directory, candidate),
            os.path.join(SCRIPT_DIR, candidate),
        ])
    for path in options:
        if os.path.exists(path):
            return os.path.abspath(path)
    return None


def resolve_mdp_template(base_directory, mdp_name_or_path, mdp_dir):
    candidates = []
    if mdp_name_or_path:
        candidates.extend([
            mdp_name_or_path,
            os.path.join(base_directory, mdp_name_or_path),
            os.path.join(SCRIPT_DIR, mdp_name_or_path),
        ])

    basename = os.path.basename(mdp_name_or_path) if mdp_name_or_path else None
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


def prepare_step1_mdp_files(directory, options):
    resolved = {
        "ions.mdp": resolve_mdp_template(directory, options.ions_mdp, options.mdp_dir),
        "em.mdp": resolve_mdp_template(directory, options.em_mdp, options.mdp_dir),
    }
    missing = [name for name, path in resolved.items() if not path]
    if missing:
        raise FileNotFoundError(
            "Missing required setup MDP template(s): "
            + ", ".join(missing)
            + f". Searched local files and '{options.mdp_dir}'."
        )

    local_map = {}
    for local_name, source_path in resolved.items():
        dest_path = os.path.join(directory, local_name)
        if os.path.abspath(source_path) != os.path.abspath(dest_path):
            shutil.copyfile(source_path, dest_path)
            print(f"📄 Copied {local_name} template from {source_path}")
        else:
            print(f"📄 Using local {local_name}: {dest_path}")
        append_setup_log(directory, f"MDP {local_name}: source={source_path} local={dest_path}")
        local_map[local_name] = dest_path
    return local_map


def detect_covalent_ligand_hints(pdb_path, ligand_code):
    if not ligand_code:
        return []

    ligand_serials = set()
    serial_is_ligand = {}
    hints = []

    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
        except ValueError:
            continue
        resname = line[17:20].strip().upper()
        is_ligand = resname == ligand_code.upper()
        serial_is_ligand[serial] = is_ligand
        if is_ligand:
            ligand_serials.add(serial)

    for line in lines:
        if line.startswith("LINK"):
            left_res = line[17:20].strip().upper()
            right_res = line[47:50].strip().upper()
            if ligand_code.upper() in {left_res, right_res}:
                hints.append(f"LINK record suggests covalent attachment: {line.strip()}")
        elif line.startswith("CONECT") and ligand_serials:
            tokens = line.split()[1:]
            try:
                serials = [int(token) for token in tokens]
            except ValueError:
                continue
            if not serials:
                continue
            anchor = serials[0]
            if anchor in ligand_serials:
                partners = [s for s in serials[1:] if not serial_is_ligand.get(s, False)]
                if partners:
                    hints.append(
                        f"CONECT record links ligand atom {anchor} to non-ligand atom(s) {partners}."
                    )
            elif serial_is_ligand.get(anchor, False) is False:
                partners = [s for s in serials[1:] if s in ligand_serials]
                if partners:
                    hints.append(
                        f"CONECT record links protein atom {anchor} to ligand atom(s) {partners}."
                    )
    return hints


def validate_input_pdb(pdb_path, ligand_code=None, strict=False, allow_covalent=False):
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    atom_records = []
    hetatm_resnames = set()
    warnings = []
    chain_ids = set()
    parsed_coords = 0

    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_records.append(line)
        resname = line[17:20].strip().upper()
        chain_id = line[21].strip()
        if resname:
            if line.startswith("HETATM"):
                hetatm_resnames.add(resname)
        else:
            raise ValueError(
                "A PDB atom record is missing a residue name. "
                "Please provide standard ATOM/HETATM residue fields before running PyMACS."
            )
        if chain_id:
            chain_ids.add(chain_id)
        try:
            float(line[30:38])
            float(line[38:46])
            float(line[46:54])
            parsed_coords += 1
        except ValueError:
            raise ValueError(
                "Failed to parse one or more PDB coordinates from columns 31-54. "
                "Please verify standard ATOM/HETATM coordinate formatting."
            )

    if not atom_records:
        raise ValueError(
            "No ATOM/HETATM records were found in the input PDB. "
            "PyMACS expects a coordinate PDB with standard atom records."
        )
    if parsed_coords == 0:
        raise ValueError("No parseable coordinates were found in the input PDB.")

    if not chain_ids:
        warnings.append(
            "No chain IDs were found in the PDB. PyMACS can continue, but chain-aware naming is strongly recommended."
        )

    detected_ligands = autodetect_ligands(pdb_path)
    if len(detected_ligands) > 1:
        warnings.append(
            "Multiple ligand-like residues were detected: "
            + ", ".join(detected_ligands)
            + ". Provide --ligand explicitly if the intended ligand is ambiguous."
        )

    if ligand_code:
        ligand_code = ligand_code.upper()
        if ligand_code not in {line[17:20].strip().upper() for line in atom_records}:
            raise ValueError(
                f"Ligand residue '{ligand_code}' was requested but was not found in {os.path.basename(pdb_path)}. "
                "Use --ligand with a residue name that exists in the PDB or rerun in protein mode."
            )
        covalent_hints = detect_covalent_ligand_hints(pdb_path, ligand_code)
        if covalent_hints:
            message = (
                "Potential covalent ligand attachment detected. The default PyMACS ligand workflow is designed "
                "for non-covalent ligands; covalent systems usually need custom topology handling.\n"
                + "\n".join(f"  - {hint}" for hint in covalent_hints)
            )
            if allow_covalent:
                warnings.append(message + "\nProceeding because --allow-covalent-ligand was supplied.")
            else:
                raise ValueError(message + "\nRe-run with --allow-covalent-ligand only if you have prepared the topology manually.")

    if strict and warnings:
        raise ValueError(
            "Strict PDB validation failed due to the following warning-level issues:\n"
            + "\n".join(f"  - {warning}" for warning in warnings)
        )

    print(f"✅ PDB preflight passed for {os.path.basename(pdb_path)}")
    for warning in warnings:
        print(f"⚠️ PDB preflight warning: {warning}")
    return {
        "warnings": warnings,
        "detected_ligands": detected_ligands,
        "chain_ids": sorted(chain_ids),
        "hetatm_resnames": sorted(hetatm_resnames),
        "atom_record_count": len(atom_records),
    }



def ensure_cgenff_prm_included(topol_path, prm_file, itp_file):
    with open(topol_path) as f:
        lines = f.readlines()

    new_lines = []
    inserted = False

    for line in lines:
        new_lines.append(line)

        # Insert PRM immediately after forcefield include
        if "forcefield.itp" in line and not inserted:
            new_lines.append(f'#include "{prm_file}"\n')
            inserted = True

    if not inserted:
        raise RuntimeError("❌ Could not locate forcefield.itp include in topol.top")

    # Ensure ligand ITP is included AFTER prm
    if itp_file not in "".join(new_lines):
        new_lines.append(f'#include "{itp_file}"\n')

    with open(topol_path, "w") as f:
        f.writelines(new_lines)

def regenerate_protein_topology(forcefield="charmm36-mar2019"):
    """Deletes existing topol.top and reruns pdb2gmx with a specific forcefield."""
    print(f"🔄 Regenerating protein topology using {forcefield}...")
    # Clean up old files to avoid conflicts
    for f in ["topol.top", "posre.itp", "protein_processed.gro"]:
        if os.path.exists(f):
            os.remove(f)
    
    # This logic assumes the forcefield folder exists in your directory
    # or is in the GROMACS path.
    return automate_pdb2gmx(os.getcwd())
    



def enforce_gromacs_prm_order(prm_path):
    """
    Reorders a CGenFF .prm file into GROMACS-friendly order WITHOUT
    requiring every possible section, and WITHOUT dropping duplicate sections
    like repeated [ dihedraltypes ] blocks.
    """
    ORDER = ["atomtypes", "bondtypes", "angletypes", "dihedraltypes", "impropertypes", "nonbond_params"]

    blocks = []  # list of (section_name, [lines])
    current = None
    current_lines = []

    with open(prm_path) as f:
        for line in f:
            m = re.match(r"\[\s*([^\]]+?)\s*\]", line)
            if m:
                # flush previous block
                if current is not None:
                    blocks.append((current, current_lines))
                current = m.group(1).strip()
                current_lines = [line]
            else:
                if current is None:
                    # preamble / comments before first section
                    blocks.append(("__preamble__", [line]))
                else:
                    current_lines.append(line)

    if current is not None:
        blocks.append((current, current_lines))

    # sanity: ensure prm has at least something useful
    present = {sec.lower() for sec, _ in blocks}
    bonded_ok = any(s in present for s in ("bondtypes", "angletypes", "dihedraltypes"))
    if not bonded_ok:
        raise RuntimeError(f"❌ PRM has no bonded sections (bond/angle/dihedral). File: {prm_path}")

    # reorder by ORDER but keep duplicates in original relative order
    out = []
    # keep preamble first
    for sec, lines in blocks:
        if sec == "__preamble__":
            out.extend(lines)

    for want in ORDER:
        for sec, lines in blocks:
            if sec.lower() == want:
                out.extend(lines)

    # append any unknown sections at end (rare, but safe)
    for sec, lines in blocks:
        if sec == "__preamble__":
            continue
        if sec.lower() not in ORDER:
            out.extend(lines)

    with open(prm_path, "w") as f:
        f.writelines(out)

    print(f"✅ PRM reordered safely (duplicates preserved): {os.path.basename(prm_path)}")





def validate_cgenff_outputs(ligand_code, directory):
    base = ligand_code.lower()
    itp = os.path.join(directory, f"{base}.itp")
    prm = os.path.join(directory, f"{base}.prm")
    ini = os.path.join(directory, f"{base}_ini.pdb")

    for path in (itp, prm, ini):
        if not os.path.exists(path):
            raise RuntimeError(f"❌ Missing required CGenFF output: {path}")
        if os.path.getsize(path) == 0:
            raise RuntimeError(f"❌ Empty CGenFF output file: {path}")

    # reorder but do NOT require atomtypes/nonbond_params
    enforce_gromacs_prm_order(prm)

    print(f"✅ CGenFF outputs validated for ligand {ligand_code}")




def build_ligand_mol2_rdkit(smiles, ligand_code, out_mol2):
    """
    Canonical RDKit → MOL2 path with explicit hydrogens and sane aromaticity.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError("❌ RDKit failed to parse SMILES")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=2025)
    AllChem.UFFOptimizeMolecule(mol)

    mol.SetProp("_Name", ligand_code)

    Chem.MolToMol2File(mol, out_mol2)
    print(f"✅ RDKit-generated MOL2 written to {out_mol2}")



def run_and_log(cmd, cwd, stdout_path=None, stderr_path=None, shell=False):
    """
    Run a command and optionally write stdout/stderr to files.
    Returns (returncode, stdout_text, stderr_text).
    """
    stdout_f = open(stdout_path, "w") if stdout_path else subprocess.PIPE
    stderr_f = open(stderr_path, "w") if stderr_path else subprocess.PIPE
    try:
        r = subprocess.run(cmd, cwd=cwd, text=True, shell=shell,
                           stdout=stdout_f, stderr=stderr_f)
    finally:
        if stdout_path: stdout_f.close()
        if stderr_path: stderr_f.close()

    # If we wrote to files, read them back for quick checks
    out_txt = ""
    err_txt = ""
    if stdout_path and os.path.exists(stdout_path):
        out_txt = open(stdout_path, "r").read()
    if stderr_path and os.path.exists(stderr_path):
        err_txt = open(stderr_path, "r").read()
    return r.returncode, out_txt, err_txt




def fix_missing_atoms(input_pdb, output_pdb, ph=7.4):
    """Uses PDBFixer to fill missing atoms and add hydrogens to the PDB file."""
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        fixer = PDBFixer(filename=input_pdb)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=ph)
        PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))
        print(f"✅ Fixed missing atoms and saved cleaned file to: {output_pdb}")
        return True
    except Exception as e:
        print(f"❌ Failed to fix atoms in {input_pdb}: {str(e)}")
        return False

def modify_topology_file(topol_path, ligand_code):
    """Updates the topology file by inserting required lines at the correct positions."""
    # Modified: ensure ligand parameters are included dynamically if ligand is present
    with open(topol_path, "r") as f:
        lines = f.readlines()
    # If ligand is present, include its parameter and topology files
    if ligand_code:
        # Find insertion point after forcefield includes
        insertion_index = None
        for i, line in enumerate(lines):
            if "; Include forcefield parameters" in line:
                insertion_index = i + 2
                break
        if insertion_index is not None:
            prm_line = f'#include "{ligand_code.lower()}.prm"\n'
            itp_line = f'; Include ligand topology\n#include "{ligand_code.lower()}.itp"\n'
            if prm_line not in lines[insertion_index:]:
                lines.insert(insertion_index, prm_line)
                lines.insert(insertion_index + 1, itp_line)
        # Append ligand molecule entry at end of file if not present
        ligand_entry = f"{ligand_code:<20}1\n"
        if not any(line.strip().startswith(ligand_code) for line in lines):
            lines.append(ligand_entry)
        print(f"✅ Added ligand includes and entry for {ligand_code} in topol.top")
    else:
        print("ℹ️ No ligand present; skipping ligand-specific topology modifications.")
    # Write changes back to topology file
    with open(topol_path, "w") as f:
        f.writelines(lines)

def merge_gro_files(protein_gro, ligand_gro, complex_gro):
    """Merges protein_processed.gro and ligand .gro into complex.gro while preserving the last row."""
    if not os.path.exists(protein_gro) or not os.path.exists(ligand_gro):
        print(f"❌ ERROR: Missing input files: {protein_gro} or {ligand_gro}")
        return
    # Read input .gro files
    with open(protein_gro, "r") as f:
        protein_lines = f.readlines()
    with open(ligand_gro, "r") as f:
        ligand_lines = f.readlines()
    correct_box = protein_lines[-1].strip()  # Box from protein
    ligand_body = ligand_lines[2:-1]         # Ligand atoms (skip header/count/box)
    total_atoms = int(protein_lines[1].strip()) + int(ligand_lines[1].strip())
    merged_lines = []
    merged_lines.append(protein_lines[0])                # Title
    merged_lines.append(f"{total_atoms}\n")              # Updated atom count
    merged_lines.extend(protein_lines[2:-1])             # Protein atoms
    merged_lines.extend(ligand_body)                     # Ligand atoms
    merged_lines.append(correct_box + "\n")              # Box dimensions
    with open(complex_gro, "w") as f:
        f.writelines(merged_lines)
    # Verify last line
    with open(complex_gro, "r") as f:
        final_lines = f.readlines()
        if final_lines[-1].strip() != correct_box:
            print("❌ ERROR: Final row of complex.gro does not match protein_processed.gro!")
            exit(1)
    print(f"✅ complex.gro merged successfully with {total_atoms} atoms.")

def run_command(command, cwd=None, input_text=None):
    """Runs a shell command in the given directory with optional input and exits on failure."""
    print(f"\n🔹 Running command: {command}")
    print(f"📂 In directory: {cwd if cwd else os.getcwd()}")

    try:
        result = subprocess.run(
            command,
            shell=True,
            cwd=cwd,
            input=input_text,
            capture_output=True,
            encoding="utf-8",     # ← force UTF-8 decoding
            errors="replace"      # ← prevent crash if weird byte appears
        )

        print(f"📜 STDOUT:\n{result.stdout}")

        if result.returncode != 0:
            print(f"❌ Error running command: {command}")
            print(f"🔺 STDERR: {result.stderr}")
            print("⚠️ Exiting script due to failure.")
            exit(1)

    except Exception as e:
        print(f"❌ Exception occurred while running command: {command}\n{str(e)}")
        print("⚠️ Exiting script due to failure.")
        exit(1)

    return True


def detect_forcefield(directory):
    """Detects if a CHARMM forcefield directory exists in the working directory."""
    for item in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, item)) and item.startswith("charmm36"):
            forcefield_path = os.path.join(directory, item)
            print(f"✅ CHARMM forcefield detected at {forcefield_path}. Selecting automatically.")
            return forcefield_path
    print("❗ CHARMM forcefield not found in the current directory. User must select manually.")
    return None

def select_pdb_file(directory):
    """Prompts the user to select a PDB file from available options."""
    pdb_files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
    if not pdb_files:
        print("❌ No PDB files found in the current directory.")
        exit(1)
    print("\n📁 Available PDB files:")
    for idx, file in enumerate(pdb_files, 1):
        print(f"{idx}. {file}")
    while True:
        selection = input("\n🔍 Select a PDB file by number (e.g., 1): ").strip()
        if selection.isdigit() and 1 <= int(selection) <= len(pdb_files):
            return pdb_files[int(selection) - 1]
        else:
            print("❗ Invalid selection. Please enter a valid number.")

def get_ligand_code():
    """Prompts user to enter ligand 3-letter code or use default."""
    choice = input("\n💊 Enter ligand 3-letter code or press ENTER to use default [PTC]: ").strip()
    return choice.upper() if choice else "PTC"

def automate_pdb2gmx(directory):
    """Automates GROMACS pdb2gmx selection for force field and termini settings."""
    omp_threads = os.environ.get("OMP_NUM_THREADS") or str(multiprocessing.cpu_count())
    os.environ["OMP_NUM_THREADS"] = omp_threads
    print(f"🧠 Using OMP_NUM_THREADS = {omp_threads}")
    protein_pdb = os.path.join(directory, "protein.pdb")
    print("\n🔍 Checking for protein.pdb in the working directory...")
    if not os.path.exists(protein_pdb):
        print(f"❌ Skipping {directory}: protein.pdb not found.")
        return False
    # Detect force field
    forcefield_path = detect_forcefield(directory)
    if forcefield_path is None:
        print("\n❌ CHARMM forcefield not found. Manual selection required.")
        return False
    print(f"✅ Using CHARMM force field from: {forcefield_path}")
    forcefield_option = "1"  # CHARMM
    chain_count = count_chains(protein_pdb)
    print(f"🔍 Detected {chain_count} chains in the protein.")
    if chain_count == 0:
        print("❌ No chains detected. Check protein.pdb formatting.")
        return False
    input_selections = [forcefield_option]
    for i in range(chain_count):
        input_selections.append("1" if i == 0 else "0")  # Start terminus: first chain NH3+, others PRO-NH2+
        input_selections.append("0")                    # End terminus: COO-
    input_text = "\n".join(input_selections) + "\n"
    print("\n✅ Auto-selecting force field and termini for pdb2gmx.")
    return run_command(
        gmx_cmd("pdb2gmx", "-f protein.pdb -o protein_processed.gro -p topol.top -ignh -water spc -ter"),
        cwd=directory,
        input_text=input_text,
    )

def count_chains(protein_pdb):
    """Counts the number of unique chains in a PDB file."""
    chains = set()
    with open(protein_pdb, "r") as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]  # Chain ID column
                chains.add(chain_id)
    return len(chains)

import os, json



def generate_atom_index_file(directory, chain_names_arg=None, chain_map_arg=None, source_pdb_path=None):
    """
    Generates atomIndex.txt with chain names, while also displaying:
      - current detected chain IDs
      - atom index ranges
      - residue counts
      - original source chain IDs / residue spans (if source_pdb_path provided)

    Priority order for naming:
      1) Explicit CLI chain_map_arg / chain_names_arg
      2) Cached .chainmap.json (if COMPLETE)
      3) Interactive prompt (fallback only)
    """

    import os, json
    from collections import OrderedDict

    protein_pdb = os.path.join(directory, "protein.pdb")
    protein_gro = os.path.join(directory, "protein_processed.gro")
    output_file = os.path.join(directory, "atomIndex.txt")
    cache_file  = os.path.join(directory, ".chainmap.json")

    if not os.path.exists(protein_pdb) or not os.path.exists(protein_gro):
        print(f"❌ Missing protein.pdb or protein_processed.gro in {directory}")
        return False

    def parse_pdb_chain_summary(pdb_path):
        chains = OrderedDict()

        with open(pdb_path) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue

                chain_id = line[21].strip() or "A"
                resname = line[17:20].strip()

                try:
                    resid = int(line[22:26].strip())
                except ValueError:
                    continue

                key = (chain_id, resid)
                if chain_id not in chains:
                    chains[chain_id] = {
                        "chain_id": chain_id,
                        "resids": [],
                        "residue_keys": set(),
                        "first_resname": resname,
                        "last_resname": resname,
                    }

                if key not in chains[chain_id]["residue_keys"]:
                    chains[chain_id]["residue_keys"].add(key)
                    chains[chain_id]["resids"].append(resid)

                chains[chain_id]["last_resname"] = resname

        out = []
        for cid, info in chains.items():
            resids = info["resids"]
            if resids:
                out.append({
                    "chain_id": cid,
                    "resids": resids,
                    "residue_count": len(resids),
                    "resid_start": min(resids),
                    "resid_end": max(resids),
                    "first_resname": info["first_resname"],
                    "last_resname": info["last_resname"],
                })
        return out

    def parse_gro_chain_atom_ranges(protein_gro, chain_ids):
        with open(protein_gro) as f:
            lines = f.readlines()

        atom_ranges = []
        atom_counter = 1
        last_resid = None
        current_start = 1
        chain_index = 0
        gro_chain_start_resid = None
        current_resid = None

        gro_lines = lines[2:-1]

        for i, line in enumerate(gro_lines):
            try:
                resid = int(line[:5].strip())
            except ValueError:
                continue

            atom_counter = i + 1
            current_resid = resid

            if gro_chain_start_resid is None:
                gro_chain_start_resid = resid

            if last_resid is not None and resid < last_resid:
                atom_ranges.append((
                    chain_ids[chain_index],
                    current_start,
                    atom_counter,
                    gro_chain_start_resid,
                    last_resid
                ))
                current_start = atom_counter + 1
                gro_chain_start_resid = resid
                chain_index = min(chain_index + 1, len(chain_ids) - 1)

            last_resid = resid

        atom_ranges.append((
            chain_ids[chain_index],
            current_start,
            atom_counter + 1,
            gro_chain_start_resid,
            current_resid
        ))

        return atom_ranges

    # Detect current chains
    current_chain_summaries = parse_pdb_chain_summary(protein_pdb)
    chain_ids = [c["chain_id"] for c in current_chain_summaries]

    if not chain_ids:
        chain_ids = ["A"]
        current_chain_summaries = [{
            "chain_id": "A",
            "resids": [],
            "residue_count": 0,
            "resid_start": None,
            "resid_end": None,
            "first_resname": None,
            "last_resname": None,
        }]

    print(f"🔍 Detected {len(chain_ids)} chain(s): {', '.join(chain_ids)}")

    # Parse GRO + source provenance EARLY
    atom_ranges = parse_gro_chain_atom_ranges(protein_gro, chain_ids)

    source_chain_summaries = []
    if source_pdb_path and os.path.exists(source_pdb_path):
        source_chain_summaries = parse_pdb_chain_summary(source_pdb_path)

    provenance_by_current_chain = {}
    for i, cid in enumerate(chain_ids):
        current_info = current_chain_summaries[i] if i < len(current_chain_summaries) else None
        source_info = source_chain_summaries[i] if i < len(source_chain_summaries) else None
        provenance_by_current_chain[cid] = {
            "current": current_info,
            "source": source_info,
        }

    # PRINT SUMMARY BEFORE ASKING FOR NAMES
    print("\n🧬 Chain provenance preview")
    print("-" * 90)
    for cid, atom_start, atom_end, gro_resid_start, gro_resid_end in atom_ranges:
        prov = provenance_by_current_chain.get(cid, {})
        cur = prov.get("current")
        src = prov.get("source")

        cur_rescount = cur["residue_count"] if cur else "?"
        cur_span = f"{cur['resid_start']}-{cur['resid_end']}" if cur and cur["resid_start"] is not None else "?"
        src_chain = src["chain_id"] if src else "?"
        src_span = f"{src['resid_start']}-{src['resid_end']}" if src and src["resid_start"] is not None else "?"
        src_rescount = src["residue_count"] if src else "?"

        print(
            f"  chain {cid} | "
            f"GRO atoms {atom_start}-{atom_end} | "
            f"current residues {cur_span} [{cur_rescount}] | "
            f"source chain {src_chain} | "
            f"source residues {src_span} [{src_rescount}]"
        )

    chain_map = {}

    # Priority 1: explicit CLI mapping
    if chain_map_arg:
        try:
            for pair in chain_map_arg.split(","):
                cid, name = pair.split(":", 1)
                chain_map[cid.strip()] = name.strip()
            print(f"🧾 Using chain map from CLI: {chain_map}")
        except Exception as e:
            print(f"❌ Failed to parse --chain-map: {e}")
            return False

    elif chain_names_arg:
        names = [n.strip() for n in chain_names_arg.split(",")]
        for i, cid in enumerate(chain_ids):
            if i < len(names):
                chain_map[cid] = names[i]
        print(f"🧾 Using chain names from CLI: {chain_map}")

    # Priority 2: cached mapping
    elif os.path.exists(cache_file):
        with open(cache_file) as f:
            saved = json.load(f)

        for cid in chain_ids:
            if cid in saved:
                chain_map[cid] = saved[cid]

        if all(cid in chain_map for cid in chain_ids):
            print(f"💾 Loaded complete chain map from cache: {chain_map}")
        else:
            print("⚠️ Cached chain map incomplete — ignoring cache.")
            chain_map = {}

    # Priority 3: interactive prompt
    if not chain_map:
        print("\n💬 Entering interactive chain naming mode.")

        for cid in chain_ids:
            while True:
                name = input(f"🧩 Enter descriptive name for chain {cid}: ").strip()
                if name:
                    chain_map[cid] = name
                    break
                print("⚠️ Name cannot be empty.")

        print("\n✅ Chain map entered:")
        for cid, name in chain_map.items():
            print(f"   {cid} → {name}")

        confirm = input("Confirm these names? [Y/n]: ").strip().lower()
        if confirm == "n":
            print("🔁 Restarting chain naming...")
            return generate_atom_index_file(
                directory,
                chain_names_arg=chain_names_arg,
                chain_map_arg=chain_map_arg,
                source_pdb_path=source_pdb_path
            )

    with open(cache_file, "w") as f:
        json.dump(chain_map, f, indent=2)

    print("💾 Chain map cached.")

    print("\n🧬 Final chain provenance summary")
    print("-" * 90)
    for cid, atom_start, atom_end, gro_resid_start, gro_resid_end in atom_ranges:
        label = chain_map.get(cid, cid)
        prov = provenance_by_current_chain.get(cid, {})
        cur = prov.get("current")
        src = prov.get("source")

        cur_rescount = cur["residue_count"] if cur else "?"
        cur_span = f"{cur['resid_start']}-{cur['resid_end']}" if cur and cur["resid_start"] is not None else "?"
        src_chain = src["chain_id"] if src else "?"
        src_span = f"{src['resid_start']}-{src['resid_end']}" if src and src["resid_start"] is not None else "?"
        src_rescount = src["residue_count"] if src else "?"

        print(
            f"  current {cid} ({label}) | "
            f"GRO atoms {atom_start}-{atom_end} | "
            f"current residues {cur_span} [{cur_rescount}] | "
            f"source chain {src_chain} | "
            f"source residues {src_span} [{src_rescount}]"
        )

    with open(output_file, "w") as f:
        for cid, atom_start, atom_end, gro_resid_start, gro_resid_end in atom_ranges:
            label = chain_map[cid]
            prov = provenance_by_current_chain.get(cid, {})
            src = prov.get("source")

            if src:
                f.write(
                    f"{label} {atom_start}-{atom_end} "
                    f"# current_chain={cid} current_residues={gro_resid_start}-{gro_resid_end} "
                    f"source_chain={src['chain_id']} source_residues={src['resid_start']}-{src['resid_end']} "
                    f"source_rescount={src['residue_count']}\n"
                )
            else:
                f.write(
                    f"{label} {atom_start}-{atom_end} "
                    f"# current_chain={cid} current_residues={gro_resid_start}-{gro_resid_end}\n"
                )

    print(f"✅ Generated atomIndex.txt in {directory}")
    return True


    

# ======== NEW HELPERS FOR REALTIME CGenFF PREP ========

def ensure_silcsbio_env():
    """Verify SILCSBio is on this machine and basic files exist."""
    home = os.environ.get("SILCSBIO_HOME", "")
    if not home or not os.path.isdir(home):
        print("❌ SILCSBIO_HOME is not set or directory missing. Make sure your ~/.bashrc SILCSBio setup is loaded.")
        return None, None, None
    cgenff_exec = os.path.join(home, "cgenff", "cgenff")
    prm36      = os.path.join(home, "cgenff", "par_all36_cgenff.prm")
    if not os.path.isfile(cgenff_exec):
        print("❌ Missing SILCSBio cgenff executable:", cgenff_exec); return None, None, None
    if not os.path.isfile(prm36):
        print("❌ Missing par_all36_cgenff.prm:", prm36); return None, None, None
    return cgenff_exec, prm36, home


def extract_ligand_from_pdb(pdb_path, ligand_code, out_pdb):
    """
    Extract ligand HETATM lines AND any CONECT records that reference ligand atom serials.
    This preserves bonding so OpenBabel doesn't guess wrong.
    """
    ligand_code = ligand_code.upper()
    ligand_serials = set()
    het_lines = []
    conect_lines = []

    with open(pdb_path, "r") as fin:
        for line in fin:
            if line.startswith("HETATM") and line[17:20].strip().upper() == ligand_code:
                het_lines.append(line)
                try:
                    serial = int(line[6:11].strip())
                    ligand_serials.add(serial)
                except:
                    pass

    if not het_lines:
        print(f"❌ Could not find ligand {ligand_code} in {pdb_path}")
        return False

    # second pass: collect CONECT lines that touch ligand atoms
    with open(pdb_path, "r") as fin:
        for line in fin:
            if line.startswith("CONECT"):
                nums = []
                for i in (6, 11, 16, 21, 26):  # PDB CONECT columns
                    chunk = line[i:i+5].strip()
                    if chunk:
                        try: nums.append(int(chunk))
                        except: pass
                if any(n in ligand_serials for n in nums):
                    conect_lines.append(line)

    with open(out_pdb, "w") as fout:
        fout.writelines(het_lines)
        if conect_lines:
            fout.writelines(conect_lines)
        fout.write("END\n")

    print(f"✅ Extracted ligand {ligand_code} (+CONECT={bool(conect_lines)}) to {out_pdb}")
    return True




def ligand_has_hydrogens(pdb_file):
    """Detect if ligand PDB has hydrogens."""
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                aname = line[12:16].strip()
                elem  = line[76:78].strip() if len(line) >= 78 else ""
                if aname.upper().startswith("H") or elem.upper() == "H":
                    return True
    return False



def autodetect_ligands(pdb_path):
    """
    Returns a sorted list of unique ligand 3-letter codes in the PDB
    (excluding water, ions, and amino acids).
    """
    # ignore common non-standard atoms
    ignore = {
        "HOH", "WAT", "SOL",  # water
        "NA", "K", "CL", "CA", "MG", "ZN", "FE", "MN", "CU", "CO", "BR", "IOD",
        "SO4", "PO4", "PI", "ACE"  # ions/buffers
    }

    # amino acids, nucleotides, common cofactors
    ignore.update({
        "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
        "DA","DT","DG","DC","A","T","G","C","U"
    })

    ligs = set()

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname.upper() not in ignore:
                    ligs.add(resname.upper())

    return sorted(ligs)




def reorder_topol_for_charmm(topol_path):
    print("🔧 Reordering topology includes for CHARMM/CGenFF safety...")

    with open(topol_path) as f:
        lines = f.readlines()

    ff = []
    protein = []
    lig_prm = []
    lig_itp = []
    water = []
    ions = []
    rest = []

    for line in lines:
        l = line.strip()

        if l.startswith("#include"):
            if "forcefield.itp" in l:
                ff.append(line)
            elif "topol_Protein_chain" in l:
                protein.append(line)
            elif l.endswith(".prm\"") or l.endswith(".prm"):
                lig_prm.append(line)
            elif l.endswith(".itp\"") or l.endswith(".itp"):
                if "spc.itp" in l:
                    water.append(line)
                elif "ions.itp" in l:
                    ions.append(line)
                elif "dr7.itp" in l or "DR7.itp" in l:
                    lig_itp.append(line)
                else:
                    rest.append(line)
            else:
                rest.append(line)
        else:
            rest.append(line)

    new_lines = (
        ff +
        lig_prm +     # ✅ PRM must come before any molecule .itp blocks
        protein +
        lig_itp +
        water +
        ions +
        rest
    )

    with open(topol_path, "w") as f:
        f.writelines(new_lines)

    print("✅ topol.top reordered safely.")




def graft_coords_onto_ini(ini_pdb, pose_pdb, out_pdb):
    """
    Preserve the atom ORDER from ini_pdb (matches .itp) but replace XYZ with coords
    from pose_pdb (extracted from complex). We match records by atom name within the only residue.
    """
    def parse_pdb_atoms(path):
        atoms = []
        with open(path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    name = line[12:16].strip()
                    resn = line[17:20].strip()
                    resi = line[22:26].strip()
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    atoms.append({"raw": line, "name": name, "resn": resn, "resi": resi, "x": x, "y": y, "z": z})
        return atoms

    ini_atoms  = parse_pdb_atoms(ini_pdb)
    pose_atoms = parse_pdb_atoms(pose_pdb)

    # Build simple name->list mapping for pose (to handle duplicates like H11/H12 etc.)
    from collections import defaultdict, deque
    by_name = defaultdict(deque)
    for a in pose_atoms:
        by_name[a["name"]].append(a)

    out_lines = []
    with open(ini_pdb) as f:
        ini_lines = f.readlines()

    for line in ini_lines:
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            if not by_name[name]:
                # fallback: if name missing, keep original ini coords (warn)
                # You could also relax to element-based matching here.
                out_lines.append(line)
                continue
            src = by_name[name].popleft()
            # write XYZ into fixed-column PDB fields (30-54)
            new = (line[:30] + f"{src['x']:8.3f}{src['y']:8.3f}{src['z']:8.3f}" + line[54:])
            out_lines.append(new)
        else:
            out_lines.append(line)

    with open(out_pdb, "w") as f:
        f.writelines(out_lines)
    return True





import os, re, shutil, subprocess, tempfile

def _which(cmd):
    from shutil import which
    return which(cmd)

def _find_pymol_exe():
    """
    Prefer a real system PyMOL binary over conda stubs.
    Return executable path or None.
    """
    # Explicit common locations first (system install)
    for p in ("/usr/bin/pymol", "/usr/local/bin/pymol"):
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p

    # Then PATH
    p = _which("pymol")
    if p and os.path.isfile(p) and os.access(p, os.X_OK):
        return p

    # Sometimes packaged as pymol2
    p2 = _which("pymol2")
    if p2 and os.path.isfile(p2) and os.access(p2, os.X_OK):
        return p2

    return None

def count_h_in_pdb(pdb_path):
    n = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                aname = line[12:16].strip().upper()
                elem  = (line[76:78].strip().upper() if len(line) >= 78 else "")
                if elem == "H" or aname.startswith("H"):
                    n += 1
    return n

def count_h_in_mol2(mol2_path):
    n = 0
    in_atoms = False
    with open(mol2_path) as f:
        for line in f:
            s = line.strip()
            if s.upper() == "@<TRIPOS>ATOM":
                in_atoms = True
                continue
            if s.upper().startswith("@<TRIPOS>") and s.upper() != "@<TRIPOS>ATOM":
                in_atoms = False
            if in_atoms and s:
                parts = s.split()
                if len(parts) >= 6:
                    name = parts[1].upper()
                    atype = parts[5].upper()
                    if name.startswith("H") or atype.startswith("H"):
                        n += 1
    return n

def ensure_mol2_has_substructure(mol2_path, ligand_code):
    """
    Some writers omit SUBSTRUCTURE. CGenFF/convert scripts can be picky.
    If missing, append a minimal SUBSTRUCTURE block.
    """
    with open(mol2_path, "r") as f:
        txt = f.read()
    if "@<TRIPOS>SUBSTRUCTURE" in txt.upper():
        return

    # Find first atom index as root_atom
    root_atom = 1
    m = re.search(r"@<TRIPOS>ATOM\s+(\d+)\s+", txt, flags=re.IGNORECASE | re.MULTILINE)
    if m:
        try:
            root_atom = int(m.group(1))
        except:
            root_atom = 1

    with open(mol2_path, "a") as f:
        f.write("\n@<TRIPOS>SUBSTRUCTURE\n")
        # subst_id subst_name root_atom subst_type ...
        f.write(f"     1 {ligand_code} {root_atom} {ligand_code} 0\n")



def normalize_mol2_atom_resnames(mol2_path, ligand_code):
    """
    Force residue name column in @<TRIPOS>ATOM records to ligand_code.
    Fixes UNK1 vs UNK mismatches that break CGenFF → GROMACS conversion.
    """
    ligand_code = ligand_code.strip()

    with open(mol2_path) as f:
        lines = f.readlines()

    out = []
    in_atoms = False

    for line in lines:
        s = line.strip().upper()

        if s == "@<TRIPOS>ATOM":
            in_atoms = True
            out.append(line)
            continue

        if s.startswith("@<TRIPOS>") and s != "@<TRIPOS>ATOM":
            in_atoms = False

        if in_atoms and line.strip():
            parts = line.split()
            if len(parts) >= 8:
                # parts layout:
                # id name x y z type resid resname charge
                parts[7] = ligand_code   # ← THIS IS THE CRITICAL FIX
                line = " ".join(parts) + "\n"

        out.append(line)

    with open(mol2_path, "w") as f:
        f.writelines(out)


def split_paramchem_top(top_path, ligand_code, out_dir):
    """
    Split ParamChem GROMACS .top into:
      - ligand .itp (moleculetype + bonded terms)
      - ligand .prm (atomtypes / parameters)
    """
    lc = ligand_code.lower()

    itp_path = os.path.join(out_dir, f"{lc}.itp")
    prm_path = os.path.join(out_dir, f"{lc}.prm")

    with open(top_path) as f:
        lines = f.readlines()

    prm_blocks = []
    itp_blocks = []

    mode = None
    for line in lines:
        if line.strip().startswith("["):
            header = line.strip().lower()

            if header in (
                "[ atomtypes ]",
                "[ nonbond_params ]",
                "[ pairtypes ]"
            ):
                mode = "prm"
            elif header in (
                "[ moleculetype ]",
                "[ atoms ]",
                "[ bonds ]",
                "[ angles ]",
                "[ dihedrals ]",
                "[ impropers ]",
                "[ cmap ]"
            ):
                mode = "itp"
            else:
                mode = None

        if mode == "prm":
            prm_blocks.append(line)
        elif mode == "itp":
            itp_blocks.append(line)

    if not itp_blocks:
        raise RuntimeError("No ligand topology blocks found in ParamChem .top")

    with open(itp_path, "w") as f:
        f.writelines(itp_blocks)

    with open(prm_path, "w") as f:
        f.writelines(prm_blocks)

    return itp_path, prm_path


def normalize_paramchem_pdb(pdb_path, ligand_code, out_dir):
    lc = ligand_code.lower()
    out_pdb = os.path.join(out_dir, f"{lc}_ini.pdb")

    with open(pdb_path) as fin, open(out_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                # Force residue name consistency
                fout.write(
                    line[:17] + ligand_code.ljust(3) + line[20:]
                )
            else:
                fout.write(line)

    return out_pdb


def import_paramchem_gromacs(directory, ligand_code):
    top = os.path.join(directory, f"{ligand_code}_gmx.top")
    pdb = os.path.join(directory, f"{ligand_code}_gmx.pdb")

    if not (os.path.exists(top) and os.path.exists(pdb)):
        return None

    print("🔒 Importing ParamChem GROMACS ligand topology")

    ini_pdb = normalize_paramchem_pdb(pdb, ligand_code, directory)
    itp, prm = split_paramchem_top(top, ligand_code, directory)

    return {
        "ini_pdb": ini_pdb,
        "itp": itp,
        "prm": prm
    }




def normalize_mol2_names(mol2_path, ligand_code):
    """
    Ensure both:
      - @<TRIPOS>MOLECULE name line
      - @<TRIPOS>SUBSTRUCTURE residue/substructure name
    match ligand_code.
    """
    ligand_code = ligand_code.strip()

    with open(mol2_path, "r") as f:
        lines = f.readlines()

    out = []
    in_sub = False
    i = 0
    while i < len(lines):
        line = lines[i]

        # Fix MOLECULE name
        if line.strip().upper() == "@<TRIPOS>MOLECULE" and i + 1 < len(lines):
            out.append(line)
            out.append(f"{ligand_code}\n")
            i += 2
            continue

        # Track SUBSTRUCTURE section
        if line.strip().upper() == "@<TRIPOS>SUBSTRUCTURE":
            in_sub = True
            out.append(line)
            i += 1
            continue

        if line.strip().upper().startswith("@<TRIPOS>") and line.strip().upper() != "@<TRIPOS>SUBSTRUCTURE":
            in_sub = False

        # Rewrite SUBSTRUCTURE lines
        if in_sub and line.strip() and not line.lstrip().startswith("#"):
            parts = line.split()
            if len(parts) >= 2:
                parts[1] = ligand_code
            if len(parts) >= 4:
                parts[3] = ligand_code
            out.append(" ".join(parts) + "\n")
            i += 1
            continue

        out.append(line)
        i += 1

    with open(mol2_path, "w") as f:
        f.writelines(out)

def normalize_str_resi_name(str_path, ligand_code):
    ligand_code = ligand_code.strip()
    with open(str_path, "r") as f:
        lines = f.readlines()

    out = []
    fixed = False
    for line in lines:
        if (not fixed) and line.startswith("RESI"):
            parts = line.split()
            if len(parts) >= 3:
                charge = parts[2]
                out.append(f"RESI {ligand_code} {charge}\n")
                fixed = True
                continue
        out.append(line)

    with open(str_path, "w") as f:
        f.writelines(out)

def assert_str_has_atoms(str_path):
    n_resi = 0
    n_atom = 0
    with open(str_path, "r") as f:
        for line in f:
            if line.startswith("RESI"):
                n_resi += 1
            elif line.startswith("ATOM"):
                n_atom += 1
    if n_resi == 0 or n_atom == 0:
        raise RuntimeError(
            f"❌ Invalid CGenFF .str: RESI={n_resi}, ATOM={n_atom}\n"
            f"   This usually happens when CGenFF prints 'skipped molecule' due to bad valence/bond orders."
        )

def run_and_log(cmd, cwd, stdout_path=None, stderr_path=None):
    stdout_f = open(stdout_path, "w") if stdout_path else subprocess.PIPE
    stderr_f = open(stderr_path, "w") if stderr_path else subprocess.PIPE
    try:
        r = subprocess.run(cmd, cwd=cwd, text=True, stdout=stdout_f, stderr=stderr_f)
    finally:
        if stdout_path: stdout_f.close()
        if stderr_path: stderr_f.close()

    out_txt = open(stdout_path, "r").read() if stdout_path and os.path.exists(stdout_path) else ""
    err_txt = open(stderr_path, "r").read() if stderr_path and os.path.exists(stderr_path) else ""
    return r.returncode, out_txt, err_txt

def run_command_ok(command, cwd=None):
    """
    Non-fatal command runner (unlike your run_command which exits).
    Returns True/False.
    """
    r = subprocess.run(command, shell=True, cwd=cwd, text=True, capture_output=True)
    if r.returncode != 0:
        print(f"❌ Command failed: {command}\n🔺 STDERR:\n{r.stderr}")
        return False
    return True


def build_cgenff_inputs_realtime(directory, ligand_code, source_complex_pdb):
    """
    Updated workflow (PyMOL-headless first, OpenBabel fallback):

    1) PyMOL headless:
       - load complex
       - select ligand by resn ligand_code (non-polymer, non-solvent)
       - create object
       - save {code}_raw.pdb
       - add H (h_add)
       - save {code}_h.pdb
       - save {code}.cgenff.mol2

    2) Normalize MOL2 (MOLECULE + SUBSTRUCTURE names), ensure SUBSTRUCTURE exists.

    3) Run SILCSBio cgenff: stdout -> {code}.str, stderr -> {code}.cgenff.stderr.log
       - If skipped/unfulfilled/rc!=0 -> run verbose (-v) and fail.

    Returns (mol2_file, str_file) or (None, None).
    """
    cgenff_exec, prm36, silcs_home = ensure_silcsbio_env()
    if not cgenff_exec:
        return None, None

    ligand_code = ligand_code.strip().upper()
    base = ligand_code.lower()

    complex_pdb = os.path.join(directory, source_complex_pdb)
    raw_pdb     = os.path.join(directory, f"{base}_raw.pdb")
    h_pdb       = os.path.join(directory, f"{base}_h.pdb")
    mol2        = os.path.join(directory, f"{ligand_code}.cgenff.mol2")
    strout      = os.path.join(directory, f"{ligand_code}.str")

    if not os.path.exists(complex_pdb):
        print(f"❌ Complex PDB not found: {complex_pdb}")
        return None, None

    # ------------------------------------------------------------------
    # 1) PyMOL headless path (preferred)
    # ------------------------------------------------------------------
    pymol_exe = _find_pymol_exe()
    used_pymol = False

    if pymol_exe:
        # Build a temporary pml script
        pml_path = os.path.join(directory, f"{ligand_code}.extract_hadd.pml")
        pymol_stdout = os.path.join(directory, f"{ligand_code}.pymol.stdout.log")
        pymol_stderr = os.path.join(directory, f"{ligand_code}.pymol.stderr.log")

        # Selection: non-polymer, non-solvent, resn matches ligand_code
        # This catches GDP/ATP etc without grabbing protein residues.
        sel = f"(resn {ligand_code} and not polymer and not solvent)"

        pml = f"""
reinitialize
load {os.path.basename(complex_pdb)}, complex
select LIGSEL, {sel}
create {ligand_code}, LIGSEL
# Ensure residue name is exactly ligand_code
alter {ligand_code}, resn='{ligand_code}'
sort
save {os.path.basename(raw_pdb)}, {ligand_code}
h_add {ligand_code}
save {os.path.basename(h_pdb)}, {ligand_code}
set retain_order, 1
save {os.path.basename(mol2)}, {ligand_code}
quit
""".lstrip()

        with open(pml_path, "w") as f:
            f.write(pml)

        print(f"🧬 PyMOL headless extract + h_add → MOL2 using: {pymol_exe}")
        rc, out_txt, err_txt = run_and_log(
            [pymol_exe, "-cq", os.path.basename(pml_path)],
            cwd=directory,
            stdout_path=pymol_stdout,
            stderr_path=pymol_stderr
        )

        if rc == 0 and os.path.exists(mol2) and os.path.getsize(mol2) > 0:
            nH_pdb = count_h_in_pdb(h_pdb) if os.path.exists(h_pdb) else 0
            nH_m2  = count_h_in_mol2(mol2)
            print(f"✅ PyMOL produced MOL2. H count: PDB={nH_pdb}, MOL2={nH_m2}")
            if nH_m2 == 0:
                print("⚠️ PyMOL MOL2 has 0 hydrogens — will try OpenBabel fallback.")
            else:
                used_pymol = True
        else:
            print("⚠️ PyMOL path failed; falling back to OpenBabel.")
            print(f"   See logs: {os.path.basename(pymol_stdout)}, {os.path.basename(pymol_stderr)}")
    else:
        print("ℹ️ No usable PyMOL executable found; using OpenBabel fallback.")

    # ------------------------------------------------------------------
    # 1b) OpenBabel fallback if PyMOL not used/successful
    # ------------------------------------------------------------------
    if not used_pymol:
        # Use your existing extractor that preserves CONECT if possible
        if not extract_ligand_from_pdb(complex_pdb, ligand_code, raw_pdb):
            return None, None

        # Always attempt to add H, and VERIFY it happened
        if not run_command_ok(f"obabel -ipdb {os.path.basename(raw_pdb)} -opdb -O {os.path.basename(h_pdb)} -h -p 7.4", cwd=directory):
            return None, None

        nH = count_h_in_pdb(h_pdb)
        print(f"🧪 OpenBabel H count after -h: {nH}")
        if nH == 0:
            # stronger fallback: rebuild chemistry in SDF then back to PDB
            tmp_sdf = os.path.join(directory, f"{base}_tmp.sdf")
            if not run_command_ok(f"obabel -ipdb {os.path.basename(raw_pdb)} -osdf -O {os.path.basename(tmp_sdf)} --gen3d -h -p 7.4", cwd=directory):
                return None, None
            if not run_command_ok(f"obabel -isdf {os.path.basename(tmp_sdf)} -opdb -O {os.path.basename(h_pdb)}", cwd=directory):
                return None, None
            nH = count_h_in_pdb(h_pdb)
            print(f"🧪 OpenBabel H count after SDF fallback: {nH}")
            if nH == 0:
                print("❌ Still 0 H after fallback — cannot build a chemically valid ligand from this PDB fragment.")
                return None, None

        # MOL2 from hydrogenated PDB
        if not run_command_ok(f"obabel -ipdb {os.path.basename(h_pdb)} -omol2 -O {os.path.basename(mol2)}", cwd=directory):
            return None, None

        nH_m2 = count_h_in_mol2(mol2)
        print(f"🧪 OpenBabel MOL2 H count: {nH_m2}")
        if nH_m2 == 0:
            print("❌ MOL2 has 0 H — aborting.")
            return None, None

    # ------------------------------------------------------------------
    # 2) Normalize MOL2 naming and ensure SUBSTRUCTURE exists
    # ------------------------------------------------------------------
    try:
        ensure_mol2_has_substructure(mol2, ligand_code)
        normalize_mol2_names(mol2, ligand_code)
        normalize_mol2_atom_resnames(mol2, ligand_code)  # ← ADD THIS

    except Exception as e:
        print(f"❌ MOL2 normalization failed: {e}")
        return None, None

    # Optional: record perceived SMILES for debugging
    ob_smiles = os.path.join(directory, f"{ligand_code}.perceived.smi")
    ob_err    = os.path.join(directory, f"{ligand_code}.perceived.stderr.log")
    run_and_log(["obabel", os.path.basename(mol2), "-osmi"], cwd=directory,
                stdout_path=ob_smiles, stderr_path=ob_err)
    print(f"🧾 Perceived SMILES written: {os.path.basename(ob_smiles)}")

    # ------------------------------------------------------------------
    # 3) Run SILCSBio cgenff (stdout -> .str)
    # ------------------------------------------------------------------
    cgenff_stderr  = os.path.join(directory, f"{ligand_code}.cgenff.stderr.log")
    cgenff_verbose = os.path.join(directory, f"{ligand_code}.cgenff.verbose.log")

    print("\n🧪 Running SILCSBio CGenFF parameterization...")
    print(f"⚙️ Command: {cgenff_exec} {os.path.basename(mol2)}")
    print(f"🧾 Output .str: {os.path.basename(strout)}")

    rc, out_txt, err_txt = run_and_log(
        [cgenff_exec, os.path.basename(mol2)],
        cwd=directory,
        stdout_path=strout,
        stderr_path=cgenff_stderr
    )

    if ("skipped molecule" in err_txt.lower()) or ("unfulfilled valence" in err_txt.lower()) or rc != 0:
        print("⚠️ CGenFF reported a chemistry/connectivity problem.")
        print(f"   See: {os.path.basename(cgenff_stderr)}")
        print("🔎 Running verbose diagnostics (-v)...")
        run_and_log(
            [cgenff_exec, "-v", os.path.basename(mol2)],
            cwd=directory,
            stdout_path=cgenff_verbose,
            stderr_path=cgenff_stderr
        )
        print(f"🧾 Verbose log: {os.path.basename(cgenff_verbose)}")
        return None, None

    if not os.path.exists(strout) or os.path.getsize(strout) == 0:
        print("❌ CGenFF produced an EMPTY .str file.")
        print(f"   See: {os.path.basename(cgenff_stderr)}")
        return None, None

    # Normalize RESI name and hard-check ATOM presence
    try:
        normalize_str_resi_name(strout, ligand_code)
        assert_str_has_atoms(strout)
    except Exception as e:
        print(str(e))
        return None, None

    print(f"✅ CGenFF successfully generated {os.path.basename(strout)}")
    return mol2, strout





def parse_make_ndx_for_water(stdout_text):
    """
    Robustly parse `gmx make_ndx` listing to find the group number for SOL/Water/WAT.
    Returns the group number as a string if found, else None.
    Lines look like:
      15 SOL                 : 55782 atoms
    or
      14 Water               : 55782 atoms
    """
    candidates_exact = {"SOL", "Water", "WAT"}
    # regex: capture group number and group name up to the colon
    pat = re.compile(r'^\s*(\d+)\s+([A-Za-z0-9_\-\+ ]+?)\s*:', re.MULTILINE)
    matches = pat.findall(stdout_text or "")
    if not matches:
        return None

    # 1) exact priority: SOL > Water > WAT
    for idx, name in matches:
        name_clean = name.strip()
        if name_clean == "SOL":
            return idx
    for idx, name in matches:
        if name.strip() == "Water":
            return idx
    for idx, name in matches:
        if name.strip() == "WAT":
            return idx

    # 2) contains-based fallback
    for idx, name in matches:
        name_u = name.strip().upper()
        if any(tag in name_u for tag in ("SOL", "WATER", "WAT", "TIP3", "TIP4", "TIP5", "HOH")):
            return idx

    # 3) last resort: pick the largest water-like group by atom count (if parsable)
    #   (parse "<idx> <name> : <count> atoms")
    pat_count = re.compile(r'^\s*(\d+)\s+([A-Za-z0-9_\-\+ ]+?)\s*:\s*(\d+)\s+atoms', re.MULTILINE)
    best = None
    best_count = -1
    for idx, name, cnt in pat_count.findall(stdout_text or ""):
        name_u = name.strip().upper()
        if any(tag in name_u for tag in ("SOL", "WATER", "WAT", "TIP3", "TIP4", "TIP5", "HOH")):
            c = int(cnt)
            if c > best_count:
                best_count = c
                best = idx
    return best


def build_water_index_with_gmx_select(directory):
    """
    Fallback: build a dedicated water-only index file using gmx select.
    Returns (ndx_path, group_number) or (None, None) on failure.
    We use ions.tpr (exists after the ions grompp) and solv.gro.
    """
    ndx_path = os.path.join(directory, "_water.ndx")
    # select resnames commonly used for water:
    selection = 'resname SOL or resname WAT or resname HOH or resname TIP3 or resname TIP4 or resname TIP5'
    ok = run_command(
        gmx_cmd("select", f'-s ions.tpr -f solv.gro -select "{selection}" -on {os.path.basename(ndx_path)}'),
        cwd=directory
    )
    if not ok:
        return None, None

    # gmx select typically writes one group named "Selection" -> group 0
    # (If multiple groups were made, the first is still a safe default)
    try:
        with open(ndx_path, "r") as f:
            content = f.read()
        # sanity check it actually contains a group
        if "[ " in content and "]" in content:
            return ndx_path, "0"
    except Exception:
        pass
    return None, None


def _mean(vals):
    return sum(vals)/len(vals) if vals else None

def centroid_from_gro(path, first=2, last=-1):
    xs = []; ys = []; zs = []
    try:
        with open(path) as f:
            lines = f.readlines()[first:last]
    except FileNotFoundError:
        return (None, None, None)
    for ln in lines:
        try:
            xs.append(float(ln[20:28])); ys.append(float(ln[28:36])); zs.append(float(ln[36:44]))
        except:
            pass
    return (_mean(xs), _mean(ys), _mean(zs))




def is_cgenff_ljpme_conflict(stderr):
    signatures = [
        "No default Bond types",
        "No default Proper Dih",
        "No default U-B types",
        "Protein_chain"
    ]
    return all(s in stderr for s in signatures)



# ======== END NEW HELPERS ========



def process_directory(directory, pdb_filename, ligand_code):
    """Processes a directory for MD setup using a selected PDB file and ligand code."""
    print(f"\n🔄 Processing directory: {directory}")
    print_gmx_preflight(directory)
    pdb_path = os.path.join(directory, pdb_filename)
    print(f"\n🔍 Checking for {pdb_filename} in the working directory...")
    if not os.path.exists(pdb_path):
        print(f"⚠️ Skipping {directory}: {pdb_filename} not found.")
        return False

    print(
        f"📦 Box configuration: type={args.box_type}, distance={args.box_distance:.3f} nm"
    )
    append_setup_log(
        directory,
        f"Box configuration: type={args.box_type}, distance={args.box_distance:.3f} nm",
    )

    try:
        mdp_files = prepare_step1_mdp_files(directory, args)
    except FileNotFoundError as exc:
        print(f"❌ {exc}")
        return False

    # Step 1: Remove ligand or copy original PDB based on presence of ligand
    if ligand_code:
        print(f"🛠️ Removing ligand ({ligand_code}) to create protein.pdb...")
        if not run_command(f"grep -v '{ligand_code}' {pdb_filename} > protein.pdb", cwd=directory):
            return False
    else:
        print("🛠️ No ligand provided. Using the original PDB as protein.pdb")
        try:
            shutil.copy(pdb_path, os.path.join(directory, "protein.pdb"))
            print("✅ Copied original PDB to protein.pdb")
        except Exception as e:
            print(f"❌ Failed to copy {pdb_filename} to protein.pdb: {e}")
            return False

    # Step 1.5: Clean up missing atoms and sidechains
    protein_pdb_path = os.path.join(directory, "protein.pdb")
    fixed_pdb_path = os.path.join(directory, "protein_fixed.pdb")
    print(f"🧼 Fixing missing atoms and adding hydrogens with PDBFixer...")
    if not fix_missing_atoms(protein_pdb_path, fixed_pdb_path):
        return False
    os.replace(fixed_pdb_path, protein_pdb_path)

    # Step 2: Generate protein topology
    if not automate_pdb2gmx(directory):
        return False

    # Step 3: Generate atom index file
    if not generate_atom_index_file(
                                    directory,
                                    chain_names_arg=args.chain_names,
                                    chain_map_arg=args.chain_map,
                                    source_pdb_path=pdb_path
                                ):
        print("⚠️ Failed to generate atom index file.")

        return False
    
    # ============================================================
    # Step 4: Generate ligand topology (if applicable)
    # ============================================================
    if ligand_code:
        print("\n🛠️ Generating ligand topology...")
        forcefield_path = detect_forcefield(directory)
        if forcefield_path is None:
            return False



        # ============================================================
        # Step 4: Generate ligand topology (if applicable)
        # ============================================================
        if ligand_code:
            print("\n🛠️ Generating ligand topology...")
            forcefield_path = detect_forcefield(directory)
            if forcefield_path is None:
                return False

            ligand_gro_file = None
            ini_pdb_file    = None

            # ------------------------------------------------------------
            # 4a0) ParamChem GROMACS import (HIGHEST PRIORITY)
            # ------------------------------------------------------------
            paramchem = import_paramchem_gromacs(directory, ligand_code)

            if paramchem:
                print("🔒 Using ParamChem GROMACS ligand topology")

                ini_pdb_file = paramchem["ini_pdb"]
                ligand_gro_file = os.path.join(directory, f"{ligand_code.lower()}.gro")

                # Convert ParamChem PDB → GRO
                if not run_command(
                    gmx_cmd(
                        "editconf",
                        f"-f {os.path.basename(ini_pdb_file)} -o {os.path.basename(ligand_gro_file)}",
                    ),
                    cwd=directory
                ):
                    return False

                # Optional sanity check
                prot_ctr = centroid_from_gro(os.path.join(directory, "protein_processed.gro"))
                lig_ctr  = centroid_from_gro(ligand_gro_file)
                print(f"🔎 Centroids (nm) — protein: {prot_ctr}, ligand: {lig_ctr}")

            else:
                # ------------------------------------------------------------
                # 4a) Ensure CGenFF inputs exist
                # ------------------------------------------------------------
                candidate_mol2 = None
                for fname in (f"{ligand_code}.cgenff.mol2", f"{ligand_code}.mol2"):
                    p = os.path.join(directory, fname)
                    if os.path.exists(p):
                        candidate_mol2 = p
                        break

                str_file = os.path.join(directory, f"{ligand_code}.str")

                if candidate_mol2 and os.path.exists(str_file):
                    mol2_file = candidate_mol2
                    print(f"✅ Using user-provided CGenFF inputs: "
                        f"{os.path.basename(mol2_file)}, {os.path.basename(str_file)}")
                else:
                    print("ℹ️ Missing CGenFF inputs — attempting SILCSBio generation...")

                    cgenff_exec, _, _ = ensure_silcsbio_env()
                    if not cgenff_exec:
                        print(
                            "❌ No SILCSBio environment AND no pre-generated CGenFF inputs.\n"
                            "   → Generate ligand parameters externally (ParamChem)\n"
                            "   → Place them in this directory and rerun."
                        )
                        return False

                    mol2_file, str_file = build_cgenff_inputs_realtime(
                        directory, ligand_code, pdb_filename
                    )
                    if not (mol2_file and str_file):
                        print("❌ Local CGenFF generation failed.")
                        return False

                # ------------------------------------------------------------
                # 4a.5) Normalize + HARD FAIL if .str is invalid
                #       (ALWAYS for CGenFF paths)
                # ------------------------------------------------------------
                try:
                    normalize_mol2_names(mol2_file, ligand_code)
                    normalize_str_resi_name(str_file, ligand_code)
                    assert_str_has_atoms(str_file)
                except Exception as e:
                    print(str(e))
                    return False

                # ------------------------------------------------------------
                # 4b) Convert CGenFF → GROMACS
                # ------------------------------------------------------------
                if not run_command(
                    f"python cgenff_charmm2gmx_py3_nx2.py {ligand_code} "
                    f"{os.path.basename(mol2_file)} {os.path.basename(str_file)} {forcefield_path}",
                    cwd=directory
                ):
                    print("❌ CGenFF → GROMACS conversion failed.")
                    return False

                try:
                    validate_cgenff_outputs(ligand_code, directory)
                except RuntimeError as e:
                    print(str(e))
                    return False

                # ------------------------------------------------------------
                # 4c) Pose-matched coordinates
                # ------------------------------------------------------------
                ini_pdb_file = os.path.join(directory, f"{ligand_code.lower()}_ini.pdb")
                pose_pdb    = os.path.join(directory, f"{ligand_code.lower()}_pose_match.pdb")
                ligand_gro_file = os.path.join(directory, f"{ligand_code.lower()}.gro")

                pose_source_pdb = os.path.join(directory, f"{ligand_code.lower()}_h.pdb")
                if not os.path.exists(pose_source_pdb):
                    extracted_pose = os.path.join(directory, f"{ligand_code.lower()}_from_pdb.pdb")
                    if not run_command(
                        f"grep '{ligand_code}' {pdb_filename} > {os.path.basename(extracted_pose)}",
                        cwd=directory
                    ):
                        return False
                    pose_source_pdb = extracted_pose

                if not graft_coords_onto_ini(ini_pdb_file, pose_source_pdb, pose_pdb):
                    print("❌ Pose grafting failed.")
                    return False

                # ------------------------------------------------------------
                # 4d) Ensure hydrogens
                # ------------------------------------------------------------
                pose_final = pose_pdb
                if not ligand_has_hydrogens(pose_pdb):
                    pose_h_pdb = os.path.join(directory, f"{ligand_code.lower()}_pose_h.pdb")
                    if not run_command(
                        f"obabel {os.path.basename(pose_pdb)} "
                        f"-O {os.path.basename(pose_h_pdb)} -h",
                        cwd=directory
                    ):
                        return False
                    pose_final = pose_h_pdb

                # ------------------------------------------------------------
                # 4e) Generate ligand .gro
                # ------------------------------------------------------------
                if not run_command(
                    gmx_cmd(
                        "editconf",
                        f"-f {os.path.basename(pose_final)} -o {os.path.basename(ligand_gro_file)}",
                    ),
                    cwd=directory
                ):
                    return False

                prot_ctr = centroid_from_gro(os.path.join(directory, "protein_processed.gro"))
                lig_ctr  = centroid_from_gro(ligand_gro_file)
                print(f"🔎 Centroids (nm) — protein: {prot_ctr}, ligand: {lig_ctr}")

        else:
            print("ℹ️ No ligand present, skipping ligand topology generation.")



    # Step 5: Merge protein and ligand structures
    if ligand_code:
        print("\n🔄 Merging protein and ligand structures...")
        protein_gro = os.path.join(directory, "protein_processed.gro")
        ligand_gro = os.path.join(directory, ligand_gro_file)
        complex_gro = os.path.join(directory, "complex.gro")
        merge_gro_files(protein_gro, ligand_gro, complex_gro)
    else:
        print("ℹ️ Skipping merge step (protein-only system). Using protein structure as complex.")
        protein_gro = os.path.join(directory, "protein_processed.gro")
        complex_gro = os.path.join(directory, "complex.gro")
        shutil.copy(protein_gro, complex_gro)
        print("✅ Single protein structure saved as complex.gro")

    # Step 6: Modify topology
    if ligand_code:
        print("\n🛠️ Modifying topology file (topol.top)...")
        topol_path = os.path.join(directory, "topol.top")
        modify_topology_file(topol_path, ligand_code)
    else:
        print("ℹ️ Skipping topology file modifications for ligand.")

    # ============================================================
    # Step 7–11: GROMACS setup steps
    # ============================================================

    print("\n⚙️ Running GROMACS setup steps...")

    # ------------------------------------------------------------
    # Step 7: Define box
    # ------------------------------------------------------------
    if not run_command(
        gmx_cmd("editconf", f"-f complex.gro -o newbox.gro -bt {args.box_type} -d {args.box_distance}"),
        cwd=directory
    ):
        return False

    # ------------------------------------------------------------
    # Step 8: Solvate
    # ------------------------------------------------------------
    if not run_command(
        gmx_cmd("solvate", "-cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro"),
        cwd=directory
    ):
        return False

    # ------------------------------------------------------------
    # Step 9: grompp for ions (with auto-patching)
    # ------------------------------------------------------------
    print("⚙️ Step 9: Running grompp for ion preparation...")

    grompp_result = subprocess.run(
        gmx_cmd("grompp", f"-f {os.path.basename(mdp_files['ions.mdp'])} -c solv.gro -p topol.top -o ions.tpr -maxwarn 2"),
        shell=True,
        cwd=directory,
        text=True,
        capture_output=True
    )

    if grompp_result.returncode != 0:
        stderr = grompp_result.stderr
        print("❌ grompp failed during ion preparation")
        print(stderr)

        # ---- CASE 1: CHARMM36-LJPME + CGenFF incompatibility ----
        if is_cgenff_ljpme_conflict(stderr):
            print("💥 Detected CHARMM36-LJPME + CGenFF incompatibility")
            print("🧹 Removing incompatible force field and restarting with classic CHARMM36")

            ljpme_ff = os.path.join(directory, "charmm36_ljpme-jul2022.ff")
            if os.path.isdir(ljpme_ff):
                shutil.rmtree(ljpme_ff)
                print("🗑️ Deleted charmm36_ljpme-jul2022.ff")

            print("🔁 Re-running pipeline with classic charmm36.ff")
            return regenerate_protein_topology(forcefield="charmm36.ff")

        # ---- CASE 2: Legit include-order issue (rare) ----
        elif (
            "No default Bond types" in stderr and
            "Protein_chain" in stderr
        ):
            print("🧠 Detected CHARMM include-order issue")
            reorder_topol_for_charmm(os.path.join(directory, "topol.top"))

            print("🔁 Retrying grompp after reorder...")
            return run_command(
                gmx_cmd("grompp", f"-f {os.path.basename(mdp_files['ions.mdp'])} -c solv.gro -p topol.top -o ions.tpr -maxwarn 2"),
                cwd=directory
            )

        # ---- CASE 3: Anything else is fatal ----
        else:
            print("❌ Unrecoverable grompp error")
            return False



 
    # ------------------------------------------------------------
    # Step 10: Detect solvent group for genion
    # ------------------------------------------------------------
    make_ndx_result = subprocess.run(
        gmx_cmd("make_ndx", "-f solv.gro"),
        shell=True,
        cwd=directory,
        text=True,
        input="q\n",
        capture_output=True
    )

    sol_group = None
    ndx_arg = ""

    if make_ndx_result.returncode == 0:
        sol_group = parse_make_ndx_for_water(make_ndx_result.stdout)

    if sol_group is None:
        print("❗ Could not detect SOL group via make_ndx. Trying gmx select fallback...")
        ndx_path, group0 = build_water_index_with_gmx_select(directory)
        if ndx_path and group0 is not None:
            ndx_arg = f"-n {os.path.basename(ndx_path)}"
            sol_group = group0

    if sol_group is None:
        print("❌ Failed to auto-detect water group.")
        print(make_ndx_result.stdout)
        sol_group = input("🔢 Enter SOL/Water group number manually: ").strip()

    # ------------------------------------------------------------
    # Step 11: Add ions
    # ------------------------------------------------------------
    salt_m = os.environ.get("GMX_SALT_M", "0.15")

    if not run_command(
        gmx_cmd(
            "genion",
            f"-s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc {salt_m} -seed 2025 {ndx_arg}",
        ),
        cwd=directory,
        input_text=f"{sol_group}\n"
    ):
        return False

    # ------------------------------------------------------------
    # Final grompp (energy minimization)
    # ------------------------------------------------------------
    if not run_command(
        gmx_cmd("grompp", f"-f {os.path.basename(mdp_files['em.mdp'])} -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2"),
        cwd=directory
    ):
        return False

    print("\n✅ All setup steps completed successfully.")
    return True





if __name__ == "__main__":
    if len(os.sys.argv) == 1:
        print("ℹ️ Running interactively — you can also provide:")
        print("   --chain-names 'Rap1B,Rap1GAP' or --chain-map 'A:Rap1B,B:Rap1GAP'\n")

    base_directory = os.getcwd()
    print("\n📂 Base Directory:", base_directory)
    if args.pdb:
        pdb_filename = args.pdb
        print(f"✅ Using specified PDB file: {pdb_filename}")
    else:
        pdb_filename = select_pdb_file(base_directory)

    pdb_path = os.path.join(base_directory, pdb_filename)
    try:
        validate_input_pdb(
            pdb_path,
            ligand_code=args.ligand.upper() if args.ligand else None,
            strict=args.strict_pdb_validation,
            allow_covalent=args.allow_covalent_ligand,
        )
    except Exception as exc:
        print(f"❌ PDB preflight failed: {exc}")
        exit(1)

    # Ask user whether a ligand is present in the system
    
    # --- AUTO-DETECT LIGANDS FROM THE PDB ---
    detected_ligs = autodetect_ligands(pdb_path)

    ligand_code = None  # default

    if args.ligand:
        ligand_code = args.ligand.upper()
        print(f"\n✅ Ligand mode detected from CLI: {ligand_code}")

    elif args.mode == "protein":
        ligand_code = None
        print("\n✅ Protein-only mode forced by CLI. Skipping ligand detection.")

    else:
        # If script is interactive OR mode not pre-selected
        if detected_ligs:
            print("\n🔍 Auto-detected ligand(s) in the PDB:")
            for lig in detected_ligs:
                print(f"   • {lig}")

            if len(detected_ligs) == 1:
                choice = input(f"\nUse detected ligand '{detected_ligs[0]}'? [Y/n]: ").strip().lower()
                if choice in ("", "y", "yes"):
                    ligand_code = detected_ligs[0]
                else:
                    ligand_code = input("Enter 3-letter ligand code manually: ").strip().upper()
            else:
                print("\nMultiple ligands detected.")
                print("Enter one of these codes, or type a different ligand manually:")
                ligand_code = input(f"Choose ligand [{', '.join(detected_ligs)}]: ").strip().upper()

                if ligand_code == "":
                    ligand_code = detected_ligs[0]  # default first
        else:
            print("\nℹ️ No ligand detected in the PDB.")
            use_lig = input("Run as protein-only? [Y/n]: ").strip().lower()

            if use_lig in ("n", "no"):
                ligand_code = input("Enter 3-letter ligand code manually: ").strip().upper()
            else:
                ligand_code = None

    # Final report
    if ligand_code:
        print(f"\n✅ Using ligand: {ligand_code}")
    else:
        print("\n✅ No ligand selected for this simulation.")

    try:
        validate_input_pdb(
            pdb_path,
            ligand_code=ligand_code,
            strict=args.strict_pdb_validation,
            allow_covalent=args.allow_covalent_ligand,
        )
    except Exception as exc:
        print(f"❌ PDB validation failed for the selected workflow: {exc}")
        exit(1)



    # Display chosen options
    print(f"\n✅ Selected PDB file: {pdb_filename}")
    if ligand_code:
        print(f"✅ Using ligand code: {ligand_code}")
    else:
        print("✅ No ligand selected for this simulation.")
    print(f"✅ Box type: {args.box_type}")
    print(f"✅ Box distance: {args.box_distance:.3f} nm")
    if not process_directory(base_directory, pdb_filename, ligand_code):
        print("⚠️ Error detected. Stopping execution.")
        exit(1)
