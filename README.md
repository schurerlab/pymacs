# 🧬 PyMACS

<p align="center">
  <strong>PyMACS: A Python-Based Automation Suite for GROMACS Molecular Dynamics Setup, Simulation, Analysis, and Figurebook Generation</strong>
</p>

<p align="center">
  <em>A modular molecular-dynamics workflow toolkit built around GROMACS, CHARMM36, CGenFF ligand support, trajectory analysis, and publication-ready reporting.</em>
</p>

<p align="center">
  <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6584500">
    <img src="https://img.shields.io/badge/Read%20the%20Preprint-SSRN-success?style=for-the-badge&logo=readthedocs" alt="Read the PyMACS preprint">
  </a>
  <a href="https://github.com/schurerlab/Pymacs/tree/main/Example_Choices">
    <img src="https://img.shields.io/badge/View%20Example-GitHub%20Folder-orange?style=for-the-badge&logo=github" alt="View the PyMACS example folder">
  </a>
  <a href="https://github.com/schurerlab/Pymacs/blob/main/Example_Choices/Example1/completed_run/MD_ANALYSIS_FIGUREBOOK.pdf">
    <img src="https://img.shields.io/badge/View%20Example-Figurebook%20PDF-blueviolet?style=for-the-badge&logo=adobeacrobatreader" alt="View the example MD analysis figurebook PDF">
  </a>
</p>

<p align="center">
  <a href="mailto:jmschulz@med.miami.edu?subject=PyMACS%20Question%20%2F%20Collaboration">
    <img src="https://img.shields.io/badge/Contact%20the%20Author-Joseph%20M.%20Schulz-blue?style=for-the-badge&logo=gmail" alt="Contact Joseph M. Schulz">
  </a>
  <a href="mailto:sschurer@miami.edu?subject=PyMACS%20Correspondence">
    <img src="https://img.shields.io/badge/Corresponding%20Author-Stephan%20C.%20Sch%C3%BCrer-lightgrey?style=for-the-badge&logo=gmail" alt="Contact Stephan C. Schürer">
  </a>
</p>

---

<a id="authors"></a>

## 👥 Authors

**Joseph M. Schulz¹, Robert C. Reynolds², Stephan C. Schürer\*¹˒³˒⁴**

¹ Department of Molecular and Cellular Pharmacology, University of Miami, Miami, FL, USA  
² O’Neal Comprehensive Cancer Center, University of Alabama at Birmingham, Birmingham, AL 35205, USA  
³ Sylvester Comprehensive Cancer Center, University of Miami Miller School of Medicine, Miami, FL, USA  
⁴ Frost Institute for Data Science & Computing, University of Miami, Miami, FL, USA  

\*Corresponding author: **sschurer@miami.edu**

**Repository:** `schurerlab/Pymacs`  
**Status:** Active development; preprint available on SSRN.

---

<a id="preprint"></a>

## 📚 Preprint

**PyMACS: A Python-Based Automation Suite for GROMACS Molecular Dynamics Setup, Simulation, and Analysis**  
Joseph M. Schulz, Robert C. Reynolds, Stephan C. Schürer

<p align="center">
  <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6584500">
    <img src="https://img.shields.io/badge/Read%20the%20Preprint-SSRN-success?style=for-the-badge&logo=readthedocs" alt="Read the PyMACS preprint">
  </a>
</p>

---

<a id="contact-and-collaboration"></a>

## 🤝 Contact and Collaboration

For questions about PyMACS, molecular dynamics workflow support, CGenFF setup, GROMACS automation, figurebook generation, or related collaborations, contact the authors below.

<p align="center">
  <a href="mailto:jmschulz@med.miami.edu?subject=PyMACS%20Question%20%2F%20Collaboration">
    <img src="https://img.shields.io/badge/Questions%20%2F%20Collaboration-Joseph%20M.%20Schulz-blue?style=for-the-badge" alt="Contact Joseph M. Schulz">
  </a>
  <a href="mailto:jmschulz@med.miami.edu?subject=PyMACS%20Workflow%20Support">
    <img src="https://img.shields.io/badge/Workflow%20Support-GROMACS%20%2B%20CGenFF-orange?style=for-the-badge" alt="PyMACS workflow support">
  </a>
  <a href="mailto:sschurer@miami.edu?subject=PyMACS%20Correspondence">
    <img src="https://img.shields.io/badge/Correspondence-Sch%C3%BCrer%20Lab-lightgrey?style=for-the-badge" alt="Schürer Lab correspondence">
  </a>
</p>

---

<p align="center">
  <img 
    src="images/Pymacs_Logo.png" 
    alt="PyMACS logo" 
    width="420"
  >
</p>

<p align="center">
  <strong>Automated molecular dynamics. Clean analysis. Publication-ready figurebooks.</strong>
</p>

<p align="center">
  <em>From structure preparation to GROMACS execution, trajectory analysis, and polished PDF reporting.</em>
</p>

---

<a id="super-quick-start-install"></a>

## ⚡ Super Quick Start Install (lightweight)

Same idea as the full Quick Start below, except this version skips the heavy example assets and avoids Git LFS downloads. Use it when you want the core PyMACS scripts, templates, force fields, docs, and environment files without copying large trajectories or completed example outputs.

<details>
<summary><strong>Show lightweight installer command</strong></summary>

```bash
bash <<'EOF'
set -e

REPO_URL="https://github.com/schurerlab/pymacs.git"
TMP_DIR="$(mktemp -d)"
TARGET_DIR="$(pwd)"

cleanup() {
  rm -rf "$TMP_DIR"
}

trap cleanup EXIT

echo "======================================"
echo " PyMACS Super Quick Start Install"
echo "======================================"
echo
echo "This will copy a lightweight PyMACS instance into:"
echo "  $TARGET_DIR"
echo
echo "This version skips Example_Choices and common heavy MD output files."
echo

if ! command -v git >/dev/null 2>&1; then
  echo "ERROR: git is not installed or not available in PATH."
  echo "Please install git first, then run this command again."
  exit 1
fi

if [ "$(find "$TARGET_DIR" -mindepth 1 -maxdepth 1 | wc -l)" -gt 0 ]; then
  echo "WARNING: This directory is not empty."
  echo "Files with the same names as PyMACS files may be overwritten."
  echo
  printf "Continue copying lightweight PyMACS files into this directory? [y/N]: "
  read -r answer
  case "$answer" in
    y|Y|yes|YES) ;;
    *)
      echo "Install cancelled."
      exit 0
      ;;
  esac
fi

echo
echo "Cloning PyMACS into a temporary folder without Git LFS payloads..."
GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 "$REPO_URL" "$TMP_DIR/pymacs"

cd "$TMP_DIR/pymacs"

echo "Removing large example/output assets from the temporary copy..."
rm -rf Example_Choices
find . -type f \( \
  -name "*.xtc" -o \
  -name "*.trr" -o \
  -name "*.tpr" -o \
  -name "*.cpt" -o \
  -name "*.edr" -o \
  -name "*.xvg" -o \
  -name "*.log" -o \
  -name "*.pdf" -o \
  -name "*.zip" -o \
  -name "*.tar" -o \
  -name "*.tar.gz" \
\) -delete

echo
echo "Copying lightweight PyMACS files into your current directory..."

shopt -s dotglob nullglob
for item in "$TMP_DIR/pymacs"/*; do
  base="$(basename "$item")"
  case "$base" in
    .git|Example_Choices)
      continue
      ;;
  esac

  if [ -e "$TARGET_DIR/$base" ]; then
    printf "Overwrite existing %s? [y/N]: " "$base"
    read -r overwrite
    case "$overwrite" in
      y|Y|yes|YES)
        rm -rf "$TARGET_DIR/$base"
        ;;
      *)
        echo "Skipping $base"
        continue
        ;;
    esac
  fi

  cp -R "$item" "$TARGET_DIR/"
done

echo "Cleaning up temporary clone..."
cleanup
trap - EXIT

cd "$TARGET_DIR"

echo
echo "Done. Lightweight PyMACS has been copied into:"
echo "  $TARGET_DIR"
echo
echo "Next recommended steps:"
echo "  conda env create -f environment_cgenff.yml"
echo "  conda env create -f environment_mdanalysis.yml"
echo
echo "Then follow the Step 1 / Step 2 / Step 3 workflow in this README."
echo "Use the full Quick Start below if you also need the curated example datasets."
EOF
```

</details>

---

<a id="quick-start-install"></a>

## ⚡ Quick Start Install

Use this when you want to create a clean PyMACS run or working folder without manually downloading the repository or copying files yourself. First, `cd` into the folder where you want the PyMACS files to appear. Your "current working directory" simply means the folder your terminal is presently pointing at. The command below temporarily clones PyMACS, copies everything into that folder, includes hidden files such as `.gitignore` and `.gitattributes`, skips the cloned repository's `.git` folder, and then removes the temporary clone when it is finished.

It also checks for `git`, warns before copying into a non-empty directory, asks for confirmation before overwriting similarly named files, and handles Git LFS when available. If Git LFS is not installed, the installer will continue and warn you that some large example assets may not be fully downloaded yet.

```bash
bash <<'EOF'
set -e

REPO_URL="https://github.com/schurerlab/pymacs.git"
TMP_DIR="$(mktemp -d)"
TARGET_DIR="$(pwd)"

cleanup() {
  rm -rf "$TMP_DIR"
}

trap cleanup EXIT

echo "======================================"
echo " PyMACS Quick Start Install"
echo "======================================"
echo
echo "This will copy a fresh PyMACS instance into:"
echo "  $TARGET_DIR"
echo

if ! command -v git >/dev/null 2>&1; then
  echo "ERROR: git is not installed or not available in PATH."
  echo "Please install git first, then run this command again."
  exit 1
fi

if [ "$(find "$TARGET_DIR" -mindepth 1 -maxdepth 1 | wc -l)" -gt 0 ]; then
  echo "WARNING: This directory is not empty."
  echo "Files with the same names as PyMACS files may be overwritten."
  echo
  printf "Continue copying PyMACS into this directory? [y/N]: "
  read -r answer
  case "$answer" in
    y|Y|yes|YES) ;;
    *)
      echo "Install cancelled."
      exit 0
      ;;
  esac
fi

echo
echo "Cloning PyMACS into a temporary folder..."
git clone "$REPO_URL" "$TMP_DIR/pymacs"

cd "$TMP_DIR/pymacs"

if command -v git-lfs >/dev/null 2>&1 || git lfs version >/dev/null 2>&1; then
  echo "Git LFS detected. Pulling large tracked files..."
  git lfs install
  git lfs pull
else
  echo "WARNING: Git LFS was not detected."
  echo "Core PyMACS files will still be copied, but large example assets may remain as LFS pointer files."
  echo "Install Git LFS later and re-run this installer if you need the full example datasets."
fi

echo
echo "Copying PyMACS files into your current directory..."

shopt -s dotglob nullglob
for item in "$TMP_DIR/pymacs"/*; do
  base="$(basename "$item")"
  if [ "$base" = ".git" ]; then
    continue
  fi

  if [ -e "$TARGET_DIR/$base" ]; then
    printf "Overwrite existing %s? [y/N]: " "$base"
    read -r overwrite
    case "$overwrite" in
      y|Y|yes|YES)
        rm -rf "$TARGET_DIR/$base"
        ;;
      *)
        echo "Skipping $base"
        continue
        ;;
    esac
  fi

  cp -R "$item" "$TARGET_DIR/"
done

echo "Cleaning up temporary clone..."
cleanup
trap - EXIT

cd "$TARGET_DIR"

echo
echo "Done. PyMACS has been copied into:"
echo "  $TARGET_DIR"
echo
echo "Next recommended steps:"
echo "  conda env create -f environment_cgenff.yml"
echo "  conda env create -f environment_mdanalysis.yml"
echo
echo "Then follow the Step 1 / Step 2 / Step 3 workflow in this README."
EOF
```

<a id="overview"></a>

## 🚀 Overview

**PyMACS** is a Python-based toolkit for automating **molecular dynamics setup**, **simulation execution**, **trajectory analysis**, and **publication-ready figure/report generation** using **GROMACS** and **CHARMM36**.

PyMACS supports:

- **Protein-only systems**
- **Protein–ligand systems**
- **Protein–protein systems**
- **PROTAC / multi-component workflows**

PyMACS Step 1 now supports flexible structure intake for real-world medicinal chemistry workflows:

- **Local PDB inputs**
- **Local CIF / mmCIF inputs**
- **Interactive download from the RCSB PDB directly into the working directory**
- **Headless fetch via PDB ID for scripted runs**
- **Interactive polymer-chain pruning to isolate a single asymmetric unit when needed**
- **Instance-level ligand and cofactor selection when multiple copies are present**

Core goals:

- reproducible MD setup and execution
- clear logging and practical failure recovery
- scriptable, headless execution suitable for `tmux`, remote servers, and HPC-style workflows
- scalable and batch-friendly trajectory analysis across many simulations
- clean downstream visualization and PDF figurebook generation

---

<a id="repository-navigation"></a>

## 🧭 Repository Navigation

<p align="center">
  <a href="#super-quick-start-install">
    <img src="https://img.shields.io/badge/Install-Super%20Quick%20Lightweight-yellow?style=for-the-badge&logo=gnubash" alt="Super Quick Start Install">
  </a>
  <a href="#quick-start-install">
    <img src="https://img.shields.io/badge/Install-Full%20Quick%20Start-yellow?style=for-the-badge&logo=gnubash" alt="Full Quick Start Install">
  </a>
  <a href="#quick-start-templates">
    <img src="https://img.shields.io/badge/Get%20Started-Quick%20Start-orange?style=for-the-badge&logo=gnubash" alt="Get started with PyMACS">
  </a>
  <a href="#guided-tutorial">
    <img src="https://img.shields.io/badge/Tutorial-Guided%20Walkthrough-blue?style=for-the-badge&logo=readthedocs" alt="Open the guided PyMACS tutorial">
  </a>
  <a href="#example-figurebook-demo">
    <img src="https://img.shields.io/badge/Run%20the%20Demo-Example%20System-9cf?style=for-the-badge&logo=python" alt="Run the PyMACS demo">
  </a>
  <a href="#citation">
    <img src="https://img.shields.io/badge/Cite-PyMACS-lightgrey?style=for-the-badge&logo=googlescholar" alt="Cite PyMACS">
  </a>
</p>

- [Super Quick Start Install (lightweight)](#super-quick-start-install)
- [Full Quick Start Install](#quick-start-install)
- [What’s in this repo](#repository-layout)
- [System requirements](#system-requirements)
- [Environment setup](#environment-setup)
- [Force fields included](#force-fields)
- [Supported structure inputs](#supported-structure-inputs)
- [Ligand parameterization and CGenFF support](#ligand-parameterization)
- [The 3 CGenFF modes](#cgenff-modes)
- [Configurable parameters and command-line flags](#configurable-parameters-and-flags)
- [Analysis modes at a glance](#analysis-modes-at-a-glance)
- [Protein-protein / protein-peptide interface RIN analysis](#interface-rin-analysis)
- [Interface RIN flags](#interface-rin-flags)
- [Interface RIN outputs](#interface-rin-outputs)
- [Quick start templates](#quick-start-templates)
- [Pipeline stages](#pipeline-stages)
- [Restarting / resuming production MD](#restart-production-md)
- [Guided tutorial / walkthrough](#guided-tutorial)
- [Example figurebook demo](#example-figurebook-demo)
- [Standard vs MPI-compatible scripts](#standard-vs-mpi-compatible-scripts)
- [Additional workflow documentation](#additional-workflow-documentation)
- [tmux quick tips](#tmux-quick-tips)
- [Troubleshooting](#troubleshooting)
- [Best practices](#best-practices)
- [Repository description](#repository-description)
- [Provenance note](#provenance-note)
- [Citation](#citation)
- [Contact](#contact)

These links use explicit README anchors so the buttons and table of contents remain stable on GitHub.

---

<a id="repository-layout"></a>

## 📦 What’s in this repo — current directory layout

```text
PyMACs/
├── 1_AutomateGromacs.py
├── 1_AutomateGromacs_MPI.py
├── 2_AutomateGromacs.py
├── 2_AutomateGromacs_MPI.py
├── 3A_AutomateGromacs.py
├── 3A_AutomateGromacs_MPI.py
├── 3B_NETWORX.py
├── 3C_Interface_RIN.py
├── 3_PROTAC_Analysis.py
├── 3_PROTAC_Analysis_MPI.py
├── 4PDF4MD.py
├── 4_MDfigs.txt
├── docs/
│   ├── README.md
│   ├── OUTPUT_GALLERY.md
│   ├── SOFTWARE_AND_CITATIONS.md
│   ├── USE_CASE_WALKTHROUGH.md
│   └── equilibration_plan_example.json
├── cgenff_charmm2gmx_py3_nx2.py
├── charmm36.ff/
├── charmm36_ljpme-jul2022.ff/
├── Example_Choices/
│   ├── Example1/
│   │   ├── input/
│   │   ├── parameters/
│   │   ├── prepared_system/
│   │   └── completed_run/
│   ├── Example2/
│   ├── Example3/
│   └── Example4/
├── em.mdp
├── MDPs/
│   ├── em.mdp
│   ├── ions.mdp
│   ├── md.mdp
│   ├── npt.mdp
│   └── nvt.mdp
├── environment_cgenff.yml
├── environment_mdanalysis.yml
├── recreate_envs.sh
└── README.md
```

The repository may include additional helper scripts or files as the project evolves.

---

<a id="system-requirements"></a>

## ⚙️ System requirements

- **OS:** Linux or WSL2 recommended  
  - macOS can be used for analysis-only workflows, but GROMACS installation and GPU support vary.
- **GROMACS:** 2022+ recommended  
  - Standard scripts expect a working `gmx` command in the active environment.
  - MPI-compatible scripts support `gmx`, `gmx_mpi`, `gmx-mpi`, or an explicit executable path.
  - MPI-compatible scripts accept `--gmx-bin`.
  - `PYMACS_GMX_BIN` can be used as an environment-variable override.
  - Auto-detection order in the MPI-compatible scripts is `gmx_mpi`, then `gmx-mpi`, then `gmx`.
- **Python:** managed through conda environments.
- **Conda / Mamba:** recommended for reproducibility.
- **tmux:** optional but strongly recommended for long simulations.
- **GPU:** optional but recommended for production MD. Step 2 can auto-select an available NVIDIA GPU and fall back to CPU when needed.

---

<a id="standard-vs-mpi-compatible-scripts"></a>

## 🔀 Standard vs MPI-compatible scripts

PyMACS now ships both the original workstation-oriented scripts and MPI-compatible variants for systems where the GROMACS executable is not exposed as plain `gmx`.

| Standard script | MPI-compatible script | Purpose |
|---|---|---|
| `1_AutomateGromacs.py` | `1_AutomateGromacs_MPI.py` | setup, `pdb2gmx`, box, solvation, ions, EM prep |
| `2_AutomateGromacs.py` | `2_AutomateGromacs_MPI.py` | EM, NVT, NPT, production MD, checkpoint resume |
| `3A_AutomateGromacs.py` | `3A_AutomateGromacs_MPI.py` | analysis, `trjconv` preprocessing, pocket/contact/RMSD/RMSF workflows |
| `3C_Interface_RIN.py` | same script | standalone protein-protein / protein-peptide interface RIN analysis; called by 3A for ligandless multi-chain systems |
| `3_PROTAC_Analysis.py` | `3_PROTAC_Analysis_MPI.py` | PROTAC-specific trajectory analysis and optional GROMACS centering |

Use the MPI-compatible scripts when:

- you are on Triton or another HPC system where GROMACS is installed as `gmx_mpi` or `gmx-mpi`
- you are running inside SLURM / MPI environments
- your environment does not provide a `gmx` alias
- you want explicit binary control through `--gmx-bin` or `PYMACS_GMX_BIN`

The standard scripts remain appropriate for local workstations where `gmx` is already the expected executable name.

`3C_Interface_RIN.py` does not require a separate MPI variant because it uses MDAnalysis-based trajectory analysis and does not call GROMACS directly.

---

<a id="analysis-modes-at-a-glance"></a>

## 🧭 Analysis modes at a glance

| System type | User mode / detection | Main script path | What PyMACS analyzes |
|---|---|---|---|
| Small-molecule ligand bound to protein | `--mode ligand` or explicit `--ligand` | `3A_AutomateGromacs.py` plus `3B_NETWORX.py` when available | ligand RMSD/RMSF, pocket contacts, ligand-residue contact maps, ligand interaction network |
| Single-chain apo protein | `--mode protein` with one polymer chain | `3A_AutomateGromacs.py` | protein RMSD/RMSF, radius of gyration, secondary-structure evolution |
| Ligandless protein-protein or protein-peptide assembly | `--mode peptide`, `--mode protein`, or multi-chain ligandless detection | `3A_AutomateGromacs.py` dispatches to `3C_Interface_RIN.py` | chain-pair contact persistence, residue-residue interface maps, interface RIN network |
| PROTAC / ternary complex | `--mode protac` | `3A_AutomateGromacs.py` dispatches to `3_PROTAC_Analysis.py` or MPI equivalent | PROTAC component RMSD/RMSF, ternary contact analysis, component-level interaction plots |
| Biological polymer system | `--mode biological` | `3A_AutomateGromacs.py` | protein/nucleic-acid-aware structural analysis and chain-level summaries where supported |

The key distinction is whether the system has a true small-molecule ligand. Ligand workflows are pocket-centered. Ligandless multi-chain workflows are interface-centered.

---

<a id="configurable-parameters-and-flags"></a>

## 🛠️ Configurable parameters and command-line flags

PyMACS can be run interactively, where the scripts ask questions one at a time, or in command-line mode, where users provide settings directly using flags.

A flag is an extra instruction added after a command. For example:

```bash
python 1_AutomateGromacs.py --box-type dodecahedron --box-distance 1.0
```

In that command, `--box-type` and `--box-distance` are flags.

For beginners, the safest choice is usually to run the scripts interactively first and press Enter to accept defaults. After that works, users can modify one setting at a time.

### 1. Structure input and setup flags

| Script / stage | Flag | What it controls | Common options or values | Beginner recommendation | What happens if you increase it / choose a larger value | What happens if you decrease it / choose a smaller value | Notes |
|---|---|---|---|---|---|---|---|
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--pdb` | Which local structure file to start from. | `.pdb`, `.cif`, `.mmcif` file name | Use this when you already have a local structure file. | Not a larger/smaller setting. Different file choices change the whole starting system. | Not a larger/smaller setting. Choosing the wrong file usually leads to wrong chains, missing ligands, or setup errors. | Safe and common. Good first explicit flag. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--fetch-pdb` | Downloads a structure directly from the PDB before setup. | PDB ID such as `9UWJ` | Very useful if you do not already have a local file. | Not numeric. Different PDB IDs give different systems. | Not numeric. If omitted, PyMACS looks for local files instead. | Fetches into the current working directory. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--fetch-format` | File format used when fetching from the PDB. | `pdb`, `cif`, `mmcif` | Leave the default unless you have a reason to prefer another format. | `cif` or `mmcif` often preserves modern structural details better than old-style PDB files. | `pdb` can be simpler to inspect manually, but sometimes carries less structural detail. | Verified default: `cif`. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--chain-names` | Simple custom names for chains in chain order. | Comma-separated names such as `Rap1B,Rap1GAP` | Leave this alone unless you want clearer labels in downstream outputs. | Not numeric. More descriptive names make reports easier to read. | Not numeric. Leaving it out keeps automatic or existing labeling. | Helpful for publication-quality naming. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--chain-map` | Explicit chain ID to name mapping. | `A:Name,B:Name` | Use this only when chain labels matter and you know the input chain IDs. | Not numeric. More explicit mapping gives you tighter control. | Not numeric. If omitted, PyMACS uses interactive/default naming behavior. | Safer than `--chain-names` when the chain order is ambiguous. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--keep-chains` | Keeps only selected polymer chains from the input structure. | Chain letters or numeric positions | Leave it alone on your first run unless you intentionally want one chain or one asymmetric unit. | Not numeric. Keeping more chains makes the system larger and more complex. | Not numeric. Keeping fewer chains simplifies the system but may remove biologically important partners. | Very useful for large deposited structures. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--drop-chains` | Removes selected polymer chains from the input structure. | Chain letters or numeric positions | Use with care. | Not numeric. Dropping more chains gives a smaller, faster system. | Not numeric. Dropping too much can remove important biology. | Easier to make mistakes here than with `--keep-chains`. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--remove-input-waters`, `--remove-input-ions`, `--remove-input-solvents`, `--remove-input-resnames` | Removes unwanted components before setup. | Boolean flags or comma-separated residue names | Keep the defaults on your first run unless you know those components should be removed. | Removing more components can make setup cleaner and simpler. | Removing less preserves the deposited context, but may leave crystallographic clutter. | Good troubleshooting tools when deposited files are messy. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--show-generated-pdbs`, `--allow-generated-inputs` | Shows generated intermediate PDB files during interactive picking. | Include the flag or omit it | Most beginners should leave this off. | Including it gives more choices, which can help advanced troubleshooting. | Omitting it keeps the file picker simpler. | Mainly a debugging convenience. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--strict-pdb-validation` | Turns selected structure warnings into hard-stop errors. | Include the flag or omit it | Useful when you want conservative setup behavior. | Including it catches possible structure problems earlier. | Omitting it is more forgiving and may let setup continue. | Safer scientifically, but more likely to stop early. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--nucleic-acid-termini`, `--protein-nterm`, `--protein-cterm`, `--chain-types` | Controls terminal handling and polymer type interpretation during setup. | Verified options include `auto`, `5ter-3ter`, `none`; `NH2`, `NH3+`, `None`; `COO-`, `COOH`, `None` | Beginners should usually leave these at default values. | More manual control can help special biomolecule cases. | Less manual control means safer defaults and fewer chances to mismatch topology choices. | These settings matter most for unusual chain chemistry and nucleic-acid systems. |

Quick example:

```bash
python 1_AutomateGromacs.py --fetch-pdb 9UWJ --fetch-format cif
```

### 2. Box and solvent-padding flags

| Script / stage | Flag | What it controls | Common options or values | Beginner recommendation | What happens if you increase it / choose a larger value | What happens if you decrease it / choose a smaller value | Notes |
|---|---|---|---|---|---|---|---|
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--box-type` | The shape of the simulation box around your molecule. | `cubic`, `dodecahedron`, `octahedron` | `cubic` is easiest to understand. `dodecahedron` is a very common compact choice once you are comfortable. | Choosing a compact shape like `dodecahedron` or `octahedron` often uses less water and can run faster. | Choosing `cubic` often uses more water and can run a bit slower, but it is conceptually simple. | `cubic`: easy to understand. `dodecahedron`: compact and often faster. `octahedron`: another compact option. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--box-distance` | The padding distance between the solute and the box wall, in nanometers. | Verified default: `1.0` | Leave the default on a first run. Increase it only for a clear reason. | A larger box distance adds more water padding. This is usually safer physically, but slower and larger on disk. | A smaller box distance speeds things up, but risks periodic-image artifacts if the molecule gets too close to its own copy. | One of the most useful and safest tuning flags. |

Quick example:

```bash
python 1_AutomateGromacs.py --box-type dodecahedron --box-distance 1.2
```

### 3. Ligand, cofactor, and system-mode flags

| Script / stage | Flag | What it controls | Common options or values | Beginner recommendation | What happens if you increase it / choose a larger value | What happens if you decrease it / choose a smaller value | Notes |
|---|---|---|---|---|---|---|---|
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py`, `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py`, `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py`, `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--ligand` | The residue name for the small molecule or PROTAC-like component of interest. | 3-letter or short residue code such as `A1D`, `A1E`, `PTC` | Use this when you know the ligand code and want to avoid interactive detection. | Not numeric. The main effect is choosing a different target molecule for setup or analysis. | Not numeric. Leaving it blank may trigger interactive or automatic detection. | One of the most important flags for ligand-centered workflows. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py`, `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py`, `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--cofactors` | Additional retained non-protein components kept alongside the main ligand. | Comma-separated names such as `CLR,HEM` | Leave empty for your first simple ligand run. Use it when your system has a cofactor, cholesterol, heme, or another context molecule that should stay. | More listed cofactors preserve more biological context, but make setup and analysis more complex. | Fewer listed cofactors simplify the system, but may remove important context. | Good for receptor systems with important bound helpers. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--auto-keep-hetero`, `--component-mode`, `--component-clash-check`, `--no-component-clash-check`, `--drop-clashing-cofactors` | Fine control over how non-protein components are kept, checked, or discarded during setup. | Verified `--component-mode` options: `selected`, `all`, `none` | Beginners should usually leave the defaults alone. | `all` keeps more detected non-protein pieces and can preserve more context. | `none` or aggressive dropping makes the system simpler, but may remove something important. | Useful when deposited structures contain multiple ligand-like molecules. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py`, `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--mode` | The kind of system you are running or analyzing. | Verified options: `ligand`, `protein`, `peptide`, `protac`, `biological` | If you are unsure, use interactive mode first. On later runs, set the mode explicitly so the script knows what kind of workflow you want. | Not numeric. Choosing more complex modes such as `protac` or `biological` activates more specialized behavior. | Not numeric. Choosing a simpler mode may skip important logic for complex systems. | Very important for getting the right workflow branch. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--allow-covalent-ligand` | Allows ligand workflows to continue even if the input suggests a covalent attachment. | Include the flag or omit it | Beginners should usually leave this off unless they know they are working with a covalent system and understand the consequences. | Including it is more permissive and may allow unusual systems to continue. | Omitting it is safer and more conservative. | High-risk override; use carefully. |
| `1_AutomateGromacs.py`, `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py`, `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py`, `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--headless` or explicit mode-selection flags | Skips interactive prompts and uses only command-line values. | Include the flag or omit it | Beginners should run interactively first, then use headless mode after the workflow is understood. | Including it is faster and reproducible for repeated runs. | Omitting it gives guided prompts and is usually easier for first runs. | Excellent once you know the correct settings. |

Quick example:

```bash
python 2_AutomateGromacs.py --mode ligand --ligand A1D --cofactors CLR
```

### 4. MDP and equilibration customization flags

| Script / stage | Flag | What it controls | Common options or values | Beginner recommendation | What happens if you increase it / choose a larger value | What happens if you decrease it / choose a smaller value | Notes |
|---|---|---|---|---|---|---|---|
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py`, `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--mdp-dir` | Where PyMACS looks for MDP parameter templates. | Verified default: `MDPs` | Leave the default unless your templates live somewhere else. | Not numeric. A different directory lets you use a different parameter set. | Not numeric. A wrong directory leads to missing-file errors. | MDP files control simulation settings such as time step and thermostat behavior. |
| `1_AutomateGromacs.py`, `1_AutomateGromacs_MPI.py` | `--ions-mdp`, `--em-mdp` | Which setup-stage MDP files are used for ion preparation and energy minimization. | Verified defaults: `ions.mdp`, `em.mdp` | Leave these alone on a first run. | Not numeric. Different files can make setup stricter or looser depending on their contents. | Not numeric. Wrong files can break setup or produce unexpected behavior. | Good for advanced method tuning. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--em-mdp`, `--nvt-mdp`, `--npt-mdp`, `--md-mdp` | Which MDP files are used for each MD stage. | Verified defaults: `em.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp` | Leave these at default until the standard workflow is working. | Not numeric. More customized MDP files can improve a specific scientific workflow, but also add troubleshooting complexity. | Not numeric. Simpler defaults are easier to troubleshoot. | A powerful advanced feature. Change one template at a time. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--equilibration-plan` | Replaces the standard NVT then NPT sequence with a custom JSON equilibration plan. | JSON file path | Beginners should leave this alone. | Not a larger/smaller setting. A more detailed plan gives more protocol control. | Not a larger/smaller setting. Leaving it out keeps the standard default sequence. | Best for advanced users who want a custom multi-stage protocol. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--index-file`, `--polymer-group` | Controls which atom groups or custom index file are used in MD setup. | Custom `index.ndx`; verified `--polymer-group` options include `auto`, `Protein`, `Protein_RNA`, `Protein_DNA`, `Protein_Nucleic`, `non-Water` | Beginners should usually let PyMACS generate and choose these automatically. | Not numeric. More specific grouping gives tighter control for unusual systems. | Not numeric. Automatic grouping is safer and easier. | These settings mainly help biomolecule systems or advanced coupling-group troubleshooting. |

Quick example:

```bash
python 2_AutomateGromacs.py --mdp-dir MDPs --md-mdp md.mdp
```

### 5. CPU, GPU, and MPI resource flags

Step 2 is designed to use the best available hardware automatically. In the current GPU-adaptive workflow, `--compute auto` is the intended default behavior: PyMACS checks for NVIDIA GPUs, selects GPU 0 automatically when only one GPU is available, and falls back to CPU-only execution if no usable GPU is detected. Use `--no-gpu` or `--compute CPU` when you intentionally want CPU-only behavior.

| Script / stage | Flag | What it controls | Common options or values | Beginner recommendation | What happens if you increase it / choose a larger value | What happens if you decrease it / choose a smaller value | Notes |
|---|---|---|---|---|---|---|---|
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--compute` | Whether MD uses automatic hardware selection, CPU-only execution, or forced GPU mode. | `auto`, `CPU`, `GPU` | Leave this at `auto` for normal use. | `GPU` forces a GPU attempt when a GPU is available. | `CPU` disables GPU acceleration and is slower, but simpler for troubleshooting. | `auto` chooses GPU when possible and CPU when needed. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--gpu`, `--gpu-id`, `--gpu-ids`, `--no-gpu` | Which GPU to use, or whether to disable GPU use entirely. | Single GPU ID like `0`, exact ID string like `01`, or `--no-gpu` | Leave this alone unless you have multiple GPUs or a GPU problem. | More explicit GPU selection helps on multi-GPU workstations and cluster nodes. | `--no-gpu` forces CPU-only mode and avoids GPU troubleshooting. | `--gpu` and `--gpu-id` are aliases in Step 2. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py`, `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--threads`, `--ntomp` | CPU thread count used for MD or analysis work. | Positive integer | Beginners should usually let PyMACS auto-detect or use a modest thread count. | More threads can speed things up until the computer is saturated. Too many can make performance worse. | Fewer threads may be slower, but usually safer for laptops and shared systems. | `--threads` in Step 2 is a legacy alias for `--ntomp`. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--ntmpi` | Number of GROMACS thread-MPI ranks for local MD commands. | Default usually `1` | Most laptop and workstation users should leave this alone. | More ranks can help on the right hardware and build, but can also complicate performance tuning. | Keeping it at `1` is simpler and often best locally. | Beginners should usually avoid MPI tuning unless instructed. |
| `2_AutomateGromacs_MPI.py` | `--mpi-ranks`, `--mpi-launcher`, `--external-mpi` | Controls external MPI launching in the MPI-compatible Step 2 script. | Rank counts; launcher strings such as `mpirun -np 4` or `srun -n 4` | Most beginners should skip these entirely. | More ranks may improve throughput on HPC systems when configured correctly. | Fewer ranks simplify debugging and avoid overcomplicating local runs. | Advanced/HPC-only area. |
| `1_AutomateGromacs_MPI.py`, `2_AutomateGromacs_MPI.py`, `3A_AutomateGromacs_MPI.py`, `3_PROTAC_Analysis_MPI.py` | `--gmx-bin` | Selects the GROMACS executable explicitly. | `auto`, `gmx`, `gmx_mpi`, `gmx-mpi`, or full path | Leave `auto` unless your environment needs a specific binary. | Not numeric. A more explicit choice gives more control on clusters or unusual installs. | Not numeric. `auto` is simpler and often works. | Related environment variable: `PYMACS_GMX_BIN`. |
| MPI-compatible scripts | `PYMACS_GMX_BIN` | Environment-variable override for the GROMACS executable. | Shell variable such as `export PYMACS_GMX_BIN=gmx_mpi` | Only use this if your environment consistently needs a non-default binary. | Not numeric. A different binary choice can make MPI-compatible scripts work on HPC systems. | Not numeric. Leaving it unset preserves auto-detection. | Useful when you want the same binary choice across many commands. |

#### Step 2 automatic GPU behavior

When `--compute auto` is active, PyMACS tries to maximize available hardware without requiring beginners to remember GPU flags:

- **one NVIDIA GPU detected:** PyMACS selects GPU 0 automatically
- **multiple NVIDIA GPUs detected:** interactive mode asks which GPU to use; headless mode defaults to GPU 0 unless `--gpu` or `--gpu-ids` is supplied
- **no NVIDIA GPU detected:** PyMACS runs CPU-only
- **GPU memory is low:** PyMACS warns, still tries GPU, and falls back if needed
- **energy minimization:** PyMACS avoids PME-on-GPU because EM uses a non-dynamical minimizer such as steepest descent
- **NVT, NPT, and production:** PyMACS tries GPU-accelerated profiles first, then safer mixed CPU/GPU profiles, then CPU fallback if needed

This means a normal run can usually be started as:

```bash
python 2_AutomateGromacs.py --ns 0.25
```

For an explicit CPU-only run:

```bash
python 2_AutomateGromacs.py --compute CPU --ns 0.25
```

For an explicit GPU run on GPU 0:

```bash
python 2_AutomateGromacs.py --compute GPU --gpu-id 0 --ns 0.25
```

### 6. Restart and checkpoint flags

| Script / stage | Flag | What it controls | Common options or values | Beginner recommendation | What happens if you increase it / choose a larger value | What happens if you decrease it / choose a smaller value | Notes |
|---|---|---|---|---|---|---|---|
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--ns` | Production simulation length in nanoseconds. | Floating-point values such as `0.25`, `1`, `25`, `50` | Start with `0.25` or `1` for a test run. | Higher values give more sampling and can answer harder scientific questions, but take longer and create larger files. | Lower values are faster and good for testing, but may be too short for real scientific conclusions. | Production-quality run length depends on the scientific question. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--production_only` | Skips EM/NVT/NPT and runs only production MD. | Include the flag or omit it | Do not use this on your very first run. | Including it is faster when your system is already prepared and equilibrated. | Omitting it runs the full standard workflow. | Best for restarts or repeated production extensions. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--resume` | Resumes from an existing checkpoint if one is present. | Include the flag or omit it | Very useful after an interrupted run. | Including it helps continue long jobs instead of starting over. | Omitting it avoids automatic resume behavior. | Best paired with `--production_only` for interrupted long runs. |
| `2_AutomateGromacs.py`, `2_AutomateGromacs_MPI.py` | `--force_restart` | Prevents resume behavior even if checkpoints exist. | Include the flag or omit it | Use only when you intentionally want a fresh restart path. | Including it forces a clean rerun path. | Omitting it allows normal resume logic when available. | Helpful when old checkpoint state is suspected to be bad. |

Quick example:

```bash
python 2_AutomateGromacs.py --resume --production_only --ns 50
```

### 7. Analysis cutoff and interaction flags

| Script / stage | Flag | What it controls | Common options or values | Beginner recommendation | What happens if you increase it / choose a larger value | What happens if you decrease it / choose a smaller value | Notes |
|---|---|---|---|---|---|---|---|
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--topo`, `--gro`, `--traj`, `--pocket_pdb`, `--pocket_xtc`, `--outdir` | Which trajectory and structure files are analyzed, and where outputs go. | File paths | Leave defaults alone unless your files are named differently. | Not numeric. Different file choices change what gets analyzed. | Not numeric. Wrong file paths cause missing-file errors. | `--gro` is an alias for `--topo`. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--contact_cutoff` | Distance cutoff for counting ligand contacts. | Verified default: `4.0 Å` | Leave the default on a first run. | A larger contact cutoff counts more interactions, including weaker or more distant ones. | A smaller cutoff is stricter, but may miss flexible or transient contacts. | Good for tuning how broad contact detection feels. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--pocket-cutoff` | Radius used to define the binding pocket around the ligand. | Verified default: `5.0 Å` | Leave the default on a first run. | A larger pocket cutoff includes more residues, which can capture broader environment effects but also adds noise. | A smaller pocket cutoff focuses on the nearest residues but may miss relevant interactions. | One of the most useful analysis-tuning flags. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--hbond-distance-cutoff` | Distance threshold for hydrogen-bond-like calls. | Verified default: `3.5 Å` | Leave the default unless you have a strong reason to tune it. | A larger value is more permissive and may count weaker contacts. | A smaller value is stricter and may miss flexible hydrogen bonds. | Helpful when comparing strict versus broad interaction definitions. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--hbond-angle-cutoff` | Angular strictness for hydrogen-bond-like geometry. | Verified default: `135°` | Leave the default on a first run. | A higher angle cutoff is stricter and favors more linear hydrogen bonds. | A lower angle cutoff is more permissive and counts less ideal geometry. | Higher usually means fewer but cleaner hydrogen-bond calls. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--min_contact_frac`, `--minfrac` | Minimum fraction of frames a residue must contact the ligand to appear in outputs. | Verified default for `--min_contact_frac`: `0.10` | Leave the default on a first run. | A larger fraction shows only the most persistent contacts, giving cleaner but narrower plots. | A smaller fraction includes more transient contacts, which can be informative but noisier. | `--minfrac` is mainly a NETWORX-style compatibility flag. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--compound-name`, `--numbering-template` | Makes outputs easier to read by changing display labels and residue numbering source. | Display string or PDB file path | Optional. Use only when you want prettier labels or reference numbering. | Not numeric. More customized naming improves presentation. | Not numeric. Leaving defaults keeps raw internal naming. | Helpful for manuscripts and reviewer-friendly figures. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--lig-atom-start`, `--lig-atom-end` | Manual ligand atom-range fallback for difficult detection cases. | Integer atom indices | Beginners should avoid this unless ligand auto-detection fails and they understand atom numbering. | Larger ranges include more atoms and risk selecting the wrong region. | Smaller ranges may miss part of the ligand. | High-risk troubleshooting feature. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--net-out`, `--ellipse-rx`, `--ellipse-ry`, `--no-edges` | Fine control over NETWORX-style network figure appearance. | Output name, ellipse radii, or boolean | Most beginners should leave these alone. | Larger ellipse radii make figure layouts broader. | Smaller radii make layouts tighter and can crowd labels. | Presentation-focused rather than science-focused. |

Quick example:

```bash
python 3A_AutomateGromacs.py --contact_cutoff 3.5 --pocket-cutoff 6.0
```

<a id="interface-rin-flags"></a>

### 8. Protein-protein / protein-peptide interface RIN flags

| Script / stage | Flag | What it controls | Common options or values | Beginner recommendation | What happens if you increase it / choose a larger value | What happens if you decrease it / choose a smaller value | Notes |
|---|---|---|---|---|---|---|---|
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py` | `--skip-interface-rin` | Prevents 3A from launching 3C even when a ligandless multi-chain interface is detected. | Include the flag or omit it | Leave this off if you want automatic interface analysis. | Not numeric. Including it skips interface RIN outputs. | Omitting it allows automatic dispatch for multi-chain ligandless systems. | Useful if you only want RMSD/RMSF/SSE outputs. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py`, `3C_Interface_RIN.py` | `--interface-contact-cutoff` / `--contact-cutoff` | Distance cutoff in Å for residue-residue chain-interface contacts. | Default usually `4.0 Å` | Leave default first. | Larger cutoffs detect broader, weaker, or more distant contacts. | Smaller cutoffs focus on tighter contacts but may miss flexible interface events. | 3A passes this to 3C as `--contact-cutoff`. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py`, `3C_Interface_RIN.py` | `--interface-min-contact-frac` / `--min-contact-frac` | Minimum fraction of analyzed frames required for residue-pair edges to appear in the network. | Default `0.10` | Good first value for readable networks. | Higher values show only persistent/core interface contacts. | Lower values include transient contacts but can make networks crowded. | 3A passes this to 3C as `--min-contact-frac`. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py`, `3C_Interface_RIN.py` | `--interface-frame-step` / `--frame-step` | Frame stride for interface analysis. | `1`, `5`, `10` | Use `1` for complete analysis; use `5` or `10` for faster screening. | Larger values analyze fewer frames and run faster. | Smaller values analyze more frames and provide better temporal resolution. | 3A passes this to 3C as `--frame-step`. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py`, `3C_Interface_RIN.py` | `--interface-max-edges` / `--max-edges` | Maximum residue-residue edges shown in the final network figure. | Default `100` | Leave default unless the figure is too crowded. | Larger values show more interface detail but can clutter the graph. | Smaller values make cleaner summary figures but may hide weaker contacts. | CSV outputs still preserve full contact tables. |
| `3A_AutomateGromacs.py`, `3A_AutomateGromacs_MPI.py`, `3C_Interface_RIN.py` | `--interface-chain-pairs` / `--chain-pairs` | Restricts analysis to specific chain pairs. | `A:B,A:C` | Omit this to analyze all chain pairs. | Not numeric. More listed pairs gives broader interface coverage. | Fewer listed pairs focuses analysis on selected interfaces. | Useful for multi-chain assemblies. |

Quick examples:

```bash
python 3A_AutomateGromacs.py --mode protein --headless --interface-frame-step 10
python 3A_AutomateGromacs.py --mode protein --headless --interface-chain-pairs A:B
python 3A_AutomateGromacs.py --mode protein --headless --skip-interface-rin
```

<a id="interface-rin-outputs"></a>

### Interface RIN outputs

All interface-RIN outputs are written to:

```text
Analysis_Results/Interface_RIN/
```

CSV outputs:

- `interface_chain_pair_summary.csv`
- `interface_contact_timeseries.csv`
- `interface_residue_pair_contacts.csv`
- `interface_node_summary.csv`
- `interface_runtime_summary.txt`

PNG outputs:

- `Interface_ChainPair_ContactFractions.png`
- `Interface_ContactCount_Timeseries.png`
- `Interface_MinDistance_Timeseries.png`
- `Interface_ResiduePair_ContactHeatmap.png`
- `Interface_RIN_Network.png`
- `Interface_ChainLevel_Network.png`
- `Interface_TopResiduePairs_Barplot.png`

Quick interpretation guide:

- Chain-pair summary tells whether two chains remain associated over time.
- Contact time series shows how many residue-residue contacts exist at each frame.
- Minimum-distance time series shows whether chains separate or remain close.
- Residue-pair heatmap shows which residue contacts persist or disappear over time.
- Residue-level RIN network shows persistent interface contacts as edges.
- Chain-level network summarizes which chains interact most persistently.

<a id="protac-specific-analysis-flags"></a>

### 9. PROTAC-specific analysis flags

| Script / stage | Flag | What it controls | Common options or values | Beginner recommendation | What happens if you increase it / choose a larger value | What happens if you decrease it / choose a smaller value | Notes |
|---|---|---|---|---|---|---|---|
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--topo`, `--traj`, `--pdb-reference`, `--atomindex`, `--outdir` | Input files and output folder for PROTAC analysis. | File paths | Leave defaults alone if your run directory uses the standard file names. | Not numeric. Different file choices change what gets analyzed. | Not numeric. Wrong paths cause immediate file errors. | `atomIndex.txt` is especially important in the PROTAC workflow. |
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--distance-cutoff` | Main protein-ligand contact cutoff in the ternary-complex analysis. | Verified default: `4.0 Å` | Leave the default on a first run. | A larger value counts more contacts and may make the ternary interface look broader. | A smaller value is stricter and may miss weak or transient contacts. | In PROTAC work, this affects contacts to both protein partners. |
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--start-frame`, `--end-frame`, `--frame-step`, `--quick-test`, `--dry-run` | How much of the trajectory is scanned, and whether the script does a smoke test or validation-only run. | Frame numbers; verified quick-test behavior defaults to `start=0`, `end=500`, `step=10` unless overridden | `--dry-run` and `--quick-test` are very beginner-friendly. | Larger frame windows and smaller frame steps analyze more data but take longer. | Smaller frame windows and larger frame steps are faster but less complete. | Excellent for checking a new PROTAC setup before a full analysis. |
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--skip-pca`, `--skip-network`, `--skip-contact-map`, `--skip-rmsf`, `--contacts-only`, `--basic-only` | Turns off expensive or advanced parts of the PROTAC analysis. | Boolean flags | Beginners can safely use `--basic-only` or `--quick-test` when they just want to confirm the workflow runs. | More outputs give a richer picture of the ternary complex, but take longer. | Fewer outputs are faster and easier to inspect. | Very useful for staged troubleshooting. |
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--residue-numbering`, `--network-min-contact-fraction`, `--network-max-residues`, `--contact-map-min-contact-fraction`, `--contact-map-max-residues` | Controls labeling style and how crowded the network and contact map become. | Verified numbering options: `source`, `current` | Leave defaults alone on a first run. | Larger max-residue values or lower fractions show more residues, which can reveal more detail but add clutter. | Smaller max-residue values or higher fractions make figures cleaner but may hide weaker contributors. | Useful when preparing figures for a paper or presentation. |
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--skip-preprocess`, `--preprocess-only`, `--write-subset-trajectories`, `--no-write-subset-trajectories`, `--write-final-trajectory`, `--no-write-final-trajectory`, `--use-final-trajectory`, `--no-use-final-trajectory`, `--overwrite-final-trajectory`, `--overwrite-structures`, `--final-trajectory-name` | Controls how the script writes centered subset structures and whether those files are reused downstream. | Boolean flags and output base name | Beginners should usually accept defaults and avoid overwrite-related flags. | Writing more helper files makes downstream inspection easier, but uses more disk space. | Writing fewer helper files saves space, but gives you less to inspect when something goes wrong. | The centered `Final_Trajectory` compatibility files are especially useful for later analysis and figurebook building. |
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--pocket-cutoff`, `--pocket-mode`, `--pocket-min-fraction` | Controls how the PROTAC binding pocket is defined. | Verified `--pocket-mode` options: `first_frame`, `any_sampled_frame`, `persistent` | Keep defaults first. | Larger pocket cutoffs and more permissive modes include more residues around the ternary interface. | Smaller cutoffs and stricter persistence focus on the most consistent pocket residues. | PROTAC pockets are often more complex than simple ligand pockets because two protein partners may contribute contacts. |
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--interaction-mode`, `--plot-allprotein` | Controls how PROTAC contacts are classified and whether full-protein aggregate traces are shown. | Verified interaction modes: `residue_chemistry`, `geometry` | Leave defaults unless you have a clear analysis reason. | `geometry`-style choices may be more granular or stricter depending on the dataset. | Default chemistry-style summaries are usually easier to interpret first. | Good for advanced comparison work. |
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--rmsf-contact-threshold`, `--rmsf-label-top-contacts`, `--rmsf-contact-marker-mode` | Controls RMSF contact overlays and labeling. | Verified marker modes: `vlines`, `scatter`, `both` | Leave defaults unless figures are too crowded. | Larger thresholds or fewer labels make plots cleaner. Larger label counts add detail but can clutter. | Lower thresholds or too many labels can create noisy figures. | Presentation-oriented tuning. |
| `3_PROTAC_Analysis.py`, `3_PROTAC_Analysis_MPI.py` | `--analyze-water-bridges`, `--water-bridge-cutoff`, `--water-pocket-cutoff`, `--write-water-pocket-trajectory`, `--water-names` | Optional water-mediated contact analysis. | Verified defaults include `3.5 Å` and `6.0 Å` | Beginners should skip water-bridge analysis until the main workflow works. | Larger cutoffs count more possible water-mediated events. | Smaller cutoffs are stricter and may miss weak or flexible bridges. | PROTAC systems often need more careful checking than simple protein-ligand systems because the linker, ligand, ligase, and target can all contribute to the contact picture. |

<a id="figurebook-reporting-flags"></a>

### 10. Figurebook/reporting flags

`4PDF4MD.py` does not expose verified command-line flags in the current script. In this repository state, the figurebook step is controlled interactively and by input files rather than by `argparse` flags.

Useful reporting controls that were verified in the script:

- report type prompt: ligand-bound, protein-centric, or biological-system report
- ligand auto-detection from local `.cgenff.mol2` files
- output PDF name prompt
- `4_MDfigs.txt`: standard figure order and inclusion list
- `4_GraphNotes.txt`: optional notes/captions helper file
- `Analysis_Results/PROTAC/QC/protac_figure_manifest.csv`: when present, the script switches into PROTAC manifest mode automatically

If a plot is missing from the PDF:

- first confirm the plot was generated in `Analysis_Results/`
- then check whether `4_MDfigs.txt` points to the expected figure
- if you are in PROTAC mode, check the manifest file under `Analysis_Results/PROTAC/QC/`

### Most useful beginner flags

| Goal | Use this flag | Example |
|---|---|---|
| Run a short test | `--ns` | `python 2_AutomateGromacs.py --ns 0.25` |
| Use a compact water box | `--box-type` | `python 1_AutomateGromacs.py --box-type dodecahedron` |
| Add more water padding | `--box-distance` | `python 1_AutomateGromacs.py --box-distance 1.2` |
| Force CPU-only Step 2 | `--compute CPU` or `--no-gpu` | `python 2_AutomateGromacs.py --compute CPU --ns 0.25` |
| Use GPU automatically | default `--compute auto` | `python 2_AutomateGromacs.py --ns 0.25` |
| Resume an interrupted run | `--resume --production_only` | `python 2_AutomateGromacs.py --resume --production_only --ns 50` |
| Make analysis stricter | lower contact or pocket cutoffs | `python 3A_AutomateGromacs.py --contact_cutoff 3.5` |
| Make analysis broader | higher contact or pocket cutoffs | `python 3A_AutomateGromacs.py --pocket-cutoff 6.0` |
| Run interface RIN automatically | `--mode protein` | `python 3A_AutomateGromacs.py --mode protein` |
| Skip interface RIN | `--skip-interface-rin` | `python 3A_AutomateGromacs.py --mode protein --skip-interface-rin` |
| Analyze fewer frames faster | `--interface-frame-step` | `python 3A_AutomateGromacs.py --mode protein --interface-frame-step 10` |
| Focus on one chain pair | `--interface-chain-pairs` | `python 3A_AutomateGromacs.py --mode protein --interface-chain-pairs A:B` |

### Beginner safety rule

If this is your first PyMACS run, do not change many settings at once.

Recommended first test:

```bash
python 2_AutomateGromacs.py --ns 0.25
```

On a GPU workstation, this now attempts GPU acceleration automatically. Add `--no-gpu` only when you intentionally want CPU-only testing.

Once the short run works, repeat the workflow with a longer production length and only one or two changes at a time.

---

<a id="supported-structure-inputs"></a>

## 🧾 Supported structure inputs

Step 1 can now begin from several different structure sources:

- local `.pdb` files already present in the working directory
- local `.cif` files from the PDB or other structural workflows
- local `.mmcif` files
- structures fetched directly from the **RCSB PDB** during the interactive Step 1 structure picker
- structures fetched in headless mode with `--fetch-pdb`

Important Step 1 intake behavior:

- `.cif` and `.mmcif` inputs are converted automatically into a setup-ready PDB before the workflow continues
- fetched structures are downloaded into the current working directory and then processed exactly like local files
- optional crystallographic waters, ions, and common buffer/solvent components can be removed before setup
- polymer chains can be retained or dropped before ligand processing, which is useful for selecting one asymmetric unit from a deposited structure
- ligand-like non-protein residues can be selected by exact instance, not only by residue name

Example commands:

```bash
# Use a local mmCIF file
python 1_AutomateGromacs.py --pdb model.mmcif

# Fetch directly from the PDB in CIF format (default)
python 1_AutomateGromacs.py --fetch-pdb 9UWJ

# Fetch directly from the PDB in PDB format
python 1_AutomateGromacs.py --fetch-pdb 9UWJ --fetch-format pdb

# MPI-compatible equivalent
python 1_AutomateGromacs_MPI.py --gmx-bin gmx_mpi --fetch-pdb 9UWJ
```

---

<a id="environment-setup"></a>

## 🧪 Environment setup

Create the environments from the YAML files in the repository root:

```bash
conda env create -f environment_cgenff.yml
conda env create -f environment_mdanalysis.yml
```

Or use the helper script:

```bash
bash recreate_envs.sh
```

---

<a id="force-fields"></a>

## 🧬 Force fields included and why BOTH are shipped

This repository ships with **two CHARMM36 variants**:

| Force field folder | Purpose |
|---|---|
| `charmm36.ff` | Classic CHARMM36. Recommended default for many CGenFF-ligand workflows. |
| `charmm36_ljpme-jul2022.ff` | CHARMM36 LJ-PME variant. Useful for workflows validated with LJ-PME assumptions. |

### ✅ Important: PyMACS may automatically try BOTH force fields

PyMACS is designed to be practical for real laboratory inputs. For some systems, especially **ligand / CGenFF** systems, the pipeline may:

1. attempt setup using one force field variant, often **LJ-PME**, then
2. fall back to the other, often **classic CHARMM36**, if the first attempt fails due to known parameter or compatibility issues.

Therefore:

- The repository root should contain **both** force field folders so fallback logic is available.
- For per-system run folders, there are two typical options:
  - **Auto-fallback enabled:** copy **both** force field folders into the run directory.
  - **Manual / controlled setup:** copy only the force field intended for the run and avoid fallback behavior in that workflow.

---

<a id="ligand-parameterization"></a>

## 🧪 Ligand parameterization and CGenFF support

PyMACS supports CGenFF ligands, but **not all CGenFF output formats are interchangeable** across conversion styles and force-field assumptions.

### ✅ Canonical path used by PyMACS

PyMACS is designed around **CGenFF CHARMM-format outputs**:

- `<LIG>.cgenff.mol2`
- `<LIG>.str`

and a local conversion step using:

- `cgenff_charmm2gmx_py3_nx2.py`

This conversion produces GROMACS-compatible ligand include files such as:

- `.itp`
- `.prm`

### ✅ CGenFF version requirement

PyMACS workflows are currently validated against **CGenFF 4.6**.

Use **CGenFF 4.6** for generating `.str` and `.cgenff.mol2` inputs.

Newer CGenFF releases may produce outputs that are not currently compatible with the conversion and include conventions used here.

If unexpected `grompp`, include, or atomtype issues occur while using a newer CGenFF release, regenerate ligand files using **CGenFF 4.6**.

### ⚠️ Common failure mode

Pre-converted GROMACS ligand bundles from servers or third-party converters can differ in:

- include ordering
- atomtype declarations
- duplication or overlap with force-field parameters
- assumptions about LJ-PME vs classic CHARMM36

PyMACS includes multiple modes and workarounds to accommodate common laboratory input formats.

---

<a id="cgenff-modes"></a>

## ✅ The 3 CGenFF modes — PyMACS-supported

### Mode 1 — Import pre-converted GROMACS ligand files

Use this when GROMACS-formatted ligand files already exist, such as:

- `.itp`
- `.prm`
- supporting include files

This is an **advanced / bring-your-own-parameters** mode.

Rules to avoid breakage:

1. Make sure the run directory has a force-field context consistent with the ligand parameters.
2. Keep include duplication under control:
   - exactly one ligand `.itp` include
   - any required ligand parameter include files, if used
   - one `[ molecules ]` entry for the ligand

This mode is powerful, but it assumes the user understands the parameter topology ecosystem being supplied.

---

### Mode 2 — CGenFF server outputs only — recommended

Use this when the CGenFF server or local CGenFF provides:

- `<LIG>.cgenff.mol2`
- `<LIG>.str`

Place them in the run directory before running PyMACS. PyMACS will convert locally using the included converter.

Naming must match the PDB residue name exactly, typically a 3-letter residue name.

Example:

```text
PDB resname: A1D
Required files:
A1D.cgenff.mol2
A1D.str
```

---

### Mode 3 — Fully automated ligand processing with local SILCSBio / CGenFF

Use this when SILCSBio / CGenFF is installed locally and available in `PATH`.

PyMACS will attempt to:

- extract ligand coordinates from the input PDB
- prepare and convert the ligand
- generate CGenFF `.str` and `.cgenff.mol2`
- convert those files into GROMACS includes
- integrate the ligand into `topol.top`

If using SILCSBio, a typical setup looks like:

```bash
export SILCSBIO_HOME=$HOME/silcsbio.2024.1
export PATH="$SILCSBIO_HOME/programs:$SILCSBIO_HOME:$PATH"
```

---

<a id="output-behavior"></a>

## ✅ Output behavior: overwrite-by-default

PyMACS is designed for iterative workflows. Re-running steps **overwrites outputs by default**, including regenerated `.tpr` files, refreshed plots, and updated PDF figurebooks.

To preserve a run state, copy or rename the run folder before re-running.

Recommended structure:

```text
RUNS/SystemA_run01/
RUNS/SystemA_run02/
RUNS/SystemA_run03/
```

---

<a id="quick-start-templates"></a>

## 🚀 Quick start templates

Below are public-facing quick-start templates for both standard local runs and MPI/HPC deployments.

---

### Template A — Standard local workstation

Recommended: do not run in the repository root. Create a per-system run folder.

```bash
mkdir -p RUNS/MySystem_01
cd RUNS/MySystem_01

# Copy core scripts
cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../3C_Interface_RIN.py .
cp ../../3_PROTAC_Analysis.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

# Copy MDP templates
cp ../../em.mdp .
cp ../../MDPs/ions.mdp .
cp ../../MDPs/nvt.mdp .
cp ../../MDPs/npt.mdp .
cp ../../MDPs/md.mdp .

# Copy BOTH force fields so auto-fallback works
cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

# Copy your input structure
cp /path/to/your/input.pdb .
```

Then run:

```bash
# Step 1: setup and ligand integration if present
conda activate cgenff
python 1_AutomateGromacs.py

# Step 2: equilibration and production
# By default, Step 2 auto-selects a usable GPU when available.
conda activate mdanalysis
python 2_AutomateGromacs.py

# Step 3: analysis
conda activate mdanalysis
python 3A_AutomateGromacs.py

# Step 4: figurebook PDF
python 4PDF4MD.py
```

---

### Template B — MPI / HPC run folder

Use this when the cluster provides `gmx_mpi` or `gmx-mpi`, or when you want explicit GROMACS binary control.

```bash
mkdir -p RUNS/MySystem_MPI_01
cd RUNS/MySystem_MPI_01

# Copy MPI-compatible scripts
cp ../../1_AutomateGromacs_MPI.py .
cp ../../2_AutomateGromacs_MPI.py .
cp ../../3A_AutomateGromacs_MPI.py .
cp ../../3_PROTAC_Analysis_MPI.py .
cp ../../3B_NETWORX.py .
cp ../../3C_Interface_RIN.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../4_GraphNotes.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

# Copy MDP templates
cp ../../em.mdp .
cp ../../MDPs/ions.mdp .
cp ../../MDPs/nvt.mdp .
cp ../../MDPs/npt.mdp .
cp ../../MDPs/md.mdp .

# Copy BOTH force fields so auto-fallback works
cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

# Copy your input structure
cp /path/to/your/input.pdb .
```

Typical MPI/HPC usage:

```bash
export PYMACS_GMX_BIN=gmx_mpi

conda activate cgenff
python 1_AutomateGromacs_MPI.py --gmx-bin gmx_mpi

conda activate mdanalysis
python 2_AutomateGromacs_MPI.py --gmx-bin gmx_mpi --external-mpi --mpi-ranks 4 --mpi-launcher "mpirun -np 4"

conda activate mdanalysis
python 3A_AutomateGromacs_MPI.py --gmx-bin gmx_mpi
```

SLURM-style Step 2 example:

```bash
python 2_AutomateGromacs_MPI.py \
  --gmx-bin gmx_mpi \
  --external-mpi \
  --mpi-ranks 4 \
  --mpi-launcher "srun -n 4"
```

---

### Template C — Analysis-only

Use this when MD is already complete and only plots and PDF reporting are needed.

```bash
conda activate mdanalysis
python 3A_AutomateGromacs.py
python 4PDF4MD.py
```

---

### Standard local and MPI-compatible command examples

Standard local workstation:

```bash
python 1_AutomateGromacs.py
python 2_AutomateGromacs.py
python 3A_AutomateGromacs.py
```

MPI / Triton-style:

```bash
export PYMACS_GMX_BIN=gmx_mpi

python 1_AutomateGromacs_MPI.py --gmx-bin gmx_mpi
python 2_AutomateGromacs_MPI.py --gmx-bin gmx_mpi --external-mpi --mpi-ranks 4 --mpi-launcher "mpirun -np 4"
python 3A_AutomateGromacs_MPI.py --gmx-bin gmx_mpi
```

Step 2 MPI-specific notes:

- `--external-mpi` omits thread-MPI `-ntmpi`
- `--ntomp` still controls OpenMP threads
- `--mpi-ranks` controls external MPI ranks
- `--mpi-launcher` can be `mpirun -np N` or `srun -n N`
- GPU fallback profiles are preserved

Step 3A / PROTAC MPI-specific notes:

- the analysis scripts do not run production `mdrun`
- the MPI editions only need binary resolution for GROMACS utilities such as `trjconv` and `make_ndx`
- `3_PROTAC_Analysis_MPI.py` can fall back to MDAnalysis behavior if GROMACS utility centering is unavailable

---

<a id="pipeline-stages"></a>

## 🧩 Pipeline stages — scripts

### 🧱 Step 1 — Setup

Script:

```text
1_AutomateGromacs.py
```

Typical actions:

- reads a local PDB, CIF, or mmCIF file, or fetches a structure directly from the RCSB PDB
- converts CIF/mmCIF inputs into a setup-ready PDB automatically
- runs a novice-friendly PDB preflight check for coordinate formatting, residue names, chain IDs, and ligand ambiguity
- detects polymer chains and can optionally prune the system to selected chains before setup
- runs `pdb2gmx` for protein topology generation
- detects non-water `HETATM` residues as ligand candidates
- supports instance-level selection of the primary ligand and any retained cofactors/context molecules
- integrates ligand parameters using one of the supported CGenFF modes
- solvates and ionizes the system
- writes helper maps such as `atomIndex.txt`, workflow-dependent
- supports configurable periodic box setup through `--box-type` and `--box-distance`
- warns when LINK/CONECT records suggest a covalent ligand and keeps the default workflow focused on non-covalent ligands

---

### ⚙️ Step 2 — Equilibration and production

Script:

```text
2_AutomateGromacs.py
```

Typical actions:

- energy minimization
- NVT equilibration
- NPT equilibration
- production MD
- automatic GPU selection when a usable NVIDIA GPU is available
- CPU fallback when no GPU is available or when a GPU profile fails
- energy-minimization handling that avoids PME-on-GPU for non-dynamical EM integrators
- GPU-accelerated NVT, NPT, and production stages with safer mixed CPU/GPU fallback profiles
- thread tuning through environment variables and script prompts
- supports custom MDP templates through `--mdp-dir`, `--em-mdp`, `--nvt-mdp`, `--npt-mdp`, and `--md-mdp`
- supports user-supplied `index.ndx` files through `--index-file`
- supports explicit resource controls through `--compute`, `--ntomp`, `--ntmpi`, `--gpu-id`, `--gpu-ids`, and `--no-gpu`
- supports optional multi-stage equilibration through `--equilibration-plan`

<a id="step-2-gpu-behavior"></a>

#### Step 2 GPU / CPU behavior in plain language

PyMACS treats energy minimization differently from true MD stages.

- **Energy minimization (`em`)** relaxes bad contacts using a non-dynamical minimizer. Because GROMACS may reject PME-on-GPU for this stage, PyMACS avoids PME-on-GPU during EM and uses a CPU-safe or mixed CPU/GPU path.
- **NVT, NPT, and production MD** are normal time-integration stages. PyMACS tries GPU acceleration first and falls back only if needed.
- **PME** is the long-range electrostatics method used by many solvated biomolecular simulations. It is not triggered by a single special atom; it comes from the MDP electrostatics settings and the periodic solvated system.

For most users, the recommended command remains simple:

```bash
python 2_AutomateGromacs.py --ns 0.25
```

Use CPU-only mode only when troubleshooting or when running on a machine without a usable GPU:

```bash
python 2_AutomateGromacs.py --compute CPU --ns 0.25
```

---

<a id="restart-production-md"></a>

## 🔁 Restarting / Resuming Production MD — checkpoint-safe mode

PyMACS supports **checkpoint-safe resumption** of production molecular dynamics runs using standard GROMACS restart mechanics. This is intended for runs interrupted by walltime limits, node failures, or manual termination.

**Key idea:** the production length is stored in `md_0_1.tpr`. On restart, GROMACS uses the existing `.tpr` and `.cpt`, and does **not** read `md.mdp`.

---

### ✅ Required files

To resume production, the run folder must contain:

- `md_0_1.tpr` — production input created by `gmx grompp`
- `md_0_1.cpt` — production checkpoint written during `gmx mdrun`

Optional but common:

- `md_0_1.xtc`
- `md_0_1.edr`
- `md_0_1.log`

---

### 🚀 Resume production only — GPU / auto GPU

Use this to skip EM / NVT / NPT and continue production only. On a GPU workstation, `--compute auto` selects a usable GPU automatically. You can still add `--compute GPU --gpu 0` when you want to force a specific GPU.

```bash
conda activate mdanalysis

python 2_AutomateGromacs.py \
  --mode ligand \
  --ligand UNK \
  --production_only \
  --resume \
  --headless \
  --ns 250
```

#### Do I need `--ns` for restarts?

Conceptually: **no**. A restart uses `md_0_1.tpr`, so the simulation continues with whatever length is encoded in that `.tpr`.

Practically, in current headless CLI behavior: **yes**. The `--ns` argument avoids an interactive prompt.

To ensure PyMACS does not alter the existing `.tpr`, set `--ns` to the original runtime, for example `250`, or any value that is less than or equal to the existing `.tpr` length.

---

### 🧠 Resume production only — CPU

```bash
conda activate mdanalysis

python 2_AutomateGromacs.py \
  --mode ligand \
  --ligand UNK \
  --compute CPU \
  --production_only \
  --resume \
  --headless \
  --ns 250
```

---

### 🧯 Force a clean restart

Use this when the run should not resume, even if a checkpoint exists.

```bash
conda activate mdanalysis

python 2_AutomateGromacs.py \
  --mode ligand \
  --ligand UNK \
  --ns 250 \
  --compute GPU \
  --gpu 0 \
  --headless \
  --force_restart
```

This runs the normal pipeline stages again:

```text
EM → NVT → NPT → Production
```

---

### 🔎 Verify planned production length in the current `.tpr`

```bash
gmx dump -s md_0_1.tpr | grep -E "delta_t|nsteps"
```

Total production time is calculated as:

```text
time (ns) = nsteps × delta_t(ps) / 1000
```

---

### 📊 Step 3A — Analysis

Script:

```text
3A_AutomateGromacs.py
```

Typical actions, workflow-dependent:

- RMSD
- RMSF
- Rg (radius of gyration)
- contacts
- pocket extraction
- additional structure and interaction analyses depending on installed tools
- configurable binding-pocket extraction through `--pocket-cutoff`
- configurable interaction reporting through `--contact_cutoff`, `--hbond-distance-cutoff`, and `--hbond-angle-cutoff`
- in **PROTAC mode**, `3A_AutomateGromacs.py` automatically dispatches to `3_PROTAC_Analysis.py`
- in **ligand mode** and explicit-ligand workflows (`--ligand`), `3A_AutomateGromacs.py` can call `3B_NETWORX.py` automatically when that helper script is present
- in **ligandless multi-chain protein-protein/protein-peptide** workflows, `3A_AutomateGromacs.py` dispatches to `3C_Interface_RIN.py` unless `--skip-interface-rin` is set
- in **single-chain apo protein** workflows, standard protein-centric analysis is run without interface-RIN dispatch
- in **protein-peptide mode without an explicit ligand**, PyMACS treats the system as an interface-analysis workflow rather than ligand-style NETWORX

---

### 🕸 Step 3B — Network analysis

Script:

```text
3B_NETWORX.py
```

Helper network renderer used by ligand- and peptide-style `3A` analyses.

Important behavior:

- for small-molecule ligand analyses, `3A_AutomateGromacs.py` can invoke `3B_NETWORX.py` automatically to render the ligand-interaction network figure
- for protein-only jobs, `3B_NETWORX.py` is not part of the normal workflow
- for PROTAC jobs, network generation is handled inside `3_PROTAC_Analysis.py` or `3_PROTAC_Analysis_MPI.py` rather than through `3B_NETWORX.py`

---

<a id="step-3c-interface-rin"></a>

### 🕸 Step 3C — Interface RIN analysis

Script:

```text
3C_Interface_RIN.py
```

This script is normally launched by 3A for ligandless multi-chain systems. You can also run it directly if `Final_Trajectory.pdb`, `Final_Trajectory.xtc`, and `atomIndex.txt` already exist.

Step 3C characteristics:

- uses MDAnalysis and NetworkX
- does not require GROMACS directly
- analyzes residue-residue contacts between polymer chains over time
- writes outputs to `Analysis_Results/Interface_RIN/`

---

<a id="interface-rin-analysis"></a>

## 🕸️ Protein-protein / protein-peptide interface RIN analysis

PyMACS now includes a dedicated standalone interface analysis script:

```text
3C_Interface_RIN.py
```

This script is used for ligandless multi-chain systems, including:

- protein-protein complexes
- protein-peptide complexes
- multi-chain assemblies
- ligandless biological interfaces where residue-residue persistence matters

3A detects this workflow automatically when:

- no ligand-centered analysis mode is active
- the run is not in PROTAC mode
- two or more polymer chains are detected from `atomIndex.txt`

In interactive mode, 3A asks:

```text
Detected a ligandless multi-chain system with N chains. Run protein-protein / protein-peptide interface RIN analysis? [Y/n]:
```

In headless mode, 3A dispatches automatically unless `--skip-interface-rin` is passed.

For chain definitions, 3C uses `atomIndex.txt` as the source of truth when available rather than relying on chain IDs alone.

Interactive ligandless multichain analysis:

```bash
conda activate mdanalysis
python 3A_AutomateGromacs.py --mode protein
```

Headless interface analysis through 3A:

```bash
conda activate mdanalysis
python 3A_AutomateGromacs.py \
  --mode protein \
  --headless \
  --interface-contact-cutoff 4.0 \
  --interface-min-contact-frac 0.10 \
  --interface-frame-step 1
```

Restrict to selected chain pairs:

```bash
python 3A_AutomateGromacs.py \
  --mode protein \
  --headless \
  --interface-chain-pairs A:B,A:C
```

Run 3C directly if Final_Trajectory files already exist:

```bash
python 3C_Interface_RIN.py \
  --topo Final_Trajectory.pdb \
  --traj Final_Trajectory.xtc \
  --atomindex atomIndex.txt \
  --outdir Analysis_Results/Interface_RIN \
  --contact-cutoff 4.0 \
  --min-contact-frac 0.10
```

---

### 🧾 Step 4 — PDF figurebook

Script:

```text
4PDF4MD.py
```

`4PDF4MD.py` has two figure-selection modes:

- **Standard ligand / protein / peptide / biological mode:** compiles figures into a PDF report ordered through `4_MDfigs.txt`
- **PROTAC manifest mode:** if `Analysis_Results/PROTAC/QC/protac_figure_manifest.csv` exists, the PDF builder switches automatically to the PROTAC manifest and uses that file list instead of `4_MDfigs.txt`

Standard figure ordering is controlled by:

```text
4_MDfigs.txt
```

Optional custom figure notes can also be supplied through:

```text
4_GraphNotes.txt
```

Interface RIN note:

Interface RIN figures are written to `Analysis_Results/Interface_RIN/`. If you want those figures included in a final PDF figurebook, add the desired PNGs to your figure list/reporting workflow or include them manually until dedicated interface-RIN PDF manifest support is added.

### 🔧 Common customization examples

```bash
# Step 1: change the box geometry and buffer distance
python 1_AutomateGromacs.py --pdb complex.pdb --box-type dodecahedron --box-distance 1.0

# Step 2: run with custom MDP templates and explicit resource controls
# GPU is selected automatically when available; --gpu-id is optional on single-GPU systems.
python 2_AutomateGromacs.py \
  --gpu-id 0 \
  --ntomp 8 \
  --ntmpi 1 \
  --mdp-dir MDPs \
  --nvt-mdp custom_nvt.mdp \
  --npt-mdp custom_npt.mdp \
  --md-mdp custom_md.mdp

# Step 2: use a user-supplied index file and a custom multi-stage equilibration plan
python 2_AutomateGromacs.py \
  --index-file index.ndx \
  --equilibration-plan docs/equilibration_plan_example.json

# Step 3A: widen the ligand pocket or tune interaction reporting
python 3A_AutomateGromacs.py \
  --pocket-cutoff 6.0 \
  --contact_cutoff 4.0 \
  --hbond-distance-cutoff 3.5 \
  --hbond-angle-cutoff 135
```

Advanced protocols such as simulated annealing, alternate integrators, or different thermostat/barostat choices are supported by editing the MDP templates or supplying your own files through these CLI options.

---

<a id="guided-tutorial"></a>

## 🧭 Guided walkthrough template — interactive

This section mirrors the practical run experience, including typical prompts.

---

### 📁 Step 0 — Prepare a clean run directory

```bash
mkdir -p RUNS/MySystem_01
cd RUNS/MySystem_01

# Copy core scripts
cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../3C_Interface_RIN.py .
cp ../../3_PROTAC_Analysis.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../4_GraphNotes.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

# Copy MDP templates
cp ../../em.mdp .
cp ../../MDPs/ions.mdp .
cp ../../MDPs/nvt.mdp .
cp ../../MDPs/npt.mdp .
cp ../../MDPs/md.mdp .

# Copy BOTH force fields so auto-fallback works
cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

# Copy your input structure
cp /path/to/your/input.pdb .
```

---

### 🧱 Step 1 — System preparation

```bash
conda activate cgenff
python 1_AutomateGromacs.py
```

Typical prompt flow:

```text
Available structure files:
1. input.pdb
2. model.cif
3. Fetch a structure directly from the RCSB PDB
Select a structure file by number, or choose the fetch option: 3
Enter a 4-character PDB ID to download (for example 9UWJ): 9UWJ
```

If multiple polymer chains are present:

```text
Selectable polymer chains:
  [1] chain A: protein
  [2] chain B: protein
Keep all detected polymer chains for this setup? [Y/n]: n
Which polymer chains should be retained? [numbers/IDs like 1 or A, comma-separated, default: 1]: 1
```

Ligand/cofactor detection, if present:

```text
Detected non-protein component species: ...
Keep all detected small-molecule residues for this setup? [Y/n]: y
All detected ligands/cofactors will be kept for now.
Which one should be treated as the PRIMARY ligand for analysis? [species number/name or instance label like 1A, default: 1A]: 1A
Optional: choose which ADDITIONAL cofactors/context molecules to keep besides the primary ligand.
Press ENTER to keep all additional cofactors/context molecules, or enter labels like 2A, 2B, "all", or "none":
```

Chain naming / mapping:

```text
Detected chains: A B
Enter descriptive name for chain A: ProteinA
Enter descriptive name for chain B: ProteinB
```

---

### ⚙️ Step 2 — EM → NVT → NPT → Production

```bash
conda activate mdanalysis
python 2_AutomateGromacs.py
```

Typical flow:

- energy minimization using EM-safe settings
- NVT equilibration with GPU acceleration when available
- NPT equilibration with GPU acceleration when available
- production MD with GPU acceleration when available
- CPU fallback if no usable GPU is detected
- duration selected interactively or through CLI arguments

For interrupted production runs, use the checkpoint-safe restart commands above.

---

### 📊 Step 3A — Analysis

```bash
conda activate mdanalysis
python 3A_AutomateGromacs.py
```

Typical prompts:

```text
Detected ligand candidate: A1D
Use A1D as ligand? [Y/n]: Y
Enter compound display name (or press ENTER to use resname): CPD32
Detected 24 CPU threads available
Enter number of threads to use [ENTER = auto]:
```

Important analysis-mode behavior:

- if you select **PROTAC** inside `3A_AutomateGromacs.py`, the script hands off automatically to `3_PROTAC_Analysis.py`
- if you run a **ligand** analysis and `3B_NETWORX.py` is present, `3A_AutomateGromacs.py` can call it automatically to build the ligand interaction network figure
- if you run a **ligandless multi-chain protein-protein/protein-peptide** analysis, `3A_AutomateGromacs.py` dispatches to `3C_Interface_RIN.py` unless `--skip-interface-rin` is used
- if you run **protein-peptide mode without an explicit ligand**, it is treated as interface analysis rather than ligand-style NETWORX
- if you provide `--ligand` explicitly, 3A preserves ligand-centered behavior
- if you run a **single-chain protein-only** or **biological-system** analysis, ligand-style NETWORX and interface-RIN dispatch are skipped

---

### 📄 Step 4 — PDF figurebook

```bash
conda activate mdanalysis
python 4PDF4MD.py
```

For standard ligand / protein / peptide / biological jobs, figure ordering is controlled by:

```text
4_MDfigs.txt
```

Optional custom notes for individual figures can be supplied with:

```text
4_GraphNotes.txt
```

For PROTAC jobs, `4PDF4MD.py` detects the generated manifest automatically:

```text
Analysis_Results/PROTAC/QC/protac_figure_manifest.csv
```

In that case, the PDF builder uses the manifest rather than `4_MDfigs.txt`.

---

<a id="example-figurebook-demo"></a>

## 🧪 Reproduce the included Example figurebook — demo

The repository now includes four curated public examples under `Example_Choices/`, with bundled inputs, prepared systems, and selected completed outputs.

<p align="center">
  <a href="https://github.com/schurerlab/Pymacs/tree/main/Example_Choices">
    <img src="https://img.shields.io/badge/Browse%20Example-GitHub%20Folder-orange?style=for-the-badge&logo=github" alt="Browse the example folder">
  </a>
  <a href="https://github.com/schurerlab/Pymacs/blob/main/Example_Choices/Example1/completed_run/MD_ANALYSIS_FIGUREBOOK.pdf">
    <img src="https://img.shields.io/badge/Open%20Figurebook-PDF-blueviolet?style=for-the-badge&logo=adobeacrobatreader" alt="Open the example figurebook PDF">
  </a>
</p>

### Curated examples

- `Example_Choices/Example1`: `CPD32_9G94` + `A1D` conventional protein-ligand workflow
- `Example_Choices/Example2`: `9UWJ` / `AVPR1A` with `A1E` and `CLR` ligand/cofactor-aware workflow
- `Example_Choices/Example3`: `1URN` RNA/protein biological assembly workflow
- `Example_Choices/Example4`: `5T35` PROTAC ternary-complex workflow with `PTC`

### Git LFS

Large curated example trajectories, GROMACS files, parameter bundles, and figurebooks are tracked with Git LFS.

```bash
git lfs install
git clone <repo>
cd Pymacs
git lfs pull
```

Included files:

- `Example_Choices/Example1/input/CPD32_9G94.pdb` — demo input
- `Example_Choices/Example1/parameters/A1D/A1D.str` — demo ligand stream file
- `Example_Choices/Example1/parameters/A1D/A1D.cgenff.mol2` — demo CGenFF MOL2 file
- `Example_Choices/Example1/completed_run/MD_ANALYSIS_FIGUREBOOK.pdf` — reference figurebook output

### Goal

Run the pipeline using the included example inputs and generate a figurebook PDF in the run directory that matches the structure and content of the included reference output.

Exact byte-for-byte identity can vary across operating systems and matplotlib versions, but the figure order and content types should match when using the same settings.

---

### Demo setup

```bash
mkdir -p RUNS/Example_CPD32
cd RUNS/Example_CPD32

# Copy scripts and templates
cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../3C_Interface_RIN.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../4_GraphNotes.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .
cp ../../em.mdp .
cp ../../MDPs/ions.mdp .
cp ../../MDPs/nvt.mdp .
cp ../../MDPs/md.mdp .
cp ../../MDPs/npt.mdp .

# Copy BOTH force fields, required for auto-fallback demo
cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

# Copy Example inputs
cp ../../Example_Choices/Example1/input/CPD32_9G94.pdb .
mkdir -p A1D
cp ../../Example_Choices/Example1/parameters/A1D/A1D.str ./A1D/
cp ../../Example_Choices/Example1/parameters/A1D/A1D.cgenff.mol2 ./A1D/
```

Run the pipeline:

```bash
conda activate cgenff
python 1_AutomateGromacs.py

conda activate mdanalysis
python 2_AutomateGromacs.py
python 3A_AutomateGromacs.py
python 4PDF4MD.py
```

### What to compare

Compare the generated figurebook against:

```text
../../Example_Choices/Example1/completed_run/MD_ANALYSIS_FIGUREBOOK.pdf
```

If the PDF is missing sections:

- confirm the analysis stage completed successfully
- confirm `4_MDfigs.txt` points to plots that exist in the run folder
- confirm CGenFF **4.6**-generated ligand files are being used, or use the included example ligand files

---

<a id="additional-workflow-documentation"></a>

## 📚 Additional workflow documentation

Additional repository documents that complement the README:

- [`docs/README.md`](docs/README.md) — index of repository documentation in the `docs/` folder
- [`docs/USE_CASE_WALKTHROUGH.md`](docs/USE_CASE_WALKTHROUGH.md) — step-by-step commands, expected inputs, and expected outputs based on the files currently present in this repository
- [`docs/SOFTWARE_AND_CITATIONS.md`](docs/SOFTWARE_AND_CITATIONS.md) — software roles, script locations, and citation placeholders to verify
- [`docs/OUTPUT_GALLERY.md`](docs/OUTPUT_GALLERY.md) — current example outputs plus expected analysis figure types and where they are generated
- [`docs/equilibration_plan_example.json`](docs/equilibration_plan_example.json) — minimal example of a multi-stage equilibration plan accepted by Step 2

These files complement the README with walkthroughs, output examples, software-role summaries, and configuration examples.

---

<a id="tmux-quick-tips"></a>

## 🧠 tmux quick tips

`tmux` is recommended for long local, remote, or server-side runs.

| Action | Command |
|---|---|
| Start a named session | `tmux new -s MDRUN` |
| Detach | `Ctrl+b`, then `d` |
| List sessions | `tmux ls` |
| Re-attach | `tmux attach -t MDRUN` |
| Kill a session | `tmux kill-session -t MDRUN` |

---

<a id="troubleshooting"></a>

## 🧯 Troubleshooting notes

Common causes of `grompp` failures:

- force-field mismatch, especially classic CHARMM36 vs LJ-PME
- duplicated includes or atomtypes when importing pre-converted ligand files
- unexpected CGenFF output conventions, often from newer CGenFF versions
- ligand residue names that do not match the `.str` / `.cgenff.mol2` filenames
- missing force-field folders inside the run directory
- attempting to resume production without both `md_0_1.tpr` and `md_0_1.cpt`

If a run fails early:

- keep both force fields present and allow fallback logic to run
- validate that ligand input files are present and named correctly
- ensure **CGenFF 4.6** is used when generating new `.str` / `.mol2` files
- inspect `topol.top` for duplicated ligand includes
- confirm `gmx` is available in the active environment

Step 2 GPU / EM troubleshooting:

- if Step 2 reports that PME cannot run on GPU during energy minimization, update to the GPU-adaptive Step 2 scripts or run EM in CPU-safe mode
- energy minimization is allowed to run without PME-on-GPU; this does not mean NVT, NPT, or production cannot use the GPU
- `nvidia-smi` may show low GPU use during `grompp`, index generation, or CPU-safe EM; check again during NVT, NPT, or production `mdrun`
- add `--no-gpu` or `--compute CPU` when you intentionally want CPU-only troubleshooting
- use `--gpu-id 0` or another explicit ID when multiple GPUs are present

MPI / GROMACS-binary troubleshooting:

- if you see `gmx: command not found`, use the MPI-compatible scripts plus `--gmx-bin` or `PYMACS_GMX_BIN`
- if `gmx_mpi` exists but `gmx` does not, prefer `1_AutomateGromacs_MPI.py`, `2_AutomateGromacs_MPI.py`, `3A_AutomateGromacs_MPI.py`, and `3_PROTAC_Analysis_MPI.py`
- external MPI builds should not use `-ntmpi`; use `--external-mpi` in `2_AutomateGromacs_MPI.py`
- `--ntomp` still controls OpenMP threading in the MPI-compatible Step 2 script
- `--mpi-ranks` and `--mpi-launcher` control external rank layout for MPI builds
- verify the selected binary with `<binary> --version`
- the analysis MPI editions only require a usable GROMACS binary for helper utilities such as `trjconv` and `make_ndx`
- `3_PROTAC_Analysis_MPI.py` can continue with MDAnalysis-based preprocessing when automatic GROMACS centering utilities are unavailable

---

## 🧬 Input structure notes

- PyMACS expects standard PDB `ATOM` / `HETATM` coordinate records with parseable x/y/z columns.
- Residue names must be present on atom records.
- Chain IDs are recommended for chain-aware naming and atomIndex provenance; the script now warns if they are missing.
- Ligand workflows work best when the ligand residue name is unique and provided explicitly with `--ligand` when multiple non-protein residues are present.
- RNA and other nucleic-acid-containing systems can be prepared when CHARMM/GROMACS residue naming and topology generation are compatible, but these cases should still be validated by the user because the repository does not yet include a dedicated RNA example workflow.
- The default ligand workflow is intended for non-covalent ligands. If `LINK` or ligand-related `CONECT` records suggest a covalent attachment, PyMACS now warns and asks the user to acknowledge that custom topology handling may be required.

---

## 📐 Interaction definitions used in analysis

- Contacts: broad distance-based proximity events using the selected contact cutoff, configurable in `3A_AutomateGromacs.py` through `--contact_cutoff`.
- Hydrogen bonds: reported separately from contacts because they represent a narrower donor/acceptor interaction class rather than generic proximity. The current Script 3A workflow exposes `--hbond-distance-cutoff` and documents `--hbond-angle-cutoff` so these assumptions are explicit in reviewer-facing reports.
- Hydrophobic interactions: residue-chemistry-guided contacts for hydrophobic side chains within the relevant distance threshold used by the analysis script.
- Ionic interactions: residue-chemistry-guided contacts involving charged side chains within the relevant distance threshold used by the analysis script.
- Water-mediated interactions: optional in the newer PROTAC analysis workflow and only computed when the corresponding water-bridge analysis mode is enabled.

These categories are intentionally separated in the reports because hydrogen bonding is a specific geometric or donor/acceptor interaction subset, while generic contacts are broader proximity measures that capture interface persistence even when strict hydrogen-bond criteria are not met.

---

<a id="best-practices"></a>

## 📌 Best practices

- Use one system per folder.
- Keep run folders clean.
- Let the pipeline overwrite outputs inside an active run folder; this is intended behavior.
- Clone or copy the folder if an archive snapshot is needed before re-running.
- For CGenFF ligands, prefer classic `charmm36.ff` unless LJ-PME has been validated end-to-end for the system.
- Keep both force-field folders present for auto-fallback workflows.
- Archive the following for reproducibility:
  - input PDB
  - ligand `.str` / `.cgenff.mol2`
  - generated `.itp` / `.prm`, if produced
  - imported ligand parameter files, if using Mode 1
  - script logs
  - final `.tpr`
  - final `.cpt`
  - analysis plots
  - generated PDF figurebook

---

<a id="repository-description"></a>

## 🧾 Repository description

> PyMACS is a Python-based automation suite for reproducible GROMACS molecular dynamics setup, CGenFF ligand integration, simulation execution, trajectory analysis, and PDF figurebook generation.

---

<a id="provenance-note"></a>

## 📌 Provenance note

This repository contains workflow utilities developed for GROMACS-based molecular dynamics simulations, with a focus on reproducible setup, ligand-aware parameterization, and downstream analysis reporting.

PyMACS is intended for:

- academic research workflows
- structure-based molecular modeling
- protein–ligand MD preparation
- protein–protein and multi-component simulations
- PROTAC / induced-proximity modeling workflows where MD analysis is needed downstream

PyMACS does not replace careful evaluation of molecular dynamics assumptions, force-field compatibility, ligand parameter quality, or system-specific simulation design.

---

<a id="citation"></a>

## 🧬 Citation

If you use PyMACS in academic work, please cite the preprint and/or the software release:

> Schulz, J. M., Reynolds, R. C., & Schürer, S. C. **PyMACS: A Python-Based Automation Suite for GROMACS Molecular Dynamics Setup, Simulation, and Analysis.** SSRN preprint.

<p align="center">
  <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6584500">
    <img src="https://img.shields.io/badge/Cite%20the%20Preprint-SSRN-success?style=for-the-badge" alt="Cite the PyMACS preprint">
  </a>
</p>

---

<a id="contact"></a>

## 📬 Contact

For questions, workflow support, bug reports, or collaboration inquiries:

<p align="center">
  <a href="mailto:jmschulz@med.miami.edu?subject=PyMACS%20Question%20%2F%20Collaboration">
    <img src="https://img.shields.io/badge/Joseph%20M.%20Schulz-jmschulz%40med.miami.edu-blue?style=for-the-badge&logo=gmail" alt="Email Joseph M. Schulz">
  </a>
  <a href="mailto:sschurer@miami.edu?subject=PyMACS%20Correspondence">
    <img src="https://img.shields.io/badge/Stephan%20C.%20Sch%C3%BCrer-sschurer%40miami.edu-lightgrey?style=for-the-badge&logo=gmail" alt="Email Stephan C. Schürer">
  </a>
</p>

- **Joseph M. Schulz** — University of Miami  
- **Stephan C. Schürer** — Corresponding author, University of Miami  

---

<a id="practical-takeaway"></a>

## 🙌 Practical takeaway

Use PyMACS to build reproducible, scriptable GROMACS workflows for:

- molecular dynamics setup,
- CHARMM36 / CGenFF ligand integration,
- MD execution,
- interrupted-run recovery,
- trajectory analysis, and
- PDF figurebook generation for downstream interpretation and reporting.
