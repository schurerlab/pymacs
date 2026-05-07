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
  <a href="https://github.com/schurerlab/Pymacs/tree/main/Example">
    <img src="https://img.shields.io/badge/View%20Example-GitHub%20Folder-orange?style=for-the-badge&logo=github" alt="View the PyMACS example folder">
  </a>
  <a href="https://github.com/schurerlab/Pymacs/blob/main/Example/MD_ANALYSIS_FIGUREBOOK.pdf">
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

<a id="overview"></a>

## 🚀 Overview

**PyMACS** is a Python-based toolkit for automating **molecular dynamics setup**, **simulation execution**, **trajectory analysis**, and **publication-ready figure/report generation** using **GROMACS** and **CHARMM36**.

PyMACS supports:

- **Protein-only systems**
- **Protein–ligand systems**
- **Protein–protein systems**
- **PROTAC / multi-component workflows**

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

- [What’s in this repo](#repository-layout)
- [System requirements](#system-requirements)
- [Environment setup](#environment-setup)
- [Force fields included](#force-fields)
- [Ligand parameterization and CGenFF support](#ligand-parameterization)
- [The 3 CGenFF modes](#cgenff-modes)
- [Quick start templates](#quick-start-templates)
- [Pipeline stages](#pipeline-stages)
- [Restarting / resuming production MD](#restart-production-md)
- [Guided tutorial / walkthrough](#guided-tutorial)
- [Example figurebook demo](#example-figurebook-demo)
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
├── 2_AutomateGromacs.py
├── 3A_AutomateGromacs.py
├── 3B_NETWORX.py
├── 4PDF4MD.py
├── 4_MDfigs.txt
├── cgenff_charmm2gmx_py3_nx2.py
├── charmm36.ff/
├── charmm36_ljpme-jul2022.ff/
├── Example/
│   ├── A1D/
│   │   ├── A1D.cgenff.mol2
│   │   ├── A1D.err
│   │   ├── A1D.mol2
│   │   └── A1D.str
│   ├── CPD32_9G94.pdb
│   └── MD_ANALYSIS_FIGUREBOOK.pdf
├── em.mdp
├── ions.mdp
├── md.mdp
├── npt.mdp
├── nvt.mdp
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
  - Non-MPI binary expected by default.
  - Must be callable as `gmx`.
  - If your system uses `gmx_mpi`, adapt the script calls accordingly.
- **Python:** managed through conda environments.
- **Conda / Mamba:** recommended for reproducibility.
- **tmux:** optional but strongly recommended for long simulations.
- **GPU:** optional but recommended for production MD.

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

Below are two short templates for common use cases.

---

### Template A — Run a new system end-to-end

Recommended: do not run in the repository root. Create a per-system run folder.

```bash
mkdir -p RUNS/MySystem_01
cd RUNS/MySystem_01

# Copy core scripts
cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

# Copy MDP templates
cp ../../em.mdp .
cp ../../ions.mdp .
cp ../../nvt.mdp .
cp ../../npt.mdp .
cp ../../md.mdp .

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
conda activate mdanalysis
python 2_AutomateGromacs.py

# Step 3: analysis
conda activate mdanalysis
python 3A_AutomateGromacs.py

# Optional network analysis
python 3B_NETWORX.py

# Step 4: figurebook PDF
python 4PDF4MD.py
```

---

### Template B — Analysis-only

Use this when MD is already complete and only plots and PDF reporting are needed.

```bash
conda activate mdanalysis
python 3A_AutomateGromacs.py
python 4PDF4MD.py
```

---

<a id="pipeline-stages"></a>

## 🧩 Pipeline stages — scripts

### 🧱 Step 1 — Setup

Script:

```text
1_AutomateGromacs.py
```

Typical actions:

- reads the input PDB
- detects protein chains
- runs `pdb2gmx` for protein topology generation
- detects non-water `HETATM` residues as ligand candidates
- integrates ligand parameters using one of the supported CGenFF modes
- solvates and ionizes the system
- writes helper maps such as `atomIndex.txt`, workflow-dependent

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
- GPU-ready execution if available
- thread tuning through environment variables and script prompts

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

### 🚀 Resume production only — GPU

Use this to skip EM / NVT / NPT and continue production only.

```bash
conda activate mdanalysis

python 2_AutomateGromacs.py \
  --mode ligand \
  --ligand UNK \
  --compute GPU \
  --gpu 0 \
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
- contacts
- pocket extraction
- additional structure and interaction analyses depending on installed tools

---

### 🕸 Step 3B — Network analysis

Script:

```text
3B_NETWORX.py
```

Optional residue / contact network utilities.

---

### 🧾 Step 4 — PDF figurebook

Script:

```text
4PDF4MD.py
```

Compiles analysis outputs into a PDF report ordered through:

```text
4_MDfigs.txt
```

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
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

# Copy MDP templates
cp ../../em.mdp .
cp ../../ions.mdp .
cp ../../nvt.mdp .
cp ../../npt.mdp .
cp ../../md.mdp .

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
Available PDB files:
1) input.pdb
Select a PDB file by number: 1
```

Ligand detection, if present:

```text
Detected ligand candidates: A1D
Use A1D as ligand? [Y/n]: Y
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

- energy minimization
- NVT equilibration
- NPT equilibration
- production MD
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

---

### 📄 Step 4 — PDF figurebook

```bash
conda activate mdanalysis
python 4PDF4MD.py
```

Figure ordering is controlled by:

```text
4_MDfigs.txt
```

---

<a id="example-figurebook-demo"></a>

## 🧪 Reproduce the included Example figurebook — demo

The repository includes a complete example system and a reference figurebook.

<p align="center">
  <a href="https://github.com/schurerlab/Pymacs/tree/main/Example">
    <img src="https://img.shields.io/badge/Browse%20Example-GitHub%20Folder-orange?style=for-the-badge&logo=github" alt="Browse the example folder">
  </a>
  <a href="https://github.com/schurerlab/Pymacs/blob/main/Example/MD_ANALYSIS_FIGUREBOOK.pdf">
    <img src="https://img.shields.io/badge/Open%20Figurebook-PDF-blueviolet?style=for-the-badge&logo=adobeacrobatreader" alt="Open the example figurebook PDF">
  </a>
</p>

Included files:

- `Example/CPD32_9G94.pdb` — demo input
- `Example/A1D/A1D.str` — demo ligand stream file
- `Example/A1D/A1D.cgenff.mol2` — demo CGenFF MOL2 file
- `Example/MD_ANALYSIS_FIGUREBOOK.pdf` — reference figurebook output

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
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .
cp ../../em.mdp .
cp ../../ions.mdp .
cp ../../nvt.mdp .
cp ../../md.mdp .
cp ../../npt.mdp .

# Copy BOTH force fields, required for auto-fallback demo
cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

# Copy Example inputs
cp ../../Example/CPD32_9G94.pdb .
mkdir -p A1D
cp ../../Example/A1D/A1D.str ./A1D/
cp ../../Example/A1D/A1D.cgenff.mol2 ./A1D/
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
../../Example/MD_ANALYSIS_FIGUREBOOK.pdf
```

If the PDF is missing sections:

- confirm the analysis stage completed successfully
- confirm `4_MDfigs.txt` points to plots that exist in the run folder
- confirm CGenFF **4.6**-generated ligand files are being used, or use the included example ligand files

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
