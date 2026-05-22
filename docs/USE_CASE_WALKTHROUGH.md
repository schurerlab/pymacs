# Use Case Walkthrough

This walkthrough documents a reviewer-visible PyMACS workflow using only files that are currently present in this repository.

## Available example inputs in this repo

- `Example/CPD32_9G94.pdb`
- `Example/A1D/A1D.cgenff.mol2`
- `Example/A1D/A1D.str`
- `Example/MD_ANALYSIS_FIGUREBOOK.pdf`

## Scope and limitations

- The repository currently contains one example PDB plus example ligand parameter inputs and a reference figurebook PDF.
- A fully packaged receptor-cofactor-dimer example requested by reviewers is not currently present in the repository.
- TODO for a future expanded example:
  - add a complete run folder with prepared MD outputs
  - add a dedicated RNA-containing example
  - add a covalent-ligand expert example if that workflow is explicitly implemented

## Recommended run-directory setup

```bash
mkdir -p RUNS/Example_CPD32
cd RUNS/Example_CPD32

cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

cp ../../em.mdp .
cp ../../MDPs/ions.mdp .
cp ../../MDPs/nvt.mdp .
cp ../../MDPs/npt.mdp .
cp ../../MDPs/md.mdp .

cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

cp ../../Example/CPD32_9G94.pdb .
mkdir -p A1D
cp ../../Example/A1D/A1D.str ./A1D/
cp ../../Example/A1D/A1D.cgenff.mol2 ./A1D/
```

## Step 1: system preparation

```bash
conda activate cgenff
python 1_AutomateGromacs.py --pdb CPD32_9G94.pdb
```

Useful reviewer-facing variants:

```bash
python 1_AutomateGromacs.py \
  --pdb CPD32_9G94.pdb \
  --box-type dodecahedron \
  --box-distance 1.0

python 1_AutomateGromacs.py \
  --pdb CPD32_9G94.pdb \
  --strict-pdb-validation
```

Expected Step 1 outputs include:

- `protein.pdb`
- `protein_processed.gro`
- `complex.gro`
- `newbox.gro`
- `solv.gro`
- `solv_ions.gro`
- `topol.top`
- `atomIndex.txt`
- `em.tpr`
- `mdrun.log`

## Step 2: equilibration and production

```bash
conda activate mdanalysis
python 2_AutomateGromacs.py
```

Useful reviewer-facing variants:

```bash
python 2_AutomateGromacs.py \
  --compute GPU \
  --gpu-id 0 \
  --ntomp 8 \
  --ntmpi 1

python 2_AutomateGromacs.py \
  --mdp-dir ../../MDPs \
  --em-mdp em.mdp \
  --nvt-mdp nvt.mdp \
  --npt-mdp npt.mdp \
  --md-mdp md.mdp

python 2_AutomateGromacs.py \
  --index-file index.ndx \
  --equilibration-plan ../../docs/equilibration_plan_example.json
```

Expected Step 2 outputs include:

- `em.gro`
- `nvt.tpr`, `nvt.gro`, `nvt.cpt`
- `npt.tpr`, `npt.gro`, `npt.cpt`
- `md_0_1.tpr`, `md_0_1.xtc`, `md_0_1.cpt`, `md_0_1.log`
- updated `mdrun.log`

## Step 3A: trajectory analysis

```bash
conda activate mdanalysis
python 3A_AutomateGromacs.py
```

Useful reviewer-facing variants:

```bash
python 3A_AutomateGromacs.py --pocket-cutoff 6.0

python 3A_AutomateGromacs.py \
  --contact_cutoff 4.0 \
  --hbond-distance-cutoff 3.5 \
  --hbond-angle-cutoff 135
```

Expected Step 3A outputs depend on mode and installed libraries, but commonly include:

- `Final_Trajectory.xtc`
- `Final_Trajectory.pdb`
- `binding_pocket_only.xtc`
- `binding_pocket_only.pdb`
- RMSD plots
- RMSF plots
- radius of gyration (Rg) plots
- contact frequency outputs
- NETWORX-ready tables and figures

## Step 3B: network rendering

```bash
conda activate mdanalysis
python 3B_NETWORX.py
```

Expected Step 3B outputs commonly include:

- residue interaction network figures
- expanded panel figures
- combined PDF panels from the NETWORX stage

## Step 4: PDF figurebook generation

```bash
conda activate mdanalysis
python 4PDF4MD.py
```

Expected Step 4 output:

- `MD_ANALYSIS_FIGUREBOOK.pdf`

Compare the result against:

- `../../Example/MD_ANALYSIS_FIGUREBOOK.pdf`

## Where logs are written

- `mdrun.log` in the run directory
- GROMACS `.log` files such as `md_0_1.log`
- script stdout/stderr in the active terminal or scheduler log

## Configuration points reviewers asked about

- Box geometry and buffer: `1_AutomateGromacs.py --box-type --box-distance`
- PDB validation strictness: `1_AutomateGromacs.py --strict-pdb-validation`
- Covalent-ligand acknowledgment: `1_AutomateGromacs.py --allow-covalent-ligand`
- MDP template location and overrides: `2_AutomateGromacs.py --mdp-dir --em-mdp --nvt-mdp --npt-mdp --md-mdp`
- User-supplied index files: `2_AutomateGromacs.py --index-file`
- CPU/GPU controls: `2_AutomateGromacs.py --ntomp --ntmpi --gpu-id --gpu-ids --no-gpu`
- Multi-stage equilibration: `2_AutomateGromacs.py --equilibration-plan`
- Pocket and interaction cutoffs: `3A_AutomateGromacs.py --pocket-cutoff --contact_cutoff --hbond-distance-cutoff --hbond-angle-cutoff`
