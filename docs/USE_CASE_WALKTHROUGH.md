# Use Case Walkthrough

This walkthrough documents reviewer-visible PyMACS workflows using the four curated examples now included in this repository under `Example_Choices/`.

## Curated example sets in this repo

- `Example_Choices/Example1`
  Purpose: conventional protein-ligand workflow for `CPD32_9G94` with ligand `A1D`.
  Main input file: `Example_Choices/Example1/input/CPD32_9G94.pdb`
  Parameter files: `Example_Choices/Example1/parameters/A1D/A1D.str`, `Example_Choices/Example1/parameters/A1D/A1D.cgenff.mol2`, plus converted `a1d.itp`, `a1d.prm`, and `posre_a1d.itp`.
  `prepared_system/`: solvated and indexed setup outputs such as `solv_ions.gro`, `topol.top`, `index.ndx`, `protein_processed.gro`, and setup registries.
  `completed_run/`: production trajectory files, analysis-ready trajectories, and `MD_ANALYSIS_FIGUREBOOK.pdf`.

- `Example_Choices/Example2`
  Purpose: cofactor-aware receptor workflow for `9UWJ` / `AVPR1A` with retained `A1E` ligand and `CLR` cofactor.
  Main input file: `Example_Choices/Example2/input/9UWJ.cif` with curated PDB exports in the same `input/` folder.
  Parameter files: `Example_Choices/Example2/parameters/A1E.str`, `Example_Choices/Example2/parameters/A1E.cgenff.mol2`, `Example_Choices/Example2/parameters/CLR.str`, `Example_Choices/Example2/parameters/CLR.cgenff.mol2`, plus converted topology files.
  `prepared_system/`: multi-chain topology outputs, indexes, component registries, and solvated system files.
  `completed_run/`: representative production MD files including `md_0_1.tpr`, `md_0_1.gro`, `md_0_1.xtc`, and `md_0_1.edr`.

- `Example_Choices/Example3`
  Purpose: RNA/protein biological assembly workflow for `1URN`.
  Main input file: `Example_Choices/Example3/input/1URN.pdb`
  Parameter files: none required beyond the bundled force-field-derived topology outputs in the example tree.
  `prepared_system/`: RNA and protein chain topology/position-restraint files, indexes, registries, and solvated setup outputs.
  `completed_run/`: representative production files, analysis trajectories, pocket exports, and `MD_ANALYSIS_FIGUREBOOK.pdf`.

- `Example_Choices/Example4`
  Purpose: PROTAC ternary-complex workflow for `5T35` with `PTC`.
  Main input file: `Example_Choices/Example4/input/5T35.pdb`
  Parameter files: `Example_Choices/Example4/parameters/PTC.str`, `Example_Choices/Example4/parameters/PTC.cgenff.mol2`, plus converted `ptc.itp`, `ptc.prm`, and `posre_ptc.itp`.
  `prepared_system/`: ternary-complex topology outputs, indexes, component registries, and solvated system files.
  `completed_run/`: representative production MD files including `md_0_1.tpr`, `md_0_1.gro`, `md_0_1.xtc`, and `md_0_1.edr`.

## Scope and limitations

- The repository now packages four curated example trees that cover protein-ligand, cofactor-aware receptor, RNA/protein assembly, and PROTAC ternary-complex workflows.
- The command walkthrough below still uses `Example1` as the shortest end-to-end reproducible demo.
- The repository also supports additional workflow branches beyond these curated examples, but the examples above are the public reference layouts included directly in the repo.

## Recommended run-directory setup

```bash
mkdir -p RUNS/Example_CPD32
cd RUNS/Example_CPD32

cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../3_PROTAC_Analysis.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../4_GraphNotes.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

cp ../../em.mdp .
cp ../../MDPs/ions.mdp .
cp ../../MDPs/nvt.mdp .
cp ../../MDPs/npt.mdp .
cp ../../MDPs/md.mdp .

cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

cp ../../Example_Choices/Example1/input/CPD32_9G94.pdb .
mkdir -p A1D
cp ../../Example_Choices/Example1/parameters/A1D/A1D.str ./A1D/
cp ../../Example_Choices/Example1/parameters/A1D/A1D.cgenff.mol2 ./A1D/
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

For this included protein-ligand walkthrough, `3A_AutomateGromacs.py` is the main analysis entry point. When the required inputs are present, it can automatically dispatch to `3B_NETWORX.py` for ligand/peptide-style interaction-network rendering. In PROTAC mode, it dispatches to the dedicated PROTAC analysis script instead.

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

- `../../Example_Choices/Example1/completed_run/MD_ANALYSIS_FIGUREBOOK.pdf`

If a run is analyzed in PROTAC mode and `Analysis_Results/PROTAC/QC/protac_figure_manifest.csv` exists, `4PDF4MD.py` switches automatically into its PROTAC figurebook path and writes `PROTAC_MD_ANALYSIS_FIGUREBOOK.pdf` instead of the standard report.

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
