# Example 1: CPD32_9G94 Protein-Ligand Workflow with A1D

## What this example shows

This example shows a standard protein-ligand PyMACS workflow using the starting structure `CPD32_9G94.pdb` and the ligand parameter set `A1D`. It is the easiest curated example in this repository and is the best place for most first-time users to start. The package includes the original input structure, ligand parameter files, setup-stage outputs in `prepared_system/`, and a finished reference run in `completed_run/`. The finished run includes `Final_Trajectory.pdb`, `Final_Trajectory.xtc`, `MD_ANALYSIS_FIGUREBOOK.pdf`, and a large `Analysis_Results/` folder with plots and tables. You can use this example both as a learning tool and as a checklist for what a successful PyMACS run should look like.

## Who this example is for

This example is for beginners who want to learn the normal protein-plus-small-molecule workflow before trying more complex systems. It is also useful for teaching, demos, reviewer checks, and quick environment testing.

## What is already included

| Folder | What it contains | Should beginners edit it? |
|---|---|---|
| `input/` | The starting structure file `CPD32_9G94.pdb`. This is usually the first place to check if a script says a structure file is missing. | Usually no |
| `parameters/` | Ligand force-field files for `A1D`, including `.str`, `.mol2`, `.itp`, `.prm`, and restraint files. Do not rename these unless you also update script inputs and topology references. | No |
| `prepared_system/` | Intermediate setup outputs such as `topol.top`, `index.ndx`, `solv_ions.gro`, `protein_processed.gro`, and registry files. These are helpful for inspection and troubleshooting. | No |
| `completed_run/` | A reference example of what a successful run can produce. You do not need to create this folder yourself. Compare your own outputs against it. | No |

## Before you start

Open a terminal. A terminal is the window where you type commands.

You can copy and paste the commands below exactly. Press Enter after each command.

For this example, the usual workflow is:

1. Use the `cgenff` environment for setup with `1_AutomateGromacs.py`.
2. Use the `mdanalysis` environment for MD, analysis, and figurebook generation.
3. Make sure GROMACS is available.
4. Run one command at a time.

This packaged example folder is best treated as a reference folder. The safest way to run the workflow is to create a new working folder so you do not overwrite the curated reference outputs.

## Step 1 - Go to this example folder

```bash
cd Example_Choices/Example1
```

## Step 2 - Check what files are here

```bash
ls
```

You should see `README.md`, `input`, `parameters`, `prepared_system`, and `completed_run`.

## Step 3 - Activate your environment

Start with the setup environment:

```bash
conda activate cgenff
```

Then check Python and GROMACS:

```bash
python --version
```

```bash
gmx --version
```

If `gmx --version` fails, read the troubleshooting section below before going further.

## Step 4 - Run the preparation script

The commands below make a clean working folder, copy in the real example input files, and then run the preparation script from the repository root.

```bash
mkdir -p ../../RUNS/example1_beginner_test
cd ../../RUNS/example1_beginner_test
cp ../../Example_Choices/Example1/input/CPD32_9G94.pdb .
mkdir -p A1D
cp ../../Example_Choices/Example1/parameters/A1D/A1D.str ./A1D/
cp ../../Example_Choices/Example1/parameters/A1D/A1D.cgenff.mol2 ./A1D/
python ../../1_AutomateGromacs.py --pdb CPD32_9G94.pdb --ligand A1D
```

This script prepares the system for MD. In plain English, it reads the starting structure, handles the protein and ligand setup, builds the topology, prepares the solvated system, and writes the setup files needed for the next stage.

## Step 5 - Run the main MD automation script

Switch to the MD environment first:

```bash
conda activate mdanalysis
```

Then run a short beginner-safe test:

```bash
python ../../2_AutomateGromacs.py --mode ligand --ligand A1D --ns 0.25 --compute CPU --headless
```

This stage runs energy minimization, equilibration, and production MD. It can take time even for a short test.

## Step 6 - Answer the script questions

If you use the commands above exactly, most prompts are skipped because the mode and ligand are already supplied. If you run the scripts without flags, these are the most important choices:

- `ligand`: protein plus small molecule
- `protein`: protein-only system
- `peptide`: peptide-containing system
- `protac`: PROTAC or ternary complex system
- `biological`: biological macromolecule system such as DNA/RNA/protein complexes, if supported by the script

For this example, the correct beginner choice is `ligand`.

If a script asks for a ligand name, use `A1D`. That name comes from the ligand files in `parameters/` and `parameters/A1D/`.

If a script asks about simulation length:

- Press Enter for the default.
- For a quick test, type a short value like `0.25`.
- Do not start with a long run if this is your first time.

For your first test, use a short run such as `0.25 ns`. Once that works, repeat the example with the full recommended runtime.

If a script asks about GPU or CPU, the safest beginner choice is CPU mode unless you already know your computer has a working GROMACS-compatible GPU setup.

## Step 7 - Optional MPI or HPC version

The MPI analysis script is for advanced users, clusters, or multi-core/HPC systems. Beginners should usually skip this unless their instructor or system administrator told them to use it.

```bash
python ../../3A_AutomateGromacs_MPI.py --mode ligand --ligand A1D --headless
```

## Step 8 - Run the analysis script

```bash
python ../../3A_AutomateGromacs.py --mode ligand --ligand A1D --headless
```

This script creates plots, CSV tables, interaction summaries, contact maps, RMSD/RMSF graphs, secondary structure summaries, and network-style outputs in `Analysis_Results/`.

Optional figurebook step:

```bash
python ../../4PDF4MD.py
```

This optional report-building step assembles the finished figures into a PDF figurebook.

## Quick test run

- Start with `0.25 ns` in Step 5.
- Confirm the script finishes without stopping on an error.
- Confirm new output files appear in your run folder.
- Then run a longer simulation later if you want the full workflow.

## Expected output

The packaged `completed_run/` folder is your reference for a successful workflow. Your own run may differ in exact file sizes or timing, but the general output pattern should be similar.

Successful runs may produce:

- structure files such as `.gro`, `.pdb`, and `.tpr`
- trajectory files such as `.xtc`
- energy files such as `.edr`
- analysis plots such as `.png`
- analysis tables such as `.csv`
- PDF figurebooks such as `MD_ANALYSIS_FIGUREBOOK.pdf`

For this example specifically, the reference `completed_run/` folder includes:

- `Final_Trajectory.pdb`
- `Final_Trajectory.xtc`
- `MD_ANALYSIS_FIGUREBOOK.pdf`
- `Analysis_Results/`
- RMSD plots
- RMSF plots
- DSSP and secondary-structure outputs
- contact maps
- ligand interaction summaries
- network figures

## How to know it worked

- The script finished without an error message.
- New output files appeared in your run folder.
- A trajectory file such as `md_0_1.xtc` was created.
- `Analysis_Results/` contains PNG and CSV files.
- Your results broadly resemble the packaged `completed_run/` folder.

## Common problems and easy fixes

### I get "python: command not found"

Python is not available in the active terminal session. Activate the correct Conda environment first and try again.

### I get "conda: command not found"

Conda is either not installed or not initialized in this terminal. Open a terminal where Conda works, or initialize Conda according to your local installation instructions.

### I get "gmx: command not found"

GROMACS is not installed, not loaded, or not visible in this terminal. Activate the environment or module that provides `gmx`, then run `gmx --version` again.

### I am in the wrong folder

```bash
pwd
```

```bash
ls
```

Inside the packaged example folder, you should see `input`, `parameters`, `prepared_system`, `completed_run`, and `README.md`.

### The script cannot find my input file

Check the `input/` folder. For this example, the starting structure is `input/CPD32_9G94.pdb`.

### The script cannot find ligand parameters

Check the `parameters/` folder and the `parameters/A1D/` subfolder. Do not rename the ligand files.

### The simulation fails during grompp

This usually means a topology, parameter, atom name, or environment problem. Copy the last 20-30 lines from the terminal output before asking for help.

### The GPU does not work

Not all computers have a compatible GPU setup for GROMACS. Try CPU mode if the script allows it.

### I do not see plots

Analysis must finish before plots appear. Check `Analysis_Results/` in your run folder.

## Tips for first-time command-line users

- Run one command at a time.
- Copy and paste exactly.
- Do not rename files.
- Do not move files out of the example folder unless you are intentionally copying them into a new run folder.
- If something fails, copy the last 20-30 lines of the error.
- Start with a short test run.
- Compare your output with `completed_run/`.

## Script guide

| Script | What it does | Beginner note |
| ------ | ------------ | ------------- |
| `../../1_AutomateGromacs.py` | Prepares the protein-ligand system from the starting structure and ligand files. | Use this first. |
| `../../2_AutomateGromacs.py` | Runs minimization, equilibration, and production MD. | Start with `--ns 0.25` for a first test. |
| `../../3A_AutomateGromacs.py` | Runs standard post-MD analysis and writes `Analysis_Results/`. | Use this after Step 5 finishes. |
| `../../4PDF4MD.py` | Builds a PDF figurebook from analysis outputs. | Optional, but useful for reports and reviewer-facing output. |

## What not to change

- Do not edit topology files unless you know what they are.
- Do not rename parameter files.
- Do not delete `completed_run/` unless you intentionally want to save space.
- Do not edit `prepared_system/` manually unless you are troubleshooting.

## Next step

Once this example works, try a longer simulation, compare your output with `completed_run/`, then move on to Example2, Example3, or Example4 depending on the kind of system you want to study.
