# Example 4: 5T35 PROTAC Ternary-Complex Workflow with PTC

## What this example shows

This example shows a PROTAC-style ternary-complex workflow using the starting structure `5T35.pdb` and the small-molecule component `PTC`. The packaged `prepared_system/` folder shows two protein partners, `VHL` and `BRD4`, which is typical of ternary-complex systems. The `parameters/` folder contains the packaged `PTC` parameter files. The `completed_run/` folder includes representative finished MD files such as `md_0_1.xtc`, `md_0_1.tpr`, and `md_0_1.edr`. This example is the best place to look if you want to learn how the repository handles PROTAC or ternary-complex analysis.

## Who this example is for

This example is for users working on PROTAC systems, ternary complexes, or other multi-component targeted-degradation style simulations. It is also useful for advanced students and reviewers who want to inspect the dedicated PROTAC analysis path.

## What is already included

| Folder | What it contains | Should beginners edit it? |
|---|---|---|
| `input/` | The starting structure file `5T35.pdb`. This is usually the first folder to check if a script says a structure is missing. | Usually no |
| `parameters/` | Force-field files for `PTC`, including `.str`, `.mol2`, `.itp`, `.prm`, and restraint files. Do not rename these unless you also update script inputs and topology references. | No |
| `prepared_system/` | Intermediate setup outputs such as `topol.top`, `index.ndx`, `solv_ions.gro`, `protein_processed.gro`, and registry files. These are setup-stage files and should usually be treated as read-only. | No |
| `completed_run/` | A reference example of successful MD-stage outputs. You do not need to create this folder yourself. Use it to compare your own results. | No |

## Before you start

Open a terminal. A terminal is the window where you type commands.

You can copy and paste the commands below exactly. Press Enter after each command.

For this example, use the `cgenff` environment for setup and the `mdanalysis` environment for MD and PROTAC analysis. Make sure GROMACS is available. Run one command at a time.

This packaged example is safest as a reference folder. The recommended beginner workflow is to create a clean run folder and copy the needed input files into it.

## Step 1 - Go to this example folder

```bash
cd Example_Choices/Example4
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

## Step 4 - Run the preparation script

These commands create a clean run folder, copy the packaged PTC files into it, and run the preparation script.

```bash
mkdir -p ../../RUNS/example4_beginner_test
cd ../../RUNS/example4_beginner_test
cp ../../Example_Choices/Example4/input/5T35.pdb .
cp ../../Example_Choices/Example4/parameters/PTC.str .
cp ../../Example_Choices/Example4/parameters/PTC.cgenff.mol2 .
python ../../1_AutomateGromacs.py --pdb 5T35.pdb --ligand PTC
```

This script prepares the ternary-complex system, writes the setup files, and creates the inputs needed for MD.

## Step 5 - Run the main MD automation script

Switch to the MD environment:

```bash
conda activate mdanalysis
```

Then run a short PROTAC test:

```bash
python ../../2_AutomateGromacs.py --mode protac --ligand PTC --ns 0.25 --compute CPU --headless
```

This gives you a short first-pass trajectory before you attempt a longer production run.

## Step 6 - Answer the script questions

If you use the exact command above, most prompts are skipped. If you run interactively, these are the key choices:

- `ligand`: protein plus small molecule
- `protein`: protein-only system
- `peptide`: peptide-containing system
- `protac`: PROTAC or ternary complex system
- `biological`: biological macromolecule system such as DNA/RNA/protein complexes, if supported by the script

For this example, the correct choice is `protac`.

If a script asks for the ligand or PROTAC name, use `PTC`. That name comes from the files in `parameters/`.

If a script asks about simulation length:

- Press Enter for the default.
- For a quick test, type a short value like `0.25`.
- Do not start with a long run if this is your first time.

For your first test, use a short run such as `0.25 ns`. Once that works, repeat the example with the full recommended runtime.

If a script asks about GPU or CPU, CPU is the safest beginner choice unless you already know your GROMACS GPU setup works.

## Step 7 - Optional MPI or HPC version

The MPI PROTAC analysis script is for advanced users, clusters, or multi-core/HPC systems. Beginners should usually skip this unless their instructor or system administrator told them to use it.

```bash
python ../../3_PROTAC_Analysis_MPI.py --ligand PTC --quick-test --headless
```

## Step 8 - Run the analysis script

```bash
python ../../3_PROTAC_Analysis.py --ligand PTC --quick-test --headless
```

This script creates PROTAC-focused plots, CSV tables, contact summaries, RMSD/RMSF outputs, QC files, interaction networks, and other ternary-complex analysis products. When you are ready for a fuller analysis, rerun it without `--quick-test`.

## Quick test run

- Use `0.25 ns` in Step 5.
- Use `--quick-test` in Step 8.
- Confirm the scripts finish.
- Confirm output files appear.
- Then run a longer simulation and a fuller analysis later.

## Expected output

Use the packaged `completed_run/` folder as the reference for what the MD stage should produce.

Successful runs may produce:

- structure files such as `.gro`, `.pdb`, and `.tpr`
- trajectory files such as `.xtc`
- energy files such as `.edr`
- analysis plots such as `.png`
- analysis tables such as `.csv`
- additional PROTAC-specific QC and network outputs if you run the dedicated analysis script

For this packaged example, the included reference outputs are:

- `md_0_1.gro`
- `md_0_1.tpr`
- `md_0_1.xtc`
- `md_0_1.edr`

## How to know it worked

- The script finished without an error message.
- New output files appeared in your run folder.
- A trajectory file such as `md_0_1.xtc` was created.
- If you ran the PROTAC analysis stage, new QC, PNG, and CSV outputs appeared.
- The main MD files resemble the packaged `completed_run/` reference.

## Common problems and easy fixes

### I get "python: command not found"

Python is not available in the active terminal. Activate the correct Conda environment first.

### I get "conda: command not found"

Conda is not installed or not initialized in this terminal.

### I get "gmx: command not found"

GROMACS is not installed, not loaded, or not visible in this terminal.

### I am in the wrong folder

```bash
pwd
```

```bash
ls
```

Inside the packaged example folder, you should see `input`, `parameters`, `prepared_system`, `completed_run`, and `README.md`.

### The script cannot find my input file

Check `input/`. For this example, the starting structure is `input/5T35.pdb`.

### The script cannot find ligand parameters

Check `parameters/` and make sure the `PTC` files were copied correctly. Do not rename them.

### The simulation fails during grompp

This usually means a topology, parameter, atom name, or environment issue. Copy the last 20-30 lines of terminal output before asking for help.

### The GPU does not work

Not all computers have a compatible GPU. Try CPU mode if needed.

### I do not see plots

Analysis must finish before plots appear. For PROTAC analysis, look in the output directory created by `3_PROTAC_Analysis.py`.

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
| `../../1_AutomateGromacs.py` | Prepares the ternary-complex starting structure for MD. | Run this first. |
| `../../2_AutomateGromacs.py` | Runs minimization, equilibration, and production MD. | Use `--mode protac`. |
| `../../3_PROTAC_Analysis.py` | Runs the dedicated PROTAC post-processing workflow. | This is the main analysis script for this example. |
| `../../3_PROTAC_Analysis_MPI.py` | Advanced MPI/HPC-flavored PROTAC analysis entry point. | Usually skip this unless you were told to use it. |

## What not to change

- Do not edit topology files unless you know what they are.
- Do not rename parameter files.
- Do not delete `completed_run/` unless you intentionally want to save space.
- Do not edit `prepared_system/` manually unless troubleshooting.

## Next step

Once this example works, try a longer simulation, rerun the dedicated PROTAC analysis without `--quick-test`, compare the outputs with your own ternary-complex system, and then adapt the workflow to your project.
