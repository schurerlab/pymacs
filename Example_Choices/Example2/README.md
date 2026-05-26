# Example 2: 9UWJ / AVPR1A Ligand-and-Cofactor Workflow

## What this example shows

This example shows a more advanced PyMACS workflow built around the `9UWJ` structure with ligand `A1E` and cofactor `CLR`. It is useful when your system has more than one important non-protein component and you want to keep that extra context during setup and simulation. The `input/` folder contains both the original `9UWJ.cif` and several curated structure files produced during cleanup and chain selection. The `parameters/` folder contains files for both `A1E` and `CLR`. The packaged `completed_run/` folder includes representative finished MD files, but this example does not currently include a packaged `Analysis_Results/` tree.

## Who this example is for

This example is for users who already understand the basic protein-ligand workflow and want to see how PyMACS handles a retained cofactor or membrane-context small molecule. It is also useful for reviewers checking that the repository includes a cofactor-aware example.

## What is already included

| Folder | What it contains | Should beginners edit it? |
|---|---|---|
| `input/` | The original `9UWJ.cif`, a PDB copy, and curated structure files such as `selected_input_source.pdb` and `selected_input_curated.pdb`. This is the first place to check if the script cannot find your structure. | Usually no |
| `parameters/` | Force-field files for `A1E` and `CLR`, including `.str`, `.mol2`, `.itp`, `.prm`, and restraint files. Do not rename these unless you also update script inputs and topology references. | No |
| `prepared_system/` | Setup-stage outputs such as `topol.top`, `index.ndx`, `solv_ions.gro`, processed structure files, and registry JSON files. These are intermediate files and should be treated as read-only unless you know exactly what you are doing. | No |
| `completed_run/` | A reference example of successful MD-stage outputs such as `md_0_1.xtc`, `md_0_1.tpr`, and `md_0_1.edr`. You do not need to create this folder yourself. Use it to compare your own results. | No |

## Before you start

Open a terminal. A terminal is the window where you type commands.

You can copy and paste the commands below exactly. Press Enter after each command.

For this example, use the `cgenff` environment for setup and the `mdanalysis` environment for MD and analysis. Make sure GROMACS is available before starting. Run one command at a time.

This packaged example folder is best used as a reference and file source. The safest way to test it is to make a clean run folder somewhere else in the repository.

## Step 1 - Go to this example folder

```bash
cd Example_Choices/Example2
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

This example includes both a main ligand and a retained cofactor. The commands below create a clean run folder, copy the real example files into it, and run the preparation script with both components named explicitly.

```bash
mkdir -p ../../RUNS/example2_beginner_test
cd ../../RUNS/example2_beginner_test
cp ../../Example_Choices/Example2/input/9UWJ.cif .
cp ../../Example_Choices/Example2/parameters/A1E.str .
cp ../../Example_Choices/Example2/parameters/A1E.cgenff.mol2 .
cp ../../Example_Choices/Example2/parameters/CLR.str .
cp ../../Example_Choices/Example2/parameters/CLR.cgenff.mol2 .
python ../../1_AutomateGromacs.py --pdb 9UWJ.cif --ligand A1E --cofactors CLR
```

This script prepares the receptor system, keeps the requested non-protein components, builds the solvated setup, and writes the files needed for the MD stage.

## Step 5 - Run the main MD automation script

Switch to the MD environment:

```bash
conda activate mdanalysis
```

Then run a short test:

```bash
python ../../2_AutomateGromacs.py --mode ligand --ligand A1E --cofactors CLR --ns 0.25 --compute CPU --headless
```

This short run is much safer for a first test than a long production simulation.

## Step 6 - Answer the script questions

If you use the exact command above, most prompts are skipped. If you run interactively, the most important choices are:

- `ligand`: protein plus small molecule
- `protein`: protein-only system
- `peptide`: peptide-containing system
- `protac`: PROTAC or ternary complex system
- `biological`: biological macromolecule system such as DNA/RNA/protein complexes, if supported by the script

For this example, the correct beginner choice is `ligand`.

If a script asks for the ligand name, use `A1E`. The retained cofactor is `CLR`.

If a script asks about simulation length:

- Press Enter for the default.
- For a quick test, type a short value like `0.25`.
- Do not start with a long run if this is your first time.

For your first test, use a short run such as `0.25 ns`. Once that works, repeat the example with the full recommended runtime.

If a script asks about GPU or CPU, CPU is the safest beginner choice unless you already know your GPU setup works with GROMACS.

## Step 7 - Optional MPI or HPC version

The MPI analysis script is for advanced users, clusters, or multi-core/HPC systems. Beginners should usually skip this unless their instructor or system administrator told them to use it.

```bash
python ../../3A_AutomateGromacs_MPI.py --mode ligand --ligand A1E --cofactors CLR --headless
```

## Step 8 - Run the analysis script

```bash
python ../../3A_AutomateGromacs.py --mode ligand --ligand A1E --cofactors CLR --headless
```

This script creates plots, CSV tables, interaction summaries, contact maps, RMSD/RMSF graphs, and other analysis outputs. Your own run may create a new `Analysis_Results/` folder even though this packaged example currently stops at the MD-stage reference outputs in `completed_run/`.

## Quick test run

- Use `0.25 ns` in Step 5.
- Confirm the script finishes.
- Confirm that MD output files appear.
- Run a longer simulation later if you want a fuller test.

## Expected output

Use the packaged `completed_run/` folder as a reference for what the MD stage should produce.

Successful runs may produce:

- structure files such as `.gro`, `.pdb`, and `.tpr`
- trajectory files such as `.xtc`
- energy files such as `.edr`
- analysis plots such as `.png`
- analysis tables such as `.csv`
- PDF figurebooks if you continue through the later analysis/report steps

For this packaged example, the included reference outputs are:

- `md_0_1.gro`
- `md_0_1.tpr`
- `md_0_1.xtc`
- `md_0_1.edr`

## How to know it worked

- The script finished without an error message.
- New output files appeared in your run folder.
- A trajectory file such as `md_0_1.xtc` was created.
- If you ran the analysis stage, new PNG and CSV files appeared in `Analysis_Results/`.
- The finished files resemble the packaged `completed_run/` reference.

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

Check `input/`. The main starting structure for this example is `9UWJ.cif`.

### The script cannot find ligand parameters

Check `parameters/` and make sure the `A1E` and `CLR` files were copied correctly into your working folder. Do not rename them.

### The simulation fails during grompp

This usually means a topology, parameter, atom name, or environment issue. Copy the last 20-30 lines of terminal output when asking for help.

### The GPU does not work

Not all computers have a compatible GPU. Try CPU mode if needed.

### I do not see plots

Analysis must finish before plots appear. Check `Analysis_Results/` if you ran the analysis script.

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
| `../../1_AutomateGromacs.py` | Prepares the receptor system and keeps the requested ligand/cofactor context. | Use this first. |
| `../../2_AutomateGromacs.py` | Runs minimization, equilibration, and production MD. | Start with `--ns 0.25`. |
| `../../3A_AutomateGromacs.py` | Runs post-MD analysis and can create `Analysis_Results/`. | Optional if you only want to test setup and MD first. |

## What not to change

- Do not edit topology files unless you know what they are.
- Do not rename parameter files.
- Do not delete `completed_run/` unless you intentionally want to save space.
- Do not edit `prepared_system/` manually unless troubleshooting.

## Next step

Once this example works, compare it with Example1 to see the difference between a simple ligand workflow and a ligand-plus-cofactor workflow, then adapt the same idea to your own system if you need retained context molecules.
