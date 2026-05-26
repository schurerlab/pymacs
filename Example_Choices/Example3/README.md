# Example 3: 1URN Biological RNA-Protein Assembly Workflow

## What this example shows

This example shows a biological macromolecule workflow built around `1URN`. It includes both RNA chains and protein chains, so it is different from the simpler small-molecule examples. The packaged `prepared_system/` folder shows that this example contains RNA chains `RRNA-A`, `RNA-B`, and `RNA-C`, plus protein chains `UOE-A`, `UOE-B`, and `UOE-C`. The `completed_run/` folder includes a finished `MD_ANALYSIS_FIGUREBOOK.pdf`, `Final_Trajectory` files, and an `Analysis_Results/` folder with chain-vs-chain biomolecule outputs. This is a good example for learning how PyMACS handles biological assemblies rather than ligand binding.

## Who this example is for

This example is for users working on RNA-protein or other biological macromolecule complexes. It is also helpful for reviewers and students who want to see that PyMACS supports non-ligand biological assemblies.

## What is already included

| Folder | What it contains | Should beginners edit it? |
|---|---|---|
| `input/` | The starting structure file `1URN.pdb`. This is the first place to check if a script cannot find the structure. | Usually no |
| `parameters/` | This folder is present, but no separate ligand parameter files are packaged here. | No |
| `prepared_system/` | Intermediate setup outputs such as `topol.top`, `index.ndx`, `solv_ions.gro`, processed structures, and chain-specific topology/restraint files for the RNA and protein components. Do not edit these unless you know exactly why you are doing it. | No |
| `completed_run/` | A reference example of what a successful run can produce. You do not need to create this folder yourself. Compare your own outputs against it. | No |

## Before you start

Open a terminal. A terminal is the window where you type commands.

You can copy and paste the commands below exactly. Press Enter after each command.

For this example, use the `cgenff` environment for the preparation stage and the `mdanalysis` environment for MD and analysis. Make sure GROMACS is available. Run one command at a time.

This packaged example is best treated as a reference folder. The safest beginner workflow is to copy the starting input into a clean run folder rather than writing directly into the curated example package.

## Step 1 - Go to this example folder

```bash
cd Example_Choices/Example3
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

Because this is a biological assembly, the preparation step is best run interactively so you can confirm the chain and component handling.

```bash
mkdir -p ../../RUNS/example3_beginner_test
cd ../../RUNS/example3_beginner_test
cp ../../Example_Choices/Example3/input/1URN.pdb .
python ../../1_AutomateGromacs.py --pdb 1URN.pdb
```

This script reads the starting structure, prepares the macromolecular system, writes the processed topology and index files, and builds the solvated setup for MD.

## Step 5 - Run the main MD automation script

Switch to the MD environment:

```bash
conda activate mdanalysis
```

Then run a short biological-system test:

```bash
python ../../2_AutomateGromacs.py --mode biological --ns 0.25 --compute CPU --headless
```

This runs minimization, equilibration, and a short production trajectory for a first-pass test.

## Step 6 - Answer the script questions

If you use the exact command above, most prompts are skipped. If you run interactively, the most important choices are:

- `ligand`: protein plus small molecule
- `protein`: protein-only system
- `peptide`: peptide-containing system
- `protac`: PROTAC or ternary complex system
- `biological`: biological macromolecule system such as DNA/RNA/protein complexes, if supported by the script

For this example, the correct choice is `biological`.

If a script asks about simulation length:

- Press Enter for the default.
- For a quick test, type a short value like `0.25`.
- Do not start with a long run if this is your first time.

For your first test, use a short run such as `0.25 ns`. Once that works, repeat the example with the full recommended runtime.

If a script asks about GPU or CPU, CPU is the safest beginner choice unless you already know your GPU setup works correctly.

## Step 7 - Optional MPI or HPC version

The MPI analysis script is for advanced users, clusters, or multi-core/HPC systems. Beginners should usually skip this unless their instructor or system administrator told them to use it.

```bash
python ../../3A_AutomateGromacs_MPI.py --mode biological --headless
```

## Step 8 - Run the analysis script

```bash
python ../../3A_AutomateGromacs.py --mode biological --headless
```

This script creates plots, CSV tables, secondary-structure summaries, RMSD/RMSF outputs, biomolecule chain-interaction summaries, and other analysis products in `Analysis_Results/`.

Optional figurebook step:

```bash
python ../../4PDF4MD.py
```

This optional step can assemble the analysis figures into a PDF report.

## Quick test run

- Use `0.25 ns` in Step 5.
- Confirm the script finishes.
- Confirm output files appear.
- Then run a longer simulation later if you need a fuller result.

## Expected output

Use the packaged `completed_run/` folder as your reference.

Successful runs may produce:

- structure files such as `.gro`, `.pdb`, and `.tpr`
- trajectory files such as `.xtc`
- energy files such as `.edr`
- analysis plots such as `.png`
- analysis tables such as `.csv`
- PDF figurebooks if present

For this example, the packaged `completed_run/` folder includes:

- `Final_Trajectory.pdb`
- `Final_Trajectory.xtc`
- `MD_ANALYSIS_FIGUREBOOK.pdf`
- `Analysis_Results/`
- chain-interaction heatmaps and timeseries files
- DSSP and secondary-structure outputs
- RMSD and RMSF outputs for both RNA and protein chains

## How to know it worked

- The script finished without an error message.
- New output files appeared in your run folder.
- A trajectory file was created.
- `Analysis_Results/` contains PNG and CSV files.
- The results broadly resemble the packaged `completed_run/` folder.

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

Check `input/`. For this example, the starting structure is `input/1URN.pdb`.

### The script cannot find ligand parameters

This example does not include a separate ligand parameter set. If a script is expecting ligand files, double-check that you selected the biological workflow rather than a ligand workflow.

### The simulation fails during grompp

This usually means a topology, parameter, atom name, or environment issue. Copy the last 20-30 lines of terminal output before asking for help.

### The GPU does not work

Not all computers have a compatible GPU. Try CPU mode if needed.

### I do not see plots

Analysis must finish before plots appear. Check `Analysis_Results/`.

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
| `../../1_AutomateGromacs.py` | Prepares the RNA-protein starting structure for MD. | Run this first. |
| `../../2_AutomateGromacs.py` | Runs minimization, equilibration, and production MD. | Use `--mode biological`. |
| `../../3A_AutomateGromacs.py` | Runs post-MD analysis for biomolecule systems. | This is the main analysis entry point for this example. |
| `../../4PDF4MD.py` | Builds a PDF figurebook from analysis outputs. | Optional but useful for reporting. |

## What not to change

- Do not edit topology files unless you know what they are.
- Do not rename packaged files.
- Do not delete `completed_run/` unless you intentionally want to save space.
- Do not edit `prepared_system/` manually unless troubleshooting.

## Next step

Once this example works, compare its outputs with Example1 so you can see the difference between ligand-centered analysis and biological assembly analysis, then adapt the workflow to your own RNA-protein or nucleic-acid system.
