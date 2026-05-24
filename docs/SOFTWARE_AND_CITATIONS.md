# Software and Citations

## Overview

PyMACS is a workflow layer around established molecular simulation, parameterization, structure-preparation, analysis, and reporting tools. Users should cite both the PyMACS manuscript and the underlying software that materially contributed to their workflow.

This document is intentionally conservative. It lists software that is clearly used, bundled, referenced, or configured by the current repository. Where an exact citation is not already documented in the repository, a verification note is included rather than guessing.

## Core simulation engine

### GROMACS

PyMACS is built around GROMACS for topology generation support, box setup, solvation, ion placement, energy minimization, equilibration, and production MD.

Repository usage:

- `1_AutomateGromacs.py`
- `1_AutomateGromacs_MPI.py`
- `2_AutomateGromacs.py`
- `2_AutomateGromacs_MPI.py`
- `3A_AutomateGromacs.py`
- `3A_AutomateGromacs_MPI.py`
- `3_PROTAC_Analysis.py`
- `3_PROTAC_Analysis_MPI.py`

Citation status:

- cite the GROMACS papers appropriate to the version used
- verify the exact recommended citation for the local GROMACS version before manuscript submission

## Force fields and parameterization

### CHARMM36

The repository ships with `charmm36.ff`, and the workflow supports CHARMM36-based protein/small-molecule preparation.

Repository usage:

- topology generation in Step 1
- system setup and all downstream MD that depends on the resulting topology

Citation status:

- cite the appropriate CHARMM36 force-field paper for the force-field variant used
- verify exact reference details before final submission

### CHARMM36/LJ-PME

The repository also ships with `charmm36_ljpme-jul2022.ff`. Some workflow branches can auto-select or fall back between classic CHARMM36 and the LJ-PME variant, depending on compatibility and available files.

Repository usage:

- Step 1 setup when the LJ-PME force-field folder is selected or auto-detected

Citation status:

- verify the exact LJ-PME CHARMM36 reference used by the manuscript
- do not cite this variant unless it was actually used in the reported simulations

### CGenFF

PyMACS supports CGenFF-style ligand parameter workflows, especially the canonical `<LIG>.cgenff.mol2` plus `<LIG>.str` input path documented in the README.

Repository usage:

- ligand and cofactor parameterization in Step 1
- primary input route for many small-molecule and cofactor workflows

Citation status:

- cite the CGenFF method and the server/software source actually used
- verify the exact CGenFF citation and version wording before publication

### `cgenff_charmm2gmx_py3_nx2.py` conversion workflow

The repository includes `cgenff_charmm2gmx_py3_nx2.py`, which converts CGenFF outputs into GROMACS-ready `.itp` and `.prm` files.

Repository usage:

- invoked by Step 1 when `.str` plus `.cgenff.mol2` files are supplied

Citation status:

- cite the original `cgenff_charmm2gmx` workflow if it is materially used
- verify the exact citation or repository URL before final release notes or manuscript text

### SILCSBio / server-assisted CGenFF workflows

The README documents an optional automated mode in which local SILCSBio/CGenFF tooling is available in `PATH`.

Repository usage:

- optional advanced ligand-processing path in Step 1

Citation status:

- only cite SILCSBio or related infrastructure if it was actually used
- verify exact software reference and version wording

## Structure preparation

### PDBFixer

PyMACS uses PDBFixer during structure cleanup and repair, including missing-atom and hydrogen-handling tasks.

Repository usage:

- `1_AutomateGromacs.py`
- `1_AutomateGromacs_MPI.py`

Citation status:

- cite PDBFixer if the structure-repair path was used
- verify the preferred reference from the project documentation

### OpenMM

PyMACS imports OpenMM I/O components during structure preparation, including PDB/PDBx read/write handling around repaired structures and CIF/mmCIF processing.

Repository usage:

- `1_AutomateGromacs.py`
- `1_AutomateGromacs_MPI.py`

Citation status:

- cite OpenMM if its preparation path is part of the reported workflow
- verify exact preferred citation for the installed version

## Trajectory analysis

### MDAnalysis

MDAnalysis is a core analysis dependency in the repository and supports trajectory reading, atom selections, contact analysis, coordinate export, and many downstream calculations.

Repository usage:

- `2_AutomateGromacs.py`
- `2_AutomateGromacs_MPI.py`
- `3A_AutomateGromacs.py`
- `3A_AutomateGromacs_MPI.py`
- `3_PROTAC_Analysis.py`
- `3_PROTAC_Analysis_MPI.py`

Citation status:

- cite MDAnalysis when its trajectory-analysis workflow is used

### MDTraj

MDTraj is used for several analysis routines, including RMSD/RMSF-related calculations and DSSP/secondary-structure analysis paths in the current codebase.

Repository usage:

- `3A_AutomateGromacs.py`
- `3A_AutomateGromacs_MPI.py`

Citation status:

- cite MDTraj when those analysis features are used

### NumPy

NumPy underpins numerical calculations across the workflow.

Repository usage:

- all major analysis scripts
- plotting and table-generation code

Citation status:

- generally optional to cite in biomolecular-methods papers unless journal or project policy requires it

### pandas

pandas is used widely for CSV export, tabular summaries, manifest creation, and downstream analysis tables.

Repository usage:

- `3A_AutomateGromacs.py`
- `3_PROTAC_Analysis.py`
- `00_Pymacs_X_analysis.py`

Citation status:

- generally optional, but can be acknowledged in software sections if desired

### matplotlib

matplotlib is the main plotting library for standard PyMACS figures.

Repository usage:

- `3A_AutomateGromacs.py`
- `3A_AutomateGromacs_MPI.py`
- `3B_NETWORX.py`
- `3_PROTAC_Analysis.py`
- `3_PROTAC_Analysis_MPI.py`

Citation status:

- generally optional to cite unless journal or project norms require it

### SciPy

SciPy is used in the analysis stack for statistics and numerical support.

Repository usage:

- `3A_AutomateGromacs.py`
- `3A_AutomateGromacs_MPI.py`
- other analysis helpers depending on environment

Citation status:

- generally optional unless the manuscript relies on specific SciPy statistical implementations

### CuPy

CuPy is present in `environment_mdanalysis.yml` and is imported optionally in the standard analysis script. It appears to provide optional GPU acceleration for supported calculations rather than being a strict requirement.

Repository usage:

- optional acceleration path in `3A_AutomateGromacs.py`
- optional acceleration path in `3A_AutomateGromacs_MPI.py`

Citation status:

- cite only if the accelerated path materially contributed to reported production analyses

### DSSP

The analysis environment includes DSSP, and the current code produces DSSP-based secondary-structure outputs.

Repository usage:

- `3A_AutomateGromacs.py`
- `3A_AutomateGromacs_MPI.py`

Citation status:

- cite DSSP or the relevant secondary-structure assignment method if those results appear in the manuscript
- verify exact reference used by the team

## Chemistry and visualization

### RDKit

RDKit is used in multiple places for chemistry-aware molecule handling and drawing.

Repository usage:

- Step 1 helper logic in structure/component workflows
- `3B_NETWORX.py` for ligand depiction
- `4PDF4MD.py` for ligand graphics in the report
- `00_Pymacs_X_analysis.py` for comparative dashboards

Citation status:

- cite RDKit if chemical depictions or RDKit-based processing are central to the workflow outputs shown

### PyMOL

The repository contains PROTAC visualization helpers that write PyMOL-related outputs such as `load_protac_visualization.pml`. PyMOL is therefore relevant as an optional downstream visualization target, even though it is not the core simulation engine.

Repository usage:

- `3_PROTAC_Analysis.py`
- `3_PROTAC_Analysis_MPI.py`

Citation status:

- cite PyMOL only if PyMOL-rendered views are included in a manuscript or workflow description

### Open Babel

`4PDF4MD.py` and parts of the workflow can call `obabel` for chemistry-format support.

Repository usage:

- molecule conversion/SMILES support in `4PDF4MD.py`
- related ligand-handling paths in setup workflows

Citation status:

- cite Open Babel if those conversion steps are essential to the reported workflow

## Report generation

### ReportLab

ReportLab is used to build the PDF figurebook.

Repository usage:

- `4PDF4MD.py`

Citation status:

- usually optional, but may be acknowledged in a software inventory

### Pillow

Pillow is used for image preparation before PDF assembly.

Repository usage:

- `4PDF4MD.py`

Citation status:

- usually optional

### svglib

`4PDF4MD.py` imports `svglib` to support vector graphics handling in the PDF workflow.

Repository usage:

- `4PDF4MD.py`

Citation status:

- usually optional

## Optional tools and specialized branches

### PROTAC-specific scripts

PyMACS contains dedicated PROTAC analysis scripts:

- `3_PROTAC_Analysis.py`
- `3_PROTAC_Analysis_MPI.py`

These generate PROTAC-specific contacts, RMSD/RMSF, network, geometry, water-bridge, QC, and figure-manifest outputs.

Citation status:

- cite the PyMACS manuscript for the workflow itself
- cite any additional methods explicitly used in the PROTAC manuscript figures, such as PCA, water-bridge analyses, or visualization tools, as appropriate

### Plotly

`00_Pymacs_X_analysis.py` uses Plotly for optional comparative dashboards.

Repository usage:

- optional comparative-analysis helper, not core numbered workflow

Citation status:

- optional unless the dashboard itself is part of distributed supporting material

## Repository citation guidance

If you use PyMACS, cite:

- the PyMACS manuscript
- GROMACS
- the force field actually used
- CGenFF and associated conversion tooling if ligand parameterization was used
- structure-preparation tools such as PDBFixer/OpenMM if that branch was used
- analysis tools such as MDAnalysis, MDTraj, DSSP, and RDKit when those outputs contribute materially to the study

If a workflow branch was not used, it does not need to be cited simply because it is supported by the repository.

## Software and citation table

| Software / method | Role in PyMACS | Required or optional | Used by script(s) | Suggested citation / DOI / URL |
|---|---|---|---|---|
| GROMACS | MD engine, topology preparation support, box/solvent/ion setup, minimization, equilibration, production, trajectory post-processing helpers | Required for most workflows | `1_AutomateGromacs.py`, `2_AutomateGromacs.py`, `3A_AutomateGromacs.py`, MPI counterparts, PROTAC scripts | Verify exact GROMACS citation for the version used |
| CHARMM36 | Biomolecular force field | Required when that force-field branch is used | Step 1 through downstream MD | Verify exact CHARMM36 reference used |
| CHARMM36/LJ-PME | Alternate CHARMM36 force-field variant | Optional, mode-dependent | Step 1 and downstream MD when selected | Verify exact LJ-PME CHARMM36 reference used |
| CGenFF | Ligand/cofactor parameter source | Optional but central for small-molecule workflows | Step 1 scripts | Verify exact CGenFF citation and version wording |
| `cgenff_charmm2gmx` workflow | Converts CGenFF outputs to GROMACS includes | Optional, used in canonical ligand path | Step 1 scripts | Verify exact converter citation or repository URL |
| SILCSBio / local CGenFF workflow | Optional automated ligand-generation branch | Optional | Step 1 scripts | Verify exact software citation if used |
| PDBFixer | Structure repair and cleanup | Optional but commonly used in setup | Step 1 scripts | Verify preferred PDBFixer citation |
| OpenMM | Structure I/O and preparation support | Optional but commonly used in setup | Step 1 scripts | Verify preferred OpenMM citation |
| MDAnalysis | Core trajectory analysis toolkit | Required for analysis workflows | Step 2, Step 3A, PROTAC scripts | Cite MDAnalysis papers |
| MDTraj | RMSD/RMSF and DSSP-related trajectory analysis support | Optional but commonly used in analysis | Step 3A scripts | Cite MDTraj |
| NumPy | Numerical array operations | Practical dependency | Many scripts | Citation optional unless required by journal/project policy |
| pandas | Tables, CSV export, manifests | Practical dependency | Analysis, NETWORX, figurebook, comparative dashboard | Citation optional |
| matplotlib | Plot generation | Practical dependency | Step 3A, NETWORX, PROTAC analysis | Citation optional |
| SciPy | Statistical/numerical helper functions | Optional but present in analysis | Step 3A and helpers | Citation optional; verify if a specific method is central |
| CuPy | Optional GPU-accelerated analysis support | Optional | Step 3A analysis environment | Cite only if used materially |
| DSSP | Secondary-structure assignment | Optional, feature-dependent | Step 3A analysis workflows | Verify exact DSSP citation |
| RDKit | Molecule depiction and cheminformatics support | Optional but important for some outputs | Step 1 helpers, NETWORX, figurebook, comparative dashboard | Cite RDKit if chemically aware depictions/processing are used |
| Open Babel | Format conversion / chemistry support | Optional | Figurebook and some setup/helper paths | Cite Open Babel if used |
| PyMOL | Optional downstream structure visualization target | Optional | PROTAC visualization outputs | Cite only if PyMOL figures are used |
| ReportLab | PDF figurebook assembly | Optional | `4PDF4MD.py` | Citation optional |
| Pillow | Image preparation for PDF assembly | Optional | `4PDF4MD.py` | Citation optional |
| svglib | SVG handling during PDF generation | Optional | `4PDF4MD.py` | Citation optional |
| Plotly | Optional comparative dashboard visualization | Optional | `00_Pymacs_X_analysis.py` | Cite only if dashboard outputs are distributed or central |

## Final verification note

Before manuscript resubmission, verify:

- the exact PyMACS manuscript citation text
- the exact GROMACS citation for the version used in the reported work
- the exact CHARMM36 or CHARMM36/LJ-PME reference actually used
- the exact CGenFF and converter references
- the exact PDBFixer, OpenMM, MDAnalysis, MDTraj, DSSP, and RDKit references if those branches contributed directly to manuscript figures or tables
