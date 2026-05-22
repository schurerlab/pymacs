# Software and Citations

This file inventories the major external software and Python libraries used by PyMACS, the role each one plays, and where it appears in the repository. Citation text should be verified before manuscript submission.

| Tool / library | Role in PyMACS | Script / file location | Citation status |
|---|---|---|---|
| GROMACS | Topology generation, solvation, ionization, minimization, equilibration, production MD, trajectory processing | `1_AutomateGromacs.py`, `2_AutomateGromacs.py`, `3A_AutomateGromacs.py` | citation to verify |
| CHARMM36 force field | Protein and biomolecular force field used during setup | `charmm36.ff/`, `charmm36_ljpme-jul2022.ff/`, setup scripts | citation to verify |
| CGenFF | Small-molecule parameter generation workflow | `1_AutomateGromacs.py`, `cgenff_charmm2gmx_py3_nx2.py`, `Example/A1D/*` | citation to verify |
| PDBFixer | Missing-atom and hydrogen repair during PDB preprocessing | `1_AutomateGromacs.py` | citation to verify |
| Open Babel (`obabel`) | Adds hydrogens and format conversion for ligand coordinates when needed | `1_AutomateGromacs.py` | citation to verify |
| SILCSBio CGenFF environment | Optional local ligand parameter generation path | `1_AutomateGromacs.py`, `environment_cgenff.yml` | citation to verify |
| MDAnalysis | Trajectory analysis and structure selection | `2_AutomateGromacs.py`, `3A_AutomateGromacs.py`, `3_PROTAC_Analysis.py` | citation to verify |
| MDTraj | Pocket extraction, ligand-neighbor detection, and trajectory utilities | `3A_AutomateGromacs.py` | citation to verify |
| matplotlib | Plot generation for analysis figures and figurebook inputs | `3A_AutomateGromacs.py`, `3B_NETWORX.py`, `4PDF4MD.py` | citation to verify |
| seaborn | Statistical plotting and heatmaps | `3A_AutomateGromacs.py` | citation to verify |
| networkx | Residue-contact network graph construction | `3B_NETWORX.py`, `3_PROTAC_Analysis.py` | citation to verify |
| reportlab | PDF figurebook assembly | `4PDF4MD.py` | citation to verify |
| Pillow / PIL | Image handling during PDF figurebook generation | `4PDF4MD.py` | citation to verify |
| PyMOL | Downstream visualization target supported by exported structures and helper scripts | documented usage, structure outputs from analysis scripts | citation to verify |
| RCSB Protein Data Bank | Source of user-provided structural inputs | user-supplied input PDBs, `Example/CPD32_9G94.pdb` | citation to verify |
| RDKit | not directly confirmed in this repository inspection | no direct usage confirmed in current repo scan | citation to verify / confirm necessity |

## Notes for manuscript revision

- If a tool is only optional or only used in one branch of the workflow, state that clearly in the manuscript.
- PyMACS currently relies on external GROMACS and associated force-field files rather than bundling those engines internally.
- RNA and nucleic-acid-containing systems should be described carefully: support depends on compatible residue naming and topology generation through CHARMM/GROMACS, and a dedicated RNA example is not yet bundled in this repository.
