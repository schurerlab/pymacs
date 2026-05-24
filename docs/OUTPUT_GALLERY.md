# Output Gallery

## Overview

PyMACS is designed as a staged workflow. Different scripts generate different classes of outputs:

- `1_AutomateGromacs.py` and `1_AutomateGromacs_MPI.py` prepare the simulation system and build the GROMACS topology.
- `2_AutomateGromacs.py` and `2_AutomateGromacs_MPI.py` perform energy minimization, equilibration, and production MD.
- `3A_AutomateGromacs.py` and `3A_AutomateGromacs_MPI.py` generate analysis-ready trajectories, plots, and summary tables.
- `3B_NETWORX.py` renders ligand-residue or peptide-contact network figures when the required contact tables exist.
- `3_PROTAC_Analysis.py` and `3_PROTAC_Analysis_MPI.py` generate a dedicated PROTAC analysis tree.
- `4PDF4MD.py` assembles figurebook PDFs from the generated outputs.
- `00_Pymacs_X_analysis.py` is an optional comparative dashboard builder for multi-run comparisons.

This guide explains the major output types that are implemented in the current repository. File names can differ by mode, ligand name, chain label, or analysis branch, so patterns are shown where appropriate.

## Setup-stage outputs from `1_AutomateGromacs.py`

Step 1 prepares an input structure, removes or retains components according to the user choices, runs `pdb2gmx`, and assembles the solvated starting system.

Common outputs include:

- `selected_input_source.pdb`
  The setup-ready PDB written after converting a `.cif` or `.mmcif` input, or after normalizing a selected PDB input.
- `selected_input_chain_curated.pdb`
  Optional intermediate structure after polymer-chain pruning.
- `selected_input_curated.pdb`
  Optional curated structure after ligand/cofactor instance filtering.
- `protein.pdb`
  Protein-only structure used for `pdb2gmx`.
- `protein_fixed.pdb`
  PDBFixer/OpenMM-repaired structure before topology generation.
- `protein_processed.gro`
  Coordinate file produced by `gmx pdb2gmx`.
- `topol.top`
  System topology updated with protein, solvent, ions, and retained small molecules.
- `atomIndex.txt`
  PyMACS registry used downstream for chain naming, PROTAC bookkeeping, and figure/report generation.
- `pymacs_system_registry.json`
  Internal registry of chain provenance and setup metadata.
- `pymacs_components.json`
  Internal registry of retained and dropped non-protein components.
- `<LIG>.itp`, `<LIG>.prm`
  GROMACS topology/parameter files produced for ligands or cofactors through the CGenFF conversion workflow.
- `<ligand>_pose_match.pdb`, `<ligand>_all.pdb`, `<ligand>.gro`, `<ligand>_all.gro`
  Small-molecule placement intermediates and coordinate files used when reinserting retained components into the complex.
- `complex.gro`
  Merged protein-plus-retained-component structure.
- `newbox.gro`
  Boxed system after `gmx editconf`.
- `solv.gro`
  Solvated system after `gmx solvate`.
- `ions.tpr`
  Ion-placement input generated before `gmx genion`.
- `solv_ions.gro`
  Final solvated and ionized structure passed to Step 2.
- `em.tpr`
  Energy-minimization input file prepared at the end of setup.

These files are usually written in the run directory.

## Simulation-stage outputs from `2_AutomateGromacs.py`

Step 2 runs the actual MD protocol. Output names follow GROMACS conventions.

Common outputs include:

- `em.gro`, `em.log`, `em.edr`
  Energy-minimized coordinates, log, and energy file.
- `nvt.tpr`, `nvt.gro`, `nvt.cpt`, `nvt.log`, `nvt.edr`
  Constant-volume equilibration outputs.
- `npt.tpr`, `npt.gro`, `npt.cpt`, `npt.log`, `npt.edr`
  Constant-pressure equilibration outputs.
- `md_0_1.tpr`, `md_0_1.xtc`, `md_0_1.trr`, `md_0_1.cpt`, `md_0_1.log`, `md_0_1.edr`
  Production-trajectory files, checkpoints, logs, and energies.
- `mdrun.log`
  Consolidated PyMACS/GROMACS run log in the working directory.
- `index.ndx`, `index_default.ndx`
  Index files used for group selection and later analysis.
- `posre_<LIG>.itp`
  Optional ligand restraint file, depending on the chosen workflow.
- `binding_pocket_only.pdb`, `binding_pocket_only.xtc`
  Optional pocket-focused exports if that extraction path is requested or reused by later analysis.

These files are also usually written in the run directory.

## Analysis-stage outputs from `3A_AutomateGromacs.py`

Step 3A writes most outputs into `Analysis_Results/`.

Common outputs across ligand-bound workflows include:

- `Protein_RMSD.csv`, `Protein_RMSD.png`
- `<LIG>_Ligand_RMSD.csv`, `<LIG>_Ligand_RMSD.png`
- `<LIG>_Ligand_RMSF.csv`, `<LIG>_Ligand_RMSF.png`
- `RMSD_<CHAIN>.csv`, `RMSD_<CHAIN>.png`
- `RMSF_<CHAIN>.csv`, `RMSF_<CHAIN>.png`
- `RMSF_<CHAIN>_per_residue.csv`
- `FullComplex_RMSD.csv`, `FullComplex_RMSD.png`
- `Radius_of_Gyration_Protein.csv`, `Radius_of_Gyration_Protein.png`
- `Radius_of_Gyration_Protein_Ligand_Complex.csv`, `Radius_of_Gyration_Overlay.png`
- `DSSP_raw.csv`, `DSSP_encoded_numeric.csv`
- `DSSP_Heatmap.png`
- `DSSP_chain_<CHAIN>.csv`, `DSSP_chain_<CHAIN>_heatmap.png`
- `SSE_Fractions.csv`, `SSE_Fractions.png`
- `Residue_SSE_Persistence.csv`
- `Residue_Helix_Persistence.png`
- `Residue_Sheet_Persistence.png`
- `Residue_Flexibility.png`
- `Final_Trajectory.pdb`, `Final_Trajectory.xtc`
- `binding_pocket_only.pdb`, `binding_pocket_only.xtc`

Ligand-contact and interaction outputs commonly include:

- `AllContacts_Framewise.csv`
- `FilteredContacts_Framewise.csv`
- `InteractionTypes_Framewise.csv`
- `InteractionTypes_Summary.csv`
- `interaction_grouped_normalized.png`
- `interaction_stacked_normalized.png`
- `interaction_type_violin.png`
- `interaction_heatmap_normalized.png`
- `binding_importance_ranking.png`
- `Composite_Distance_Distributions.png`
- `<LIG>_interaction_timeline.png`
- `<LIG>_interaction_event_timeline.png`
- `<LIG>_persistence_vs_distance.png`
- `<LIG>_SSE_vs_LigandContacts_Dynamics.png`
- `<LIG>_SSE_Contact_Coupling_Timeseries.csv`

Biological macromolecule mode can additionally write chain-interface outputs such as:

- `Biomolecule_Chain_Interaction_Heatmap.png`
- `Biomolecule_Chain_Interaction_Distances.png`

Some outputs are mode-dependent. Protein-only or biological-system analyses intentionally skip ligand-only products when no ligand-style component is present.

## NETWORX / interaction-network outputs from `3B_NETWORX.py`

When `Analysis_Results/InteractionTypes_Summary.csv`, `Analysis_Results/FilteredContacts_Framewise.csv`, and a suitable ligand structure file are available, NETWORX writes outputs into:

- `Analysis_Results/NETWORX/`

Common outputs include:

- `<LIG>_bar.png`
  Contact-frequency or interaction-summary bar panel.
- `<LIG>_network.png`
  Network diagram showing retained ligand-residue contacts.
- `<LIG>_expanded.png`
  Expanded panel version of the ligand-residue network.
- `<LIG>_expanded_clean.png`
  Cleaned expanded panel for presentation use.
- `<LIG>_combined.pdf`
  Combined NETWORX multi-panel PDF.
- `<name>_peptide_network.png`
  Peptide-network variant used in peptide contact mode.

These files are usually generated automatically when `3A_AutomateGromacs.py` dispatches to `3B_NETWORX.py` for ligand/peptide workflows. Users can also run `3B_NETWORX.py` directly.

## PROTAC analysis outputs

If the trajectory is analyzed in PROTAC mode, `3A_AutomateGromacs.py` dispatches to `3_PROTAC_Analysis.py` or `3_PROTAC_Analysis_MPI.py`. Those scripts write a dedicated tree under:

- `Analysis_Results/PROTAC/`

Major subdirectories include:

- `RMSD_RMSF/`
- `Contacts/`
- `Networks/`
- `Geometry/`
- `Structures/`
- `Water_Bridges/` if water-bridge analysis is enabled
- `QC/`

Important PROTAC outputs include:

- `QC/protac_figure_manifest.csv`
- `QC/protac_figure_manifest.json`
- `QC/protac_system_roles.csv`
- `QC/protac_contact_summary.csv`
- `QC/protac_contact_summary_by_component.csv`
- `QC/protac_interaction_analysis_summary.csv`
- `QC/protein_pca.csv`
- `RMSD_RMSF/ligand_rmsd.csv`
- `RMSD_RMSF/protein_rmsd_by_component.csv`
- `RMSD_RMSF/ligand_rmsf.csv`
- `RMSD_RMSF/protein_rmsf_by_component.csv`
- `RMSD_RMSF/radius_of_gyration_by_component.csv`
- `Contacts/ligand_contacts_detailed.csv`
- `Contacts/ligand_contacts_by_component_timeseries.csv`
- `Contacts/ligand_contact_residue_summary.csv`
- `Contacts/protac_contacts_long.csv`
- `Contacts/protac_contact_frequency_by_residue.csv`
- `Networks/protac_interaction_network_all_components.png`
- `Structures/load_protac_visualization.pml`
- `Structures/README_visualization.txt`

The PROTAC branch also creates figure aliases such as `Protein_RMSD.png` and `Radius_of_Gyration_Overlay.png` so downstream figurebook generation can reuse standard names where possible.

## Figurebook/report outputs from `4PDF4MD.py`

`4PDF4MD.py` assembles figures into a PDF report.

Standard output:

- `MD_ANALYSIS_FIGUREBOOK.pdf`

PROTAC-manifest mode output:

- `PROTAC_MD_ANALYSIS_FIGUREBOOK.pdf`

The standard mode uses:

- `4_MDfigs.txt`
- optionally `4_GraphNotes.txt`

The PROTAC mode uses:

- `Analysis_Results/PROTAC/QC/protac_figure_manifest.csv`

The figurebook does not create new scientific measurements. It packages existing plots into a reviewer-friendly report with titles, captions, and optional notes.

## Optional comparative/dashboard outputs

The repository also contains:

- `00_Pymacs_X_analysis.py`

This is an optional comparative dashboard builder for multi-run comparisons such as apo-versus-ligand or ligand-versus-ligand studies. It is not part of the four main numbered stages, but it can generate:

- comparative CSV summaries
- overlay plots
- ligand structure cards
- an HTML dashboard for side-by-side run comparison

Because this script works from previously generated `Analysis_Results/` folders, its outputs are project-dependent and should be treated as optional downstream analyses.

## Common CSV/table outputs

Frequently encountered CSV outputs in PyMACS include:

- RMSD tables
- RMSF tables
- radius-of-gyration tables
- DSSP raw and encoded matrices
- framewise contact tables
- interaction-type summaries
- residue-contact frequency summaries
- PROTAC QC and manifest tables

These tables are the machine-readable counterparts to the plots and are useful for manuscript statistics, downstream plotting, and reproducibility.

## Common figure outputs

Frequently encountered figure outputs include:

- RMSD line plots
- RMSF residue plots
- radius-of-gyration overlays
- DSSP heatmaps
- secondary-structure composition charts
- contact heatmaps
- interaction-type summary plots
- persistence-versus-distance plots
- NETWORX network figures
- PROTAC contact and component-summary figures
- final PDF figurebooks

## How to interpret key analysis outputs

### Protein RMSD

Protein RMSD reports how far the protein backbone or selected reference atoms move relative to a starting or reference structure over time. A plateau often suggests that the simulation has reached a relatively stable conformational ensemble, whereas sustained drift or abrupt shifts can indicate slower conformational rearrangements, domain motion, or incomplete equilibration.

### Ligand RMSD

Ligand RMSD tracks how much the ligand pose changes during the trajectory. Low or bounded ligand RMSD is often consistent with a stable binding mode. Larger excursions can reflect ligand reorientation, partial pocket escape, or true changes in the dominant binding pose.

### RMSF

RMSF measures per-residue or per-atom flexibility. High RMSF values often identify mobile loops, termini, solvent-exposed regions, or flexible ligand atoms, while low RMSF values usually mark rigid cores or persistent contact regions.

### Radius of gyration

Radius of gyration measures how compact the structure remains over time. Stable values suggest that the protein or complex maintains a similar overall size and compactness, whereas a rising trend can indicate expansion, partial unfolding, or changes in domain packing.

### Residue-ligand contact maps

Contact maps show which residues approach the ligand within the chosen cutoff and how persistently those contacts appear across frames. They help identify the binding pocket, distinguish persistent anchor residues from transient contacts, and compare different complexes or trajectories.

### Contact-frequency plots

Contact-frequency plots summarize how often each residue contacts the ligand. High-frequency residues are often strong candidates for binding-pocket definition, mutagenesis prioritization, and medicinal-chemistry interpretation.

### Binding-importance ranking

`binding_importance_ranking.png` is intended as a quick ranking of residues that contribute most strongly to persistent ligand engagement. It is useful for identifying binding “hot spots,” but it should be interpreted together with contact chemistry, distance distributions, and structural context rather than as a standalone free-energy metric.

### Interaction-type summaries

Interaction-type summary tables and plots break contact events into broad classes such as hydrogen-bond, hydrophobic, ionic, aromatic, or other contact categories. These outputs help users understand not only which residues contact the ligand, but also how they contact it chemically.

### NETWORX ligand-residue networks

NETWORX figures convert residue-contact persistence into an interpretable graph. Nodes typically represent the ligand and contact residues, while edges emphasize retained interactions. These figures are especially useful for presentations and manuscript figures because they make recurring contact patterns easier to see than raw time series alone.

### `MD_ANALYSIS_FIGUREBOOK.pdf`

The figurebook is a narrative summary of the run, not a new analysis method. It collects the key plots produced by PyMACS into a single PDF so that users, collaborators, and reviewers can quickly evaluate structural stability, flexibility, contact behavior, and interface organization without manually browsing the output directory.

## Plain-language output table

| Output / file pattern | Generated by | What it shows | Why it matters | Notes |
|---|---|---|---|---|
| `protein_processed.gro` | `1_AutomateGromacs.py` | Protein coordinates after `pdb2gmx` | Confirms the force field and termini were applied successfully | Written in the run directory |
| `topol.top` | `1_AutomateGromacs.py` | System topology for GROMACS | Central topology file used by all later MD steps | Updated as solvent, ions, and ligands are added |
| `atomIndex.txt` | `1_AutomateGromacs.py` | PyMACS chain/component registry | Supports chain naming, PROTAC bookkeeping, and report generation | Required by some downstream branches |
| `complex.gro` | `1_AutomateGromacs.py` | Protein plus retained ligand/cofactor coordinates | Confirms the curated complex before boxing/solvation | Especially useful after ligand/cofactor selection |
| `solv_ions.gro` | `1_AutomateGromacs.py` | Final solvated, ionized system | Starting structure for production pipeline | Usually the last major Step 1 coordinate file |
| `em.gro`, `nvt.gro`, `npt.gro` | `2_AutomateGromacs.py` | Structures after minimization and equilibration stages | Lets users verify that each stage completed and remained stable | Standard GROMACS progression |
| `md_0_1.xtc` | `2_AutomateGromacs.py` | Production trajectory | Primary time-resolved coordinate file for downstream analysis | Naming follows GROMACS defaults |
| `Protein_RMSD.csv`, `Protein_RMSD.png` | `3A_AutomateGromacs.py` | Global protein RMSD over time | Fast assessment of stability and equilibration | Commonly highlighted in reports |
| `<LIG>_Ligand_RMSD.csv`, `<LIG>_Ligand_RMSD.png` | `3A_AutomateGromacs.py` | Ligand positional RMSD | Helps judge whether a binding pose remains stable | Only present in ligand-style workflows |
| `RMSF_<CHAIN>.csv`, `RMSF_<CHAIN>.png` | `3A_AutomateGromacs.py` | Chain-resolved flexibility | Helps localize mobile domains and loops | Chain labels depend on the processed system |
| `Radius_of_Gyration_Protein.png` | `3A_AutomateGromacs.py` | Protein compactness over time | Tracks breathing motions or expansion/contraction | Protein-only and general workflows |
| `Radius_of_Gyration_Overlay.png` | `3A_AutomateGromacs.py` | Protein/complex compactness comparison | Useful for ligand-bound or multi-component comparisons | Mode-dependent |
| `DSSP_Heatmap.png` | `3A_AutomateGromacs.py` | Secondary-structure evolution across time | Reveals structural persistence or unfolding | Requires DSSP-capable analysis environment |
| `FilteredContacts_Framewise.csv` | `3A_AutomateGromacs.py` | Filtered residue-ligand contacts by frame | Powers higher-level contact statistics and NETWORX | Important intermediate table |
| `InteractionTypes_Summary.csv` | `3A_AutomateGromacs.py` | Per-residue interaction-type summary | Supports contact chemistry interpretation | Used directly by NETWORX |
| `binding_importance_ranking.png` | `3A_AutomateGromacs.py` | Ranked view of persistent binding residues | Helps identify likely binding hot spots | Best interpreted alongside contact summaries |
| `Analysis_Results/NETWORX/<LIG>_network.png` | `3B_NETWORX.py` | Ligand-residue interaction network | Presentation-friendly view of persistent contacts | Requires Step 3A contact outputs |
| `Analysis_Results/NETWORX/<LIG>_combined.pdf` | `3B_NETWORX.py` | Combined NETWORX panels | Convenient shareable summary of network panels | Ligand/peptide mode |
| `Analysis_Results/PROTAC/QC/protac_figure_manifest.csv` | `3_PROTAC_Analysis.py` | Inventory of PROTAC figures and captions | Drives automated PROTAC figurebook generation | PROTAC mode only |
| `Analysis_Results/PROTAC/Contacts/protac_contact_frequency_by_residue.csv` | `3_PROTAC_Analysis.py` | Residue-level PROTAC contact frequencies | Helps dissect ternary-complex contact architecture | PROTAC mode only |
| `MD_ANALYSIS_FIGUREBOOK.pdf` | `4PDF4MD.py` | Standard PyMACS analysis PDF | Consolidates key plots for users and reviewers | Built from existing figures |
| `PROTAC_MD_ANALYSIS_FIGUREBOOK.pdf` | `4PDF4MD.py` | PROTAC figurebook PDF | Consolidates PROTAC-specific figures through manifest mode | Only generated when PROTAC manifest exists |
