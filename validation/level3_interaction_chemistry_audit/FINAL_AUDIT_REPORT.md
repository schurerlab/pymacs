# PyMACS Level 3 interaction-chemistry refactor audit

Date: 2026-10-22

## What is now authoritative

`pymacs_interaction_chemistry.py` emits schema `pymacs.interaction-events.v1` in `ValidatedInteractionEvents.csv`. A chemical label is accepted only when its evidence is present:

- **Ionic:** explicitly opposite formal-charge groups within 4.5 Å. A partial charge or absent charge annotation is not silently promoted to a formal ionic state.
- **Hydrophobic:** two nonpolar heavy atoms within 4.5 Å. Hydrogen and virtual sites are excluded.
- **H-bond:** donor, attached hydrogen, acceptor, 3.5 Å heavy-atom distance, and donor-H-acceptor angle >=135°.
- **Water bridge:** both ligand-water and water-protein legs satisfy that H-bond geometry in the same frame. Solvent-free input produces no validated bridge.
- **Other Contact:** retained proximity diagnostic, never used as a chemically specific label.

Occupancy is unique observed frame count divided by the number of analyzed frames. Events are deduplicated at residue/canonical ligand atom/interaction class/frame before aggregation.

## Modified pathways

| Pathway | Validated output |
|---|---|
| `3A_AutomateGromacs.py` | `Analysis_Results/ValidatedInteractionEvents.csv` plus legacy-compatible framewise CSV |
| `3A_AutomateGromacs_MPI.py` | same schema/logic |
| `3A_AutomateGromacs_chunks.py` | same schema/logic |
| `3A_AutomateGromacs_chunks_MPI.py` | same schema/logic |
| `3_PROTAC_Analysis.py` | `Contacts/ValidatedInteractionEvents.csv`, occupancy; water screen distinguished from `ValidatedWaterBridgeEvents.csv` |
| `3_PROTAC_Analysis_MPI.py` | same schema/logic |
| `3C_Interface_RIN.py` | `ValidatedInterfaceInteractionEvents.csv`, atom-level interface chemistry |
| `3B_NETWORX.py` | prefers validated event inputs and warns for historical legacy input |

The interface RIN still retains its geometric contact network. It no longer assigns salt bridge, hydrophobic, or H-bond merely from two residue names.

## Test results actually executed

```text
python -m py_compile pymacs_interaction_chemistry.py 3A_AutomateGromacs.py 3A_AutomateGromacs_MPI.py 3A_AutomateGromacs_chunks.py 3A_AutomateGromacs_chunks_MPI.py 3_PROTAC_Analysis.py 3_PROTAC_Analysis_MPI.py 3C_Interface_RIN.py 3B_NETWORX.py  PASS
tests/test_interaction_chemistry.py   6 passed
tests/test_networx_layout.py          2 passed
git diff --check                      PASS
```

The current desktop Python environment has no `MDAnalysis` installation; full trajectory execution, serial/MPI/chunk equivalence, PBC integration, and actual example reruns were therefore not executed here. This is a blocking environment limitation, not a pass.

## UMF607-Z stress-test audit

Source trajectory/topology are present at the user-supplied UMF607-Z directory, including `Final_Trajectory.pdb`, `Final_Trajectory.xtc`, `binding_pocket_only.pdb`, `binding_pocket_only.xtc`, `UNK.cgenff.mol2`, and historical Level-3 CSVs. Historical data was never modified.

The non-destructive endpoint screen found:

| Class | Historical rows | Rows after excluding H/virtual ligand endpoints | Historical residues | Residues after screen |
|---|---:|---:|---:|---:|
| Ionic | 250,704 | 28,295 | 17 | 17 |
| Hydrophobic | 234,919 | 31,862 | 29 | 20 |
| H-bond | 74,755 | 74,755 | 24 | 24 |

- Hydrogen-anchored ionic artifacts removed by this endpoint rule: **222,409**.
- Hydrogen-anchored hydrophobic artifacts removed: **203,057**.
- Virtual-site ionic/hydrophobic artifacts in this historical CSV: **0**.
- Exact duplicate residue/atom/type/time rows: **0**. The visual crowding instead comes principally from large numbers of distinct but chemically invalid hydrogen-propagated rows, not duplicate CSV rows.

These are diagnostic counts only. The endpoint-screened data is not called “validated” because historical CSVs cannot establish donor-H-acceptor geometry or formal-charge provenance. The resulting artifacts are in `umf607_z/`:

- `historical_sid.png` — original historical SID.
- `umf607_z_endpoint_screen_before_after.png` — labelled endpoint-screen comparison.
- CSVs with the supporting counts.

## Example-regression matrix

| Dataset | Mode | Data availability | Executed here | Status |
|---|---|---|---|---|
| UMF607-Z | ligand / large-molecule stress test | PDB + XTC + historical outputs present | Historical audit only | Blocked: MDAnalysis absent |
| Example_Choices/Example1 | standard ligand | repository reference inputs present | Not rerun | Blocked: MDAnalysis absent |
| Example_Choices/Example2 | standard ligand | reference inputs inspected, completeness not established | Not rerun | Blocked: MDAnalysis absent |
| Example_Choices/Example3 | interface/peptide candidate | reference inputs present | Not rerun | Blocked: MDAnalysis absent |
| Example_Choices/Example4 | PROTAC candidate | reference inputs present | Not rerun | Blocked: MDAnalysis absent |

## Remaining limitations and next execution

1. Run all trajectory commands in the project MDAnalysis environment, including serial/MPI/chunk comparisons at multiple chunk sizes, then replace this audit’s blocked rows with actual metrics.
2. A topology adapter should expose verified ligand formal charges from a source that distinguishes them from force-field partial charges; unknown charge deliberately suppresses ionic labels.
3. Normal 3A pocket trajectories generally omit water. Water bridges require solvent-containing coordinates and are optional, not inferred.
4. This audit preserves legacy read compatibility; legacy files are explicitly warned/labelled unvalidated by NETWORX.
