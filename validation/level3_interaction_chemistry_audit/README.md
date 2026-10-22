# Level 3 interaction-chemistry audit

The shared `pymacs_interaction_chemistry.py` engine is the authoritative source for new chemically labelled events. Its versioned `ValidatedInteractionEvents.csv` distinguishes `validated` events from `proximity_only` contacts. Historical `InteractionTypes_Framewise.csv` files remain readable as diagnostics but are never silently promoted to validated chemistry.

Rules: hydrogen and virtual sites are excluded as hydrophobic or ionic endpoints; ionic events require explicit opposite formal-charge groups; hydrophobic events require nonpolar heavy atoms; hydrogen bonds require donor--attached-H--acceptor geometry (3.5 Å heavy-atom cutoff, 135° DHA angle). Occupancy is unique frame count divided by analyzed eligible frames, grouped by residue, canonical ligand atom, and interaction class.

`run_historical_umf_audit.py` produces a non-destructive screen of UMF607-Z historical output. A full geometry rerun must be performed in the PyMACS MDAnalysis environment; the desktop shell used for this audit has no MDAnalysis installation, so a static CSV cannot establish hydrogen-bond geometry or ligand formal-charge provenance.

Known limitations: ligand formal charge is accepted only when topology exposes an explicit `formalcharge` attribute; force-field partial charge alone is deliberately reported as insufficient to assert a formal ionic state. Solvent-free trajectories cannot produce validated water bridges.
