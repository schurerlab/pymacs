# Interface explorer coordinate and timeline data

`3C_Interface_RIN.py` remains the scientific source of truth.  Its existing
CSV, NPZ, and figure outputs are unchanged.  The generated HTML payload now
also contains:

- `stable_id` for each residue (`display chain`, MDAnalysis residue index, and
  residue number), while retaining the existing display label for compatibility;
- `position`, taken from C-alpha or (only when C-alpha is absent) the residue
  centre of geometry in the final analysed trajectory frame;
- a base64-encoded row from `interface_residue_pair_presence.npz` for each
  residue pair plus the analysed time axis; and
- an interface-only PDB representing that same final analysed frame.

The presence rows are binary geometric-contact evidence.  They may be used to
recompute geometric occupancy for a selected time interval without rerunning
distance calculations.  They must not be relabelled as hydrogen bonds, ionic
contacts, or other validated chemistry.  Chemistry remains in the existing
validated event outputs.

## Regenerating a completed run

From the completed simulation directory, run the normal 3C command with the
original topology, trajectory, chain selections, and output directory.  If the
saved interface contact tables and `interface_residue_pair_presence.npz` are
still valid, PyMACS' normal cache path regenerates presentation outputs without
altering the trajectory or raw MD results.  Then regenerate the dashboard with
`5_PyMACS_Result_Dashboard.py` using the existing results directory.

Historical runs without `interface_residue_pair_presence.npz` can still show
aggregate network data; time-window filtering must be presented as unavailable.

## Frontend dependency

The interactive 2D graph uses Cytoscape.js 3.34.4, vendored as the MIT-licensed
`pymacs_visualization_assets/cytoscape-3.34.4.min.js` with its license retained
alongside it.  The generator embeds the pinned source directly into every
explorer; npm is only needed by a maintainer when refreshing that asset, never
by the HPC analysis or dashboard viewer.  If a browser refuses the embedded
renderer, the generated SVG RIN remains present as the dependency-free fallback.
