# PyMACS Curated Examples

This folder contains four curated example packages for PyMACS.

Each example has its own README file with beginner-friendly instructions.

## Open an example

| Example | Best for | Click to open |
| ------- | -------- | ------------- |
| Example 1 | Protein-ligand beginner walkthrough | [Open Example 1](./Example1/README.md) |
| Example 2 | Ligand plus cofactor workflow | [Open Example 2](./Example2/README.md) |
| Example 3 | RNA-protein biological assembly workflow | [Open Example 3](./Example3/README.md) |
| Example 4 | PROTAC ternary-complex workflow | [Open Example 4](./Example4/README.md) |

## Quick links

- [Example 1: Protein-ligand walkthrough](./Example1/README.md)
- [Example 2: Ligand plus cofactor workflow](./Example2/README.md)
- [Example 3: RNA-protein biological assembly workflow](./Example3/README.md)
- [Example 4: PROTAC ternary-complex workflow](./Example4/README.md)

If you are new to PyMACS, start with `Example1` unless another example matches your scientific system better.

The packaged `completed_run/` folders are reference outputs. You do not need to create them yourself. They are there so you can compare your own results with a known successful layout.

In general:

- `input/` contains the starting structure files.
- `parameters/` contains ligand or molecule parameter files.
- `prepared_system/` contains generated setup-stage files.
- `completed_run/` contains example outputs from a finished workflow.

## Which example should I start with?

| Example | Best for | Starting file(s) | Main thing it teaches |
| ------- | -------- | ---------------- | --------------------- |
| `Example1` | First-time users, standard protein-ligand systems | `Example1/input/CPD32_9G94.pdb` | Basic protein-ligand setup, MD, analysis, and figurebook workflow with `A1D` |
| `Example2` | Users with a ligand plus a retained cofactor or context molecule | `Example2/input/9UWJ.cif`, `Example2/input/9UWJ.pdb` | Cofactor-aware setup using ligand `A1E` and cofactor `CLR` |
| `Example3` | RNA-protein or other biological macromolecule complexes | `Example3/input/1URN.pdb` | Biological-system workflow with RNA and protein chains |
| `Example4` | PROTAC or ternary-complex users | `Example4/input/5T35.pdb` | Dedicated PROTAC-style workflow with `PTC`, `VHL`, and `BRD4` context |

## Beginner advice

- Start with `Example1` if you are unsure.
- Read the README inside the example folder before running anything.
- Start with a short test run instead of a long simulation.
- Compare your results with the packaged `completed_run/` folder.
- Do not rename parameter files unless you know exactly what else must change.

## Example README files

- [Example1/README.md](./Example1/README.md)
- [Example2/README.md](./Example2/README.md)
- [Example3/README.md](./Example3/README.md)
- [Example4/README.md](./Example4/README.md)
