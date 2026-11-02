# Cheminformatics lab environments

This directory is the source-controlled companion to the shared environments
installed at `/projectnb/cheminfo/share/pymacs-lab/lab-tools`.

The shared, *golden* environments are copied from a known working environment
and are deliberately read-only to ordinary lab members:

| Cluster | Environment | Purpose |
| --- | --- | --- |
| Triton (POWER9) | `cgenff` | PyMACS setup and CGenFF/OpenMM work |
| Triton (POWER9) | `mdanalysis` | trajectory analysis and figurebooks |
| Triton (POWER9) | `vina_env` | Triton Vina workflows |
| Pegasus (x86-64) | `vina_env` | Pegasus docking workflows |

Triton and Pegasus have different CPU architectures, so the environments are
kept in architecture-specific directories. Do not copy an environment folder
between the two platforms.

## New-user setup

After their accounts and SSH access are approved, a user runs this once on
each cluster they will use:

```bash
bash /projectnb/cheminfo/share/pymacs-lab/lab-tools/bootstrap-cluster-envs.sh
```

The command is safe to repeat. It writes a small shell helper in the user's
`~/.config/pymacs/` directory and adds one marked line to `~/.bashrc`. It does
not modify the shared environments.

Open a new shell (or run `source ~/.bashrc`), then use:

```bash
pymacs-env cgenff       # Triton only
pymacs-env mdanalysis   # Triton only
pymacs-env vina_env     # Triton or Pegasus
```

For non-interactive scripts, use an absolute executable path rather than
relying on shell activation:

```bash
/projectnb/cheminfo/share/pymacs-lab/lab-tools/triton/envs/mdanalysis/bin/python analysis.py
```

## Administration and updates

The bootstrap uses the shared copies by default: they are quick to use,
require no package solver, and all users run the same tested software. Updates
must be installed as a new dated release directory, smoke-tested, and then made
the `current` release; never run `conda install` inside a released environment.

The initial golden copies are made with `conda create --clone` from the known
working accounts. The copy operation is platform-local: Triton sources are
cloned on Triton and the Pegasus Vina source is cloned on Pegasus.
