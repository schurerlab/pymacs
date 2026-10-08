# Triton → Kyle CGenFF Funnel

This opt-in workflow allows a Triton setup to send an already hydrogenated,
normalized ligand MOL2 to Kyle over HTTPS. Kyle runs the x86-64 SILCSBio CGenFF
binary and returns only the generated CHARMM stream (`.str`). Triton then uses
the existing PyMACS `cgenff_charmm2gmx_py3_nx2.py` conversion and continues its
normal GROMACS setup.

The normal **Protein–small-molecule complex** selection remains unchanged. In
`python pymacs_run.py configure`, choose option **5**:

```text
Protein–small-molecule complex — remote CGenFF via Kyle Funnel
```

The Funnel URL is public, but the service accepts only a token-authenticated
parameterization request. The CGenFF executable and Kyle shell remain private.

## One-time Kyle service setup

Run these commands on Kyle. Do not put the token in a Git repository, a run
folder, or a shell command history.

```bash
# Kyle's system Python does not include venv, so use the existing personal Conda installation.
~/miniconda3/bin/conda create -y -p ~/.local/share/pymacs-cgenff-funnel/conda-env python=3.13 flask waitress

install -d -m 700 ~/.config/pymacs
umask 077
openssl rand -hex 32 > ~/.config/pymacs/cgenff-funnel.token
```

Copy the same token to Triton through a trusted, interactive method. On both
hosts it must be readable only by its owner:

```bash
chmod 600 ~/.config/pymacs/cgenff-funnel.token
```

Start the local-only service on Kyle. Select a maintained SILCSBio release;
the known installed location is shown below as an example.

```bash
~/.local/share/pymacs-cgenff-funnel/conda-env/bin/python /path/to/Pymacs/cgenff_funnel_service.py \
  --token-file ~/.config/pymacs/cgenff-funnel.token \
  --cgenff-exec ~/silcsbio.2025.1/cgenff/cgenff \
  --port 8787
```

In another Kyle terminal, expose *only* that localhost service through
Tailscale Funnel:

```bash
tailscale funnel --bg http://127.0.0.1:8787
tailscale funnel status
```

Record the resulting `https://…ts.net` URL. Confirm that institutional policy
allows this public Funnel before enabling it. Funnel exposes an internet-facing
HTTPS route; the token and strict input validation are still required.

For persistence, run the service and Funnel through a user-level systemd unit
only after the manual end-to-end test passes.

## Triton use

Create the local token file once:

```bash
install -d -m 700 ~/.config/pymacs
umask 077
# Paste the token interactively; do not use a command that records it in history.
read -r -s -p 'CGenFF Funnel token: ' token; echo
printf '%s\n' "$token" > ~/.config/pymacs/cgenff-funnel.token
unset token
chmod 600 ~/.config/pymacs/cgenff-funnel.token
```

During `python pymacs_run.py configure`, choose option 5 and provide:

- the Funnel HTTPS URL;
- `~/.config/pymacs/cgenff-funnel.token`.

The path is saved in `pymacs_run.json`; the token contents are not. Then use
the normal commands:

```bash
python pymacs_run.py validate
python pymacs_run.py setup
```

During Step 1, PyMACS extracts the selected ligand, adds hydrogens, writes and
validates `LIG.cgenff.mol2`, calls the Funnel, receives `LIG.str`, and resumes
the existing topology conversion. A remote failure writes
`LIG.cgenff.stderr.log` and stops setup without submitting an MD job.

## Service contract and limits

`POST /v1/parameterize` accepts only a ligand code and a base64-encoded MOL2.
It accepts no remote paths, arbitrary commands, or executable names. The
request must include `Authorization: Bearer <token>`. Input is limited to 2 MiB
and CGenFF execution to 120 seconds. Each request uses a temporary directory
that is removed immediately after the response.
