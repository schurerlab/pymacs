# Lab Job Portal

Internal Flask portal that **actually stages and submits** jobs through the
service account. It calls `scp` and `ssh`, then runs fixed `bsub` commands on
Pegasus or Triton. It does not provide a terminal, accept commands, or trust
paths supplied by a browser.

## Inputs

- **Docking:** a ZIP made from an AutoDock Vina LSF package. It must contain
  one `run_confgen_job.lsf` and one `run_vina_job.lsf`. The portal submits the
  Vina job dependent on the conformer job succeeding.
- **MD:** a configured PyMACS data ZIP. It must contain `pymacs_run.json` and
  the input structure, but no Python scripts. The portal copies the reviewed
  PyMACS checkout already installed on Triton into the new project, then runs
  `validate`, `setup`, and `submit --with-analysis`.

The MD input is intentionally a configured PyMACS folder, not a bare PDB: a
reproducible simulation requires choices such as ligand residue, chains, box,
and run length. This is the safe first portal release. A structure-only wizard
can be added after the lab confirms its default choices.

## Deploy on Kyle

1. Ensure Kyle can SSH non-interactively to both `pegasus` and `triton` as
   `jxs794`, and that `bsub` works on both.
2. Copy `job_portal/` to `/home/jxs794/lab-job-portal/`, configure `.env`, then
   run `bash /home/jxs794/lab-job-portal/install-kyle-service.sh` interactively.
3. Copy `lab-job-portal.env.example` to `/etc/lab-job-portal.env`, set a long
   secret and the correct staging directory for each HPC system, and restrict it to root.
4. Restart the service and verify locally:
   `curl -f http://127.0.0.1:8788/`
5. Expose it only to Tailnet users, for example:
   `tailscale serve --https=443 http://127.0.0.1:8788`

Use `tailscale funnel` only if the lab explicitly wants internet exposure;
Tailnet Serve is the VPN-only choice.
