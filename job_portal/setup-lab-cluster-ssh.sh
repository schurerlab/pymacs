#!/usr/bin/env bash
# Run interactively as jxs794 on Kyle.  It gives each lab workstation its own
# SSH key and authorizes that public key on Acorn, Triton, and Pegasus.  Token
# is deliberately not in WORKSTATIONS: its SSH setup is already managed.
set -euo pipefail

die() { printf 'Error: %s\n' "$*" >&2; exit 1; }
note() { printf '\n==> %s\n' "$*"; }
warn() { printf 'Warning: %s\n' "$*" >&2; }

[[ -t 0 && -t 1 ]] || die 'Run this from an interactive terminal on Kyle.'
[[ "$(hostname -s | tr '[:lower:]' '[:upper:]')" == 'KYLE' ]] || die 'This script must run on Kyle.'
command -v ssh >/dev/null || die 'ssh is not installed.'
command -v ssh-keygen >/dev/null || die 'ssh-keygen is not installed.'
command -v iconv >/dev/null || die 'iconv is not installed (needed for the Windows workstation).'
command -v base64 >/dev/null || die 'base64 is not installed (needed for the Windows workstation).'

# These are Kyle SSH aliases.  They are intentionally explicit so Token cannot
# be changed accidentally.  Offline machines are reported and can be retried.
WORKSTATIONS=(jimmy timmy mrhankey randy stan butters kenny towelie)
KEY_FILE='.ssh/id_ed25519_pymacs_clusters'
START_MARKER='# >>> PyMACS cluster access >>>'
END_MARKER='# <<< PyMACS cluster access <<<'

hostname_for() {
  ssh -G "$1" 2>/dev/null | awk '$1 == "hostname" {print $2; exit}'
}

cluster_user="$(id -un)"
read -r -p "Cluster username [$cluster_user]: " supplied_user
cluster_user="${supplied_user:-$cluster_user}"
[[ "$cluster_user" =~ ^[A-Za-z0-9._-]+$ ]] || die 'Cluster username contains unsupported characters.'

acorn_host="$(hostname_for acorn)"
triton_host="$(hostname_for triton)"
pegasus_host="$(hostname_for pegasus)"
[[ -n "$acorn_host" && -n "$triton_host" && -n "$pegasus_host" ]] || die "Kyle's acorn, triton, and pegasus aliases must be configured first."

reachable=()
platforms=()
note 'Checking workstation connectivity from Kyle (Token is excluded)'
for workstation in "${WORKSTATIONS[@]}"; do
  if ssh -o BatchMode=yes -o ConnectTimeout=8 "$workstation" true 2>/dev/null; then
    platform="linux"
    if ssh -o BatchMode=yes -o ConnectTimeout=8 "$workstation" 'uname -s >/dev/null 2>&1' 2>/dev/null; then
      platform="linux"
    else
      platform="windows"
    fi
    printf '  %-10s reachable (%s)\n' "$workstation" "$platform"
    reachable+=("$workstation")
    platforms+=("$platform")
  else
    warn "$workstation is unreachable; it was skipped. Rerun this script later to configure it."
  fi
done
(( ${#reachable[@]} > 0 )) || die 'No eligible workstations are reachable from Kyle.'

run_windows_ps() {
  local workstation="$1" script="$2" encoded
  script="\$ProgressPreference = 'SilentlyContinue'\n$script"
  encoded="$(printf '%b' "$script" | iconv -f UTF-8 -t UTF-16LE | base64 | tr -d '\n')"
  ssh -o BatchMode=yes -o ConnectTimeout=15 "$workstation" "powershell.exe -NoProfile -NonInteractive -OutputFormat Text -EncodedCommand $encoded"
}

linux_public_key() {
  local workstation="$1"
  ssh -o BatchMode=yes -o ConnectTimeout=15 "$workstation" bash -s -- "$cluster_user" "$KEY_FILE" <<'REMOTE'
set -euo pipefail
cluster_user="$1"
key_file="$HOME/$2"
mkdir -p "$HOME/.ssh"
chmod 700 "$HOME/.ssh"
if [[ -e "$key_file" && ! -e "$key_file.pub" ]] || [[ ! -e "$key_file" && -e "$key_file.pub" ]]; then
  printf 'Only one half of %s exists; repair it before retrying.\n' "$key_file" >&2
  exit 1
fi
if [[ ! -e "$key_file" ]]; then
  ssh-keygen -q -t ed25519 -N '' -f "$key_file" -C "${cluster_user}@pymacs-$(hostname -s)"
fi
chmod 600 "$key_file"
cat "$key_file.pub"
REMOTE
}

windows_public_key() {
  local workstation="$1"
  run_windows_ps "$workstation" "\$ErrorActionPreference = 'Stop'
\$keyDir = Join-Path \$env:USERPROFILE '.ssh'
\$key = Join-Path \$keyDir 'id_ed25519_pymacs_clusters'
New-Item -ItemType Directory -Force -Path \$keyDir | Out-Null
if ((Test-Path \$key) -xor (Test-Path (\"\$key.pub\"))) { throw \"Only one half of \$key exists; repair it before retrying.\" }
if (-not (Test-Path \$key)) { & ssh-keygen.exe -q -t ed25519 -N '\"\"' -f \$key -C '${cluster_user}@pymacs-windows' }
(Get-Content -Raw \"\$key.pub\").Trim()"
}

public_keys=()
note 'Creating or reusing the dedicated key on each reachable workstation'
for index in "${!reachable[@]}"; do
  workstation="${reachable[$index]}"
  if [[ "${platforms[$index]}" == 'windows' ]]; then
    key="$(windows_public_key "$workstation")" || { warn "Could not create/read the key on $workstation; skipped."; continue; }
  else
    key="$(linux_public_key "$workstation")" || { warn "Could not create/read the key on $workstation; skipped."; continue; }
  fi
  [[ "$key" == ssh-ed25519\ * ]] || { warn "$workstation returned an invalid public key; skipped."; continue; }
  public_keys+=("$workstation|$key")
done
(( ${#public_keys[@]} > 0 )) || die 'No workstation public keys were collected.'

install_keys() {
  local label="$1" host="$2" payload
  local -a ssh_options=(
    -o PreferredAuthentications=keyboard-interactive,password
    -o PubkeyAuthentication=no
    -o NumberOfPasswordPrompts=3
    -o ConnectTimeout=20
  )
  payload="$(printf '%s\n' "${public_keys[@]#*|}")"
  note "Authorizing workstation keys on $label"
  printf '%s\n' "Enter the $label password when SSH asks. It is sent directly to $label and is never read or stored by this script."
  # Triton and Pegasus are reached from Kyle through Acorn, matching the
  # normal Kyle route and avoiding a dependency on direct cluster access.
  [[ "$label" == 'Acorn' ]] || ssh_options+=(-o ProxyJump=acorn-gw)
  ssh "${ssh_options[@]}" "${cluster_user}@${host}" \
    'umask 077; mkdir -p ~/.ssh; touch ~/.ssh/authorized_keys; chmod 700 ~/.ssh; chmod 600 ~/.ssh/authorized_keys; while IFS= read -r entry; do key=${entry#*|}; grep -Fqx "$key" ~/.ssh/authorized_keys || printf "%s\n" "$key" >> ~/.ssh/authorized_keys; done' \
    <<< "$payload"
}

# One login to each destination installs every reachable workstation key.
install_keys 'Acorn' "$acorn_host"
install_keys 'Triton' "$triton_host"
install_keys 'Pegasus' "$pegasus_host"

linux_configure() {
  local workstation="$1"
  ssh -o BatchMode=yes -o ConnectTimeout=15 "$workstation" bash -s -- \
    "$cluster_user" "$KEY_FILE" "$acorn_host" "$triton_host" "$pegasus_host" "$START_MARKER" "$END_MARKER" <<'REMOTE'
set -euo pipefail
user="$1"; key="$2"; acorn="$3"; triton="$4"; pegasus="$5"; start="$6"; end="$7"
config="$HOME/.ssh/config"
mkdir -p "$HOME/.ssh"; chmod 700 "$HOME/.ssh"; touch "$config"; chmod 600 "$config"
tmp="$(mktemp "${config}.pymacs.XXXXXX")"
awk -v start="$start" -v end="$end" '$0 == start {skip=1; next} $0 == end {skip=0; next} !skip {print}' "$config" > "$tmp"
cat >> "$tmp" <<EOF

$start
Host pymacs-acorn
    HostName $acorn
    User $user
    IdentityFile ~/$key
    IdentitiesOnly yes

Host pymacs-triton
    HostName $triton
    User $user
    IdentityFile ~/$key
    IdentitiesOnly yes
    ProxyJump pymacs-acorn

Host pymacs-pegasus
    HostName $pegasus
    User $user
    IdentityFile ~/$key
    IdentitiesOnly yes
    ProxyJump pymacs-acorn
$end
EOF
mv "$tmp" "$config"; chmod 600 "$config"
REMOTE
}

windows_configure() {
  local workstation="$1" ps
  ps="\$ErrorActionPreference = 'Stop'
\$dir = Join-Path \$env:USERPROFILE '.ssh'
\$path = Join-Path \$dir 'config'
New-Item -ItemType Directory -Force -Path \$dir | Out-Null
\$block = @'
$START_MARKER
Host pymacs-acorn
    HostName $acorn_host
    User $cluster_user
    IdentityFile ~/$KEY_FILE
    IdentitiesOnly yes

Host pymacs-triton
    HostName $triton_host
    User $cluster_user
    IdentityFile ~/$KEY_FILE
    IdentitiesOnly yes
    ProxyJump pymacs-acorn

Host pymacs-pegasus
    HostName $pegasus_host
    User $cluster_user
    IdentityFile ~/$KEY_FILE
    IdentitiesOnly yes
    ProxyJump pymacs-acorn
$END_MARKER
'@
\$old = if (Test-Path \$path) { Get-Content -Raw \$path } else { '' }
\$pattern = '(?ms)^# >>> PyMACS cluster access >>>\\r?\\n.*?^# <<< PyMACS cluster access <<<\\r?\\n?'
\$old = [regex]::Replace(\$old, \$pattern, '').TrimEnd()
[System.IO.File]::WriteAllText(\$path, \"\$old\`r\`n\`r\`n\$block\`r\`n\", [System.Text.UTF8Encoding]::new(\$false))"
  run_windows_ps "$workstation" "$ps"
}

note 'Writing passwordless Acorn, Triton, and Pegasus profiles on each workstation'
configured=()
for index in "${!reachable[@]}"; do
  workstation="${reachable[$index]}"
  if [[ ! " ${public_keys[*]} " == *"$workstation|"* ]]; then continue; fi
  if [[ "${platforms[$index]}" == 'windows' ]]; then
    windows_configure "$workstation" || { warn "Could not write profiles on $workstation."; continue; }
  else
    linux_configure "$workstation" || { warn "Could not write profiles on $workstation."; continue; }
  fi
  configured+=("$workstation")
done

note 'Verifying passwordless cluster access from each configured workstation'
failed=()
for index in "${!reachable[@]}"; do
  workstation="${reachable[$index]}"
  [[ " ${configured[*]} " == *" $workstation "* ]] || continue
  if [[ "${platforms[$index]}" == 'windows' ]]; then
    test_ps="& ssh.exe -o BatchMode=yes -o NumberOfPasswordPrompts=0 -o ConnectTimeout=15 pymacs-triton exit; if (\$LASTEXITCODE -ne 0) { exit \$LASTEXITCODE }; & ssh.exe -o BatchMode=yes -o NumberOfPasswordPrompts=0 -o ConnectTimeout=15 pymacs-pegasus exit"
    run_windows_ps "$workstation" "$test_ps" || failed+=("$workstation")
  else
    ssh -o BatchMode=yes -o ConnectTimeout=15 "$workstation" 'ssh -o BatchMode=yes -o NumberOfPasswordPrompts=0 -o ConnectTimeout=15 pymacs-triton true && ssh -o BatchMode=yes -o NumberOfPasswordPrompts=0 -o ConnectTimeout=15 pymacs-pegasus true' || failed+=("$workstation")
  fi
done

printf '\nConfigured: %s\n' "${configured[*]:-none}"
printf 'Skipped/offline: %s\n' "$(comm -23 <(printf '%s\n' "${WORKSTATIONS[@]}" | sort) <(printf '%s\n' "${configured[@]}" | sort) | tr '\n' ' ')"
if (( ${#failed[@]} )); then
  warn "Profiles were written, but verification failed on: ${failed[*]}"
  exit 1
fi
printf 'Success. On each configured machine, use: ssh pymacs-triton or ssh pymacs-pegasus\n'
