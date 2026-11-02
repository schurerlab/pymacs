#!/usr/bin/env bash
# Guided passwordless-SSH setup for the PyMACS portal on Kyle.
# Run from an interactive Kyle shell. Password prompts are handled directly by
# SSH; this script never reads, stores, or prints a password.
set -euo pipefail

die() { printf 'Error: %s\n' "$*" >&2; exit 1; }
note() { printf '\n==> %s\n' "$*"; }

[[ -t 0 && -t 1 ]] || die "Run this from an interactive terminal on Kyle."
command -v ssh >/dev/null || die "ssh is not installed."
command -v ssh-keygen >/dev/null || die "ssh-keygen is not installed."

read -r -p "First name for the portal (for example Mainak): " first_name
[[ "$first_name" =~ ^[A-Za-z][A-Za-z-]{0,39}$ ]] || die "Use letters or hyphens only for the first name."
read -r -p "Miami email (for example abc123@miami.edu): " miami_email
[[ "$miami_email" =~ ^([A-Za-z0-9._-]+)@miami\.edu$ ]] || die "Enter a valid @miami.edu address."

username="${BASH_REMATCH[1]}"
name_slug="$(printf '%s' "$first_name" | tr '[:upper:]' '[:lower:]' | tr -cd 'a-z-')"
triton_alias="${first_name}-Triton"
pegasus_alias="${first_name}-Pegasus"
key_file="$HOME/.ssh/id_ed25519_${name_slug}_clusters"
config_file="$HOME/.ssh/config"

mkdir -p "$HOME/.ssh"
chmod 700 "$HOME/.ssh"
touch "$config_file"
chmod 600 "$config_file"

# Use the reviewed hostnames already configured on Kyle rather than guessing
# which campus DNS alias is currently canonical.
triton_host="$(ssh -G triton 2>/dev/null | awk '$1 == "hostname" {print $2; exit}')"
pegasus_host="$(ssh -G pegasus 2>/dev/null | awk '$1 == "hostname" {print $2; exit}')"
[[ -n "$triton_host" && -n "$pegasus_host" ]] || die "Kyle's reviewed Triton/Pegasus SSH hosts are not configured."

if [[ -e "$key_file" && -e "$key_file.pub" ]]; then
  note "Reusing the dedicated portal SSH key from an earlier setup attempt"
elif [[ -e "$key_file" || -e "$key_file.pub" ]]; then
  die "Only part of $key_file exists. Ask the administrator to review it before retrying."
else
  note "Creating a dedicated portal SSH key on Kyle"
  # The portal needs to authenticate unattended after setup, so this key has no
  # passphrase. It is protected by Kyle's filesystem permissions.
  ssh-keygen -q -t ed25519 -N '' -f "$key_file" -C "${username}@pymacs-kyle"
  chmod 600 "$key_file"
fi

install_key() {
  local label="$1" host="$2"
  note "Installing the public key on $label"
  printf '%s\n' "Enter your $label cluster password when SSH asks. It is sent directly to $label and is not handled by this script."
  # Force the initial login to use password or keyboard-interactive auth so a
  # stale local key cannot make the setup appear to succeed for the wrong user.
  cat "$key_file.pub" | ssh \
    -o ProxyJump=acorn-gw \
    -o PreferredAuthentications=keyboard-interactive,password \
    -o PubkeyAuthentication=no \
    "${username}@${host}" \
    'umask 077; mkdir -p ~/.ssh; touch ~/.ssh/authorized_keys; chmod 700 ~/.ssh; chmod 600 ~/.ssh/authorized_keys; cat >> ~/.ssh/authorized_keys'
}

cleanup_key_on_failure() {
  printf '\nSetup did not finish. The new key remains on Kyle but was not added to both aliases.\n' >&2
  printf 'If Pegasus reported “No space left on device”, run: bash ~/cleanup-pegasus-space.sh\n' >&2
  printf 'Then rerun this setup script with the same name and email.\n' >&2
}
trap cleanup_key_on_failure ERR

key_works() {
  local host="$1"
  ssh -o ProxyJump=acorn-gw -o BatchMode=yes -o NumberOfPasswordPrompts=0 \
    -o ConnectTimeout=15 -o IdentitiesOnly=yes -i "$key_file" "${username}@${host}" true >/dev/null 2>&1
}

if key_works "$triton_host"; then
  note "Triton is already configured; skipping its password prompt"
else
  install_key "Triton" "$triton_host"
fi
if key_works "$pegasus_host"; then
  note "Pegasus is already configured; skipping its password prompt"
else
  install_key "Pegasus" "$pegasus_host"
fi

note "Adding the two named SSH profiles"
add_profile() {
  local alias="$1" host="$2"
  if grep -Fqx "Host $alias" "$config_file"; then
    return
  fi
  {
    printf '\nHost %s\n' "$alias"
    printf '    HostName %s\n    User %s\n    IdentityFile %s\n    IdentitiesOnly yes\n    ProxyJump acorn-gw\n' "$host" "$username" "$key_file"
  } >> "$config_file"
}
add_profile "$triton_alias" "$triton_host"
add_profile "$pegasus_alias" "$pegasus_host"

note "Verifying passwordless access"
ssh -o BatchMode=yes -o NumberOfPasswordPrompts=0 -o ConnectTimeout=15 "$triton_alias" true
ssh -o BatchMode=yes -o NumberOfPasswordPrompts=0 -o ConnectTimeout=15 "$pegasus_alias" true

note "Configuring the Cheminformatics shared PyMACS environments"
for cluster_alias in "$triton_alias" "$pegasus_alias"; do
  ssh -o BatchMode=yes -o NumberOfPasswordPrompts=0 -o ConnectTimeout=30 "$cluster_alias" \
    'bash /projectnb/cheminfo/share/pymacs-lab/lab-tools/bootstrap-cluster-envs.sh'
done
trap - ERR

printf '\nSuccess. %s and %s are ready for the PyMACS portal.\n' "$triton_alias" "$pegasus_alias"
printf 'Your shared lab environments are now available through pymacs-env after you open a new cluster shell.\n'
printf 'Return to the portal and choose “I need an account” → “Check compatibility.”\n'
