#!/usr/bin/env bash
# Guided Pegasus-home cleanup for a user whose SSH key cannot be installed
# because their home quota is full. Run from an interactive Kyle terminal.
set -euo pipefail

die() { printf 'Error: %s\n' "$*" >&2; exit 1; }
[[ -t 0 && -t 1 ]] || die "Run this from an interactive Kyle terminal."

read -r -p "Miami email (for example abc123@miami.edu): " miami_email
[[ "$miami_email" =~ ^([A-Za-z0-9._-]+)@miami\.edu$ ]] || die "Enter a valid @miami.edu address."
username="${BASH_REMATCH[1]}"
pegasus_host="$(ssh -G pegasus 2>/dev/null | awk '$1 == "hostname" {print $2; exit}')"
[[ -n "$pegasus_host" ]] || die "Kyle's reviewed Pegasus host is not configured."
remote=(ssh -o ProxyJump=acorn-gw -o PreferredAuthentications=keyboard-interactive,password -o PubkeyAuthentication=no "${username}@${pegasus_host}")

printf '\n==> Checking Pegasus storage for %s\n' "$username"
printf 'Enter the university password when SSH asks. It is not read or saved by this script.\n\n'
"${remote[@]}" 'df -h "$HOME"; quota -s 2>&1 || true; printf "\nLargest home-directory entries:\n"; du -xhd1 "$HOME" 2>/dev/null | sort -h | tail -25; printf "\nLargest files (up to 30):\n"; find "$HOME" -xdev -type f -printf "%s\t%p\n" 2>/dev/null | sort -n | tail -30 | numfmt --field=1 --to=iec-i --suffix=B 2>/dev/null || true'

printf '\nNothing has been deleted. Review the paths above carefully.\n'
read -r -p 'Relative path under your Pegasus home to permanently delete (blank to exit): ' relative
[[ -n "$relative" ]] || { printf 'No files deleted.\n'; exit 0; }
[[ "$relative" =~ ^[A-Za-z0-9._/-]+$ ]] || die "Use a simple relative path copied from the report."
[[ "$relative" != /* && "$relative" != .* && "$relative" != *".."* && "$relative" != .ssh* ]] || die "Refusing paths outside your home or anything in .ssh."

printf '\n==> Candidate for permanent deletion: ~/%s\n' "$relative"
"${remote[@]}" "target=\"\$HOME/$relative\"; test -e \"\$target\" || { echo 'Path does not exist.'; exit 2; }; du -sh -- \"\$target\""
read -r -p 'Type DELETE to permanently remove exactly this path, or press Enter to cancel: ' confirmation
[[ "$confirmation" == "DELETE" ]] || { printf 'Cancelled; no files deleted.\n'; exit 0; }

"${remote[@]}" "target=\"\$HOME/$relative\"; rm -rf -- \"\$target\"; df -h \"\$HOME\"; quota -s 2>&1 || true"
printf '\nDeleted ~/%s. You can now rerun: bash ~/setup-cluster-ssh.sh\n' "$relative"
