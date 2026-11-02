#!/usr/bin/env python3
"""Administrator-only portal account approval and password changes.

Run this directly on the portal host to approve a working SSH profile:

    python manage_accounts.py approve "Maddison" --email mrm431@miami.edu
"""

import argparse
import datetime as dt
import getpass
import json
from pathlib import Path

from portal import DEFAULT_PORTAL_PASSWORD, Settings, configured_profiles, password_hash, settings_from_env, submitter_slug


def load_accounts(path: Path) -> list[dict]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return []
    if not isinstance(data, list):
        raise SystemExit(f"{path} does not contain an account list.")
    return [item if isinstance(item, dict) else {"name": str(item)} for item in data]


def approve(settings: Settings, name: str, email: str) -> None:
    # Validate syntax and prove the SSH setup the same way the portal does.
    submitter_slug(name)
    profiles = configured_profiles(settings, name)
    if not profiles:
        raise SystemExit(f"{name} has no working passwordless SSH profile; account was not changed.")
    path = settings.stage_root / "submitters.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    accounts = load_accounts(path)
    record = next((item for item in accounts if str(item.get("name", "")).casefold() == name.casefold()), None)
    if record is None:
        record = {"name": name.strip(), "created_at": dt.datetime.now(dt.timezone.utc).isoformat()}
        accounts.append(record)
    record.update({"email": email.strip(), "profiles": profiles, "approved": True,
                   "password_hash": password_hash(DEFAULT_PORTAL_PASSWORD)})
    accounts.sort(key=lambda item: str(item.get("name", "")).casefold())
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(accounts, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)
    print(f"Approved {record['name']} for {', '.join(profiles)} with the default portal password.")


def set_password(settings: Settings, name: str) -> None:
    submitter_slug(name)
    first = getpass.getpass(f"New portal password for {name}: ")
    second = getpass.getpass("Confirm portal password: ")
    if not first:
        raise SystemExit("Password may not be blank.")
    if first != second:
        raise SystemExit("Passwords did not match; account was not changed.")
    path = settings.stage_root / "submitters.json"
    accounts = load_accounts(path)
    record = next((item for item in accounts if str(item.get("name", "")).casefold() == name.casefold()), None)
    if record is None or not record.get("approved"):
        raise SystemExit(f"{name} is not an approved portal account.")
    record["password_hash"] = password_hash(first)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(accounts, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)
    print(f"Updated {record['name']}'s portal password.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Approve a PyMACS portal account with working SSH.")
    subcommands = parser.add_subparsers(dest="command", required=True)
    command = subcommands.add_parser("approve", help="verify SSH and approve a user with the default password")
    command.add_argument("name", help="first name matching Name-Triton / Name-Pegasus")
    command.add_argument("--email", required=True)
    reset = subcommands.add_parser("set-password", help="set an approved user's portal password")
    reset.add_argument("name", help="approved portal account name")
    args = parser.parse_args()
    if args.command == "approve":
        approve(settings_from_env(), args.name, args.email)
    elif args.command == "set-password":
        set_password(settings_from_env(), args.name)


if __name__ == "__main__":
    main()
