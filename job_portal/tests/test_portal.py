import json
import zipfile
from pathlib import Path

import pytest

import portal


def make_zip(path: Path, files: dict[str, str]) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        for name, text in files.items():
            archive.writestr(name, text)


def test_safe_extract_rejects_path_traversal(tmp_path):
    archive = tmp_path / "bad.zip"
    make_zip(archive, {"../outside.txt": "no"})
    with pytest.raises(portal.PortalError, match="unsafe path"):
        portal.safe_extract(archive, tmp_path / "out")


def test_collaboration_acl_is_limited_to_a_portal_owned_path():
    command = portal.collaborative_permissions_command(
        "/projectnb/cheminfo/share/pymacs-lab/submissions/joey/demo/campaign-001",
        "cheminfo",
    )
    assert "setfacl -R -m g:cheminfo:rwX,m::rwX" in command
    assert "d:g:cheminfo:rwx,d:m::rwx" in command
    assert "find '/projectnb/cheminfo/share/pymacs-lab/submissions/joey/demo/campaign-001' -type d" in command
    with pytest.raises(portal.PortalError, match="group name is invalid"):
        portal.collaborative_permissions_command("/project/path", "cheminfo; rm -rf /")


def test_docking_submission_uses_lsf_dependency(tmp_path, monkeypatch):
    root = tmp_path / "dock"
    root.mkdir()
    (root / "run_confgen_job.lsf").write_text("confgen")
    (root / "run_vina_job.lsf").write_text("vina")
    commands = []
    def fake_run(command, **kwargs):
        commands.append(command)
        return "Job <12> is submitted" if "run_confgen" in command[-1] else "Job <13> is submitted"
    monkeypatch.setattr(portal, "run", fake_run)
    monkeypatch.setattr(portal, "stage_project", lambda *args: "/remote/job")
    settings = portal.Settings(
        stage_root=tmp_path, pegasus_remote_root="/peg", triton_remote_root="/tri",
        pegasus_host="peg", triton_host="tri", pymacs_root="/pymacs",
        cgenff_url="http://localhost", cgenff_token_file=tmp_path / "token",
    )
    result = portal.submit_docking(settings, root, "test", "Joey")
    assert result == {"remote": "/remote/job", "confgen": "12", "vina": "13"}
    assert "done(12)" in commands[-1][-1]


def test_campaign_board_includes_four_legacy_ie1_screens(tmp_path, monkeypatch):
    settings = portal.Settings(
        stage_root=tmp_path, pegasus_remote_root="/peg", triton_remote_root="/tri",
        pegasus_host="peg", triton_host="tri", pymacs_root="/pymacs",
        cgenff_url="http://localhost", cgenff_token_file=tmp_path / "token",
    )
    received = []

    def fake_monitor(_settings, _host, campaigns):
        received.extend(campaigns)
        return {item["key"]: {"legacy": True, "batches": [], "done": 0, "total": 15,
                              "ligand_total": 1, "attempted": 0, "successful": 0,
                              "failed": 0, "progress": 0} for item in campaigns}

    monkeypatch.setattr(portal, "remote_monitor", fake_monitor)
    board = portal.monitored_submissions(settings)
    legacy = [item for item in board if str(item["project"]).startswith("IE1 ·")]
    assert len(legacy) == 4
    assert {item["project"] for item in legacy} == {name for name, _ in portal.LEGACY_IE1_CAMPAIGNS}
    assert all(item["legacy_ie1"] for item in received)


def test_md_project_requires_complete_runner(tmp_path):
    root = tmp_path / "md"
    root.mkdir()
    (root / "pymacs_run.json").write_text(json.dumps({}))
    (root / "input.pdb").write_text("ATOM")
    with pytest.raises(portal.PortalError, match="Python scripts"):
        portal.require_md_project(root)


def test_new_submitter_can_verify_ssh_but_cannot_self_approve(tmp_path, monkeypatch):
    settings = portal.Settings(
        stage_root=tmp_path, pegasus_remote_root="/peg", triton_remote_root="/tri",
        pegasus_host="pegasus", triton_host="triton", pymacs_root="/pymacs",
        cgenff_url="http://localhost", cgenff_token_file=tmp_path / "token",
    )
    (tmp_path / "eligible_accounts.json").write_text(json.dumps({"accounts": [
        {"name": "Maddison", "email": "mrm431@miami.edu"},
    ]}))

    def fake_run(command, **kwargs):
        host = command[-2] if command[-1] == "true" else command[-1]
        configs = {
            "triton": "hostname t2.idsc.miami.edu\nidentityfile ~/.ssh/SOUTHPARK",
            "pegasus": "hostname pegasus2.idsc.miami.edu\nidentityfile ~/.ssh/SOUTHPARK",
            "Maddison-Triton": "hostname t2.idsc.miami.edu\nidentityfile ~/.ssh/id_ed25519_maddison_triton",
            "Maddison-Pegasus": "hostname pegasus2.idsc.miami.edu\nidentityfile ~/.ssh/id_ed25519_maddison_triton",
        }
        return configs.get(host, f"hostname {host}\nidentityfile ~/.ssh/id_ed25519")

    monkeypatch.setattr(portal, "run", fake_run)
    app = portal.create_app(settings)
    client = app.test_client()

    checked = client.post("/api/submitters/validate", json={"name": "Maddison", "email": "mrm431@miami.edu"})
    assert checked.status_code == 200
    assert checked.get_json()["profiles"] == ["Triton", "Pegasus"]

    created = client.post("/api/submitters", json={"name": "Maddison", "email": "mrm431@miami.edu"})
    assert created.status_code == 403
    assert "administrator must approve" in created.get_json()["error"]

    denied = client.post("/api/submitters", json={"name": "Jordan", "email": "jordan@miami.edu"})
    assert denied.status_code == 403
    assert denied.get_json()["administrator_email"] == "jmschulz@med.miami.edu"


def test_new_submitter_can_bootstrap_cluster_ssh_from_browser(tmp_path, monkeypatch):
    settings = portal.Settings(
        stage_root=tmp_path, pegasus_remote_root="/peg", triton_remote_root="/tri",
        pegasus_host="pegasus", triton_host="triton", pymacs_root="/pymacs",
        cgenff_url="http://localhost", cgenff_token_file=tmp_path / "token",
    )
    (tmp_path / "eligible_accounts.json").write_text(json.dumps({"accounts": [
        {"name": "Angela", "email": "angela@miami.edu"},
    ]}))
    received = {}

    def fake_bootstrap(actual_settings, name, email, password):
        received.update(settings=actual_settings, name=name, email=email, password=password)
        return ["Triton", "Pegasus"]

    monkeypatch.setattr(portal, "bootstrap_cluster_profiles", fake_bootstrap)
    client = portal.create_app(settings).test_client()
    response = client.post("/api/submitters/bootstrap-ssh", json={
        "name": "Angela", "email": "angela@miami.edu", "cluster_password": "not-logged",
    })
    assert response.status_code == 200
    assert response.get_json()["profiles"] == ["Triton", "Pegasus"]
    assert received == {"settings": settings, "name": "Angela", "email": "angela@miami.edu", "password": "not-logged"}


def test_bootstrap_cluster_profiles_answers_two_private_ssh_prompts(tmp_path, monkeypatch):
    script = tmp_path / "setup.sh"
    script.write_text("""#!/usr/bin/env bash
set -euo pipefail
read -r -p 'First name: ' name
read -r -p 'Email: ' email
read -r -s -p 'Password: ' first; echo
read -r -s -p 'Password: ' second; echo
[[ "$name" == Angela && "$email" == angela@miami.edu && "$first" == secret && "$second" == secret ]]
""")
    settings = portal.Settings(
        stage_root=tmp_path, pegasus_remote_root="/peg", triton_remote_root="/tri",
        pegasus_host="pegasus", triton_host="triton", pymacs_root="/pymacs",
        cgenff_url="http://localhost", cgenff_token_file=tmp_path / "token", ssh_setup_script=script,
    )
    monkeypatch.setattr(portal, "require_configured_submitter", lambda *_: ["Triton", "Pegasus"])
    assert portal.bootstrap_cluster_profiles(settings, "Angela", "angela@miami.edu", "secret") == ["Triton", "Pegasus"]


def test_only_approved_password_protected_accounts_can_sign_in(tmp_path):
    settings = portal.Settings(
        stage_root=tmp_path, pegasus_remote_root="/peg", triton_remote_root="/tri",
        pegasus_host="pegasus", triton_host="triton", pymacs_root="/pymacs",
        cgenff_url="http://localhost", cgenff_token_file=tmp_path / "token",
    )
    (tmp_path / "submitters.json").write_text(json.dumps([
        {"name": "Joey", "email": "joey@example.edu", "approved": True,
         "password_hash": portal.password_hash("test-password")},
        {"name": "Mainak", "email": "mainak@example.edu", "approved": False},
    ]))
    client = portal.create_app(settings).test_client()
    assert client.get("/api/submitters").get_json()["names"] == [{"name": "Joey", "email": "joey@example.edu", "profiles": [], "approved": True}]
    assert client.post("/api/session", json={"name": "Joey", "password": "wrong"}).status_code == 403
    signed_in = client.post("/api/session", json={"name": "Joey", "password": "test-password"})
    assert signed_in.status_code == 200
    assert signed_in.get_json()["user"]["name"] == "Joey"


def test_published_result_route_serves_only_local_published_files(tmp_path):
    settings = portal.Settings(
        stage_root=tmp_path, pegasus_remote_root="/peg", triton_remote_root="/tri",
        pegasus_host="pegasus", triton_host="triton", pymacs_root="/pymacs",
        cgenff_url="http://localhost", cgenff_token_file=tmp_path / "token",
    )
    report = tmp_path / "published-md-results" / "external" / "test" / "report.html"
    report.parent.mkdir(parents=True)
    report.write_text("<h1>Result</h1>")
    client = portal.create_app(settings).test_client()
    response = client.get("/results/files/external/test/report.html?inline=1")
    assert response.status_code == 200
    assert response.mimetype == "text/html"
    assert client.get("/results/files/../../etc/passwd").status_code == 404
