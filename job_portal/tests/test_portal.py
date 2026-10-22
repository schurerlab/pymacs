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
    result = portal.submit_docking(portal.Settings(tmp_path, "/peg", "/tri", "peg", "tri", "/pymacs"), root, "test", "Joey")
    assert result == {"remote": "/remote/job", "confgen": "12", "vina": "13"}
    assert "done(12)" in commands[-1][-1]


def test_md_project_requires_complete_runner(tmp_path):
    root = tmp_path / "md"
    root.mkdir()
    (root / "pymacs_run.json").write_text(json.dumps({}))
    (root / "input.pdb").write_text("ATOM")
    with pytest.raises(portal.PortalError, match="Python scripts"):
        portal.require_md_project(root)
