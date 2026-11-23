"""Keep the README installers aligned with the default analysis runtime."""

from pathlib import Path


README = Path(__file__).resolve().parents[1] / "README.md"

REQUIRED_RUNTIME_FILES = (
    "pymacs_run.py",
    "3A_AutomateGromacs_chunks_MPI.py",
    "pymacs_compute_native_contact_q.py",
    "pymacs_native_contacts.py",
    "pymacs_component_utils.py",
    "pymacs_analysis_cache.py",
    "pymacs_interaction_chemistry.py",
)


def _installer_blocks(readme: str) -> list[str]:
    """Return the two bash installer blocks in the README."""
    return readme.split("```bash\n")[1:3]


def test_quick_installers_copy_and_verify_the_default_runtime() -> None:
    blocks = _installer_blocks(README.read_text())

    assert len(blocks) == 2
    for block in blocks:
        assert 'cp -R "$item" "$TARGET_DIR/"' in block
        for filename in REQUIRED_RUNTIME_FILES:
            assert filename in block


def test_required_runtime_files_exist_in_the_repository() -> None:
    repo_root = README.parent
    for filename in REQUIRED_RUNTIME_FILES:
        assert (repo_root / filename).is_file()
