from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_every_3a_runner_dispatches_the_shared_q_suite():
    for name in (
        "3A_AutomateGromacs.py",
        "3A_AutomateGromacs_MPI.py",
        "3A_AutomateGromacs_chunks.py",
        "3A_AutomateGromacs_chunks_MPI.py",
    ):
        source = (ROOT / name).read_text(encoding="utf-8")
        assert "def run_native_contact_q_suite" in source
        assert '"pymacs_compute_native_contact_q.py"' in source
        assert "run_native_contact_q_suite()" in source


def test_shared_suite_publishes_all_supported_q_contracts():
    source = (ROOT / "pymacs_compute_native_contact_q.py").read_text(encoding="utf-8")
    for filename in (
        "protein_native_contact_q.csv",
        "native_contact_q.csv",
        "peptide_native_contact_q.csv",
        "protein_interface_native_contact_q.csv",
        "native_contact_q_metrics.csv",
    ):
        assert filename in source
