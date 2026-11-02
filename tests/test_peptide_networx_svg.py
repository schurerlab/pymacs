import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "4D_Peptide_NETWORX.py"
SPEC = importlib.util.spec_from_file_location("peptide_networx", SCRIPT)
peptide_networx = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(peptide_networx)


def test_peptide_node_letters_and_sequence_numbers_have_dark_outlines():
    svg = peptide_networx.svg_scene(peptide_networx.model_from_data(peptide_networx.mock_data(8)))

    assert 'paint-order="stroke fill"' in svg
    assert 'stroke="#07111e" stroke-width="1.8"' in svg
    assert 'stroke="#07111e" stroke-width="1.25"' in svg
