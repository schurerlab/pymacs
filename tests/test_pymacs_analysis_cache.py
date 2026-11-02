from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from pymacs_analysis_cache import (
    compute_signature,
    fingerprint_files,
    load_npz_dict,
    record_stage_completion,
    save_npz_atomic,
    validate_stage_cache,
)


class AnalysisCacheTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = tempfile.TemporaryDirectory()
        self.base = Path(self.tmpdir.name)
        self.inputs = {
            "topology": self._write("topology.pdb", "MODEL\nENDMDL\n"),
            "trajectory": self._write("traj.xtc", "trajectory-bytes"),
            "atomindex": self._write("atomIndex.txt", "A 1 - 10\nB 11 - 20\n"),
        }
        self.output = self._write("analysis/required.csv", "x,y\n1,2\n")
        self.checkpoint = self._write("chk_D1.txt", "done")
        self.manifest = self.base / "analysis" / "analysis_manifest.json"

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def _write(self, relative: str, content: str) -> Path:
        path = self.base / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def _signature(self, cutoff: float = 4.0, frame_step: int = 1, chain_pair: str = "A:B") -> str:
        return compute_signature(
            {
                "inputs": fingerprint_files(self.inputs),
                "parameters": {
                    "contact_cutoff_A": cutoff,
                    "frame_step": frame_step,
                    "selected_chain_pairs": [chain_pair],
                },
            }
        )

    def test_cache_valid_for_same_data_affecting_inputs(self) -> None:
        data_params = {
            "contact_cutoff_A": 4.0,
            "frame_step": 1,
            "selected_chain_pairs": ["A:B"],
        }
        record_stage_completion(
            self.manifest,
            stage_name="INTERFACE_RIN",
            analysis_script="3C_Interface_RIN.py",
            analysis_source_path=self.inputs["topology"],
            input_files=self.inputs,
            data_parameters=data_params,
            plot_parameters={"max_edges": 100},
            outputs=[self.output],
        )
        status = validate_stage_cache(
            manifest_path=self.manifest,
            stage_name="INTERFACE_RIN",
            required_outputs=[self.output],
            current_data_signature=self._signature(),
            checkpoint_path=self.checkpoint,
        )
        self.assertTrue(status["reusable"])
        self.assertFalse(status["legacy_adopt"])

    def test_contact_cutoff_change_invalidates_numerical_cache(self) -> None:
        record_stage_completion(
            self.manifest,
            stage_name="INTERFACE_RIN",
            analysis_script="3C_Interface_RIN.py",
            analysis_source_path=self.inputs["topology"],
            input_files=self.inputs,
            data_parameters={"contact_cutoff_A": 4.0, "frame_step": 1, "selected_chain_pairs": ["A:B"]},
            plot_parameters={"max_edges": 100},
            outputs=[self.output],
        )
        status = validate_stage_cache(
            manifest_path=self.manifest,
            stage_name="INTERFACE_RIN",
            required_outputs=[self.output],
            current_data_signature=self._signature(cutoff=5.0),
            checkpoint_path=self.checkpoint,
        )
        self.assertFalse(status["reusable"])
        self.assertEqual(status["reason"], "data_signature_mismatch")

    def test_frame_step_change_invalidates_numerical_cache(self) -> None:
        record_stage_completion(
            self.manifest,
            stage_name="INTERFACE_RIN",
            analysis_script="3C_Interface_RIN.py",
            analysis_source_path=self.inputs["topology"],
            input_files=self.inputs,
            data_parameters={"contact_cutoff_A": 4.0, "frame_step": 1, "selected_chain_pairs": ["A:B"]},
            plot_parameters={"max_edges": 100},
            outputs=[self.output],
        )
        status = validate_stage_cache(
            manifest_path=self.manifest,
            stage_name="INTERFACE_RIN",
            required_outputs=[self.output],
            current_data_signature=self._signature(frame_step=5),
            checkpoint_path=self.checkpoint,
        )
        self.assertFalse(status["reusable"])
        self.assertEqual(status["reason"], "data_signature_mismatch")

    def test_plot_only_change_does_not_invalidate_numerical_cache(self) -> None:
        record_stage_completion(
            self.manifest,
            stage_name="INTERFACE_RIN",
            analysis_script="3C_Interface_RIN.py",
            analysis_source_path=self.inputs["topology"],
            input_files=self.inputs,
            data_parameters={"contact_cutoff_A": 4.0, "frame_step": 1, "selected_chain_pairs": ["A:B"]},
            plot_parameters={"max_edges": 100},
            outputs=[self.output],
        )
        status = validate_stage_cache(
            manifest_path=self.manifest,
            stage_name="INTERFACE_RIN",
            required_outputs=[self.output],
            current_data_signature=self._signature(),
            checkpoint_path=self.checkpoint,
        )
        self.assertTrue(status["reusable"])

    def test_checkpoint_plus_missing_output_is_not_valid(self) -> None:
        self.output.unlink()
        status = validate_stage_cache(
            manifest_path=self.manifest,
            stage_name="D1",
            required_outputs=[self.output],
            current_data_signature=self._signature(),
            checkpoint_path=self.checkpoint,
        )
        self.assertFalse(status["reusable"])
        self.assertEqual(status["reason"], "required_outputs_missing")

    def test_checkpoint_plus_required_outputs_is_reusable_as_legacy(self) -> None:
        status = validate_stage_cache(
            manifest_path=self.manifest,
            stage_name="D1",
            required_outputs=[self.output],
            current_data_signature=self._signature(),
            checkpoint_path=self.checkpoint,
        )
        self.assertTrue(status["reusable"])
        self.assertTrue(status["legacy_adopt"])

    def test_legacy_interface_missing_temporal_dataset_requires_backfill(self) -> None:
        old_summary = self._write("analysis/interface_chain_pair_summary.csv", "a,b\n1,2\n")
        old_times = self._write("analysis/interface_contact_timeseries.csv", "a,b\n1,2\n")
        old_pairs = self._write("analysis/interface_residue_pair_contacts.csv", "a,b\n1,2\n")
        old_nodes = self._write("analysis/interface_node_summary.csv", "a,b\n1,2\n")
        missing_new = self.base / "analysis" / "interface_residue_pair_presence.npz"
        status = validate_stage_cache(
            manifest_path=self.manifest,
            stage_name="INTERFACE_RIN",
            required_outputs=[old_summary, old_times, old_pairs, old_nodes, missing_new],
            current_data_signature=self._signature(),
            checkpoint_path=self.checkpoint,
        )
        self.assertFalse(status["reusable"])
        self.assertEqual(status["reason"], "required_outputs_missing")

    def test_npz_round_trip(self) -> None:
        npz_path = self.base / "presence.npz"
        payload = {
            "contact_matrix": np.array([[1, 0, 1], [0, 1, 1]], dtype=np.uint8),
            "frame_indices": np.array([0, 5, 10], dtype=np.int64),
            "times_ns": np.array([0.0, 1.5, 3.0], dtype=np.float64),
            "pair_ids": np.array(["A__B__A:GLY1__B:ASP2", "A__B__A:LYS3__B:GLU4"]),
            "chain_a_labels": np.array(["A", "A"]),
            "chain_b_labels": np.array(["B", "B"]),
            "residue_a_labels": np.array(["A:GLY1", "A:LYS3"]),
            "residue_b_labels": np.array(["B:ASP2", "B:GLU4"]),
        }
        save_npz_atomic(npz_path, **payload)
        loaded = load_npz_dict(npz_path)
        for key, value in payload.items():
            self.assertTrue(np.array_equal(loaded[key], value), msg=f"Mismatch for {key}")

    def test_legacy_object_arrays_require_explicit_trusted_cache_mode(self) -> None:
        npz_path = self.base / "legacy-presence.npz"
        save_npz_atomic(npz_path, residue_a_labels=np.array(["A:GLY1"], dtype=object))
        with self.assertRaises(ValueError):
            load_npz_dict(npz_path)
        loaded = load_npz_dict(npz_path, allow_legacy_object_arrays=True)
        self.assertEqual(loaded["residue_a_labels"].tolist(), ["A:GLY1"])


if __name__ == "__main__":
    unittest.main()
