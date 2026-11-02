import json
import re
import unittest
from pathlib import Path

from pymacs_interface_explorer import build_cytoscape_interface_explorer, build_interface_explorer


class InterfaceMultiviewExplorerTests(unittest.TestCase):
    def test_mock_peptide_lengths_keep_all_residues_and_views(self):
        for length in (5, 10, 15, 20, 25):
            peptide = [{"label": f"Pep:ALA{i}", "stable_id": f"B|{i}", "one": "A", "number": i,
                        "position": [i, 0, 0]} for i in range(1, length + 1)]
            receptor = [{"label": f"Rec:LYS{i}", "stable_id": f"A|{i}", "one": "K", "number": i,
                         "position": [i, 1, 0]} for i in range(1, length + 1)]
            data = {
                "pairs": [{"id": "p", "left": "Rec", "right": "Pep"}],
                "chains": {"Rec": {"type": "protein", "residues": receptor},
                           "Pep": {"type": "peptide", "residues": peptide}},
                "edges": [{"a": peptide[0]["label"], "b": receptor[0]["label"],
                           "fraction": 0.5, "type": "Hydrogen bond"}],
                "default_threshold": 0.25, "presence": {"times_ns": [0, 1], "rows": {}},
            }
            page = build_interface_explorer(json.dumps(data))
            for view in ("Network", "3D Structure", "3D Nodes", "Contact Map", "Peptide Sequence"):
                self.assertIn(view, page)
            self.assertIn("Time window", page)
            self.assertIn("stable_id", page)
            self.assertIn(f"Pep:ALA{length}", page)

    def test_cytoscape_is_vendored_and_layouts_are_embedded(self):
        page = build_cytoscape_interface_explorer(json.dumps({"pairs": [], "chains": {}, "edges": {}, "presence": {}}))
        self.assertIn("Cytoscape.js 3.34.4", page)
        self.assertIn(".cytoscape=t()", page)
        for layout in ("Force-directed", "Bipartite", "Concentric", "Circular", "Chain-separated"):
            self.assertIn(layout, page)

    def test_selection_bridge_orientation_and_geometric_window_contract(self):
        peptide = [{"label": "Pep:TYR7", "stable_id": "B|7", "name": "TYR", "one": "Y", "number": 7,
                    "position": [0, 0, 0]}]
        receptor = [{"label": "Rec:LYS7", "stable_id": "A|7", "name": "LYS", "one": "K", "number": 7,
                     "position": [1, 0, 0]}]
        page = build_cytoscape_interface_explorer(json.dumps({
            "pairs": [{"id": "reversed", "left": "Pep", "right": "Rec"}],
            "chains": {"Pep": {"type": "peptide", "residues": peptide},
                       "Rec": {"type": "protein", "residues": receptor}},
            "edges": [{"a": "Pep:TYR7", "b": "Rec:LYS7", "fraction": .5, "type": "Hydrogen bond"}],
            "default_threshold": .25, "presence": {"times_ns": [0, 1], "rows": {}},
        }))
        self.assertIn("window.PyMACSSelection", page)
        self.assertIn("pymacs-select", page)
        self.assertIn("left.type==='peptide'?left:right", page)
        self.assertIn("e.window_geometric=$('#window').value!=='all'", page)
        self.assertIn("!e.window_geometric&&e.type===state.kind", page)
        self.assertNotIn("/static/3Dmol-min.js", page)

    def test_embedded_views_are_contact_focused_and_contained(self):
        page = build_cytoscape_interface_explorer(json.dumps({
            "pairs": [{"id": "p", "left": "Rec", "right": "Pep",
                       "left_chain": "A", "right_chain": "B"}],
            "chains": {"Rec": {"type": "protein", "residues": []},
                       "Pep": {"type": "peptide", "residues": []}},
            "edges": [], "presence": {"available": True, "times_ns": [0, 1], "rows": {}}
        }))
        self.assertIn("visible=new Set(es.flatMap", page)
        self.assertIn("No contacts meet the active filter", page)
        self.assertIn(".map-scroll{max-width:100%;max-height:535px;overflow:auto", page)
        self.assertIn("+'</table></div>';$('#map')", page)
        self.assertIn("centre=bounds.reduce", page)
        self.assertIn("pepChain=p.left===pep?.label?p.left_chain:p.right_chain", page)
        self.assertIn("DATA.presence?.available===true", page)

    def test_structure_selector_and_network_resize_contracts(self):
        page = build_cytoscape_interface_explorer(json.dumps({
            "pairs": [{"id": "p", "left": "Rec", "right": "Pep",
                       "left_chain": "A", "right_chain": "B"}],
            "chains": {"Rec": {"type": "protein", "residues": []},
                       "Pep": {"type": "peptide", "residues": []}},
            "edges": [], "default_threshold": .5,
            "presence": {"available": True, "times_ns": [0, 1], "rows": {}},
        }))
        self.assertIn('class="structure-grid"', page)
        self.assertIn('id="mol"', page)
        self.assertIn('id="structure-detail"', page)
        self.assertIn("function structureInspector()", page)
        self.assertIn("viewer.resize()", page)
        self.assertIn("defaultOcc=()=>", page)
        self.assertIn("layout='bipartite'", page)
        self.assertIn("label:''", page)
        self.assertIn("if(window.ResizeObserver)", page)

    def test_sandbox_network_fallback_and_selection_colours(self):
        page = build_cytoscape_interface_explorer(json.dumps({
            "pairs": [{"id": "p", "left": "Rec", "right": "Pep"}],
            "chains": {"Rec": {"type": "protein", "residues": []},
                       "Pep": {"type": "peptide", "residues": []}},
            "edges": [], "presence": {"available": True, "times_ns": [], "rows": {}},
        }))
        self.assertIn("#cy{display:none!important}#rin{display:block", page)
        self.assertIn("state.selected===id){state.selected=null", page)
        self.assertIn("#24e6ff", page)
        self.assertIn("#ff5a78", page)
        self.assertIn("$('#space').onwheel", page)
        self.assertIn("document.querySelector('#cy').offsetWidth===0", page)

    def test_reusable_3c_cache_rebuilds_the_explorer(self):
        source = Path(__file__).resolve().parents[1].joinpath("3C_Interface_RIN.py").read_text()
        match = re.search(
            r"if cached_plot_parameters == plot_parameters and plots_ready:(.*?)progress\.complete\(\"contact_scan\"\)",
            source,
            re.DOTALL,
        )
        self.assertIsNotNone(match)
        branch = match.group(1)
        self.assertIn("write_interactive_interface_rin_html(", branch)
        self.assertIn("cached_events", branch)


if __name__ == "__main__":
    unittest.main()
