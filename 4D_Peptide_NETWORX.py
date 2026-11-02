#!/usr/bin/env python3
"""Render a sequence-preserving peptide--protein interaction explorer.

This is a presentation-only successor to 3C_Interface_RIN.py.  It consumes
the self-contained DATA payload written by 3C and never recalculates contacts.
Run directly for repeatable layout review::

    python 4D_Peptide_NETWORX.py --mock --all-lengths
    python 4D_Peptide_NETWORX.py --interface-dir Analysis_Results/Interface_RIN
"""
from __future__ import annotations

import argparse
import html
import json
import math
import random
import re
from pathlib import Path


AA1 = "ACDEFGHIKLMNPQRSTVWY"
COLOURS = {
    "Hydrogen bond": "#ff5f85", "Ionic": "#bc8cff", "Hydrophobic": "#59e39a",
    "Water bridge": "#48c8ff", "Other contact": "#8ca3c0",
}


def mock_data(length: int) -> dict:
    """Deterministic, explicitly synthetic data used only for display review."""
    randomizer = random.Random(400 + length)
    peptide, receptor = "Mock peptide (peptide, chain B)", "Mock receptor (protein, chain A)"
    residues = [
        {"label": f"{peptide}:{AA1[(index * 7 + 3) % 20]}{index + 1}", "one": AA1[(index * 7 + 3) % 20],
         "number": index + 1, "name": "MOCK"}
        for index in range(length)
    ]
    protein_names = ("ARG", "ASP", "TYR", "LYS", "SER", "GLU", "TRP", "ASN", "PHE", "HIS", "THR", "LEU")
    protein = [
        {"label": f"{receptor}:{protein_names[index % len(protein_names)]}{88 + index * 11}",
         "one": "X", "number": 88 + index * 11, "name": protein_names[index % len(protein_names)]}
        for index in range(min(length + 3, 15))
    ]
    edges = []
    for index, residue in enumerate(residues):
        if index % 4 == 0 or index % 7 == 3:
            continue
        for partner_index in range(1 + int(index % 6 == 0)):
            edges.append({"a": residue["label"], "b": protein[(index * 3 + partner_index * 5) % len(protein)]["label"],
                          "type": ("Hydrogen bond", "Ionic", "Hydrophobic", "Other contact")[(index + partner_index) % 4],
                          "fraction": round(randomizer.uniform(.16, .98), 3)})
    return {"pairs": [{"id": "mock", "left": receptor, "right": peptide}],
            "chains": {peptide: {"label": peptide, "type": "peptide", "residues": residues},
                       receptor: {"label": receptor, "type": "protein", "residues": protein}},
            "edges": edges, "mock": True}


def load_data(interface_dir: Path) -> dict:
    source = interface_dir / "Interface_RIN_Explorer.html"
    contents = source.read_text(encoding="utf-8")
    match = re.search(r"\bconst\s+DATA\s*=\s*", contents)
    if not match:
        raise ValueError(f"{source} does not contain 3C's embedded DATA payload")
    data, _ = json.JSONDecoder().raw_decode(contents[match.end():])
    if not isinstance(data, dict) or not {"chains", "edges"}.issubset(data):
        raise ValueError("3C explorer DATA payload is incomplete")
    return data


def model_from_data(data: dict, requested: str = "") -> dict:
    chains = data["chains"]
    peptide_candidates = [(key, chain) for key, chain in chains.items()
                          if str(chain.get("type", "")).casefold() == "peptide"]
    if requested:
        peptide_candidates = [(key, chain) for key, chain in peptide_candidates
                              if requested == key or requested in key]
    if not peptide_candidates:
        raise ValueError("No explicitly typed peptide sequence was supplied by 3C")
    peptide_name, peptide = peptide_candidates[0]
    residues = list(peptide.get("residues") or [])
    if not residues:
        raise ValueError(f"No complete sequence was supplied for {peptide_name}")
    residue_ids = {str(item.get("label", "")) for item in residues}
    other_residues = {str(item.get("label", "")) for name, chain in chains.items() if name != peptide_name
                      for item in chain.get("residues", [])}
    edges = []
    for edge in data.get("edges", []):
        first, second = str(edge.get("a", "")), str(edge.get("b", ""))
        if first in residue_ids and second in other_residues:
            peptide_id, protein_id = first, second
        elif second in residue_ids and first in other_residues:
            peptide_id, protein_id = second, first
        else:
            continue
        try:
            occupancy = float(edge.get("fraction"))
        except (TypeError, ValueError):
            continue
        if math.isfinite(occupancy) and 0 <= occupancy <= 1:
            edges.append({"peptide": peptide_id, "protein": protein_id,
                          "kind": str(edge.get("type") or "Other contact"), "occupancy": occupancy})
    return {"name": peptide_name, "residues": residues, "edges": edges,
            "is_mock": bool(data.get("mock"))}


def svg_scene(model: dict, threshold: float = .25) -> str:
    """Produce a clean schematic: sequence geometry is fixed; contacts vary."""
    residues = model["residues"]
    count = len(residues)
    columns = count if count <= 8 else 7 if count <= 15 else 8 if count <= 20 else 9
    rows = math.ceil(count / columns)
    # Give the sequence and its satellites generous breathing room.  The
    # explorer is embedded in the dossier, but remains horizontally scrollable
    # on narrower screens rather than compressing labels into each other.
    width, height = max(1560, 400 + (columns - 1) * 170), 570 + (rows - 1) * 330
    step_x = (width - 320) / max(1, columns - 1)
    peptide_xy = []
    for index in range(count):
        row, column = divmod(index, columns)
        x = 160 + (column if row % 2 == 0 else columns - 1 - column) * step_x
        peptide_xy.append((x, 285 + row * 330))
    index_by_id = {str(item.get("label")): index for index, item in enumerate(residues)}
    visible = [edge for edge in model["edges"] if edge["occupancy"] >= threshold]
    partners: dict[str, list[dict]] = {}
    for edge in visible:
        partners.setdefault(edge["protein"], []).append(edge)
    # Place satellites in outer lanes nearest their sequence neighborhood.  The
    # alternating lanes reserve the peptide's central reading path and reduce
    # crossings in dense 20--25 residue diagrams.
    slots = [(x, y) for y in ([145] + [445 + 330 * row for row in range(rows)])
             for x in range(120, width - 40, 174)]
    protein_xy: dict[str, tuple[int, int]] = {}
    for partner, partner_edges in sorted(partners.items(), key=lambda item: (
            sum(index_by_id[edge["peptide"]] for edge in item[1]) / len(item[1]), item[0])):
        target_x = sum(peptide_xy[index_by_id[edge["peptide"]]][0] for edge in partner_edges) / len(partner_edges)
        target_y = sum(peptide_xy[index_by_id[edge["peptide"]]][1] for edge in partner_edges) / len(partner_edges)
        available = [slot for slot in slots if all(math.dist(slot, placed) >= 128 for placed in protein_xy.values())
                     and all(math.dist(slot, point) >= 82 for point in peptide_xy)] or slots
        protein_xy[partner] = min(available, key=lambda slot: (slot[0] - target_x) ** 2 + .6 * (slot[1] - target_y) ** 2)
    escape = lambda value: html.escape(str(value), quote=True)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Peptide interaction topology">',
             '<rect width="100%" height="100%" fill="#07111e"/>',
             '<style>text{font-family:Inter,Arial,sans-serif}.node:hover{filter:brightness(1.22)}</style>',
             f'<text x="42" y="42" fill="#69dfff" font-size="22" font-weight="800">PEPTIDE NETWORX · {count} AA</text>',
             '<text x="42" y="74" fill="#96aeca" font-size="15">N → C sequence topology · schematic, not 3D conformation</text>']
    if model["is_mock"]:
        parts.append(f'<text x="{width - 28}" y="44" text-anchor="end" fill="#ffd26b" font-size="18" font-weight="900">SYNTHETIC DEMO</text>')
    for index in range(1, len(peptide_xy)):
        ax, ay = peptide_xy[index - 1]
        bx, by = peptide_xy[index]
        path = f'M {ax} {ay} L {bx} {by}' if ay == by else f'M {ax} {ay} C {ax + (72 if bx > ax else -72)} {ay} {bx + (72 if bx > ax else -72)} {by} {bx} {by}'
        parts.append(f'<path d="{path}" fill="none" stroke="#627995" stroke-width="9" stroke-linecap="round"/>')
    for edge in visible:
        x1, y1 = peptide_xy[index_by_id[edge["peptide"]]]
        x2, y2 = protein_xy[edge["protein"]]
        colour = next((value for name, value in COLOURS.items() if name.casefold() in edge["kind"].casefold()), COLOURS["Other contact"])
        parts.append(f'<path d="M {x1} {y1} Q {(x1 + x2) / 2:.1f} {(y1 + y2) / 2:.1f} {x2} {y2}" fill="none" stroke="{colour}" stroke-width="{2.5 + edge["occupancy"] * 3:.1f}" stroke-dasharray="8 6" opacity=".9"><title>{escape(edge["peptide"])} ↔ {escape(edge["protein"])} · {escape(edge["kind"])} · {edge["occupancy"] * 100:.1f}%</title></path>')
    for partner, (x, y) in protein_xy.items():
        label = escape(partner.split(":")[-1])
        residue_name = label[:3].upper()
        parts.append(f'<g class="node protein" data-aa="{escape(residue_name)}"><rect x="{x - 64}" y="{y - 24}" width="128" height="48" rx="22" fill="#153e54" stroke="#52d7f6" stroke-width="2.5"/><text x="{x}" y="{y + 6}" text-anchor="middle" fill="#c0eeff" font-size="16" font-weight="800">{label}</text></g>')
    active = {edge["peptide"] for edge in visible}
    for index, (residue, (x, y)) in enumerate(zip(residues, peptide_xy), 1):
        active_node = str(residue.get("label")) in active
        fill, stroke, text = ("#ffdb80", "#ffe8ad", "#1d2430") if active_node else ("#25394f", "#7695b5", "#edf5ff")
        residue_name = str(residue.get("label", "")).split(":")[-1][:3].upper()
        # Chemistry mode recolours the text.  Draw its dark stroke first so
        # every label remains readable against the bright residue fills.
        parts.append(f'<g class="node peptide" data-aa="{escape(residue_name)}"><circle cx="{x}" cy="{y}" r="40" fill="{fill}" stroke="{stroke}" stroke-width="3.5"/><text x="{x}" y="{y + 2}" text-anchor="middle" fill="{text}" stroke="#07111e" stroke-width="1.8" stroke-linejoin="round" paint-order="stroke fill" font-size="27" font-weight="900">{escape(residue.get("one", "X"))}</text><text x="{x}" y="{y + 26}" text-anchor="middle" fill="{text}" stroke="#07111e" stroke-width="1.25" stroke-linejoin="round" paint-order="stroke fill" font-size="12" font-weight="800">{index}</text><title>{escape(residue.get("label", ""))} · sequence position {index}</title></g>')
    first, last = peptide_xy[0], peptide_xy[-1]
    parts.extend((f'<text x="{first[0] - 78}" y="{first[1] + 8}" fill="#5ce5fa" font-size="22" font-weight="900">N</text>',
                  f'<text x="{last[0] + 62}" y="{last[1] + 8}" fill="#5ce5fa" font-size="22" font-weight="900">C</text>', "</svg>"))
    return "\n".join(parts)


def write_outputs(model: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    scenes = {str(step): svg_scene(model, step / 100) for step in range(0, 101, 5)}
    (output / "Peptide_NETWORX.svg").write_text(scenes["25"], encoding="utf-8")
    payload = json.dumps(model, separators=(",", ":")).replace("</", "<\\/")
    scenes_payload = json.dumps(scenes).replace("</", "<\\/")
    html_document = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Peptide NETWORX</title>
<style>:root{{--bg:#050911;--panel:#0e1828;--text:#edf7ff;--muted:#9cb6d0;--cyan:#4bd7ff}}*{{box-sizing:border-box}}body{{margin:0;background:radial-gradient(ellipse at top right,#102c45,var(--bg) 58%);font-family:Inter,system-ui,sans-serif;color:var(--text)}}main{{max-width:1800px;margin:auto;padding:24px}}h1{{margin:5px 0;font-size:clamp(28px,4vw,48px);letter-spacing:-.05em}}p{{color:var(--muted);line-height:1.5;max-width:980px}}.eyebrow{{font-size:12px;letter-spacing:.18em;color:var(--cyan);font-weight:900}}.toolbar{{display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin:18px 0}}input{{accent-color:var(--cyan)}}button{{border:1px solid #36536d;border-radius:9px;background:#132840;color:var(--text);padding:9px 13px;font-weight:750;cursor:pointer}}button.active{{border-color:#55ddff;background:#15516a;box-shadow:0 0 0 1px #55ddff inset}}.panel{{background:#0d1725;border:1px solid #253a54;border-radius:18px;overflow:auto}}.diagram svg{{display:block;width:100%;min-width:1320px;height:auto}}.notice{{color:#ffd26b;font-weight:800}}.hint{{font-size:.88rem;margin:.2rem 0 1rem}}</style>
<main><div class="eyebrow">PYMACS · PEPTIDE NETWORX</div><h1>Peptide interaction topology</h1><p>Every peptide residue stays in N → C sequence order. Gold residues pass the selected occupancy threshold. This is a 2D sequence schematic, not a structural conformation.</p><p class="notice" id="notice"></p><div class="toolbar"><label>Minimum occupancy <b id="value">25%</b> <input id="threshold" type="range" min="0" max="100" value="25" step="5"></label><button data-colour="contact" class="active">Contact status</button><button data-colour="peptide">Color peptide by chemistry</button><button data-colour="protein">Color receptor by chemistry</button><button id="export">Export SVG</button></div><p class="hint">This sequence layout is fixed: occupancy changes only contact edges and highlighting. Chemistry colors change only a residue's outline and label; its original fill remains visible.</p><section class="panel diagram" id="diagram"></section></main>
<script>const DATA={payload},SCENES={scenes_payload};const threshold=document.querySelector('#threshold'),diagram=document.querySelector('#diagram');let colourMode='contact';const chemistry={{ASP:'#ff6b7a',GLU:'#ff6b7a',LYS:'#c58aff',ARG:'#c58aff',HIS:'#c58aff',ALA:'#6be79d',VAL:'#6be79d',LEU:'#6be79d',ILE:'#6be79d',MET:'#6be79d',PHE:'#6be79d',TRP:'#6be79d',TYR:'#6be79d',PRO:'#6be79d',SER:'#66dfff',THR:'#66dfff',ASN:'#66dfff',GLN:'#66dfff',CYS:'#66dfff'}};function applyColours(){{diagram.querySelectorAll('.node').forEach(node=>{{const apply=colourMode!=='contact'&&node.classList.contains(colourMode),colour=chemistry[node.dataset.aa]||'#dbe9f8';node.querySelectorAll('circle,rect').forEach(shape=>shape.setAttribute('stroke',apply?colour:shape.dataset.stroke));node.querySelectorAll('text').forEach(text=>text.setAttribute('fill',apply?colour:text.dataset.fill))}})}}function draw(){{document.querySelector('#value').textContent=threshold.value+'%';diagram.innerHTML=SCENES[threshold.value]||SCENES['25'];diagram.querySelectorAll('.node').forEach(node=>{{node.querySelectorAll('circle,rect').forEach(shape=>shape.dataset.stroke=shape.getAttribute('stroke'));node.querySelectorAll('text').forEach(text=>text.dataset.fill=text.getAttribute('fill'))}});applyColours()}}threshold.oninput=draw;document.querySelectorAll('[data-colour]').forEach(button=>button.onclick=()=>{{colourMode=button.dataset.colour;document.querySelectorAll('[data-colour]').forEach(item=>item.classList.toggle('active',item===button));applyColours()}});document.querySelector('#notice').textContent=DATA.is_mock?'SYNTHETIC DEMONSTRATION — contacts and occupancies are mock data.':'';document.querySelector('#export').onclick=()=>{{const blob=new Blob([diagram.innerHTML],{{type:'image/svg+xml'}}),link=document.createElement('a');link.href=URL.createObjectURL(blob);link.download='Peptide_NETWORX.svg';link.click();setTimeout(()=>URL.revokeObjectURL(link.href),500)}};draw();</script>'''
    (output / "Peptide_NETWORX_Explorer.html").write_text(html_document, encoding="utf-8")
    manifest = {"schema_version": 1, "renderer": Path(__file__).name,
                "mode": "mock" if model["is_mock"] else "3C_interface", "chain": model["name"],
                "residue_count": len(model["residues"]), "edge_count": len(model["edges"]),
                "generated": ["Peptide_NETWORX.svg", "Peptide_NETWORX_Explorer.html"],
                "notes": "2D sequence-topology schematic; noncontacting residues are always retained."}
    (output / "peptide_networx_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mock", action="store_true")
    parser.add_argument("--mock-length", type=int, default=25)
    parser.add_argument("--all-lengths", action="store_true")
    parser.add_argument("--interface-dir", type=Path, default=Path("Analysis_Results/Interface_RIN"))
    parser.add_argument("--peptide-chain", default="")
    parser.add_argument("--outdir", type=Path, default=Path("Analysis_Results/Peptide_NETWORX"))
    args = parser.parse_args()
    lengths = (5, 10, 15, 20, 25) if args.mock and args.all_lengths else (args.mock_length,)
    if args.mock:
        for length in lengths:
            if not 5 <= length <= 25:
                parser.error("Mock peptide length must be in the 5--25 range")
            write_outputs(model_from_data(mock_data(length)), args.outdir / f"Mock_{length}mer" if args.all_lengths else args.outdir)
    else:
        write_outputs(model_from_data(load_data(args.interface_dir), args.peptide_chain), args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
