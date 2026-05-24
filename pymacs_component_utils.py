import json
import math
import os
from collections import OrderedDict, defaultdict

try:
    import numpy as _np
except Exception:  # pragma: no cover - optional at runtime
    _np = None


COMPONENT_REGISTRY_FILENAME = "pymacs_components.json"
SYSTEM_REGISTRY_FILENAME = "pymacs_system_registry.json"
WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P", "SPC", "SPCE"}
ORDINARY_IONS = {
    "NA", "CL", "K", "LI", "CS", "RB", "F", "BR", "I", "IOD",
}
STRUCTURAL_METALS = {"ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI", "CD", "HG"}
COMMON_SOLVENTS_AND_BUFFERS = {
    "SO4", "PO4", "PI", "ACE", "DMS", "MES", "TRS", "EDO", "GOL", "PEG", "MPD",
}
STANDARD_PROTEIN_AND_NA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL", "ASH", "GLH", "HID", "HIE", "HIP", "LYN", "CYX",
    "DA", "DT", "DG", "DC", "DI", "A", "T", "G", "C", "U",
}
PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYX", "GLN", "GLU", "GLH", "GLY",
    "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "LYN",
    "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL",
}
DNA_ONLY_RESNAMES = {"DA", "DC", "DG", "DT", "DI", "THY"}
RNA_ONLY_RESNAMES = {"U", "URA", "URI", "RU", "RA", "RC", "RG", "RU5", "RA5", "RC5", "RG5"}
AMBIGUOUS_NUCLEIC_RESNAMES = {"A", "C", "G", "I", "ADE", "CYT", "GUA", "URI", "URA"}
RNA_MARKER_ATOMS = {"O2'", "O2*"}


def parse_component_list(value):
    if not value:
        return []
    items = []
    for item in value.split(","):
        cleaned = item.strip().upper()
        if cleaned and cleaned not in items:
            items.append(cleaned)
    return items


def ordered_unique_resnames(items):
    ordered = []
    for item in items:
        cleaned = (item or "").strip().upper()
        if cleaned and cleaned not in ordered:
            ordered.append(cleaned)
    return ordered


def get_all_retained_components(primary_ligand, cofactors=None):
    items = []
    if primary_ligand:
        items.append(primary_ligand)
    items.extend(cofactors or [])
    return ordered_unique_resnames(items)


def _record_resname(line):
    return line[17:20].strip().upper()


def _record_chain_id(line):
    return line[21].strip() or "_"


def _record_resid(line):
    return line[22:26].strip()


def _record_icode(line):
    return line[26].strip() or "_"


def _record_atom_name(line):
    return line[12:16].strip()


def _record_element(line):
    if len(line) >= 78:
        element = line[76:78].strip().upper()
        if element:
            return element
    atom_name = _record_atom_name(line).strip()
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if not letters:
        return ""
    return letters[:2].upper() if len(letters) > 1 and letters[:2].upper() in {"CL", "BR", "NA", "MG", "ZN", "FE", "MN", "CA", "CU"} else letters[:1].upper()


def residue_instance_key(line):
    return (
        _record_resname(line),
        _record_chain_id(line),
        _record_resid(line),
        _record_icode(line),
    )


def polymer_type_label(chain_type):
    normalized = (chain_type or "").strip().lower()
    return {
        "protein": "protein",
        "rna": "RNA",
        "dna": "DNA",
        "mixed_or_unknown": "mixed_or_unknown",
    }.get(normalized, normalized or "mixed_or_unknown")


def parse_chain_type_overrides(value):
    overrides = {}
    if not value:
        return overrides
    for item in value.split(","):
        token = item.strip()
        if not token:
            continue
        if ":" not in token:
            raise ValueError(f"Invalid chain type override '{token}'. Expected CHAIN:TYPE.")
        chain_id, chain_type = token.split(":", 1)
        chain_id = chain_id.strip() or "_"
        normalized = chain_type.strip().lower().replace("-", "_").replace(" ", "_")
        normalized = {
            "rna": "rna",
            "dna": "dna",
            "protein": "protein",
            "unknown": "mixed_or_unknown",
            "mixed": "mixed_or_unknown",
            "mixed_or_unknown": "mixed_or_unknown",
            "nucleic_acid_unknown": "mixed_or_unknown",
        }.get(normalized)
        if not normalized:
            raise ValueError(f"Unsupported chain type override '{chain_type}' for chain {chain_id}.")
        overrides[chain_id] = normalized
    return overrides


def _classify_polymer_residue(resname, atom_names):
    resname = (resname or "").strip().upper()
    atom_names = {name.strip().upper().replace('"', "'") for name in (atom_names or set()) if name}
    if resname in PROTEIN_RESNAMES:
        return "protein"
    if resname in DNA_ONLY_RESNAMES:
        return "dna"
    if resname in RNA_ONLY_RESNAMES:
        return "rna"
    if resname == "T":
        return "dna"
    if resname in AMBIGUOUS_NUCLEIC_RESNAMES:
        if atom_names & RNA_MARKER_ATOMS:
            return "rna"
        return "nucleic_unknown"
    return "mixed_or_unknown"


def _polymer_system_flags(chain_summaries):
    chain_types = [entry.get("chain_type", "mixed_or_unknown") for entry in (chain_summaries or [])]
    has_protein = any(chain_type == "protein" for chain_type in chain_types)
    has_rna = any(chain_type == "rna" for chain_type in chain_types)
    has_dna = any(chain_type == "dna" for chain_type in chain_types)
    has_nucleic_acid = has_rna or has_dna
    if has_protein and has_rna and has_dna:
        system_class = "protein_nucleic"
    elif has_protein and has_rna:
        system_class = "protein_rna"
    elif has_protein and has_dna:
        system_class = "protein_dna"
    elif has_rna and has_dna:
        system_class = "nucleic_acid"
    elif has_rna:
        system_class = "rna"
    elif has_dna:
        system_class = "dna"
    elif has_protein:
        system_class = "protein_only"
    else:
        system_class = "mixed_or_unknown"
    return {
        "has_protein": has_protein,
        "has_rna": has_rna,
        "has_dna": has_dna,
        "has_nucleic_acid": has_nucleic_acid,
        "system_class": system_class,
    }


def classify_polymer_chains(pdb_path):
    chains = OrderedDict()
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line.startswith("HETATM") and _record_resname(line) not in STANDARD_PROTEIN_AND_NA:
                continue
            chain_id = _record_chain_id(line)
            resid = _record_resid(line)
            icode = _record_icode(line)
            resname = _record_resname(line)
            atom_name = _record_atom_name(line).upper().replace('"', "'")
            chain_entry = chains.setdefault(
                chain_id,
                {
                    "chain_id": chain_id,
                    "residue_order": [],
                    "residues": OrderedDict(),
                },
            )
            residue_key = (resid, icode, resname)
            residue_entry = chain_entry["residues"].setdefault(
                residue_key,
                {
                    "resid": resid,
                    "insertion_code": icode,
                    "resname": resname,
                    "atom_names": set(),
                },
            )
            if residue_key not in chain_entry["residue_order"]:
                chain_entry["residue_order"].append(residue_key)
            residue_entry["atom_names"].add(atom_name)

    summaries = []
    warnings = []
    for chain_id, chain_entry in chains.items():
        residue_types = []
        residues = [chain_entry["residues"][key] for key in chain_entry["residue_order"]]
        for residue in residues:
            residue_type = _classify_polymer_residue(residue["resname"], residue["atom_names"])
            residue_types.append(residue_type)
        if "protein" in residue_types and any(kind in {"rna", "dna", "nucleic_unknown"} for kind in residue_types):
            chain_type = "mixed_or_unknown"
            warnings.append(
                f"Chain {chain_id} mixes protein and nucleic-acid residue signatures; treating it as mixed_or_unknown."
            )
        else:
            normalized_types = {kind for kind in residue_types if kind not in {"mixed_or_unknown", "nucleic_unknown"}}
            has_unknown_nucleic = "nucleic_unknown" in residue_types
            if not normalized_types:
                chain_type = "mixed_or_unknown"
            elif normalized_types == {"protein"}:
                chain_type = "protein"
            elif normalized_types == {"rna"}:
                chain_type = "rna"
                if has_unknown_nucleic:
                    warnings.append(
                        f"Chain {chain_id} contains ambiguous nucleic-acid residues; classified as RNA from chain context."
                    )
            elif normalized_types == {"dna"}:
                chain_type = "dna"
                if has_unknown_nucleic:
                    warnings.append(
                        f"Chain {chain_id} contains ambiguous nucleic-acid residues; classified as DNA from chain context."
                    )
            else:
                chain_type = "mixed_or_unknown"
                warnings.append(
                    f"Chain {chain_id} contains mixed or unrecognized polymer residue types; treating it as mixed_or_unknown."
                )
        first_residue = residues[0] if residues else {}
        last_residue = residues[-1] if residues else {}
        first_atom_names = set(first_residue.get("atom_names") or set())
        last_atom_names = set(last_residue.get("atom_names") or set())
        has_5prime_phosphate = bool(
            first_atom_names & {"P"}
            and first_atom_names & {"OP1", "OP2", "O1P", "O2P"}
        )
        summary = {
            "chain_id": chain_id,
            "chain_type": chain_type,
            "residue_count": len(residues),
            "first_resname": first_residue.get("resname"),
            "first_resid": first_residue.get("resid"),
            "last_resname": last_residue.get("resname"),
            "last_resid": last_residue.get("resid"),
            "first_atom_names": sorted(first_atom_names),
            "last_atom_names": sorted(last_atom_names),
            "has_5prime_phosphate": has_5prime_phosphate,
        }
        summaries.append(summary)

    flags = _polymer_system_flags(summaries)
    return {
        "chains": summaries,
        "chain_types": {entry["chain_id"]: entry["chain_type"] for entry in summaries},
        "warnings": warnings,
        **flags,
    }


def apply_chain_type_overrides(polymer_info, overrides):
    if not overrides:
        return polymer_info
    updated = dict(polymer_info)
    chains = []
    for entry in polymer_info.get("chains", []):
        revised = dict(entry)
        if entry["chain_id"] in overrides:
            revised["chain_type"] = overrides[entry["chain_id"]]
        chains.append(revised)
    updated["chains"] = chains
    updated["chain_types"] = {entry["chain_id"]: entry["chain_type"] for entry in chains}
    updated.update(_polymer_system_flags(chains))
    return updated


def format_polymer_chain_report(polymer_info):
    lines = []
    for entry in polymer_info.get("chains", []):
        lines.append(
            f"  chain {entry['chain_id']}: {polymer_type_label(entry['chain_type'])} | "
            f"residues={entry['residue_count']} | "
            f"first={entry.get('first_resname')}{entry.get('first_resid')} | "
            f"last={entry.get('last_resname')}{entry.get('last_resid')}"
        )
    for warning in polymer_info.get("warnings", []):
        lines.append(f"  warning: {warning}")
    return "\n".join(lines)


def polymer_chain_type_list(polymer_info):
    ordered = []
    for entry in polymer_info.get("chains", []):
        label = polymer_type_label(entry.get("chain_type"))
        if label not in ordered:
            ordered.append(label)
    return ordered


def build_polymer_group_name(system_registry=None, components=None, override="auto"):
    components = ordered_unique_resnames(components or [])
    override = (override or "auto").strip()
    if override != "auto":
        base_group = override
    else:
        registry = system_registry or {}
        has_protein = bool(registry.get("has_protein"))
        has_rna = bool(registry.get("has_rna"))
        has_dna = bool(registry.get("has_dna"))
        if has_protein and has_rna and has_dna:
            base_group = "Protein_Nucleic"
        elif has_protein and has_rna:
            base_group = "Protein_RNA"
        elif has_protein and has_dna:
            base_group = "Protein_DNA"
        elif has_rna and has_dna:
            base_group = "Nucleic_Acid"
        elif has_rna:
            base_group = "RNA"
        elif has_dna:
            base_group = "DNA"
        else:
            base_group = "Protein"
    if components and base_group != "non-Water":
        return base_group + "_" + "_".join(components)
    return base_group


def build_polymer_selection_terms(system_registry=None):
    registry = system_registry or {}
    terms = []
    if registry.get("has_protein"):
        terms.append('group "Protein"')
    if registry.get("has_rna"):
        terms.append('group "RNA"')
    if registry.get("has_dna"):
        terms.append('group "DNA"')
    return terms


def write_system_registry(directory, polymer_info):
    polymer_info = polymer_info or {}
    registry = {
        "system_class": polymer_info.get("system_class", "mixed_or_unknown"),
        "polymer_chains": {
            chain_id: polymer_type_label(chain_type)
            for chain_id, chain_type in (polymer_info.get("chain_types") or {}).items()
        },
        "polymer_chain_summaries": [
            {
                "chain_id": entry.get("chain_id"),
                "chain_type": polymer_type_label(entry.get("chain_type")),
                "residue_count": entry.get("residue_count"),
                "first_resname": entry.get("first_resname"),
                "first_resid": entry.get("first_resid"),
                "last_resname": entry.get("last_resname"),
                "last_resid": entry.get("last_resid"),
                "has_5prime_phosphate": bool(entry.get("has_5prime_phosphate")),
                "chosen_start_terminus": entry.get("chosen_start_terminus"),
                "chosen_end_terminus": entry.get("chosen_end_terminus"),
            }
            for entry in (polymer_info.get("chains") or [])
        ],
        "has_protein": bool(polymer_info.get("has_protein")),
        "has_nucleic_acid": bool(polymer_info.get("has_nucleic_acid")),
        "has_rna": bool(polymer_info.get("has_rna")),
        "has_dna": bool(polymer_info.get("has_dna")),
        "warnings": list(polymer_info.get("warnings") or []),
    }
    path = os.path.join(directory, SYSTEM_REGISTRY_FILENAME)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(registry, handle, indent=2, sort_keys=False)
        handle.write("\n")
    return path, registry


def load_system_registry(directory="."):
    path = os.path.join(directory, SYSTEM_REGISTRY_FILENAME)
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    chain_types = {}
    for chain_id, chain_type in (data.get("polymer_chains") or {}).items():
        normalized = str(chain_type).strip().lower()
        chain_types[chain_id] = {
            "rna": "rna",
            "dna": "dna",
            "protein": "protein",
        }.get(normalized, "mixed_or_unknown")
    chains = []
    for entry in data.get("polymer_chain_summaries") or []:
        chains.append(
            {
                "chain_id": entry.get("chain_id"),
                "chain_type": chain_types.get(entry.get("chain_id"), "mixed_or_unknown"),
                "residue_count": entry.get("residue_count"),
                "first_resname": entry.get("first_resname"),
                "first_resid": entry.get("first_resid"),
                "last_resname": entry.get("last_resname"),
                "last_resid": entry.get("last_resid"),
                "has_5prime_phosphate": bool(entry.get("has_5prime_phosphate")),
                "chosen_start_terminus": entry.get("chosen_start_terminus"),
                "chosen_end_terminus": entry.get("chosen_end_terminus"),
            }
        )
    normalized = {
        "chains": chains,
        "chain_types": chain_types,
        "warnings": list(data.get("warnings") or []),
    }
    normalized.update(_polymer_system_flags(chains))
    return normalized


def categorize_component(resname):
    resname = (resname or "").strip().upper()
    if resname in WATER_RESNAMES:
        return "water", "ordinary_ion", True
    if resname in ORDINARY_IONS:
        return "ordinary_ion", "ordinary_ion", True
    if resname in COMMON_SOLVENTS_AND_BUFFERS:
        return "common_solvent_or_buffer", "likely_solvent_or_buffer", True
    if resname in STRUCTURAL_METALS:
        return "structural_metal", "candidate_cofactor", False
    if resname in STANDARD_PROTEIN_AND_NA:
        return "standard_biopolymer", "standard_biopolymer", True
    return "ligand_like", "candidate_primary_ligand", False


def summarize_detected_hetero_components(pdb_path):
    summary = OrderedDict()
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("HETATM"):
                continue
            resname = _record_resname(line)
            category, recommendation, excluded = categorize_component(resname)
            serial_key = residue_instance_key(line)
            entry = summary.setdefault(
                resname,
                {
                    "resname": resname,
                    "atom_count_raw_pdb": 0,
                    "copies": 0,
                    "instances": OrderedDict(),
                    "category": category,
                    "recommended_role": recommendation,
                    "excluded_by_default": excluded,
                    "warning": None,
                },
            )
            entry["atom_count_raw_pdb"] += 1
            if serial_key not in entry["instances"]:
                entry["instances"][serial_key] = {
                    "chain_id": serial_key[1],
                    "resid": serial_key[2],
                    "insertion_code": serial_key[3],
                    "atom_count": 0,
                }
            entry["instances"][serial_key]["atom_count"] += 1

    for resname, entry in summary.items():
        entry["copies"] = len(entry["instances"])
        entry["instance_labels"] = [
            f"{data['chain_id']}:{data['resid']}{'' if data['insertion_code'] == '_' else data['insertion_code']}"
            for data in entry["instances"].values()
        ]
        if entry["category"] == "structural_metal":
            entry["warning"] = (
                f"{resname} looks like a structural/catalytic metal; topology and coordination may require manual validation."
            )
        elif entry["copies"] > 1 and entry["atom_count_raw_pdb"] / max(entry["copies"], 1) >= 30:
            entry["recommended_role"] = "candidate_cofactor"
    return summary


def selectable_component_resnames(component_summary):
    ordered = []
    for resname, entry in component_summary.items():
        if entry.get("excluded_by_default"):
            continue
        ordered.append(resname)
    return ordered


def render_detected_component_table(component_summary):
    total_species = len(component_summary)
    total_instances = sum(int(entry.get("copies", 0)) for entry in component_summary.values())
    lines = [
        f"Detected non-protein component species: {total_species} unique residue name(s), {total_instances} total instance(s)."
    ]
    selectable = selectable_component_resnames(component_summary)
    index_lookup = {}
    row_index = 1
    for resname, entry in component_summary.items():
        locations = ",".join(entry.get("instance_labels", [])) or "n/a"
        status = "excluded" if entry.get("excluded_by_default") else entry.get("recommended_role", "candidate")
        prefix = f"[{row_index}] " if resname in selectable else "[-] "
        lines.append(
            f"  {prefix}{resname} | instances: {entry.get('copies', 0)} | atoms total: {entry.get('atom_count_raw_pdb', 0)} "
            f"| chains/resids: {locations} | suggested: {status}"
        )
        for instance_key, instance_data in entry.get("instances", {}).items():
            insertion_code = instance_data["insertion_code"]
            icode_suffix = "" if insertion_code == "_" else insertion_code
            lines.append(
                f"      - instance {instance_data['chain_id']}:{instance_data['resid']}{icode_suffix} "
                f"| atoms: {instance_data['atom_count']}"
            )
        if entry.get("warning"):
            lines.append(f"      warning: {entry['warning']}")
        if resname in selectable:
            index_lookup[str(row_index)] = resname
            row_index += 1
    return "\n".join(lines), index_lookup


def resolve_component_selection(raw_value, index_lookup, allowed_resnames):
    value = (raw_value or "").strip()
    if not value:
        return []
    if value.lower() == "none":
        return []
    if value.lower() == "all":
        return list(allowed_resnames)

    chosen = []
    allowed = {code.upper() for code in allowed_resnames}
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        token_upper = token.upper()
        if token in index_lookup:
            code = index_lookup[token]
        elif token_upper in allowed:
            code = token_upper
        else:
            raise ValueError(f"Unknown component selection token: {token}")
        if code not in chosen:
            chosen.append(code)
    return chosen


def write_protein_pdb_without_components(input_pdb, output_pdb, remove_resnames):
    remove_set = set(ordered_unique_resnames(remove_resnames))
    removed_counts = OrderedDict((resname, 0) for resname in remove_set)
    written_lines = []

    with open(input_pdb, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                resname = _record_resname(line)
                if resname in remove_set:
                    removed_counts[resname] = removed_counts.get(resname, 0) + 1
                    continue
                written_lines.append(line)
            elif line.startswith("CONECT"):
                continue
            elif line.startswith(("TER", "END", "ENDMDL")):
                written_lines.append(line)
            else:
                written_lines.append(line)

    if not written_lines or not any(line.startswith(("ATOM", "HETATM")) for line in written_lines):
        raise RuntimeError(f"No coordinate records remained after removing {sorted(remove_set)} from {input_pdb}")

    with open(output_pdb, "w", encoding="utf-8") as handle:
        for line in written_lines:
            handle.write(line)
        if not written_lines[-1].startswith("END"):
            handle.write("END\n")

    return dict(removed_counts)


def count_component_instances(pdb_path, resnames):
    targets = set(ordered_unique_resnames(resnames))
    counts = OrderedDict((resname, 0) for resname in ordered_unique_resnames(resnames))
    seen = {resname: set() for resname in counts}
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = _record_resname(line)
            if resname not in targets:
                continue
            residue_key = residue_instance_key(line)
            if residue_key not in seen[resname]:
                seen[resname].add(residue_key)
                counts[resname] += 1
    return dict(counts)


def list_component_instances(pdb_path, component_code):
    component_code = component_code.strip().upper()
    residues = parse_pdb_atoms_by_residue(pdb_path, component_code)

    def _residue_sort_key(item):
        residue_key, _atoms = item
        _resname, chain_id, resid, icode = residue_key
        resid_text = str(resid).strip()
        try:
            resid_rank = (0, int(resid_text))
        except (TypeError, ValueError):
            resid_rank = (1, resid_text)
        return (chain_id or "", resid_rank, icode or "")

    instances = []
    for residue_key, atom_records in sorted(residues.items(), key=_residue_sort_key):
        instances.append(
            {
                "instance_key": residue_key,
                "resname": residue_key[0],
                "chain_id": residue_key[1],
                "resid": residue_key[2],
                "insertion_code": residue_key[3],
                "atom_count": len(atom_records),
                "atom_records": atom_records,
            }
        )
    return instances


def select_representative_component_instance(pdb_path, component_code):
    instances = list_component_instances(pdb_path, component_code)
    if not instances:
        return None
    return instances[0]


def _write_component_atom_records_to_pdb(atom_records, out_pdb):
    serials = set()
    atom_lines = []
    for record in atom_records:
        line = record["line"]
        atom_lines.append(line if line.endswith("\n") else line + "\n")
        try:
            serials.add(int(line[6:11].strip()))
        except ValueError:
            continue

    conect_lines = []
    if serials:
        for record in atom_records:
            source_path = record.get("source_pdb_path")
            if source_path:
                break
        else:
            source_path = None
        if source_path and os.path.exists(source_path):
            with open(source_path, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if not line.startswith("CONECT"):
                        continue
                    numbers = []
                    for start in (6, 11, 16, 21, 26):
                        chunk = line[start:start + 5].strip()
                        if not chunk:
                            continue
                        try:
                            numbers.append(int(chunk))
                        except ValueError:
                            continue
                    if numbers and any(number in serials for number in numbers):
                        conect_lines.append(line if line.endswith("\n") else line + "\n")

    with open(out_pdb, "w", encoding="utf-8") as handle:
        for line in atom_lines:
            handle.write(line)
        for line in conect_lines:
            handle.write(line)
        handle.write("END\n")
    return True


def extract_component_instances_from_pdb(pdb_path, component_code, out_pdb):
    component_code = component_code.strip().upper()
    instances = list_component_instances(pdb_path, component_code)
    if not instances:
        return False
    all_records = []
    for instance in instances:
        for record in instance["atom_records"]:
            record = dict(record)
            record["source_pdb_path"] = pdb_path
            all_records.append(record)
    return _write_component_atom_records_to_pdb(all_records, out_pdb)


def extract_component_instance_from_pdb(pdb_path, component_code, instance_key, output_pdb):
    component_code = component_code.strip().upper()
    instances = list_component_instances(pdb_path, component_code)
    for instance in instances:
        if instance["instance_key"] == instance_key:
            records = []
            for record in instance["atom_records"]:
                copied = dict(record)
                copied["source_pdb_path"] = pdb_path
                records.append(copied)
            return _write_component_atom_records_to_pdb(records, output_pdb)
    return False


def validate_single_instance_template(template_pdb, component_code):
    component_code = component_code.strip().upper()
    instances = list_component_instances(template_pdb, component_code)
    if len(instances) != 1:
        raise RuntimeError(
            f"ERROR: Component topology template for {component_code} contains {len(instances)} residue instances.\n"
            "Topology generation must use exactly one representative instance.\n"
            "Extracting all copies together creates an invalid aggregate template.\n"
            "Use extract_component_instance_from_pdb(...) before parameterization."
        )
    return instances[0]


def parse_pdb_atoms_by_residue(pdb_path, resname_filter=None):
    residues = OrderedDict()
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = _record_resname(line)
            if resname_filter and resname != resname_filter.upper():
                continue
            key = residue_instance_key(line)
            entry = {
                "line": line,
                "record": line[:6],
                "atom_name": _record_atom_name(line),
                "resname": resname,
                "chain_id": _record_chain_id(line),
                "resid": _record_resid(line),
                "icode": _record_icode(line),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "element": _record_element(line),
            }
            residues.setdefault(key, []).append(entry)
    return residues


def parse_pdb_atom_records(pdb_path):
    atoms = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atoms.append(
                {
                    "line": line,
                    "record": line[:6],
                    "atom_name": _record_atom_name(line),
                    "resname": _record_resname(line),
                    "chain_id": _record_chain_id(line),
                    "resid": _record_resid(line),
                    "icode": _record_icode(line),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": _record_element(line),
                }
            )
    return atoms


def parse_gro_atom_records(gro_path):
    atoms = []
    with open(gro_path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    for index, line in enumerate(lines[2:-1], start=1):
        if len(line) < 44:
            continue
        resid = line[0:5].strip()
        resname = line[5:10].strip().upper()
        atom_name = line[10:15].strip()
        atom_index = int(line[15:20].strip())
        x = float(line[20:28].strip())
        y = float(line[28:36].strip())
        z = float(line[36:44].strip())
        atoms.append(
            {
                "line": line,
                "record": "GRO",
                "atom_name": atom_name,
                "resname": resname,
                "chain_id": "_",
                "resid": resid,
                "icode": "_",
                "x": x,
                "y": y,
                "z": z,
                "element": _record_element(line),
                "atom_index": atom_index,
            }
        )
    return atoms


def load_structure_atoms(path):
    lower = path.lower()
    if lower.endswith(".gro"):
        return parse_gro_atom_records(path), "gro"
    return parse_pdb_atom_records(path), "pdb"


def _coords_array(records):
    if _np is None:
        raise RuntimeError("NumPy is required for topology-consistent component rebuilding.")
    return _np.array([[float(rec["x"]), float(rec["y"]), float(rec["z"])] for rec in records], dtype=float)


def _kabsch_transform(template_coords, target_coords):
    template_centroid = template_coords.mean(axis=0)
    target_centroid = target_coords.mean(axis=0)
    template_centered = template_coords - template_centroid
    target_centered = target_coords - target_centroid
    covariance = template_centered.T @ target_centered
    v, _, wt = _np.linalg.svd(covariance)
    determinant = _np.linalg.det(wt.T @ v.T)
    correction = _np.eye(3)
    if determinant < 0:
        correction[2, 2] = -1
    rotation = wt.T @ correction @ v.T
    translation = target_centroid - (template_centroid @ rotation.T)
    return rotation, translation


def _apply_transform(coords, rotation, translation):
    return coords @ rotation.T + translation


def _rmsd(coords_a, coords_b):
    diff = coords_a - coords_b
    return float(_np.sqrt(_np.mean(_np.sum(diff * diff, axis=1))))


def _fit_with_outlier_pruning(
    template_heavy_by_name,
    raw_heavy_by_name,
    common_names,
    min_common_heavy_atoms,
    max_postfit_rmsd_ang,
):
    active_names = list(common_names)
    dropped_names = []
    best = None

    while len(active_names) >= min_common_heavy_atoms:
        template_common = [template_heavy_by_name[name] for name in active_names]
        raw_common = [raw_heavy_by_name[name] for name in active_names]
        template_common_coords = _coords_array(template_common)
        raw_common_coords = _coords_array(raw_common)
        rotation, translation = _kabsch_transform(template_common_coords, raw_common_coords)
        fitted_common = _apply_transform(template_common_coords, rotation, translation)
        residuals = _np.sqrt(_np.sum((fitted_common - raw_common_coords) ** 2, axis=1)) * 10.0
        postfit_rmsd_ang = _rmsd(fitted_common, raw_common_coords) * 10.0
        worst_index = int(_np.argmax(residuals))
        worst_name = active_names[worst_index]
        worst_residual_ang = float(residuals[worst_index])
        candidate = {
            "active_names": list(active_names),
            "dropped_names": list(dropped_names),
            "rotation": rotation,
            "translation": translation,
            "postfit_rmsd_ang": postfit_rmsd_ang,
            "worst_residual_ang": worst_residual_ang,
            "worst_name": worst_name,
            "residuals_ang": {name: float(value) for name, value in zip(active_names, residuals)},
        }
        if best is None or candidate["postfit_rmsd_ang"] < best["postfit_rmsd_ang"]:
            best = candidate
        if postfit_rmsd_ang <= max_postfit_rmsd_ang:
            return candidate
        if len(active_names) == min_common_heavy_atoms:
            break
        dropped_names.append(worst_name)
        del active_names[worst_index]

    return best


def _minimum_distance_nm(atom_records_a, atom_records_b):
    if not atom_records_a or not atom_records_b:
        return math.inf, None, None
    coords_a = _coords_array(atom_records_a)
    coords_b = _coords_array(atom_records_b)
    best = math.inf
    best_pair = (None, None)
    for i, coord_a in enumerate(coords_a):
        deltas = coords_b - coord_a
        distances = _np.sqrt(_np.sum(deltas * deltas, axis=1))
        local_index = int(_np.argmin(distances))
        local_distance = float(distances[local_index])
        if local_distance < best:
            best = local_distance
            best_pair = (atom_records_a[i], atom_records_b[local_index])
    return best, best_pair[0], best_pair[1]


def _format_pdb_line(serial, atom_name, resname, chain_id, resid, icode, x, y, z, source_line, element):
    atom_field = atom_name.rjust(4) if len(atom_name) < 4 else atom_name[:4]
    chain_field = (chain_id or " ")[:1]
    resid_int = int(str(resid).strip())
    icode_field = " " if icode in {"", "_"} else str(icode)[0]
    line = list(source_line.rstrip("\n").ljust(80))
    line[0:6] = list("HETATM")
    line[6:11] = list(f"{serial:5d}")
    line[12:16] = list(atom_field)
    line[17:20] = list(resname.rjust(3)[:3])
    line[21] = chain_field
    line[22:26] = list(f"{resid_int:4d}")
    line[26] = icode_field
    line[30:38] = list(f"{x:8.3f}")
    line[38:46] = list(f"{y:8.3f}")
    line[46:54] = list(f"{z:8.3f}")
    if element:
        element_field = element.rjust(2)[:2]
        line[76:78] = list(element_field)
    return "".join(line).rstrip() + "\n"


def count_atoms_in_itp_moleculetype(itp_path, moleculetype_name=None):
    atoms_section = False
    count = 0
    current_molecule = None
    capture_next_noncomment = False
    with open(itp_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped or stripped.startswith(";"):
                continue
            if stripped.startswith("["):
                section = stripped.strip("[]").strip().lower()
                atoms_section = section == "atoms" and (moleculetype_name is None or current_molecule == moleculetype_name)
                capture_next_noncomment = section == "moleculetype"
                continue
            if capture_next_noncomment:
                current_molecule = stripped.split()[0]
                capture_next_noncomment = False
                continue
            if atoms_section:
                if stripped.startswith("["):
                    atoms_section = False
                else:
                    count += 1
    return count


def count_gro_atoms(gro_path):
    with open(gro_path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    return int(lines[1].strip())


def _classify_structure_atom(atom, component_resnames):
    resname = atom["resname"].upper()
    if resname in {code.upper() for code in component_resnames}:
        return "component"
    if resname in WATER_RESNAMES:
        return "water"
    if resname in ORDINARY_IONS:
        return "ion"
    if resname in STANDARD_PROTEIN_AND_NA:
        return "protein"
    return "other"


def detect_component_clashes(gro_or_pdb, component_resnames, cutoff_nm=0.08, include_water=False):
    component_resnames = ordered_unique_resnames(component_resnames or [])
    atoms, source_type = load_structure_atoms(gro_or_pdb)
    component_set = set(component_resnames)
    if not atoms or not component_set:
        return []

    clashes = []
    component_atoms = [atom for atom in atoms if atom["resname"].upper() in component_set]
    other_atoms = atoms
    for atom_a in component_atoms:
        class_a = _classify_structure_atom(atom_a, component_resnames)
        for atom_b in other_atoms:
            if atom_a is atom_b:
                continue
            if atom_a["resname"] == atom_b["resname"] and atom_a["resid"] == atom_b["resid"] and atom_a["atom_name"] == atom_b["atom_name"]:
                continue
            class_b = _classify_structure_atom(atom_b, component_resnames)
            if not include_water and "water" in {class_a, class_b}:
                continue
            if atom_a["resname"] == atom_b["resname"] and atom_a["resid"] == atom_b["resid"]:
                continue
            dx = atom_a["x"] - atom_b["x"]
            dy = atom_a["y"] - atom_b["y"]
            dz = atom_a["z"] - atom_b["z"]
            distance = math.sqrt(dx * dx + dy * dy + dz * dz)
            if source_type == "pdb":
                distance_nm = distance / 10.0
            else:
                distance_nm = distance
            if distance_nm < cutoff_nm:
                pair = sorted(
                    [
                        (atom_a.get("atom_index", 0), atom_a["resname"], atom_a["resid"], atom_a["atom_name"]),
                        (atom_b.get("atom_index", 0), atom_b["resname"], atom_b["resid"], atom_b["atom_name"]),
                    ]
                )
                clash_key = tuple(pair)
                clashes.append(
                    {
                        "clash_key": clash_key,
                        "distance_nm": distance_nm,
                        "atom_a": atom_a,
                        "atom_b": atom_b,
                        "class_a": class_a,
                        "class_b": class_b,
                        "pair_type": f"{class_a}-{class_b}",
                    }
                )
    deduped = OrderedDict()
    for clash in sorted(clashes, key=lambda item: item["distance_nm"]):
        deduped.setdefault(clash["clash_key"], clash)
    return list(deduped.values())


def format_clash_report(clashes, limit=20):
    lines = ["Detected severe retained-component clashes:"]
    for clash in clashes[:limit]:
        atom_a = clash["atom_a"]
        atom_b = clash["atom_b"]
        lines.append(
            f"  {atom_a['resname']}:{atom_a['resid']}:{atom_a['atom_name']} "
            f"<-> {atom_b['resname']}:{atom_b['resid']}:{atom_b['atom_name']} "
            f"| distance={clash['distance_nm']:.4f} nm | type={clash['pair_type']}"
        )
    if len(clashes) > limit:
        lines.append(f"  ... {len(clashes) - limit} more clashes omitted")
    return "\n".join(lines)


def describe_gro_atom(gro_path, atom_number):
    atoms = parse_gro_atom_records(gro_path)
    atom_number = int(atom_number)
    for atom in atoms:
        if atom.get("atom_index") == atom_number:
            return {
                "atom_index": atom_number,
                "resid": atom["resid"],
                "resname": atom["resname"],
                "atom_name": atom["atom_name"],
                "coordinate_nm": (atom["x"], atom["y"], atom["z"]),
                "likely_group": atom["resname"],
            }
    return None


def build_component_all_copies_from_template(
    original_pdb,
    component_code,
    complete_template_pdb,
    output_all_pdb,
    min_common_heavy_atoms=3,
    max_postfit_rmsd_ang=2.5,
    min_nonself_distance_ang=0.6,
):
    if _np is None:
        raise RuntimeError("NumPy is required for topology-consistent component rebuilding.")

    component_code = component_code.strip().upper()
    raw_residues = parse_pdb_atoms_by_residue(original_pdb, component_code)
    if not raw_residues:
        raise RuntimeError(f"No residue instances for {component_code} were found in {original_pdb}.")

    template_atoms = parse_pdb_atom_records(complete_template_pdb)
    if not template_atoms:
        raise RuntimeError(f"Template PDB {complete_template_pdb} does not contain any atoms.")

    template_heavy_by_name = OrderedDict()
    for atom in template_atoms:
        if atom["atom_name"] not in template_heavy_by_name and atom["element"] not in {"H"} and not atom["atom_name"].upper().startswith("LP"):
            template_heavy_by_name[atom["atom_name"]] = atom

    all_atoms = parse_pdb_atom_records(original_pdb)
    scaffold_atoms = [atom for atom in all_atoms if atom["resname"] != component_code]
    def _residue_sort_key(item):
        residue_key, _atoms = item
        _resname, chain_id, resid, icode = residue_key
        resid_text = str(resid).strip()
        try:
            resid_value = int(resid_text)
            resid_rank = (0, resid_value)
        except (TypeError, ValueError):
            resid_rank = (1, resid_text)
        return (chain_id or "", resid_rank, icode or "")

    serial = 1
    transformed_lines = []
    copy_reports = []
    prior_rebuilt_atoms = []
    for residue_key, raw_atoms in sorted(raw_residues.items(), key=_residue_sort_key):
        raw_heavy_by_name = {
            atom["atom_name"]: atom
            for atom in raw_atoms
            if atom["element"] not in {"H"} and not atom["atom_name"].upper().startswith("LP")
        }
        common_names = [name for name in template_heavy_by_name if name in raw_heavy_by_name]
        if len(common_names) < min_common_heavy_atoms:
            raise RuntimeError(
                f"Cannot build topology-complete coordinates for {component_code} residue {residue_key}: "
                f"only {len(common_names)} common heavy atoms between template and raw PDB."
            )
        fit_result = _fit_with_outlier_pruning(
            template_heavy_by_name,
            raw_heavy_by_name,
            common_names,
            min_common_heavy_atoms=min_common_heavy_atoms,
            max_postfit_rmsd_ang=max_postfit_rmsd_ang,
        )
        rotation = fit_result["rotation"]
        translation = fit_result["translation"]
        postfit_rmsd_ang = fit_result["postfit_rmsd_ang"]
        active_names = fit_result["active_names"]
        dropped_names = fit_result["dropped_names"]
        if postfit_rmsd_ang > max_postfit_rmsd_ang:
            raise RuntimeError(
                f"Alignment RMSD is too high for {component_code} residue {residue_key}: "
                f"{postfit_rmsd_ang:.3f} Å using {len(active_names)} of {len(common_names)} common heavy atoms. "
                f"Dropped outliers: {', '.join(dropped_names) if dropped_names else 'none'}. "
                f"Worst remaining match: {fit_result['worst_name']} ({fit_result['worst_residual_ang']:.3f} Å)."
            )
        transformed = _apply_transform(_coords_array(template_atoms), rotation, translation)
        transformed_records = []
        for atom, coords in zip(template_atoms, transformed):
            transformed_records.append(
                {
                    "atom_name": atom["atom_name"],
                    "resname": component_code,
                    "chain_id": residue_key[1],
                    "resid": residue_key[2],
                    "icode": residue_key[3],
                    "x": float(coords[0]),
                    "y": float(coords[1]),
                    "z": float(coords[2]),
                    "element": atom["element"],
                    "line": atom["line"],
                }
            )
            transformed_lines.append(
                _format_pdb_line(
                    serial=serial,
                    atom_name=atom["atom_name"],
                    resname=component_code,
                    chain_id=residue_key[1],
                    resid=residue_key[2],
                    icode=residue_key[3],
                    x=coords[0],
                    y=coords[1],
                    z=coords[2],
                    source_line=atom["line"],
                    element=atom["element"],
                )
            )
            serial += 1
        transformed_lines.append("TER\n")
        min_scaffold_nm, close_a, close_b = _minimum_distance_nm(transformed_records, scaffold_atoms)
        min_scaffold_ang = min_scaffold_nm * 10.0 if math.isfinite(min_scaffold_nm) else math.inf
        if min_scaffold_ang < min_nonself_distance_ang:
            raise RuntimeError(
                f"Rebuilt {component_code} residue {residue_key} is clashing after placement: "
                f"minimum scaffold distance {min_scaffold_ang:.3f} Å between "
                f"{close_a['atom_name']} and {close_b['resname']}:{close_b['resid']}:{close_b['atom_name']}."
            )
        min_prior_nm, close_a2, close_b2 = _minimum_distance_nm(transformed_records, prior_rebuilt_atoms)
        min_prior_ang = min_prior_nm * 10.0 if math.isfinite(min_prior_nm) else math.inf
        if min_prior_ang < min_nonself_distance_ang:
            raise RuntimeError(
                f"Rebuilt {component_code} residue {residue_key} is clashing with another rebuilt {component_code} copy: "
                f"minimum inter-copy distance {min_prior_ang:.3f} Å between "
                f"{close_a2['atom_name']} and {close_b2['resname']}:{close_b2['resid']}:{close_b2['atom_name']}."
            )
        copy_reports.append(
            {
                "component_code": component_code,
                "residue_key": residue_key,
                "common_heavy_atoms": len(common_names),
                "fit_heavy_atoms": len(active_names),
                "dropped_heavy_atom_names": list(dropped_names),
                "postfit_rmsd_ang": postfit_rmsd_ang,
                "worst_remaining_residual_ang": fit_result["worst_residual_ang"],
                "worst_remaining_atom_name": fit_result["worst_name"],
                "min_scaffold_distance_ang": min_scaffold_ang,
                "min_prior_component_distance_ang": min_prior_ang,
                "min_scaffold_partner": (
                    None
                    if not math.isfinite(min_scaffold_ang) or close_b is None
                    else f"{close_b['resname']}:{close_b['chain_id']}:{close_b['resid']}:{close_b['atom_name']}"
                ),
                "min_prior_component_partner": (
                    None
                    if not math.isfinite(min_prior_ang) or close_b2 is None
                    else f"{close_b2['resname']}:{close_b2['chain_id']}:{close_b2['resid']}:{close_b2['atom_name']}"
                ),
            }
        )
        prior_rebuilt_atoms.extend(transformed_records)

    with open(output_all_pdb, "w", encoding="utf-8") as handle:
        for line in transformed_lines:
            handle.write(line)
        handle.write("END\n")
    return {
        "output_all_pdb": output_all_pdb,
        "copy_reports": copy_reports,
    }


def validate_component_gro_against_topology(component_code, gro_path, itp_path, molecule_count):
    atoms_per_molecule = count_atoms_in_itp_moleculetype(itp_path)
    coordinate_atoms = count_gro_atoms(gro_path)
    expected_atoms = atoms_per_molecule * int(molecule_count)
    return {
        "component": component_code,
        "atoms_per_molecule": atoms_per_molecule,
        "molecule_count": int(molecule_count),
        "expected_atoms": expected_atoms,
        "coordinate_atoms": coordinate_atoms,
        "ok": expected_atoms == coordinate_atoms,
    }


def validate_complex_against_topology_components(component_infos, molecule_counts):
    results = []
    for info in component_infos:
        component_code = info["component_code"]
        molecule_count = molecule_counts.get(component_code, 0)
        results.append(
            validate_component_gro_against_topology(
                component_code,
                info["all_gro"],
                info["itp_path"],
                molecule_count,
            )
        )
    return results


def format_component_validation_report(results):
    lines = ["Component coordinate/topology validation:"]
    for result in results:
        status = "OK" if result["ok"] else "MISMATCH"
        lines.append(
            f"  {result['component']}: topology atoms/molecule={result['atoms_per_molecule']}, "
            f"copies={result['molecule_count']}, expected={result['expected_atoms']}, "
            f"coordinates={result['coordinate_atoms']} {status}"
        )
    return "\n".join(lines)


def write_component_registry(
    directory,
    primary_ligand=None,
    cofactors=None,
    molecule_counts=None,
    analysis_ligand=None,
    dropped_components=None,
    detected_components=None,
):
    detected_components = detected_components or {}
    registry = {
        "primary_ligand": (primary_ligand or None),
        "cofactors": ordered_unique_resnames(cofactors or []),
        "all_components": get_all_retained_components(primary_ligand, cofactors or []),
        "dropped_components": ordered_unique_resnames(dropped_components or []),
        "detected_components": detected_components,
        "molecule_counts": {k: int(v) for k, v in (molecule_counts or {}).items()},
        "analysis_ligand": analysis_ligand or primary_ligand or None,
    }
    path = os.path.join(directory, COMPONENT_REGISTRY_FILENAME)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(registry, handle, indent=2, sort_keys=False)
        handle.write("\n")
    return path, registry


def load_component_registry(directory="."):
    path = os.path.join(directory, COMPONENT_REGISTRY_FILENAME)
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    data["primary_ligand"] = (data.get("primary_ligand") or None)
    data["cofactors"] = ordered_unique_resnames(data.get("cofactors") or [])
    data["all_components"] = get_all_retained_components(
        data.get("primary_ligand"),
        data.get("cofactors") or [],
    )
    data["dropped_components"] = ordered_unique_resnames(data.get("dropped_components") or [])
    molecule_counts = {}
    for key, value in (data.get("molecule_counts") or {}).items():
        try:
            molecule_counts[key.upper()] = int(value)
        except (TypeError, ValueError):
            continue
    data["molecule_counts"] = molecule_counts
    if not data.get("analysis_ligand"):
        data["analysis_ligand"] = data.get("primary_ligand")
    return data


def merged_component_group_name(components, fallback_ligand=None, protein_only=False):
    ordered = ordered_unique_resnames(components or [])
    if protein_only or not ordered:
        if fallback_ligand:
            return f"Protein_{fallback_ligand.strip().upper()}"
        return "Protein"
    return "Protein_" + "_".join(ordered)


def read_index_groups(index_path):
    groups = []
    if not os.path.exists(index_path):
        return groups
    with open(index_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                groups.append(line[1:-1].strip())
    return groups


def index_has_group(index_path, group_name):
    return group_name in read_index_groups(index_path)
