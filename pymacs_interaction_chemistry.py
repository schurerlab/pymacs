"""Shared, conservative Level 3 interaction chemistry.

This module is deliberately independent of a particular execution strategy
(serial, MPI, or chunked).  It separates raw proximity from chemically
validated events: when topology cannot establish a specific interaction, the
result is ``Other Contact`` rather than an invented ionic/hydrophobic label.
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import numpy as np


SCHEMA_VERSION = "pymacs.interaction-events.v1"
EVENT_COLUMNS = [
    "SchemaVersion", "Frame", "Time_ns", "InteractionClass",
    "ValidationStatus", "ClassificationSource", "Partner", "Component",
    "LigandAtomName", "LigandAtomIndex", "LigandCanonicalAtom",
    "ProteinAtomName", "ProteinAtomIndex", "ProteinCanonicalAtom",
    "Residue", "ChainID", "ResidueID", "LigandRole", "ProteinRole",
    "LigandFunctionalGroup", "ProteinFunctionalGroup", "Distance",
    "Angle", "DonorAtom", "AcceptorAtom", "HydrogenEvidenceAtom",
    "RawLigandAtomName", "RawProteinAtomName",
]

HBOND_DISTANCE_DEFAULT = 3.5
HBOND_ANGLE_DEFAULT = 135.0
HYDROPHOBIC_DISTANCE_DEFAULT = 4.5
IONIC_DISTANCE_DEFAULT = 4.5

_POSITIVE_PROTEIN = {
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
    ("LYS", "NZ"),
}
_NEGATIVE_PROTEIN = {
    ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"),
}
_PROTEIN_ACCEPTORS = {
    ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"),
    ("ASN", "OD1"), ("GLN", "OE1"), ("HIS", "ND1"), ("HIS", "NE2"),
    ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("CYS", "SG"),
}
_PROTEIN_DONORS = {
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("LYS", "NZ"),
    ("ASN", "ND2"), ("GLN", "NE2"), ("SER", "OG"), ("THR", "OG1"),
    ("TYR", "OH"), ("CYS", "SG"), ("TRP", "NE1"),
}
_HYDROPHOBIC_PROTEIN_RESIDUES = {"ALA", "VAL", "ILE", "LEU", "MET", "PRO", "PHE", "TRP", "TYR"}


def infer_element(atom) -> str:
    """Best-effort element inference without treating atom names as chemistry."""
    element = str(getattr(atom, "element", "") or "").strip()
    if element:
        return element.title()
    name = str(getattr(atom, "name", "") or "").strip()
    while name and name[0].isdigit():
        name = name[1:]
    if name[:2].title() in {"Cl", "Br"}:
        return name[:2].title()
    return name[:1].upper()


def atom_is_hydrogen(atom) -> bool:
    return infer_element(atom) == "H"


def atom_is_virtual(atom) -> bool:
    name = str(getattr(atom, "name", "") or "").upper()
    atype = str(getattr(atom, "type", "") or "").upper()
    return name.startswith(("LP", "VS", "MWS")) or atype.startswith(("LP", "VS"))


def heavy_atoms(atomgroup):
    return [atom for atom in atomgroup.atoms if not atom_is_hydrogen(atom) and not atom_is_virtual(atom)]


def _resname(atom) -> str:
    return str(atom.residue.resname).upper()


def _atom_name(atom) -> str:
    return str(atom.name).upper()


def protein_role(atom) -> str:
    key = (_resname(atom), _atom_name(atom))
    if key in _POSITIVE_PROTEIN:
        return "cation"
    if key in _NEGATIVE_PROTEIN:
        return "anion"
    if key in _PROTEIN_DONORS:
        return "donor"
    if key in _PROTEIN_ACCEPTORS or _atom_name(atom) == "O":
        return "acceptor"
    if infer_element(atom) in {"C", "S"} and _resname(atom) in _HYDROPHOBIC_PROTEIN_RESIDUES:
        return "hydrophobe"
    return "other"


def ligand_role(atom) -> str:
    """Return only a defensible role from universally available topology data.

    Formal charge is accepted only from an explicit ``formalcharge`` topology
    attribute.  MD force-field partial charges are intentionally not used to
    invent a formal ionic state.
    """
    formal = getattr(atom, "formalcharge", None)
    try:
        if formal is not None and float(formal) > 0.5:
            return "cation"
        if formal is not None and float(formal) < -0.5:
            return "anion"
    except (TypeError, ValueError):
        pass
    element = infer_element(atom)
    if element in {"C", "S", "F", "Cl", "Br", "I"}:
        return "hydrophobe"
    # N/O/S can participate in H-bonds; whether they are donor or acceptor is
    # established below using bonded-H geometry, not assumed from proximity.
    if element in {"N", "O", "S"}:
        return "hetero"
    return "other"


def _attached_hydrogens(atom, group, maximum=1.35):
    hydrogens = [candidate for candidate in group.atoms if atom_is_hydrogen(candidate)]
    if not hydrogens:
        return []
    return [
        hydrogen for hydrogen in hydrogens
        if 0.6 <= float(np.linalg.norm(hydrogen.position - atom.position)) <= maximum
    ]


def _dha_angle(donor, hydrogen, acceptor) -> Optional[float]:
    first = donor.position - hydrogen.position
    second = acceptor.position - hydrogen.position
    denom = float(np.linalg.norm(first) * np.linalg.norm(second))
    if denom == 0:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(first, second) / denom, -1.0, 1.0))))


def _ligand_acceptor(atom) -> bool:
    # Conservative generic fallback: O/S are accepted; a neutral ligand N is
    # only accepted when explicit topology chemistry is supplied in a future
    # adapter.  This avoids calling amide/charged nitrogen an acceptor.
    return infer_element(atom) in {"O", "S"}


def _ligand_donor(atom, ligand_group) -> bool:
    return infer_element(atom) in {"N", "O", "S"} and bool(_attached_hydrogens(atom, ligand_group))


def _protein_acceptor(atom) -> bool:
    return (_resname(atom), _atom_name(atom)) in _PROTEIN_ACCEPTORS or _atom_name(atom) == "O"


def _protein_donor(atom, residue_group) -> bool:
    return ((_resname(atom), _atom_name(atom)) in _PROTEIN_DONORS or _atom_name(atom) == "N") and bool(_attached_hydrogens(atom, residue_group))


def _event(frame, time_ns, interaction_class, ligand_atom, protein_atom, distance,
           residue, *, angle=np.nan, donor=None, acceptor=None, hydrogen=None,
           ligand_role_value=None, protein_role_value=None, status="validated",
           source="shared_chemistry_v1", component="", partner="protein") -> dict:
    # Residue-level chain metadata is inconsistent across topology formats;
    # atom-level values are the most dependable fallback.
    residue_atoms = list(residue.atoms)
    first_atom = residue_atoms[0] if residue_atoms else None
    chain = str(
        getattr(residue, "segid", "") or getattr(residue, "chainID", "")
        or getattr(first_atom, "segid", "") or getattr(first_atom, "chainID", "") or ""
    )
    residue_id = f"{chain}:{residue.resname}{residue.resid}" if chain else f"{residue.resname}{residue.resid}"
    return {
        "SchemaVersion": SCHEMA_VERSION, "Frame": int(frame), "Time_ns": float(time_ns),
        "InteractionClass": interaction_class, "ValidationStatus": status,
        "ClassificationSource": source, "Partner": partner, "Component": component,
        "LigandAtomName": str(ligand_atom.name), "LigandAtomIndex": int(ligand_atom.index),
        "LigandCanonicalAtom": str(ligand_atom.name), "ProteinAtomName": str(protein_atom.name),
        "ProteinAtomIndex": int(protein_atom.index), "ProteinCanonicalAtom": str(protein_atom.name),
        "Residue": f"{residue.resname}{residue.resid}", "ChainID": chain, "ResidueID": residue_id,
        "LigandRole": ligand_role_value or ligand_role(ligand_atom),
        "ProteinRole": protein_role_value or protein_role(protein_atom),
        "LigandFunctionalGroup": ligand_role_value or ligand_role(ligand_atom),
        "ProteinFunctionalGroup": protein_role_value or protein_role(protein_atom),
        "Distance": float(distance), "Angle": float(angle) if angle is not None else np.nan,
        "DonorAtom": str(donor.name) if donor is not None else "",
        "AcceptorAtom": str(acceptor.name) if acceptor is not None else "",
        "HydrogenEvidenceAtom": str(hydrogen.name) if hydrogen is not None else "",
        "RawLigandAtomName": str(ligand_atom.name), "RawProteinAtomName": str(protein_atom.name),
    }


def classify_ligand_residue_events(ligand_group, residue, frame, time_ns, *,
                                   hbond_distance=HBOND_DISTANCE_DEFAULT,
                                   hbond_angle=HBOND_ANGLE_DEFAULT,
                                   hydrophobic_distance=HYDROPHOBIC_DISTANCE_DEFAULT,
                                   ionic_distance=IONIC_DISTANCE_DEFAULT,
                                   component="", partner="protein") -> List[dict]:
    """Return validated events for one ligand/residue/frame.

    The method records one best event per class/residue/frame.  This prevents a
    cloud of equivalent hydrogen contacts while retaining distinct classes and
    distinct ligand anchors across frames for later aggregation.
    """
    residue_group = residue.atoms
    lig_heavy, prot_heavy = heavy_atoms(ligand_group), heavy_atoms(residue_group)
    if not lig_heavy or not prot_heavy:
        return []
    events: List[dict] = []
    candidates = []
    for lig_atom in lig_heavy:
        for prot_atom in prot_heavy:
            distance = float(np.linalg.norm(lig_atom.position - prot_atom.position))
            candidates.append((distance, lig_atom, prot_atom))
    candidates.sort(key=lambda item: item[0])

    # Geometry-validated H bonds.  Both directions are considered and the
    # endpoints remain heavy atoms; H is evidence only.
    best_hbond = None
    for distance, lig_atom, prot_atom in candidates:
        if distance > hbond_distance:
            break
        directions = []
        if _ligand_donor(lig_atom, ligand_group) and _protein_acceptor(prot_atom):
            directions.append((lig_atom, ligand_group, prot_atom))
        if _protein_donor(prot_atom, residue_group) and _ligand_acceptor(lig_atom):
            directions.append((prot_atom, residue_group, lig_atom))
        for donor, donor_group, acceptor in directions:
            for hydrogen in _attached_hydrogens(donor, donor_group):
                angle = _dha_angle(donor, hydrogen, acceptor)
                if angle is not None and angle >= hbond_angle:
                    score = (angle, -distance)
                    if best_hbond is None or score > best_hbond[0]:
                        best_hbond = (score, distance, lig_atom, prot_atom, donor, acceptor, hydrogen, angle)
    if best_hbond is not None:
        _, distance, lig_atom, prot_atom, donor, acceptor, hydrogen, angle = best_hbond
        events.append(_event(frame, time_ns, "H-bond", lig_atom, prot_atom, distance, residue,
                             angle=angle, donor=donor, acceptor=acceptor, hydrogen=hydrogen,
                             source="donor_hydrogen_acceptor_geometry", component=component, partner=partner))

    # Verified salt bridges require an explicit ligand formal-charge attribute.
    # If absent, a neutral/unknown ligand cannot emit Ionic events.
    ionic = [
        (distance, lig_atom, prot_atom) for distance, lig_atom, prot_atom in candidates
        if distance <= ionic_distance
        and {ligand_role(lig_atom), protein_role(prot_atom)} == {"cation", "anion"}
    ]
    if ionic:
        distance, lig_atom, prot_atom = ionic[0]
        events.append(_event(frame, time_ns, "Ionic", lig_atom, prot_atom, distance, residue,
                             source="explicit_opposite_charge_groups", component=component, partner=partner))

    hydrophobic = [
        (distance, lig_atom, prot_atom) for distance, lig_atom, prot_atom in candidates
        if distance <= hydrophobic_distance
        and ligand_role(lig_atom) == "hydrophobe" and protein_role(prot_atom) == "hydrophobe"
    ]
    if hydrophobic:
        distance, lig_atom, prot_atom = hydrophobic[0]
        events.append(_event(frame, time_ns, "Hydrophobic", lig_atom, prot_atom, distance, residue,
                             source="nonpolar_heavy_atom_pair", component=component, partner=partner))

    if not events and candidates:
        distance, lig_atom, prot_atom = candidates[0]
        events.append(_event(frame, time_ns, "Other Contact", lig_atom, prot_atom, distance, residue,
                             status="proximity_only", source="heavy_atom_proximity", component=component, partner=partner))
    return events


def classify_protein_residue_pair_events(residue_a, residue_b, frame, time_ns, *,
                                         hbond_distance=HBOND_DISTANCE_DEFAULT,
                                         hbond_angle=HBOND_ANGLE_DEFAULT,
                                         hydrophobic_distance=HYDROPHOBIC_DISTANCE_DEFAULT,
                                         ionic_distance=IONIC_DISTANCE_DEFAULT,
                                         component="interface") -> List[dict]:
    """Classify one protein--protein residue pair using atom-level evidence.

    The event schema is deliberately reused, with ``Ligand*`` meaning the
    first residue endpoint.  Consumers should rely on ``Partner`` and the
    residue-pair CSV for interface identity rather than treating it as ligand
    chemistry.  A residue-name match alone never produces salt-bridge or
    hydrophobic labels.
    """
    atoms_a, atoms_b = heavy_atoms(residue_a.atoms), heavy_atoms(residue_b.atoms)
    if not atoms_a or not atoms_b:
        return []
    candidates = sorted(
        (float(np.linalg.norm(a.position - b.position)), a, b)
        for a in atoms_a for b in atoms_b
    )
    events = []
    best_hbond = None
    for distance, atom_a, atom_b in candidates:
        if distance > hbond_distance:
            break
        directions = []
        if _protein_donor(atom_a, residue_a.atoms) and _protein_acceptor(atom_b):
            directions.append((atom_a, residue_a.atoms, atom_b))
        if _protein_donor(atom_b, residue_b.atoms) and _protein_acceptor(atom_a):
            directions.append((atom_b, residue_b.atoms, atom_a))
        for donor, donor_group, acceptor in directions:
            for hydrogen in _attached_hydrogens(donor, donor_group):
                angle = _dha_angle(donor, hydrogen, acceptor)
                if angle is not None and angle >= hbond_angle:
                    score = (angle, -distance)
                    if best_hbond is None or score > best_hbond[0]:
                        best_hbond = (score, distance, atom_a, atom_b, donor, acceptor, hydrogen, angle)
    if best_hbond:
        _, distance, atom_a, atom_b, donor, acceptor, hydrogen, angle = best_hbond
        events.append(_event(frame, time_ns, "H-bond", atom_a, atom_b, distance, residue_b,
                             angle=angle, donor=donor, acceptor=acceptor, hydrogen=hydrogen,
                             source="protein_donor_hydrogen_acceptor_geometry",
                             component=component, partner=str(residue_a.resname)))
    ionic = [(d, a, b) for d, a, b in candidates if d <= ionic_distance
             and {protein_role(a), protein_role(b)} == {"cation", "anion"}]
    if ionic:
        distance, atom_a, atom_b = ionic[0]
        events.append(_event(frame, time_ns, "Ionic", atom_a, atom_b, distance, residue_b,
                             source="protein_opposite_charge_groups", component=component,
                             partner=str(residue_a.resname), ligand_role_value=protein_role(atom_a),
                             protein_role_value=protein_role(atom_b)))
    hydrophobic = [(d, a, b) for d, a, b in candidates if d <= hydrophobic_distance
                   and protein_role(a) == "hydrophobe" and protein_role(b) == "hydrophobe"]
    if hydrophobic:
        distance, atom_a, atom_b = hydrophobic[0]
        events.append(_event(frame, time_ns, "Hydrophobic", atom_a, atom_b, distance, residue_b,
                             source="protein_nonpolar_heavy_atom_pair", component=component,
                             partner=str(residue_a.resname), ligand_role_value=protein_role(atom_a),
                             protein_role_value=protein_role(atom_b)))
    if not events:
        distance, atom_a, atom_b = candidates[0]
        events.append(_event(frame, time_ns, "Other Contact", atom_a, atom_b, distance, residue_b,
                             status="proximity_only", source="protein_heavy_atom_proximity",
                             component=component, partner=str(residue_a.resname),
                             ligand_role_value=protein_role(atom_a), protein_role_value=protein_role(atom_b)))
    return events


def classify_water_bridge_events(ligand_group, water_residue, protein_residue, frame, time_ns, *,
                                 hbond_distance=HBOND_DISTANCE_DEFAULT,
                                 hbond_angle=HBOND_ANGLE_DEFAULT,
                                 component="", partner="protein") -> List[dict]:
    """Return a water bridge only when both legs pass H-bond geometry checks.

    A solvent-free topology yields no event.  This deliberately does not use
    oxygen--oxygen or ligand--residue distance as a water-bridge surrogate.
    """
    water_atoms = heavy_atoms(water_residue.atoms)
    if not water_atoms:
        return []
    water = next((atom for atom in water_atoms if infer_element(atom) == "O"), None)
    if water is None:
        return []

    def bridge_leg(left_atoms, left_group, left_donor, left_acceptor):
        for atom in left_atoms:
            distance = float(np.linalg.norm(atom.position - water.position))
            if distance > hbond_distance:
                continue
            # ligand/protein donor -> water acceptor
            if left_donor(atom, left_group):
                for hydrogen in _attached_hydrogens(atom, left_group):
                    angle = _dha_angle(atom, hydrogen, water)
                    if angle is not None and angle >= hbond_angle:
                        return atom, distance, hydrogen, angle
            # water donor -> ligand/protein acceptor
            if left_acceptor(atom):
                for hydrogen in _attached_hydrogens(water, water_residue.atoms):
                    angle = _dha_angle(water, hydrogen, atom)
                    if angle is not None and angle >= hbond_angle:
                        return atom, distance, hydrogen, angle
        return None

    ligand_leg = bridge_leg(heavy_atoms(ligand_group), ligand_group, _ligand_donor, _ligand_acceptor)
    protein_leg = bridge_leg(heavy_atoms(protein_residue.atoms), protein_residue.atoms, _protein_donor, _protein_acceptor)
    if ligand_leg is None or protein_leg is None:
        return []
    ligand_atom, ligand_distance, ligand_h, ligand_angle = ligand_leg
    protein_atom, protein_distance, protein_h, protein_angle = protein_leg
    event = _event(frame, time_ns, "Water bridge", ligand_atom, protein_atom,
                   max(ligand_distance, protein_distance), protein_residue,
                   angle=min(ligand_angle, protein_angle), donor=water, acceptor=protein_atom,
                   hydrogen=protein_h or ligand_h, source="two_leg_water_hbond_geometry",
                   component=component, partner=partner)
    event["WaterResidue"] = f"{water_residue.resname}{water_residue.resid}"
    event["LigandWaterDistance"] = ligand_distance
    event["WaterProteinDistance"] = protein_distance
    return [event]


def aggregate_event_occupancy(events, total_frames: int):
    """Aggregate by residue, canonical ligand atom, class, and unique frame."""
    import pandas as pd
    if events is None or len(events) == 0:
        return pd.DataFrame(columns=["ResidueID", "Residue", "LigandCanonicalAtom", "InteractionClass", "Frames", "Occupancy"])
    frame_table = events if isinstance(events, pd.DataFrame) else pd.DataFrame(events)
    keys = ["ResidueID", "Residue", "LigandCanonicalAtom", "InteractionClass"]
    unique = frame_table.drop_duplicates(keys + ["Frame"])
    out = unique.groupby(keys, as_index=False)["Frame"].nunique().rename(columns={"Frame": "Frames"})
    out["Occupancy"] = out["Frames"] / float(max(1, total_frames))
    return out.sort_values(["Occupancy", "Residue", "LigandCanonicalAtom"], ascending=[False, True, True]).reset_index(drop=True)


def events_to_legacy_framewise(events):
    """Provide old column names alongside the explicit versioned schema.

    This is a read-compatibility bridge only.  Consumers should prefer
    ``ValidatedInteractionEvents.csv`` and inspect ``ValidationStatus``.
    """
    import pandas as pd
    table = events if isinstance(events, pd.DataFrame) else pd.DataFrame(events, columns=EVENT_COLUMNS)
    if table.empty:
        for column in EVENT_COLUMNS:
            if column not in table:
                table[column] = pd.Series(dtype=object)
    table = table.copy()
    table["Type"] = table.get("InteractionClass", pd.Series(dtype=object))
    table["LigAtomName"] = table.get("LigandCanonicalAtom", pd.Series(dtype=object))
    table["LigAtomIndex"] = table.get("LigandAtomIndex", pd.Series(dtype=float))
    table["ProtAtomName"] = table.get("ProteinCanonicalAtom", pd.Series(dtype=object))
    table["ProtAtomIndex"] = table.get("ProteinAtomIndex", pd.Series(dtype=float))
    table["HbondDonorAtom"] = table.get("DonorAtom", pd.Series(dtype=object))
    table["HbondAcceptorAtom"] = table.get("AcceptorAtom", pd.Series(dtype=object))
    table["HbondAngle"] = table.get("Angle", pd.Series(dtype=float))
    return table


def legacy_summary(events, interaction_classes=("H-bond", "Hydrophobic", "Ionic")):
    """Unique-frame legacy summary; type counts never exceed sampled frames."""
    import pandas as pd
    table = events_to_legacy_framewise(events)
    if table.empty:
        return pd.DataFrame(columns=["Residue", *interaction_classes, "Total"])
    table = table[table["InteractionClass"].isin(interaction_classes)].copy()
    unique = table.drop_duplicates(["ResidueID", "Residue", "InteractionClass", "Frame"])
    summary = unique.groupby(["Residue", "InteractionClass"])["Frame"].nunique().unstack(fill_value=0)
    summary = summary.reindex(columns=list(interaction_classes), fill_value=0)
    summary["Total"] = unique.drop_duplicates(["ResidueID", "Residue", "Frame"]).groupby("Residue")["Frame"].nunique()
    return summary.fillna(0).astype(int).reset_index()
