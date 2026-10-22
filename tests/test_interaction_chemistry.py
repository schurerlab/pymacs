import numpy as np

from pymacs_interaction_chemistry import (
    aggregate_event_occupancy,
    atom_is_hydrogen,
    classify_ligand_residue_events,
    classify_protein_residue_pair_events,
    ligand_role,
)


class Atom:
    def __init__(self, name, element, position=(0, 0, 0), formalcharge=None, residue=None, index=0, charge=None):
        self.name, self.element, self.position, self.formalcharge = name, element, np.asarray(position, dtype=float), formalcharge
        self.index = index
        self.residue = residue
        self.charge = charge


class Group:
    def __init__(self, atoms):
        self.atoms = atoms

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def __getitem__(self, item):
        return self.atoms[item]


class Residue:
    def __init__(self, resname, resid, atoms, segid="A"):
        self.resname, self.resid, self.segid = resname, resid, segid
        self.atoms = Group(atoms)
        for atom in atoms:
            atom.residue = self


def test_hydrogen_is_not_a_ligand_hydrophobic_or_ionic_endpoint():
    hydrogen = Atom("H1", "H")
    assert atom_is_hydrogen(hydrogen)
    assert ligand_role(hydrogen) == "other"


def test_unknown_ligand_charge_is_not_promoted_to_ionic():
    neutral_carbon = Atom("C1", "C")
    assert ligand_role(neutral_carbon) == "hydrophobe"
    assert ligand_role(Atom("N1", "N")) != "cation"
    assert ligand_role(Atom("Npartial", "N", charge=0.9)) != "cation"


def test_occupancy_uses_unique_frames_and_is_bounded():
    events = [
        {"ResidueID": "A:ARG1", "Residue": "ARG1", "LigandCanonicalAtom": "C1", "InteractionClass": "Hydrophobic", "Frame": 0},
        {"ResidueID": "A:ARG1", "Residue": "ARG1", "LigandCanonicalAtom": "C1", "InteractionClass": "Hydrophobic", "Frame": 0},
        {"ResidueID": "A:ARG1", "Residue": "ARG1", "LigandCanonicalAtom": "C1", "InteractionClass": "Hydrophobic", "Frame": 2},
    ]
    row = aggregate_event_occupancy(events, 3).iloc[0]
    assert row.Frames == 2
    assert row.Occupancy == 2 / 3


def test_geometry_validated_hbond_uses_hydrogen_as_evidence_not_endpoint():
    ligand_n = Atom("N1", "N", (0, 0, 0), index=1)
    ligand_h = Atom("H1", "H", (1, 0, 0), index=2)
    ligand = Group([ligand_n, ligand_h])
    oxygen = Atom("OD1", "O", (2.8, 0, 0), index=3)
    residue = Residue("ASP", 10, [oxygen])
    events = classify_ligand_residue_events(ligand, residue, 3, 0.3)
    hbond = next(event for event in events if event["InteractionClass"] == "H-bond")
    assert hbond["LigandAtomName"] == "N1"
    assert hbond["HydrogenEvidenceAtom"] == "H1"
    assert hbond["Angle"] >= 135


def test_unknown_ligand_charge_cannot_form_ionic_contact():
    ligand = Group([Atom("N1", "N", (0, 0, 0), index=1)])
    lys_nz = Atom("NZ", "N", (3, 0, 0), index=2)
    residue = Residue("LYS", 11, [lys_nz])
    events = classify_ligand_residue_events(ligand, residue, 0, 0.0)
    assert all(event["InteractionClass"] != "Ionic" for event in events)


def test_protein_salt_bridge_requires_specific_charged_atoms():
    lys = Residue("LYS", 1, [Atom("NZ", "N", (0, 0, 0), index=1)])
    asp = Residue("ASP", 2, [Atom("OD1", "O", (3, 0, 0), index=2)])
    events = classify_protein_residue_pair_events(lys, asp, 0, 0.0)
    assert any(event["InteractionClass"] == "Ionic" for event in events)
    lys_carbon = Residue("LYS", 1, [Atom("CB", "C", (0, 0, 0), index=1)])
    events = classify_protein_residue_pair_events(lys_carbon, asp, 0, 0.0)
    assert all(event["InteractionClass"] != "Ionic" for event in events)
