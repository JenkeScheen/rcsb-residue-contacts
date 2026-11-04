"""Plot rendering (synthetic tables) and the interaction-mode classification."""

import pandas as pd

from rcsb_residue_contacts import plots
from rcsb_residue_contacts.pipeline import CONTACT_COLUMNS, LIGAND_COLUMNS, STACK_COLUMNS


def _stack(entry, comp, kind, polar=False, metal=False, pattern="benzene", smiles="c1ccccc1",
           year="2021", distance=3.8, angle=10.0, offset=1.0):
    return {
        "entry_id": entry, "asym_id": "B", "comp_id": comp, "is_loi": True,
        "residue_name": "HIS", "residue_chain": "A", "residue_seqid": "1",
        "residue_ring": "imidazole", "ring_pattern": pattern, "ring_type": pattern,
        "fused": False, "ring_system_smiles": smiles, "ring_size": 6, "distance": distance,
        "angle": angle, "offset": offset, "vertical": 3.5, "stack_type": kind,
        "polar_contact": polar, "metal_coordination": metal, "residue_polar_min_distance": 4.0,
        "release_date": f"{year}-01-01", "resolution": 2.0,
    }


def _ligand(entry, comp, polar=0, metal=0, status="ok"):
    row = dict.fromkeys(LIGAND_COLUMNS)
    row.update(entry_id=entry, asym_id="B", comp_id=comp, status=status,
               polar_contact_residues=polar, metal_contact_residues=metal)
    return row


def _tables():
    stacks = pd.DataFrame(
        [
            _stack("1AAA", "L1", "face-to-face"),
            _stack("1AAA", "L1", "face-to-face"),  # second chain: same pair, counted once
            _stack("2BBB", "L2", "face-to-face", polar=True, pattern="pyridine",
                   smiles="c1ccncc1", year="2015"),
            _stack("3CCC", "L3", "edge-to-face", angle=80),
            _stack("4DDD", "L4", "none", distance=6.2, angle=45, offset=3.0),
        ],
        columns=STACK_COLUMNS,
    )
    ligands = pd.DataFrame(
        [_ligand("1AAA", "L1"), _ligand("2BBB", "L2", polar=1), _ligand("3CCC", "L3"),
         _ligand("4DDD", "L4", polar=1), _ligand("5EEE", "L5"),
         _ligand("6FFF", "L6", status="unparseable")],
        columns=LIGAND_COLUMNS,
    )
    contacts = pd.DataFrame(
        [
            {"entry_id": "1AAA", "asym_id": "B", "comp_id": "L1", "residue_name": "HIS",
             "residue_chain": "A", "residue_seqid": "1", "pharmacophore": p,
             "n_contact_atoms": 3, "min_distance": 3.5}
            for p in ["benzene/phenyl", "amide"]
        ]
        + [
            {"entry_id": "4DDD", "asym_id": "B", "comp_id": "L4", "residue_name": "HIS",
             "residue_chain": "A", "residue_seqid": "1", "pharmacophore": "hydroxyl",
             "n_contact_atoms": 1, "min_distance": 2.9}
        ],
        columns=CONTACT_COLUMNS,
    )
    return ligands, stacks, contacts


def test_interaction_modes_priority_and_dedup():
    ligands, stacks, _ = _tables()
    modes = plots.interaction_modes(ligands, stacks)
    assert modes.to_dict() == {
        "face-to-face, no polar/metal contact": 1,
        "face-to-face + polar/metal contact": 1,
        "edge-to-face only": 1,
        "polar/metal contact only": 1,
        "other contacts only": 1,
    }


def test_make_all_writes_every_plot(tmp_path):
    ligands, stacks, contacts = _tables()
    out = plots.make_all(ligands, stacks, contacts, tmp_path, residue="HIS")
    names = {p.name for p in out.glob("*.png")}
    assert names == {
        "face_to_face_ring_patterns.png", "face_to_face_ring_systems.png",
        "stacking_geometry.png", "interaction_modes.png", "face_to_face_by_year.png",
        "contact_pharmacophores.png",
    }
    assert all((out / n).stat().st_size > 1000 for n in names)


def test_make_all_handles_empty_tables(tmp_path):
    empty = (
        pd.DataFrame(columns=LIGAND_COLUMNS),
        pd.DataFrame(columns=STACK_COLUMNS),
        pd.DataFrame(columns=CONTACT_COLUMNS),
    )
    out = plots.make_all(*empty, tmp_path, residue="HIS")
    assert len(list(out.glob("*.png"))) == 6
