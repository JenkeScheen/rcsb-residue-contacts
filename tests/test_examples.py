"""Example-structure rendering: full pocket + interacting molecule, with per-residue fallback."""

import numpy as np

from pdbbind_residue_contacts.contacts import ResidueInstance
from pdbbind_residue_contacts.examples import render_example_pdb
from pdbbind_residue_contacts.molecules import AtomRecord


def _residue() -> ResidueInstance:
    atoms = [AtomRecord(name="CA", element="C", x=1.0, y=1.0, z=1.0)]
    return ResidueInstance(
        chain="A", seqid="38", name="TRP", atoms=atoms,
        coords=np.array([[1.0, 1.0, 1.0]]),
    )


_MOL_ATOMS = [("C1", "C", 5.0, 5.0, 5.0), ("O1", "O", 6.0, 5.0, 5.0)]


def test_render_includes_full_pocket_and_molecule():
    pocket = [
        ("ATOM", "N", "HIS", "A", "27", "N", 0.0, 0.0, 0.0),
        ("ATOM", "CA", "HIS", "A", "27", "C", 0.5, 0.0, 0.0),
        ("ATOM", "CA", "TRP", "A", "38", "C", 1.0, 1.0, 1.0),
    ]
    text = render_example_pdb(_residue(), "LIG", _MOL_ATOMS, pocket_atoms=pocket)
    atom_lines = [l for l in text.splitlines() if l.startswith("ATOM")]
    het_lines = [l for l in text.splitlines() if l.startswith("HETATM")]
    # Whole pocket is present (both residues), not just the queried TRP.
    assert len(atom_lines) == 3
    assert any("HIS" in l for l in atom_lines) and any("TRP" in l for l in atom_lines)
    assert len(het_lines) == len(_MOL_ATOMS)
    # Continuous serial numbering across pocket + molecule.
    serials = [int(l[6:11]) for l in atom_lines + het_lines]
    assert serials == list(range(1, len(serials) + 1))


def test_render_falls_back_to_residue_when_no_pocket():
    for pocket in (None, []):
        text = render_example_pdb(_residue(), "LIG", _MOL_ATOMS, pocket_atoms=pocket)
        atom_lines = [l for l in text.splitlines() if l.startswith("ATOM")]
        assert len(atom_lines) == 1 and "TRP" in atom_lines[0]
        assert sum(l.startswith("HETATM") for l in text.splitlines()) == len(_MOL_ATOMS)
