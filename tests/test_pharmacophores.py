"""Pharmacophore assignment on molecules with known functional groups."""

from rdkit import Chem

from rcsb_residue_contacts import pharmacophores


def _labels(smiles: str) -> set[str]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    all_idx = set(range(mol.GetNumAtoms()))
    return pharmacophores.assign(mol, all_idx)


def test_benzene():
    assert "benzene/phenyl" in _labels("c1ccccc1")


def test_nitrile():
    assert "nitrile" in _labels("CC#N")


def test_carboxylic_acid():
    assert "carboxylate/acid" in _labels("CC(=O)O")


def test_sulfonamide():
    labels = _labels("O=S(=O)(N)c1ccccc1")
    assert "sulfonamide" in labels
    assert "benzene/phenyl" in labels


def test_halogen():
    assert "halogen" in _labels("c1ccc(Cl)cc1")


def test_hydroxyl():
    assert "hydroxyl" in _labels("CCO")


def test_every_atom_gets_a_label():
    mol = Chem.MolFromSmiles("CCO")
    labels = pharmacophores.assign(mol, set(range(mol.GetNumAtoms())))
    assert labels  # nothing silently dropped


def test_indices_outside_contact_are_ignored():
    mol = Chem.MolFromSmiles("CC#N")
    # Only the terminal methyl carbon (index 0) is in contact -> no nitrile.
    labels = pharmacophores.assign(mol, {0})
    assert "nitrile" not in labels
