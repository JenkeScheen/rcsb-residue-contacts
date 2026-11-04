"""Parsing ModelServer ligand SDF into RDKit molecules."""

from rdkit import Chem
from rdkit.Chem import AllChem

from conftest import fixture_text
from rcsb_residue_contacts import molecules
from rcsb_residue_contacts.stacking import ligand_rings


def test_rcsb_sdf_parses_with_bond_orders():
    mol = molecules.mol_from_sdf(fixture_text("6PML_D.sdf"))
    assert mol is not None
    assert mol.GetNumConformers() == 1
    assert all(a.GetAtomicNum() > 1 for a in mol.GetAtoms())
    assert sum(a.GetIsAromatic() for a in mol.GetAtoms()) >= 15  # benzofuran + phenyl


def test_hydrogens_are_removed(tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    AllChem.EmbedMolecule(mol, randomSeed=2)
    loaded = molecules.mol_from_sdf(Chem.MolToMolBlock(mol))
    assert loaded.GetNumAtoms() == 7


def _without_charges(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    return "\n".join(l for l in Chem.MolToMolBlock(mol).splitlines() if not l.startswith("M  CHG"))


def test_missing_formal_charges_are_repaired():
    # Zwitterionic glycine with explicit NH3 hydrogens but no charge record (as in some
    # ModelServer SDFs for free amino acids): N has valence 4 and RDKit rejects it as-is.
    block = _without_charges("[NH3+]CC(=O)[O-]")
    assert Chem.MolFromMolBlock(block, removeHs=False) is None
    mol = molecules.mol_from_sdf(block)
    assert mol is not None
    assert {a.GetSymbol(): a.GetFormalCharge() for a in mol.GetAtoms()}["N"] == 1
    boronate = molecules.mol_from_sdf(_without_charges("O[B-](O)(O)c1ccccc1"))
    assert boronate is not None and len(ligand_rings(boronate)) == 1


def test_empty_or_garbage_sdf():
    assert molecules.mol_from_sdf("") is None
    assert molecules.mol_from_sdf("not an sdf") is None


def test_element_coords_and_smiles():
    mol = molecules.mol_from_sdf(fixture_text("6PML_D.sdf"))
    hetero = molecules.element_coords(mol, {"N", "O"})
    assert hetero.shape[1] == 3 and len(hetero) > 0
    assert molecules.element_coords(mol, molecules.METALS).shape == (0, 3)
    assert Chem.MolFromSmiles(molecules.to_smiles(mol)) is not None
