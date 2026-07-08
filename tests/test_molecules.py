"""Building RDKit molecules from bare PDB het atoms + loading ligands."""

from rdkit import Chem
from rdkit.Chem import AllChem

from pdbbind_residue_contacts import molecules, pharmacophores
from pdbbind_residue_contacts.molecules import AtomRecord


def _heavy_atom_records(smiles: str) -> list[AtomRecord]:
    """Generate 3D coordinates for a molecule and return heavy-atom records (as a PDB would)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    records = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        records.append(
            AtomRecord(
                name=f"{atom.GetSymbol()}{atom.GetIdx()}",
                element=atom.GetSymbol(),
                x=pos.x, y=pos.y, z=pos.z,
            )
        )
    return records


def test_single_atom_ion():
    mol = molecules.build_mol_from_atoms([AtomRecord("ZN", "Zn", 0.0, 0.0, 0.0)])
    assert mol is not None
    assert mol.GetNumAtoms() == 1
    labels = pharmacophores.assign(mol, {0})
    assert labels == {"metal ion"}


def test_bond_perception_benzene():
    # Het groups are built connectivity-only, so an aromatic ring is labelled via ring membership.
    records = _heavy_atom_records("c1ccccc1")
    mol = molecules.build_mol_from_atoms(records)
    assert mol is not None
    assert mol.GetNumAtoms() == 6
    labels = pharmacophores.assign_connectivity(mol, set(range(mol.GetNumAtoms())))
    assert "aromatic ring (other)" in labels


def test_bond_perception_ethanol_has_bonds():
    records = _heavy_atom_records("CCO")
    mol = molecules.build_mol_from_atoms(records)
    assert mol is not None
    assert mol.GetNumBonds() >= 2  # C-C and C-O perceived from geometry


def test_atom_names_preserved():
    mol = molecules.build_mol_from_atoms(
        [AtomRecord("CX1", "C", 0, 0, 0), AtomRecord("OX2", "O", 1.4, 0, 0)]
    )
    names = {a.GetProp("atomName") for a in mol.GetAtoms()}
    assert names == {"CX1", "OX2"}


def test_load_ligand_from_sdf(tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    AllChem.EmbedMolecule(mol, randomSeed=2)
    sdf = tmp_path / "lig.sdf"
    with Chem.SDWriter(str(sdf)) as w:
        w.write(mol)
    loaded = molecules.load_ligand(sdf, None)
    assert loaded is not None
    assert loaded.GetNumConformers() == 1
