"""Build RDKit molecules (with 3D coordinates) from PDBBind ligand files and PDB het atoms.

The drug-like ligand comes from an SDF (fallback mol2), which carries a full bond graph. Solvent
het groups come straight from the protein PDB as bare atoms + coordinates, so their bonds are
perceived locally with RDKit -- no external chemical-component lookup is needed.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdDetermineBonds

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)


@dataclass
class AtomRecord:
    """A single atom read from a PDB het group."""

    name: str
    element: str
    x: float
    y: float
    z: float


def _conformer_coords(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())], dtype=float
    )


def load_ligand(sdf_path, mol2_path) -> Chem.Mol | None:
    """Load the PDBBind ligand as an RDKit molecule with a 3D conformer.

    Tries the SDF first, then the mol2. Returns ``None`` if neither can be read.
    """
    mol = None
    if sdf_path is not None and sdf_path.exists():
        try:
            supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
            mol = next((m for m in supplier if m is not None), None)
        except Exception as exc:  # pragma: no cover - malformed file
            log.debug("SDMolSupplier failed for %s: %s", sdf_path, exc)
    if mol is None and mol2_path is not None and mol2_path.exists():
        try:
            mol = Chem.MolFromMol2File(str(mol2_path), removeHs=False, sanitize=True)
        except Exception as exc:  # pragma: no cover - malformed file
            log.debug("MolFromMol2File failed for %s: %s", mol2_path, exc)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    # Drop every hydrogen so contacts are heavy-atom to heavy-atom (matching the residue side and
    # the documented cutoff semantics). RemoveAllHs also strips hydrogens that RemoveHs would keep
    # (e.g. on stereocentres). The conformer is preserved.
    try:
        mol = Chem.rdmolops.RemoveAllHs(mol)
    except Exception:  # pragma: no cover - already H-free or unsanitisable
        pass
    return mol


# Elements that we treat as monatomic ions when they appear as a lone-atom het group.
METALS = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "MN", "FE", "CO", "NI",
    "CU", "ZN", "CD", "HG", "MO", "W", "V", "CR", "PT", "PD", "AU", "AG", "AL",
    "GA", "PB", "SN", "TL", "BE", "SB",
}
HALIDES = {"F", "CL", "BR", "I"}


def build_mol_from_atoms(atoms: list[AtomRecord]) -> Chem.Mol | None:
    """Build an RDKit molecule (conformer + per-atom name + connectivity) from bare PDB het atoms.

    PDB het groups carry no hydrogens, so RDKit bond-*order* perception is unreliable (it fills
    missing valence with spurious multiple bonds). We therefore perceive **connectivity only**
    (single bonds, from geometry) plus ring membership, and label such molecules with the
    bond-order-independent graph rules in :func:`pharmacophores.assign_connectivity`.
    Single-atom groups (ions) need no perception. Returns ``None`` only if nothing can be built.
    """
    if not atoms:
        return None
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(atoms))
    for i, rec in enumerate(atoms):
        element = rec.element.strip().capitalize() or "X"
        try:
            atom = Chem.Atom(element)
        except RuntimeError:
            atom = Chem.Atom(0)
        atom.SetProp("atomName", rec.name)
        atom.SetNoImplicit(True)
        rw.AddAtom(atom)
        conf.SetAtomPosition(i, (float(rec.x), float(rec.y), float(rec.z)))
    mol = rw.GetMol()
    mol.AddConformer(conf, assignId=True)

    if mol.GetNumAtoms() == 1:
        return mol  # monatomic ion: no bonds to perceive

    try:
        rdDetermineBonds.DetermineConnectivity(mol)
        Chem.FastFindRings(mol)  # ring membership without kekulize/aromaticity
    except Exception:
        pass  # leave disconnected; element fallbacks still apply
    return mol


def coords_of(mol: Chem.Mol) -> np.ndarray:
    """Return an (N, 3) array of atom coordinates from the molecule's conformer."""
    return _conformer_coords(mol)
