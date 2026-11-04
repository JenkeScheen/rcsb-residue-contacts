"""Build RDKit molecules from ModelServer ligand SDF.

The RCSB ModelServer ``/ligand`` query writes the deposited (crystal) pose of a ligand instance
with bond orders and formal charges from the Chemical Component Dictionary, so RDKit can sanitise
it directly -- no bond-order perception from geometry is needed. Disordered atoms that were not
modelled are simply absent.
"""

from __future__ import annotations

import logging

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

# Elements treated as metal ions (pharmacophore fallback + metal-coordination flag).
METALS = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "MN", "FE", "CO", "NI",
    "CU", "ZN", "CD", "HG", "MO", "W", "V", "CR", "PT", "PD", "AU", "AG", "AL",
    "GA", "PB", "SN", "TL", "BE", "SB", "RU", "RH", "IR", "OS", "RE",
}
HALIDES = {"F", "CL", "BR", "I"}


# (element, explicit valence) -> formal charge, for hypervalent atoms whose charge is missing or
# wrong in the SDF (seen in ModelServer output: zwitterion NH3, tetrahedral boronate written as
# B+ or B0, metal-bound O). Only applied after RDKit has rejected the molecule.
_EXPECTED_CHARGE = {("N", 4): 1, ("B", 4): -1, ("O", 3): 1, ("S", 3): 1}


def _explicit_valence(atom: Chem.Atom) -> int:
    try:
        return atom.GetValence(Chem.ValenceType.EXPLICIT)
    except AttributeError:  # pragma: no cover - RDKit < 2024.09
        return atom.GetExplicitValence()


def _repair_charges(text: str) -> Chem.Mol | None:
    """Re-parse without sanitisation, fix the formal charges of hypervalent atoms, sanitise."""
    mol = Chem.MolFromMolBlock(text, removeHs=False, sanitize=False)
    if mol is None:
        return None
    mol.UpdatePropertyCache(strict=False)
    for atom in mol.GetAtoms():
        charge = _EXPECTED_CHARGE.get((atom.GetSymbol(), _explicit_valence(atom)))
        if charge is not None and atom.GetFormalCharge() != charge:
            atom.SetFormalCharge(charge)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def mol_from_sdf(text: str) -> Chem.Mol | None:
    """Parse one SDF record into a sanitised, heavy-atom-only RDKit molecule with a conformer.

    If RDKit rejects the record because a hypervalent atom lacks its formal charge, the charge is
    repaired (see ``_EXPECTED_CHARGE``). Returns ``None`` if the record is empty or unusable.
    """
    if not text or not text.strip():
        return None
    mol = Chem.MolFromMolBlock(text, removeHs=False, sanitize=True)
    if mol is None:
        mol = _repair_charges(text)
    if mol is None or mol.GetNumAtoms() == 0 or mol.GetNumConformers() == 0:
        return None
    # Contacts are heavy-atom to heavy-atom; RemoveAllHs also drops H that RemoveHs keeps.
    try:
        mol = Chem.RemoveAllHs(mol)
    except Exception:  # pragma: no cover - already H-free
        pass
    return mol


def coords_of(mol: Chem.Mol) -> np.ndarray:
    """(N, 3) array of conformer coordinates."""
    return np.asarray(mol.GetConformer().GetPositions(), dtype=float).reshape(-1, 3)


def element_coords(mol: Chem.Mol, elements: set[str]) -> np.ndarray:
    """Coordinates of atoms whose (upper-case) element symbol is in ``elements``."""
    pos = coords_of(mol)
    idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol().upper() in elements]
    return pos[idx].reshape(-1, 3)


def to_smiles(mol: Chem.Mol) -> str:
    """Canonical SMILES of the modelled ligand (may be partial if atoms were not modelled)."""
    try:
        return Chem.MolToSmiles(mol)
    except Exception:  # pragma: no cover
        return ""
