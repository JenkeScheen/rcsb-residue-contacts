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
    "GA", "PB", "SN", "TL", "BE", "SB", "RU", "RH", "IR", "OS", "RE", "IN", "BI",
    "TI", "ZR", "HF", "NB", "TA", "TC", "SC", "Y", "LA", "CE", "PR", "ND", "SM", "EU",
    "GD", "TB", "DY", "HO", "ER", "YB", "LU", "U",
}
HALIDES = {"F", "CL", "BR", "I"}


# (element, explicit valence) -> formal charge, for atoms whose charge is missing or wrong in
# the SDF. Seen in ModelServer output: zwitterion NH3 without a charge, tetrahedral boronate
# written as B+ or B0, and carboxylate/phosphate O written as O- *with* its H (heme propionates,
# cobalamin). Only applied after RDKit has rejected the molecule.
_EXPECTED_CHARGE = {("N", 4): 1, ("B", 4): -1, ("O", 2): 0, ("O", 3): 1, ("S", 3): 1}


def _explicit_valence(atom: Chem.Atom) -> int:
    try:
        return atom.GetValence(Chem.ValenceType.EXPLICIT)
    except AttributeError:  # pragma: no cover - RDKit < 2024.09
        return atom.GetExplicitValence()


def _is_metal(atom: Chem.Atom) -> bool:
    return atom.GetSymbol().upper() in METALS


def _dative_metal_bonds(mol: Chem.Mol) -> Chem.Mol:
    """Rewrite every metal--non-metal bond as a dative bond from the non-metal to the metal.

    The CCD draws coordination bonds as ordinary covalent bonds, which over-fills the valence of
    the donor (Cl- on Ru, S- on Mo, eta-6 arene carbons, ...). Dative bonds do not count towards
    the donor's valence, so the organic part of the complex sanitises as written.
    """
    rw = Chem.RWMol(mol)
    for bond in list(mol.GetBonds()):
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if _is_metal(a) == _is_metal(b):
            continue
        donor, metal = (b, a) if _is_metal(a) else (a, b)
        rw.RemoveBond(a.GetIdx(), b.GetIdx())
        rw.AddBond(donor.GetIdx(), metal.GetIdx(), Chem.BondType.DATIVE)
    return rw.GetMol()


def _repair(text: str) -> Chem.Mol | None:
    """Re-parse without sanitisation, make metal bonds dative, fix formal charges, sanitise."""
    mol = Chem.MolFromMolBlock(text, removeHs=False, sanitize=False)
    if mol is None:
        return None
    if any(_is_metal(a) for a in mol.GetAtoms()):
        mol = _dative_metal_bonds(mol)
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

    If RDKit rejects the record, it is repaired once: bonds to metals become dative bonds and
    formal charges that contradict an atom's valence are fixed (see ``_EXPECTED_CHARGE``).
    Returns ``None`` if the record is empty or still unusable (e.g. carborane boron clusters).
    """
    if not text or not text.strip():
        return None
    mol = Chem.MolFromMolBlock(text, removeHs=False, sanitize=True)
    if mol is None:
        mol = _repair(text)
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
