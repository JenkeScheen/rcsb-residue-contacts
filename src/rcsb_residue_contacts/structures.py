"""Parse ModelServer mmCIF subsets (gemmi) into residue instances with ring geometry.

The ModelServer ``/atoms`` query returns every atom of a given residue type for an entry; this
module turns that into :class:`ResidueInstance` objects carrying the in-scope heavy atoms (for
distance contacts), the aromatic ring atom sets (for stacking geometry) and the polar side-chain
atoms (to flag H-bond / metal-coordination contacts). All coordinates are in the deposited
asymmetric-unit frame, shared with the ligand SDF from the same server.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import gemmi
import numpy as np

log = logging.getLogger(__name__)

BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

# Aromatic side-chain rings per residue type: ring label -> atom names (in ring order).
RING_ATOMS: dict[str, dict[str, tuple[str, ...]]] = {
    "HIS": {"imidazole": ("CG", "ND1", "CE1", "NE2", "CD2")},
    "PHE": {"phenyl": ("CG", "CD1", "CE1", "CZ", "CE2", "CD2")},
    "TYR": {"phenol": ("CG", "CD1", "CE1", "CZ", "CE2", "CD2")},
    "TRP": {
        "indole 5-ring": ("CG", "CD1", "NE1", "CE2", "CD2"),
        "indole 6-ring": ("CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"),
    },
}
# Side-chain atoms that can H-bond or coordinate a metal.
POLAR_ATOMS: dict[str, tuple[str, ...]] = {
    "HIS": ("ND1", "NE2"),
    "PHE": (),
    "TYR": ("OH",),
    "TRP": ("NE1",),
}
SUPPORTED_RESIDUES = tuple(RING_ATOMS)


@dataclass
class AtomRecord:
    """A single heavy atom."""

    name: str
    element: str
    x: float
    y: float
    z: float


@dataclass
class ResidueInstance:
    """One instance of the queried residue."""

    chain: str
    seqid: str
    name: str
    atoms: list[AtomRecord] = field(default_factory=list)
    coords: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    rings: dict[str, np.ndarray] = field(default_factory=dict)
    polar_coords: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))

    @property
    def label(self) -> str:
        return f"{self.name}{self.chain}{self.seqid}"


def read_cif_text(text: str) -> gemmi.Structure | None:
    """Parse mmCIF text into a gemmi structure (first alt-conformer kept, hydrogens removed)."""
    try:
        doc = gemmi.cif.read_string(text)
        structure = gemmi.make_structure_from_block(doc.sole_block())
        structure.setup_entities()
        structure.remove_alternative_conformations()
        structure.remove_hydrogens()
        return structure
    except Exception as exc:  # malformed / empty response
        log.debug("Failed to parse mmCIF: %s", exc)
        return None


def _is_polymer(res: gemmi.Residue) -> bool:
    """True for residues of a polymer chain (a free amino acid bound as a ligand is excluded)."""
    if res.entity_type == gemmi.EntityType.Unknown:
        return res.het_flag != "H"
    return res.entity_type == gemmi.EntityType.Polymer


def residue_instances(
    structure: gemmi.Structure, residue: str, atom_scope: str = "sidechain"
) -> list[ResidueInstance]:
    """Every instance of ``residue`` in the first model, with contact atoms, rings, polar atoms.

    Rings with missing atoms (disordered side chains) are omitted from ``rings``.
    """
    residue = residue.upper()
    ring_defs = RING_ATOMS.get(residue, {})
    polar_names = POLAR_ATOMS.get(residue, ())
    instances: list[ResidueInstance] = []
    if structure is None or not len(structure):
        return instances
    for chain in structure[0]:
        for res in chain:
            if res.name != residue or not _is_polymer(res):
                continue
            by_name = {atom.name: atom for atom in res}
            records = [
                AtomRecord(a.name, a.element.name, a.pos.x, a.pos.y, a.pos.z)
                for a in res
                if not (atom_scope == "sidechain" and a.name in BACKBONE_ATOMS)
            ]
            rings = {}
            for label, names in ring_defs.items():
                if all(n in by_name for n in names):
                    rings[label] = np.array([by_name[n].pos.tolist() for n in names], dtype=float)
            polar = [by_name[n].pos.tolist() for n in polar_names if n in by_name]
            instances.append(
                ResidueInstance(
                    chain=chain.name,
                    seqid=str(res.seqid.num) + (res.seqid.icode.strip() or ""),
                    name=residue,
                    atoms=records,
                    coords=np.array([[r.x, r.y, r.z] for r in records], dtype=float).reshape(-1, 3),
                    rings=rings,
                    polar_coords=np.array(polar, dtype=float).reshape(-1, 3),
                )
            )
    return instances


def contacting_atoms(
    residue_coords: np.ndarray, molecule_coords: np.ndarray, cutoff: float
) -> tuple[set[int], float]:
    """Molecule atom indices within ``cutoff`` of any residue atom, and the minimum distance.

    Returns ``(set(), inf)`` when there is no contact.
    """
    if residue_coords.size == 0 or molecule_coords.size == 0:
        return set(), float("inf")
    diff = molecule_coords[:, None, :] - residue_coords[None, :, :]
    per_atom_min = np.sqrt((diff * diff).sum(axis=-1)).min(axis=1)
    hit = np.where(per_atom_min <= cutoff)[0]
    if hit.size == 0:
        return set(), float("inf")
    return {int(i) for i in hit}, float(per_atom_min[hit].min())


def min_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Smallest pairwise distance between two (N, 3) coordinate sets (``inf`` if either empty)."""
    if a.size == 0 or b.size == 0:
        return float("inf")
    diff = a[:, None, :] - b[None, :, :]
    return float(np.sqrt((diff * diff).sum(axis=-1)).min())
