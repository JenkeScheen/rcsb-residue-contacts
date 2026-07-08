"""Parse PDBBind protein PDBs and detect residue<->molecule contacts by 3D distance.

The target-residue side comes from the protein PDB (via gemmi). The molecule side is any RDKit
molecule with a conformer -- the PDBBind ligand (drug-like) or a solvent het group. Because all
PDBBind files for a complex share one coordinate frame, contacts are simple distance tests.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np

from .molecules import AtomRecord

log = logging.getLogger(__name__)

BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


@dataclass
class ResidueInstance:
    """One instance of the queried residue, with its in-scope heavy atoms."""

    chain: str
    seqid: str
    name: str
    atoms: list[AtomRecord] = field(default_factory=list)
    coords: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))


@dataclass
class HetGroup:
    """A non-water non-polymer het group (candidate solvent molecule)."""

    ccd: str
    chain: str
    seqid: str
    atoms: list[AtomRecord] = field(default_factory=list)


def read_protein(path: Path) -> gemmi.Structure | None:
    """Read a protein PDB and set up entity typing (needed to tell het from polymer)."""
    try:
        structure = gemmi.read_structure(str(path))
        structure.setup_entities()
        return structure
    except Exception as exc:  # pragma: no cover - malformed file
        log.debug("Failed to parse %s: %s", path, exc)
        return None


def _atom_record(res_atom: gemmi.Atom) -> AtomRecord:
    return AtomRecord(
        name=res_atom.name,
        element=res_atom.element.name,
        x=res_atom.pos.x,
        y=res_atom.pos.y,
        z=res_atom.pos.z,
    )


def _is_amino_acid(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return bool(info and info.is_amino_acid())


def residue_instances(
    structure: gemmi.Structure, residue: str, atom_scope: str
) -> list[ResidueInstance]:
    """Collect every instance of ``residue`` with its in-scope heavy atoms."""
    residue = residue.upper()
    instances: list[ResidueInstance] = []
    if not len(structure):
        return instances
    model = structure[0]
    for chain in model:
        for res in chain:
            if res.name != residue:
                continue
            records: list[AtomRecord] = []
            for atom in res:
                if atom.is_hydrogen():
                    continue
                if atom_scope == "sidechain" and atom.name in BACKBONE_ATOMS:
                    continue
                records.append(_atom_record(atom))
            if not records:
                continue
            coords = np.array([[r.x, r.y, r.z] for r in records], dtype=float)
            instances.append(
                ResidueInstance(
                    chain=chain.name,
                    seqid=str(res.seqid),
                    name=residue,
                    atoms=records,
                    coords=coords,
                )
            )
    return instances


def solvent_het_groups(structure: gemmi.Structure) -> list[HetGroup]:
    """Collect non-water, non-polymer het groups (ions, buffers, additives) from the protein."""
    groups: list[HetGroup] = []
    if not len(structure):
        return groups
    model = structure[0]
    for chain in model:
        for res in chain:
            if res.is_water():
                continue
            if res.entity_type != gemmi.EntityType.NonPolymer:
                continue
            if _is_amino_acid(res.name):
                # Free/modified amino acid flagged HETATM: not a solvent molecule.
                continue
            records = [_atom_record(a) for a in res if not a.is_hydrogen()]
            if not records:
                continue
            groups.append(
                HetGroup(ccd=res.name, chain=chain.name, seqid=str(res.seqid), atoms=records)
            )
    return groups


def contacting_atoms(
    residue_coords: np.ndarray, molecule_coords: np.ndarray, cutoff: float
) -> tuple[set[int], float]:
    """Return molecule atom indices within ``cutoff`` of any residue atom, and the min distance.

    Returns ``(set(), inf)`` when there is no contact.
    """
    if residue_coords.size == 0 or molecule_coords.size == 0:
        return set(), float("inf")
    diff = molecule_coords[:, None, :] - residue_coords[None, :, :]
    dists = np.sqrt((diff * diff).sum(axis=-1))
    per_atom_min = dists.min(axis=1)
    hit = np.where(per_atom_min <= cutoff)[0]
    if hit.size == 0:
        return set(), float("inf")
    return {int(i) for i in hit}, float(per_atom_min[hit].min())
