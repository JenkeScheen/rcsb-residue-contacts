"""Aromatic ring perception and pi-stacking geometry (face-to-face vs edge-to-face).

Geometric criteria follow PLIP (Salentin et al., Nucleic Acids Res. 2015) defaults:

* ring centroid distance <= 5.5 A,
* lateral **offset** <= 2.0 A -- the distance between one centroid and the projection of the
  other centroid onto its ring plane (minimum of the two projections),
* inter-plane angle <= 30 deg -> **face-to-face** (parallel / sandwich / parallel-displaced);
  angle >= 60 deg -> **edge-to-face** (T-shaped). Intermediate tilts are not counted as stacks.

Ring pairs that fail the criteria are still reported (up to ``StackCriteria.window``) with
``stack_type == "none"`` so the geometry distribution can be inspected.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from .molecules import METALS


FACE_TO_FACE = "face-to-face"
EDGE_TO_FACE = "edge-to-face"
NO_STACK = "none"


@dataclass(frozen=True)
class StackCriteria:
    max_distance: float = 5.5
    max_parallel_angle: float = 30.0
    min_perpendicular_angle: float = 60.0
    max_offset: float = 2.0
    window: float = 7.0  # ring pairs are recorded up to this centroid distance


@dataclass
class LigandRing:
    """One aromatic ring of a ligand."""

    atoms: tuple[int, ...]
    centroid: np.ndarray
    normal: np.ndarray
    label: str
    fused: bool
    system_smiles: str | None


@dataclass(frozen=True)
class StackGeometry:
    distance: float  # centroid-centroid, A
    angle: float  # between ring planes, deg in [0, 90]
    offset: float  # lateral displacement, A
    vertical: float  # centroid separation along the residue ring normal, A


# ----------------------------------------------------------------------------- geometry
def plane(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centroid and unit normal of the best-fit plane through ``coords`` (N >= 3)."""
    centroid = coords.mean(axis=0)
    _, _, vt = np.linalg.svd(coords - centroid)
    normal = vt[-1]
    return centroid, normal / np.linalg.norm(normal)


def stack_geometry(
    c_res: np.ndarray, n_res: np.ndarray, c_lig: np.ndarray, n_lig: np.ndarray
) -> StackGeometry:
    v = c_lig - c_res
    distance = float(np.linalg.norm(v))
    cos = abs(float(np.dot(n_res, n_lig)))
    angle = float(np.degrees(np.arccos(min(1.0, cos))))
    off_res = float(np.linalg.norm(v - np.dot(v, n_res) * n_res))
    off_lig = float(np.linalg.norm(v - np.dot(v, n_lig) * n_lig))
    return StackGeometry(
        distance=distance,
        angle=angle,
        offset=min(off_res, off_lig),
        vertical=abs(float(np.dot(v, n_res))),
    )


def classify(geom: StackGeometry, criteria: StackCriteria = StackCriteria()) -> str:
    """``face-to-face``, ``edge-to-face`` or ``none`` for one ring pair."""
    if geom.distance > criteria.max_distance or geom.offset > criteria.max_offset:
        return NO_STACK
    if geom.angle <= criteria.max_parallel_angle:
        return FACE_TO_FACE
    if geom.angle >= criteria.min_perpendicular_angle:
        return EDGE_TO_FACE
    return NO_STACK


# --------------------------------------------------------------------------- ring types
def _cyclic_gap(i: int, j: int, n: int) -> int:
    d = abs(i - j)
    return min(d, n - d)


def ring_label(mol: Chem.Mol, ring: tuple[int, ...]) -> str:
    """Name an aromatic ring by size and heteroatom pattern (e.g. ``pyridine``, ``thiazole``)."""
    syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
    size = len(syms)
    het = [(k, s) for k, s in enumerate(syms) if s != "C"]
    kinds = sorted({s for _, s in het})
    n_n = syms.count("N")

    def gap(a: int, b: int) -> int:
        return _cyclic_gap(het[a][0], het[b][0], size)

    if size == 6:
        if not het:
            return "benzene"
        if kinds == ["N"]:
            if n_n == 1:
                return "pyridine"
            if n_n == 2:
                return {1: "pyridazine", 2: "pyrimidine", 3: "pyrazine"}[gap(0, 1)]
            if n_n == 3:
                return "triazine"
        return "other 6-ring heteroaromatic"
    if size == 5:
        if kinds == ["S"] and len(het) == 1:
            return "thiophene"
        if kinds == ["O"] and len(het) == 1:
            return "furan"
        if kinds == ["N"]:
            if n_n == 1:
                return "pyrrole"
            if n_n == 2:
                return "pyrazole" if gap(0, 1) == 1 else "imidazole"
            if n_n == 3:
                return "triazole"
            if n_n == 4:
                return "tetrazole"
        if kinds in (["N", "S"], ["N", "O"]):
            x = "S" if "S" in kinds else "O"
            stem = {"S": ("thiazole", "isothiazole", "thiadiazole"),
                    "O": ("oxazole", "isoxazole", "oxadiazole")}[x]
            if n_n == 1 and len(het) == 2:
                return stem[1] if gap(0, 1) == 1 else stem[0]
            if n_n == 2 and len(het) == 3:
                return stem[2]
        return "other 5-ring heteroaromatic"
    return f"other {size}-ring aromatic"


def _ring_systems(rings: list[tuple[int, ...]]) -> list[set[int]]:
    """Group rings that share a bond (>= 2 atoms) into fused systems; returns ring-index sets."""
    parent = list(range(len(rings)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    sets = [set(r) for r in rings]
    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            if len(sets[i] & sets[j]) >= 2:
                parent[find(i)] = find(j)
    groups: dict[int, set[int]] = {}
    for i in range(len(rings)):
        groups.setdefault(find(i), set()).add(i)
    return list(groups.values())


def ring_system_smiles(mol: Chem.Mol, atoms: set[int]) -> str | None:
    """SMILES of an aromatic ring system with substituents replaced by hydrogen.

    Exocyclic double-bonded atoms (e.g. the carbonyl O of a pyridone) are kept so that the
    fragment stays a valid aromatic system. Every other cut bond (single, or aromatic -- e.g. the
    methine bridges of a porphyrin) is capped with hydrogen, bonds to metals are cut without
    capping, and the fragment is then neutralised (pyridinium -> pyridine, pyrrolide -> pyrrole)
    and put in RDKit's canonical tautomer, so one ring system gets one name whatever its deposited
    charge or tautomer. Returns ``None`` if the fragment cannot be sanitised.
    """
    keep = set(atoms)
    for i in atoms:
        for bond in mol.GetAtomWithIdx(i).GetBonds():
            j = bond.GetOtherAtomIdx(i)
            if j not in keep and bond.GetBondType() == Chem.BondType.DOUBLE:
                keep.add(j)
    capped = (Chem.BondType.SINGLE, Chem.BondType.AROMATIC)
    rw = Chem.RWMol(mol)
    for i in keep:
        atom = rw.GetAtomWithIdx(i)
        cut = sum(
            1
            for b in atom.GetBonds()
            if b.GetOtherAtomIdx(i) not in keep
            and b.GetBondType() in capped
            and b.GetOtherAtom(atom).GetSymbol().upper() not in METALS
        )
        atom.SetNumExplicitHs(atom.GetTotalNumHs() + cut)
        atom.SetNoImplicit(True)
    for i in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        rw.RemoveAtom(i)
    frag = rw.GetMol()
    try:
        Chem.SanitizeMol(frag)
        frag = rdMolStandardize.Uncharger().uncharge(frag)  # per call: thread-safe
        # One name per ring system regardless of deposited tautomer (e.g. 7H/9H purines).
        frag = rdMolStandardize.TautomerEnumerator().Canonicalize(frag)
        return Chem.MolToSmiles(frag)
    except Exception:
        return None


def ligand_rings(mol: Chem.Mol) -> list[LigandRing]:
    """All 5/6-membered fully aromatic rings of ``mol`` with geometry and pattern labels."""
    if mol.GetNumConformers() == 0:
        return []
    rings = [
        tuple(r)
        for r in mol.GetRingInfo().AtomRings()
        if len(r) in (5, 6) and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)
    ]
    if not rings:
        return []
    pos = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    systems = _ring_systems(rings)
    system_of = {ri: k for k, group in enumerate(systems) for ri in group}
    system_smiles: dict[int, str | None] = {}
    out: list[LigandRing] = []
    for ri, ring in enumerate(rings):
        k = system_of[ri]
        if k not in system_smiles:
            atoms = set().union(*(rings[j] for j in systems[k]))
            system_smiles[k] = ring_system_smiles(mol, atoms)
        centroid, normal = plane(pos[list(ring)])
        out.append(
            LigandRing(
                atoms=ring,
                centroid=centroid,
                normal=normal,
                label=ring_label(mol, ring),
                fused=len(systems[k]) > 1,
                system_smiles=system_smiles[k],
            )
        )
    return out
