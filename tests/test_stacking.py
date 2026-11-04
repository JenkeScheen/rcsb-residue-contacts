"""Ring geometry, stacking classification, ring labels and ring-system extraction."""

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from conftest import fixture_text
from rcsb_residue_contacts import molecules, stacking
from rcsb_residue_contacts.stacking import (
    EDGE_TO_FACE, FACE_TO_FACE, NO_STACK, StackCriteria, classify, plane, stack_geometry,
)


def hexagon(center=(0.0, 0.0, 0.0), normal_axis="z", radius=1.39):
    t = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    a, b = radius * np.cos(t), radius * np.sin(t)
    zero = np.zeros(6)
    pts = {"z": (a, b, zero), "x": (zero, a, b), "y": (a, zero, b)}[normal_axis]
    return np.stack(pts, axis=1) + np.asarray(center)


def geom(ring_a, ring_b):
    return stack_geometry(*plane(ring_a), *plane(ring_b))


def test_plane_normal_and_centroid():
    c, n = plane(hexagon(center=(1, 2, 3)))
    assert np.allclose(c, [1, 2, 3])
    assert np.isclose(abs(n[2]), 1.0)


def test_sandwich_is_face_to_face():
    g = geom(hexagon(), hexagon(center=(0, 0, 3.6)))
    assert g.distance == pytest.approx(3.6)
    assert g.angle == pytest.approx(0.0, abs=1e-6)
    assert g.offset == pytest.approx(0.0, abs=1e-6)
    assert g.vertical == pytest.approx(3.6)
    assert classify(g) == FACE_TO_FACE


def test_parallel_displaced_within_offset_is_face_to_face():
    g = geom(hexagon(), hexagon(center=(1.5, 0, 3.5)))
    assert g.offset == pytest.approx(1.5)
    assert classify(g) == FACE_TO_FACE


def test_too_much_offset_is_not_stacked():
    g = geom(hexagon(), hexagon(center=(3.0, 0, 3.5)))
    assert classify(g) == NO_STACK


def test_t_shaped_is_edge_to_face():
    # Ring B stands on edge above ring A's face (normals perpendicular).
    g = geom(hexagon(), hexagon(center=(0, 0, 4.8), normal_axis="x"))
    assert g.angle == pytest.approx(90.0)
    assert classify(g) == EDGE_TO_FACE


def test_tilted_45_is_not_stacked():
    a = hexagon()
    theta = np.radians(45)
    rot = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
    b = hexagon() @ rot.T + np.array([0, 0, 3.8])
    g = geom(a, b)
    assert g.angle == pytest.approx(45.0, abs=0.1)
    assert classify(g) == NO_STACK


def test_too_far_is_not_stacked_and_criteria_are_tunable():
    g = geom(hexagon(), hexagon(center=(0, 0, 6.0)))
    assert classify(g) == NO_STACK
    assert classify(g, StackCriteria(max_distance=6.5)) == FACE_TO_FACE


@pytest.mark.parametrize(
    "smiles,label",
    [
        ("c1ccccc1", "benzene"), ("c1ccncc1", "pyridine"), ("c1cncnc1", "pyrimidine"),
        ("c1cnccn1", "pyrazine"), ("c1ccnnc1", "pyridazine"), ("c1ncncn1", "triazine"),
        ("c1cc[nH]c1", "pyrrole"), ("c1c[nH]cn1", "imidazole"), ("c1cn[nH]c1", "pyrazole"),
        ("c1nc[nH]n1", "triazole"), ("c1nn[nH]n1", "tetrazole"), ("c1ccsc1", "thiophene"),
        ("c1ccoc1", "furan"), ("c1cscn1", "thiazole"), ("c1cnsc1", "isothiazole"),
        ("c1cocn1", "oxazole"), ("c1cnoc1", "isoxazole"), ("c1nnco1", "oxadiazole"),
        ("c1nncs1", "thiadiazole"),
    ],
)
def test_ring_labels(smiles, label):
    mol = Chem.MolFromSmiles(smiles)
    ring = mol.GetRingInfo().AtomRings()[0]
    assert stacking.ring_label(mol, ring) == label


def embedded(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=7) == 0
    return Chem.RemoveHs(mol)


def test_fused_ring_system_strips_substituents():
    rings = stacking.ligand_rings(embedded("CC(=O)Nc1ccc2[nH]c(C)cc2c1"))  # substituted indole
    assert {r.label for r in rings} == {"benzene", "pyrrole"}
    assert all(r.fused for r in rings)
    assert {r.system_smiles for r in rings} == {"c1ccc2[nH]ccc2c1"}


def test_biaryl_rings_are_separate_systems():
    rings = stacking.ligand_rings(embedded("c1ccc(-c2ccncc2)cc1"))
    assert sorted((r.label, r.fused, r.system_smiles) for r in rings) == [
        ("benzene", False, "c1ccccc1"), ("pyridine", False, "c1ccncc1"),
    ]


def test_charged_and_exocyclic_patterns_are_normalised():
    rings = stacking.ligand_rings(embedded("C[n+]1ccccc1"))  # N-methylpyridinium
    assert rings[0].system_smiles == "c1ccncc1"
    rings = stacking.ligand_rings(embedded("CCn1ccccc1=O"))  # N-alkyl pyridone keeps its C=O
    assert rings[0].system_smiles == "O=c1cccc[nH]1"


def test_tautomers_share_one_ring_system_name():
    seven_h = stacking.ligand_rings(embedded("O=c1[nH]cnc2[nH]cnc12"))  # hypoxanthine, 7H
    nine_h = stacking.ligand_rings(embedded("O=c1[nH]cnc2nc[nH]c12"))  # hypoxanthine, 9H
    assert {r.system_smiles for r in seven_h} == {r.system_smiles for r in nine_h}


def test_heme_pyrroles_from_rcsb_sdf():
    mol = molecules.mol_from_sdf(fixture_text("1XJ6_G.sdf"))
    rings = stacking.ligand_rings(mol)
    assert rings and {r.system_smiles for r in rings} == {"c1cc[nH]c1"}


def test_molecule_without_aromatic_rings():
    assert stacking.ligand_rings(embedded("C1CCCCC1O")) == []
