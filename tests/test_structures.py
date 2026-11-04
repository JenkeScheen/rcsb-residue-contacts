"""Parsing ModelServer mmCIF into residue instances with rings and polar atoms."""

import numpy as np

from conftest import fixture_text
from rcsb_residue_contacts import structures


def test_his_instances_have_ring_and_polar_atoms():
    st = structures.read_cif_text(fixture_text("6PML_HIS.cif"))
    his = structures.residue_instances(st, "HIS")
    assert his and all(r.name == "HIS" for r in his)
    a32 = next(r for r in his if r.chain == "A" and r.seqid == "32")
    assert a32.rings["imidazole"].shape == (5, 3)
    assert a32.polar_coords.shape == (2, 3)  # ND1 + NE2
    assert {a.name for a in a32.atoms}.isdisjoint(structures.BACKBONE_ATOMS)


def test_atom_scope_all_keeps_backbone():
    st = structures.read_cif_text(fixture_text("6PML_HIS.cif"))
    his = structures.residue_instances(st, "HIS", atom_scope="all")
    assert any("CA" in {a.name for a in r.atoms} for r in his)


def test_free_histidine_ligand_is_not_a_target_residue():
    # 2ABE has a free HIS bound as a ligand (auth seq 300, non-polymer); the ModelServer /atoms
    # query returns it alongside the protein histidines. It must not be treated as a residue.
    st = structures.read_cif_text(fixture_text("2ABE_HIS.cif"))
    seqids = {r.seqid for r in structures.residue_instances(st, "HIS")}
    assert "300" not in seqids
    assert len(seqids) > 5


def test_unparseable_cif_returns_none():
    assert structures.read_cif_text("this is not mmcif") is None
    assert structures.residue_instances(None, "HIS") == []


def test_contacting_atoms_and_min_distance():
    res = np.array([[0.0, 0.0, 0.0]])
    mol = np.array([[3.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    idxs, d = structures.contacting_atoms(res, mol, 4.0)
    assert idxs == {0} and d == 3.0
    assert structures.contacting_atoms(res, mol, 2.0) == (set(), float("inf"))
    assert structures.min_distance(res, mol) == 3.0
    assert structures.min_distance(np.empty((0, 3)), mol) == float("inf")


def test_supported_residues_all_have_rings():
    for name in structures.SUPPORTED_RESIDUES:
        assert structures.RING_ATOMS[name]
