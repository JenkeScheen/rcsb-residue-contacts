"""Example structures: ring-system ranking, pick order, PDB output with REMARKs."""

import pandas as pd

from conftest import fixture_text
from rcsb_residue_contacts.examples import (
    group_dir_name, ring_system_examples, slug, to_pdb, write_examples,
)

BENZ, PYR = "c1ccccc1", "c1ccncc1"


def _row(entry, smiles, distance, asym="B", kind="face-to-face", polar=False, comp="LIG"):
    return {
        "entry_id": entry, "asym_id": asym, "comp_id": comp, "ring_system_smiles": smiles,
        "stack_type": kind, "distance": distance, "angle": 5.0, "offset": 0.5,
        "polar_contact": polar, "metal_coordination": False,
        "residue_name": "HIS", "residue_chain": "A", "residue_seqid": "57",
    }


STACKS = pd.DataFrame(
    [
        _row("1AAA", BENZ, 3.9),
        _row("1AAA", BENZ, 3.7, asym="C"),  # same entry/ligand: one pair, closer copy wins
        _row("2BBB", BENZ, 3.4, polar=True),  # closest, but has a polar contact -> ranked last
        _row("3CCC", BENZ, 3.8),
        _row("4DDD", PYR, 4.0),
        _row("5EEE", "c1ccsc1", 3.5, kind="edge-to-face"),  # not face-to-face: ignored
    ]
)


def test_slug_and_group_names():
    assert slug("c1ncc2nc[nH]c2n1") == "c1ncc2nc_nH_c2n1"
    assert slug("***") == "group"
    assert group_dir_name(3, BENZ) == "03_c1ccccc1"


def test_ranking_and_pick_order():
    picks = ring_system_examples(STACKS, per_group=3)
    assert picks.groupby("rank")["ring_system_smiles"].first().to_dict() == {1: BENZ, 2: PYR}
    assert picks.groupby("rank")["n_pairs"].first().to_dict() == {1: 3, 2: 1}
    benz = picks[picks["rank"] == 1]
    # Pure stacks first (closest first), the polar-contact stack last; one pick per entry.
    assert benz["entry_id"].tolist() == ["1AAA", "3CCC", "2BBB"]
    assert benz["asym_id"].tolist() == ["C", "B", "B"]
    assert benz["pure_stacking"].tolist() == [True, True, False]
    assert ring_system_examples(STACKS, per_group=2, max_groups=1)["entry_id"].tolist() == [
        "1AAA", "3CCC"
    ]


def test_to_pdb_adds_remarks_and_atoms():
    pdb = to_pdb(fixture_text("6PML_HIS.cif"), ["REMARK 999 hello"])
    lines = pdb.splitlines()
    assert lines[0] == "REMARK 999 hello"
    assert any(l.startswith("ATOM") and " HIS " in l for l in lines)
    assert to_pdb("not a cif", []) is None


def test_write_examples(tmp_path, fake_client):
    index = write_examples(fake_client, STACKS, tmp_path, per_group=3, max_groups=2)
    assert index["entry_id"].tolist() == ["1AAA", "3CCC", "2BBB", "4DDD"]
    assert index["file"].iloc[0] == "01_c1ccccc1/1AAA_LIG_C__HISA57.pdb"
    text = (tmp_path / index["file"].iloc[0]).read_text()
    assert text.startswith("REMARK 999 face-to-face pi-stack")
    assert "stacks on HISA57: centroid 3.70 A" in text
    assert (tmp_path / "index.csv").exists()
    assert [c[0] for c in fake_client.calls] == ["pocket"] * 4


def test_write_examples_with_nothing_to_write(tmp_path, fake_client):
    index = write_examples(fake_client, STACKS.iloc[5:], tmp_path)
    assert index.empty and (tmp_path / "index.csv").exists()
