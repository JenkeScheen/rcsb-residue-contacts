"""Offline end-to-end run on real RCSB data for 6PML (served by the fake client)."""

import json

import pandas as pd
import pytest

from rcsb_residue_contacts.pipeline import Config, run
from rcsb_residue_contacts.rcsb import NotFound

EXPECTED_FILES = [
    "ligands.csv", "stacking.csv", "contacts.csv", "face_to_face_hits.csv", "run_info.json",
    "examples/index.csv",
    "plots/face_to_face_ring_patterns.png", "plots/face_to_face_ring_systems.png",
    "plots/stacking_geometry.png", "plots/interaction_modes.png",
    "plots/face_to_face_by_year.png", "plots/contact_pharmacophores.png",
]


@pytest.mark.parametrize("entry_ids", [None, ["6pml"]])
def test_end_to_end_6pml(tmp_path, fake_client, entry_ids):
    cfg = Config(entry_ids=entry_ids, out_dir=tmp_path, cache_dir=None, workers=2)
    tables = run(cfg, client=fake_client)
    for name in EXPECTED_FILES:
        assert (tmp_path / name).exists(), name

    ligands = tables["ligands"]
    assert ligands[["entry_id", "asym_id", "comp_id", "status"]].values.tolist() == [
        ["6PML", "D", "OR4", "ok"]
    ]
    hits = tables["hits"]
    assert len(hits) == 1
    hit = hits.iloc[0]
    assert hit["residue"] == "HISA32"
    assert hit["pure_stacking"]
    assert hit["ring_system_smiles"] == "c1ccc2occc2c1"  # benzofuran
    assert 3.3 < hit["distance"] < 4.2 and hit["angle"] < 10

    stacks = tables["stacks"]
    assert set(stacks["stack_type"]) <= {"face-to-face", "edge-to-face", "none"}
    assert (stacks["distance"] > 3.0).all()
    assert not tables["contacts"].empty

    index = pd.read_csv(tmp_path / "examples" / "index.csv")
    assert set(index["entry_id"]) == {"6PML"}
    assert index["file"].str.endswith(".pdb").all()
    assert (tmp_path / "examples" / index["file"].iloc[0]).read_text().startswith("REMARK 999")
    info = json.loads((tmp_path / "run_info.json").read_text())
    assert info["face_to_face_entries"] == 1
    if entry_ids is None:
        assert info["mode"] == "search" and info["matching_entries"] == 1


def test_residue_without_ring_is_rejected(tmp_path, fake_client):
    with pytest.raises(RuntimeError, match="no aromatic ring"):
        run(Config(residue="ALA", out_dir=tmp_path, cache_dir=None), client=fake_client)


def test_missing_ligand_coordinates_are_reported(tmp_path, fake_client):
    def no_coordinates(entry_id, asym):
        raise NotFound(f"{entry_id}/{asym}")

    fake_client.ligand_sdf = no_coordinates
    tables = run(Config(entry_ids=["6PML"], out_dir=tmp_path, cache_dir=None), client=fake_client)
    assert tables["ligands"]["status"].tolist() == ["ligand-download-failed"]
    assert tables["hits"].empty
