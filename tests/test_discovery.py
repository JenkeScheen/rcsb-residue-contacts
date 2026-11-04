"""Search-query construction, sampling and per-instance filtering of Data API records."""

from rcsb_residue_contacts import discovery
from rcsb_residue_contacts.rcsb import RcsbClient


def _attrs(query):
    return {n["parameters"]["attribute"]: n["parameters"] for n in query["nodes"]}


def test_default_search_query():
    q = discovery.build_search_query("his")
    attrs = _attrs(q)
    assert q["logical_operator"] == "and"
    assert attrs["rcsb_target_neighbors.target_comp_id"]["value"] == "HIS"
    assert attrs["rcsb_chem_comp_info.bond_count_aromatic"]["value"] == 5
    assert (
        attrs["rcsb_nonpolymer_instance_validation_score.is_subject_of_investigation"]["value"]
        == "Y"
    )
    assert attrs["rcsb_entry_info.resolution_combined"]["value"] == 3.0
    assert "rcsb_accession_info.initial_release_date" not in attrs


def test_optional_filters():
    q = discovery.build_search_query(
        "PHE", min_aromatic_bonds=0, loi_only=False, max_resolution=None,
        released_after="2020-01-01", released_before="2024-01-01",
    )
    nodes = q["nodes"]
    assert [n["parameters"]["attribute"] for n in nodes] == [
        "rcsb_target_neighbors.target_comp_id",
        "rcsb_accession_info.initial_release_date",
        "rcsb_accession_info.initial_release_date",
    ]
    assert [n["parameters"]["operator"] for n in nodes[1:]] == ["greater_or_equal", "less"]


def test_sample_entries_is_reproducible_and_sorted():
    ids = [f"{i}ABC" for i in range(100)]
    a = discovery.sample_entries(ids, 10, seed=3)
    assert a == discovery.sample_entries(list(reversed(ids)), 10, seed=3)
    assert a == sorted(a) and len(a) == 10
    assert discovery.sample_entries(ids, None) == sorted(ids)
    assert discovery.sample_entries(ids, 500) == sorted(ids)


def test_search_entries_dedupes_entity_ids(fake_client):
    entries, n_entities = discovery.search_entries(fake_client, {})
    assert entries == ["6PML"] and n_entities == 1


def test_ligand_instances_filters(entry_6pml):
    insts = discovery.ligand_instances(entry_6pml, "HIS")
    # PEG (no aromatic bonds) and CU are dropped; OR4 is the aromatic ligand of interest.
    assert [(i.comp_id, i.asym_id, i.is_loi) for i in insts] == [("OR4", "D", True)]
    inst = insts[0]
    assert inst.release_date == "2019-11-06" and inst.resolution == 2.0
    assert any(seq == "32" for _, seq in inst.neighbor_residues)
    assert discovery.ligand_instances(entry_6pml, "TRP") == []  # no TRP neighbour
    assert discovery.ligand_instances(entry_6pml, "HIS", exclude_comps=frozenset({"OR4"})) == []
    assert discovery.ligand_instances(entry_6pml, "HIS", min_aromatic_bonds=17) == []


def test_fetch_entry_metadata_uses_cache(tmp_path, fake_client):
    fake_client.cache_dir = tmp_path
    first = discovery.fetch_entry_metadata(fake_client, ["6pml", "0XXX"])
    assert set(first) == {"6PML"}  # unknown ids are silently absent
    n_calls = len(fake_client.calls)
    again = discovery.fetch_entry_metadata(fake_client, ["6PML"])
    assert again == first and len(fake_client.calls) == n_calls  # served from cache


def test_model_urls():
    url = RcsbClient._model_url("6PML", "ligand", label_asym_id="D", model_nums=1, encoding="sdf")
    assert url == "https://models.rcsb.org/v1/6pml/ligand?label_asym_id=D&model_nums=1&encoding=sdf"
