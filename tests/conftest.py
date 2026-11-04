"""Shared fixtures: real RCSB responses for 6PML (a benzofuran ligand stacking face-to-face on
HIS A32), served by an offline fake client so the whole pipeline is testable without network."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from rcsb_residue_contacts.rcsb import NotFound, RcsbClient

FIXTURES = Path(__file__).parent / "fixtures"


def fixture_text(name: str) -> str:
    return (FIXTURES / name).read_text()


class FakeClient(RcsbClient):
    """Serves fixture files instead of hitting RCSB; records every call."""

    def __init__(self, cache_dir=None):
        super().__init__(cache_dir=cache_dir)
        self.calls: list[tuple] = []
        self.entries = {"6PML": json.loads(fixture_text("6PML_entry.json"))}

    def _request(self, url, data=None):  # pragma: no cover - must never be reached
        raise AssertionError(f"unexpected network request: {url}")

    def search(self, query, return_type="entry"):
        self.calls.append(("search", return_type))
        return ["6PML_2"] if return_type == "non_polymer_entity" else ["6PML"]

    def graphql(self, query, variables=None):
        ids = (variables or {}).get("ids", [])
        self.calls.append(("graphql", tuple(ids)))
        return {"entries": [self.entries[i] for i in ids if i in self.entries]}

    def ligand_sdf(self, entry_id, label_asym_id):
        self.calls.append(("ligand", entry_id, label_asym_id))
        try:
            return fixture_text(f"{entry_id.upper()}_{label_asym_id}.sdf")
        except FileNotFoundError as exc:
            raise NotFound(entry_id) from exc

    def residue_atoms_cif(self, entry_id, comp_id):
        self.calls.append(("residues", entry_id, comp_id))
        try:
            return fixture_text(f"{entry_id.upper()}_{comp_id.upper()}.cif")
        except FileNotFoundError as exc:
            raise NotFound(entry_id) from exc

    def pocket_cif(self, entry_id, label_asym_id, radius=6.0):
        self.calls.append(("pocket", entry_id, label_asym_id))
        return fixture_text("6PML_HIS.cif")  # any real mmCIF will do as a stand-in pocket


@pytest.fixture
def fake_client():
    return FakeClient()


@pytest.fixture
def entry_6pml():
    return json.loads(fixture_text("6PML_entry.json"))
