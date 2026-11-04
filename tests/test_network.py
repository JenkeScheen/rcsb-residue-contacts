"""Live RCSB smoke test. Skipped unless RCSB_NETWORK_TESTS=1 (the rest of the suite is offline)."""

import os

import pytest

from rcsb_residue_contacts.pipeline import Config, run

pytestmark = [
    pytest.mark.network,
    pytest.mark.skipif(
        os.environ.get("RCSB_NETWORK_TESTS") != "1", reason="set RCSB_NETWORK_TESTS=1 to run"
    ),
]


def test_live_6pml(tmp_path):
    tables = run(Config(entry_ids=["6PML"], out_dir=tmp_path, cache_dir=None, workers=2))
    assert "HISA32" in set(tables["hits"]["residue"])
