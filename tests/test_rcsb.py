"""HTTP client behaviour: caching, empty search responses, permanent vs retryable errors."""

import io
import urllib.error

import pytest

from rcsb_residue_contacts import rcsb
from rcsb_residue_contacts.rcsb import NotFound, RcsbClient, RcsbError


def test_cache_roundtrip_skips_network(tmp_path, monkeypatch):
    client = RcsbClient(cache_dir=tmp_path)
    calls = []
    monkeypatch.setattr(client, "_request", lambda url, data=None: calls.append(url) or b"SDF")
    assert client.ligand_sdf("6pml", "D") == "SDF"
    assert client.ligand_sdf("6PML", "D") == "SDF"
    assert len(calls) == 1
    assert (tmp_path / "ligand" / "6PML_D.sdf").read_bytes() == b"SDF"


def test_no_cache_dir_always_fetches(monkeypatch):
    client = RcsbClient(cache_dir=None)
    calls = []
    monkeypatch.setattr(client, "_request", lambda url, data=None: calls.append(url) or b"x")
    client.residue_atoms_cif("6PML", "HIS")
    client.residue_atoms_cif("6PML", "HIS")
    assert len(calls) == 2


def test_search_handles_no_hits(monkeypatch):
    client = RcsbClient()
    monkeypatch.setattr(client, "_request", lambda url, data=None: b"")
    assert client.search({}) == []


def _http_error(code):
    return urllib.error.HTTPError("u", code, "msg", {}, io.BytesIO(b"body"))


def test_404_is_not_retried(monkeypatch):
    attempts = []

    def fake_urlopen(req, timeout):
        attempts.append(1)
        raise _http_error(404)

    monkeypatch.setattr(rcsb.urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(NotFound):
        RcsbClient(retries=3, backoff=0).ligand_sdf("0XXX", "A")
    assert len(attempts) == 1


def test_5xx_is_retried_then_gives_up(monkeypatch):
    attempts = []

    def fake_urlopen(req, timeout):
        attempts.append(1)
        raise _http_error(503)

    monkeypatch.setattr(rcsb.urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(RcsbError):
        RcsbClient(retries=2, backoff=0).ligand_sdf("6PML", "D")
    assert len(attempts) == 3
