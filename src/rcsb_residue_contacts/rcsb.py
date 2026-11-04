"""Minimal client for the three RCSB PDB web services this package depends on.

* **Search API** (``search.rcsb.org``) -- find candidate ligands server-side. Its index includes
  ``rcsb_target_neighbors`` (the closest polymer residue atom per ligand instance, <= 5 A), the
  "ligand of interest" flag and chemical-component descriptors, so "aromatic ligands next to a
  histidine" is a single query over the whole archive.
* **Data API** (``data.rcsb.org/graphql``) -- per-entry metadata: ligand instances, their
  neighbouring residues, release date, resolution.
* **ModelServer** (``models.rcsb.org``) -- coordinates. ``/ligand`` returns one ligand instance
  in its deposited (crystal) pose as SDF, with bond orders and formal charges taken from the
  Chemical Component Dictionary; ``/atoms`` and ``/residueSurroundings`` return mmCIF subsets.

Only the standard library is used (no ``requests``). Downloads are cached on disk so re-runs with
different geometric thresholds are fast and offline.
"""

from __future__ import annotations

import json
import logging
import os
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from . import __version__

log = logging.getLogger(__name__)

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
MODELSERVER_URL = "https://models.rcsb.org/v1"
USER_AGENT = f"rcsb-residue-contacts/{__version__} (python-urllib)"

_RETRY_STATUS = {429, 500, 502, 503, 504}


class RcsbError(RuntimeError):
    """A request to an RCSB service failed permanently."""


class NotFound(RcsbError):
    """The requested resource does not exist (HTTP 404)."""


class RcsbClient:
    """Thread-safe HTTP client with retries and an optional on-disk cache.

    ``cache_dir=None`` disables caching. Search results are never cached (the archive is
    updated weekly); Data API entry records and ModelServer coordinates are.
    """

    def __init__(
        self,
        cache_dir: Path | None = None,
        timeout: float = 60.0,
        retries: int = 4,
        backoff: float = 1.5,
    ) -> None:
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        self.timeout = timeout
        self.retries = retries
        self.backoff = backoff

    # ------------------------------------------------------------------ transport
    def _request(self, url: str, data: bytes | None = None) -> bytes:
        headers = {"User-Agent": USER_AGENT, "Accept-Encoding": "identity"}
        if data is not None:
            headers["Content-Type"] = "application/json"
        last_exc: Exception | None = None
        for attempt in range(self.retries + 1):
            req = urllib.request.Request(url, data=data, headers=headers)
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    return resp.read()
            except urllib.error.HTTPError as exc:
                if exc.code == 404:
                    raise NotFound(f"404 for {url}") from exc
                if exc.code not in _RETRY_STATUS:
                    body = exc.read()[:300].decode(errors="replace")
                    raise RcsbError(f"HTTP {exc.code} for {url}: {body}") from exc
                last_exc = exc
            except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
                last_exc = exc
            if attempt < self.retries:
                time.sleep(self.backoff * (2**attempt))
        raise RcsbError(f"giving up on {url}: {last_exc}")

    # ---------------------------------------------------------------------- cache
    def _cache_path(self, kind: str, key: str) -> Path | None:
        if self.cache_dir is None:
            return None
        return self.cache_dir / kind / key

    def cache_get(self, kind: str, key: str) -> bytes | None:
        path = self._cache_path(kind, key)
        if path is not None and path.exists():
            return path.read_bytes()
        return None

    def cache_put(self, kind: str, key: str, payload: bytes) -> None:
        path = self._cache_path(kind, key)
        if path is None:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        # Atomic write: concurrent workers never observe a half-written file.
        fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".tmp-")
        with os.fdopen(fd, "wb") as fh:
            fh.write(payload)
        os.chmod(tmp, 0o644)  # mkstemp creates 0600
        os.replace(tmp, path)

    def _cached_get(self, kind: str, key: str, url: str) -> bytes:
        hit = self.cache_get(kind, key)
        if hit is not None:
            return hit
        payload = self._request(url)
        self.cache_put(kind, key, payload)
        return payload

    # ----------------------------------------------------------------- Search API
    def search(self, query: dict, return_type: str = "entry") -> list[str]:
        """Return every identifier matching ``query`` (experimental structures only)."""
        body = {
            "query": query,
            "return_type": return_type,
            "request_options": {
                "return_all_hits": True,
                "results_content_type": ["experimental"],
            },
        }
        raw = self._request(SEARCH_URL, json.dumps(body).encode())
        if not raw:  # HTTP 204: no hits
            return []
        return [hit["identifier"] for hit in json.loads(raw)["result_set"]]

    # ------------------------------------------------------------------- Data API
    def graphql(self, query: str, variables: dict | None = None) -> dict:
        body = json.dumps({"query": query, "variables": variables or {}}).encode()
        payload = json.loads(self._request(GRAPHQL_URL, body))
        if payload.get("errors"):
            raise RcsbError(f"GraphQL errors: {payload['errors'][:2]}")
        return payload["data"]

    # ---------------------------------------------------------------- ModelServer
    @staticmethod
    def _model_url(entry_id: str, query: str, **params) -> str:
        params = {k: v for k, v in params.items() if v is not None}
        return f"{MODELSERVER_URL}/{entry_id.lower()}/{query}?{urllib.parse.urlencode(params)}"

    def ligand_sdf(self, entry_id: str, label_asym_id: str) -> str:
        """One ligand instance in its deposited pose, as SDF (CCD bond orders + charges)."""
        url = self._model_url(
            entry_id, "ligand", label_asym_id=label_asym_id, model_nums=1, encoding="sdf"
        )
        key = f"{entry_id.upper()}_{label_asym_id}.sdf"
        return self._cached_get("ligand", key, url).decode()

    def residue_atoms_cif(self, entry_id: str, comp_id: str) -> str:
        """Every atom of every ``comp_id`` residue in the entry (first model), as mmCIF."""
        url = self._model_url(
            entry_id, "atoms", label_comp_id=comp_id, model_nums=1, encoding="cif"
        )
        key = f"{entry_id.upper()}_{comp_id.upper()}.cif"
        return self._cached_get("residues", key, url).decode()

    def pocket_cif(self, entry_id: str, label_asym_id: str, radius: float = 6.0) -> str:
        """All residues within ``radius`` A of a ligand instance (incl. the ligand), as mmCIF."""
        url = self._model_url(
            entry_id, "residueSurroundings", label_asym_id=label_asym_id,
            radius=radius, model_nums=1, encoding="cif",
        )
        key = f"{entry_id.upper()}_{label_asym_id}_r{radius:g}.cif"
        return self._cached_get("pocket", key, url).decode()

