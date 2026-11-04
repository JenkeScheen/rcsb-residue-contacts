"""Find candidate ligand instances in the PDB archive via the RCSB Search and Data APIs.

Discovery is two-staged so that coordinates are only downloaded for ligands that can matter:

1. **Search API** (server-side, whole archive): non-polymer entities that are aromatic
   (``rcsb_chem_comp_info.bond_count_aromatic``), flagged as the *ligand of interest* and have a
   residue of the queried type among their ``rcsb_target_neighbors`` (closest heavy atom
   <= 5 A). Optional resolution / release-date filters.
2. **Data API** (GraphQL, batched): for the selected entries, every ligand instance with its own
   neighbour list, so per-*instance* filtering happens before any ModelServer request.
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass, field

from .rcsb import RcsbClient

ENTRY_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_accession_info { initial_release_date }
    rcsb_entry_info { resolution_combined experimental_method }
    nonpolymer_entities {
      rcsb_nonpolymer_entity_container_identifiers { entity_id }
      nonpolymer_comp {
        chem_comp { id name type formula_weight }
        rcsb_chem_comp_info { atom_count_heavy bond_count_aromatic }
      }
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers {
          asym_id auth_asym_id auth_seq_id comp_id
        }
        rcsb_nonpolymer_instance_validation_score { is_subject_of_investigation }
        rcsb_target_neighbors { target_comp_id target_asym_id target_auth_seq_id distance }
      }
    }
  }
}
"""


def _terminal(attribute: str, operator: str, value) -> dict:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {"attribute": attribute, "operator": operator, "value": value},
    }


def build_search_query(
    residue: str,
    min_aromatic_bonds: int = 5,
    loi_only: bool = True,
    max_resolution: float | None = 3.0,
    released_after: str | None = None,
    released_before: str | None = None,
) -> dict:
    """Search API query for ligand entities neighbouring ``residue`` (see module docstring)."""
    nodes = [_terminal("rcsb_target_neighbors.target_comp_id", "exact_match", residue.upper())]
    if min_aromatic_bonds > 0:
        nodes.append(
            _terminal(
                "rcsb_chem_comp_info.bond_count_aromatic", "greater_or_equal", min_aromatic_bonds
            )
        )
    if loi_only:
        nodes.append(
            _terminal(
                "rcsb_nonpolymer_instance_validation_score.is_subject_of_investigation",
                "exact_match", "Y",
            )
        )
    if max_resolution:
        nodes.append(
            _terminal("rcsb_entry_info.resolution_combined", "less_or_equal", max_resolution)
        )
    if released_after:
        nodes.append(
            _terminal("rcsb_accession_info.initial_release_date", "greater_or_equal", released_after)
        )
    if released_before:
        nodes.append(_terminal("rcsb_accession_info.initial_release_date", "less", released_before))
    return {"type": "group", "logical_operator": "and", "nodes": nodes}


def search_entries(client: RcsbClient, query: dict) -> tuple[list[str], int]:
    """Return (sorted unique entry ids, number of matching ligand entities) for ``query``."""
    entity_ids = client.search(query, return_type="non_polymer_entity")
    entries = sorted({eid.split("_")[0].upper() for eid in entity_ids})
    return entries, len(entity_ids)


def sample_entries(entries: list[str], max_entries: int | None, seed: int = 0) -> list[str]:
    """Reproducible random subset (sorted) of at most ``max_entries`` ids; all if uncapped."""
    if max_entries is None or max_entries < 0 or max_entries >= len(entries):
        return sorted(entries)
    return sorted(random.Random(seed).sample(sorted(entries), max_entries))


def fetch_entry_metadata(
    client: RcsbClient, entry_ids: list[str], batch_size: int = 50
) -> dict[str, dict]:
    """Data API records for ``entry_ids`` keyed by upper-case id (cached per entry).

    Ids unknown to the Data API (e.g. obsoleted entries) are silently absent from the result.
    """
    out: dict[str, dict] = {}
    missing: list[str] = []
    for eid in (e.upper() for e in entry_ids):
        hit = client.cache_get("entry", f"{eid}.json")
        if hit is not None:
            out[eid] = json.loads(hit)
        else:
            missing.append(eid)
    for start in range(0, len(missing), batch_size):
        batch = missing[start : start + batch_size]
        data = client.graphql(ENTRY_QUERY, {"ids": batch})
        for record in data.get("entries") or []:
            if not record:
                continue
            eid = record["rcsb_id"].upper()
            out[eid] = record
            client.cache_put("entry", f"{eid}.json", json.dumps(record).encode())
    return out


@dataclass
class LigandInstance:
    """One non-polymer instance (label_asym_id) worth fetching coordinates for."""

    entry_id: str
    asym_id: str
    auth_asym_id: str
    auth_seq_id: str
    comp_id: str
    name: str
    is_loi: bool
    aromatic_bonds: int
    heavy_atoms: int
    release_date: str
    resolution: float | None
    method: str
    neighbor_residues: list[tuple[str, str]] = field(default_factory=list)


def _first(values):
    return values[0] if values else None


def ligand_instances(
    record: dict,
    residue: str,
    min_aromatic_bonds: int = 5,
    loi_only: bool = True,
    exclude_comps: frozenset[str] = frozenset(),
) -> list[LigandInstance]:
    """Ligand instances of one entry that pass the per-instance filters.

    An instance is kept when its component has >= ``min_aromatic_bonds`` aromatic bonds, it is
    (optionally) flagged as ligand of interest, it is not in ``exclude_comps``, and at least one
    ``residue`` is among its RCSB-computed neighbours.
    """
    residue = residue.upper()
    entry_id = record["rcsb_id"].upper()
    release = ((record.get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10]
    info = record.get("rcsb_entry_info") or {}
    resolution = _first(info.get("resolution_combined") or [])
    method = info.get("experimental_method") or ""

    out: list[LigandInstance] = []
    for entity in record.get("nonpolymer_entities") or []:
        comp = entity.get("nonpolymer_comp") or {}
        chem = comp.get("chem_comp") or {}
        cinfo = comp.get("rcsb_chem_comp_info") or {}
        comp_id = (chem.get("id") or "").upper()
        aromatic = int(cinfo.get("bond_count_aromatic") or 0)
        if aromatic < min_aromatic_bonds or comp_id in exclude_comps:
            continue
        for inst in entity.get("nonpolymer_entity_instances") or []:
            ids = inst.get("rcsb_nonpolymer_entity_instance_container_identifiers") or {}
            scores = inst.get("rcsb_nonpolymer_instance_validation_score") or []
            is_loi = any((s or {}).get("is_subject_of_investigation") == "Y" for s in scores)
            if loi_only and not is_loi:
                continue
            neighbors = [
                (n.get("target_asym_id") or "", str(n.get("target_auth_seq_id") or ""))
                for n in inst.get("rcsb_target_neighbors") or []
                if (n.get("target_comp_id") or "").upper() == residue
            ]
            if not neighbors:
                continue
            out.append(
                LigandInstance(
                    entry_id=entry_id,
                    asym_id=ids.get("asym_id") or "",
                    auth_asym_id=ids.get("auth_asym_id") or "",
                    auth_seq_id=str(ids.get("auth_seq_id") or ""),
                    comp_id=comp_id,
                    name=(chem.get("name") or "").replace("\n", " "),
                    is_loi=is_loi,
                    aromatic_bonds=aromatic,
                    heavy_atoms=int(cinfo.get("atom_count_heavy") or 0),
                    release_date=release,
                    resolution=resolution,
                    method=method,
                    neighbor_residues=neighbors,
                )
            )
    return out
