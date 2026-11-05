"""End-to-end pipeline: find ligands stacking on a residue type across the PDB via RCSB APIs.

1. Search API -> candidate entries (aromatic ligand of interest next to the residue type).
2. Data API   -> per-instance metadata and neighbour lists (batched GraphQL).
3. ModelServer -> ligand instance SDF + all atoms of the residue type, per entry.
4. Local geometry -> ring-ring stacking classification, polar/metal flags, pharmacophore contacts.
5. Tables, plots and example pocket structures under ``out_dir``.
"""

from __future__ import annotations

import json
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from . import __version__, discovery, molecules, pharmacophores, plots, stacking, structures
from .examples import write_examples
from .rcsb import RcsbClient, RcsbError
from .stacking import EDGE_TO_FACE, FACE_TO_FACE, StackCriteria

log = logging.getLogger(__name__)

LIGAND_COLUMNS = [
    "entry_id", "asym_id", "auth_asym_id", "auth_seq_id", "comp_id", "name", "is_loi",
    "release_date", "resolution", "method", "heavy_atoms_ccd", "heavy_atoms_modelled",
    "aromatic_rings", "residues_in_contact", "polar_contact_residues", "metal_contact_residues",
    "face_to_face", "edge_to_face", "status", "smiles",
]
STACK_COLUMNS = [
    "entry_id", "asym_id", "comp_id", "is_loi", "residue_name", "residue_chain", "residue_seqid",
    "residue_ring", "ring_pattern", "ring_type", "fused", "ring_system_smiles", "ring_size",
    "distance", "angle", "offset", "vertical", "stack_type", "polar_contact",
    "metal_coordination", "residue_polar_min_distance", "release_date", "resolution",
]
CONTACT_COLUMNS = [
    "entry_id", "asym_id", "comp_id", "residue_name", "residue_chain", "residue_seqid",
    "pharmacophore", "n_contact_atoms", "min_distance",
]
HIT_COLUMNS = [
    "entry_id", "comp_id", "name", "asym_id", "residue", "ring_patterns", "ring_system_smiles",
    "distance", "angle", "offset", "pure_stacking", "polar_contact", "metal_coordination",
    "release_date", "resolution", "url",
]


@dataclass
class Config:
    residue: str = "HIS"
    max_entries: int | None = 300
    seed: int = 0
    entry_ids: list[str] | None = None
    released_after: str | None = None
    released_before: str | None = None
    max_resolution: float | None = 3.0
    loi_only: bool = True
    min_aromatic_bonds: int = 5
    exclude_comps: frozenset[str] = frozenset()
    criteria: StackCriteria = field(default_factory=StackCriteria)
    contact_cutoff: float = 4.0
    atom_scope: str = "sidechain"
    polar_cutoff: float = 3.5
    metal_cutoff: float = 2.8
    out_dir: Path = Path("results")
    cache_dir: Path | None = Path(".rcsb_cache")
    workers: int = 8
    examples_per_group: int = 3
    max_groups: int = 30
    pocket_radius: float = 6.0


@dataclass
class EntryResult:
    ligands: list[dict] = field(default_factory=list)
    stacks: list[dict] = field(default_factory=list)
    contacts: list[dict] = field(default_factory=list)


def analyse_instance(
    inst: discovery.LigandInstance,
    mol,
    residues: list[structures.ResidueInstance],
    cfg: Config,
) -> tuple[dict, list[dict], list[dict]]:
    """Stacking + contacts between one ligand instance and every residue instance."""
    rings = stacking.ligand_rings(mol)
    coords = molecules.coords_of(mol)
    hetero = molecules.element_coords(mol, {"N", "O"})
    metals = molecules.element_coords(mol, molecules.METALS)
    label_map = pharmacophores.atom_label_map(mol)
    stacks: list[dict] = []
    contacts: list[dict] = []
    in_contact = polar_residues = metal_residues = 0
    for res in residues:
        if structures.min_distance(res.coords, coords) > max(cfg.criteria.window, cfg.contact_cutoff):
            continue
        idxs, mind = structures.contacting_atoms(res.coords, coords, cfg.contact_cutoff)
        if idxs:
            in_contact += 1
            for label in sorted(pharmacophores.assign_with_map(mol, label_map, idxs)):
                contacts.append(
                    {
                        "entry_id": inst.entry_id, "asym_id": inst.asym_id,
                        "comp_id": inst.comp_id, "residue_name": res.name,
                        "residue_chain": res.chain, "residue_seqid": res.seqid,
                        "pharmacophore": label, "n_contact_atoms": len(idxs),
                        "min_distance": round(mind, 3),
                    }
                )
        polar_d = structures.min_distance(res.polar_coords, hetero)
        metal_d = structures.min_distance(res.polar_coords, metals)
        polar_residues += polar_d <= cfg.polar_cutoff
        metal_residues += metal_d <= cfg.metal_cutoff
        if not rings or not res.rings:
            continue
        for res_ring, ring_coords in res.rings.items():
            c_res, n_res = stacking.plane(ring_coords)
            for ring in rings:
                geom = stacking.stack_geometry(c_res, n_res, ring.centroid, ring.normal)
                if geom.distance > cfg.criteria.window:
                    continue
                stacks.append(
                    {
                        "entry_id": inst.entry_id, "asym_id": inst.asym_id,
                        "comp_id": inst.comp_id, "is_loi": inst.is_loi,
                        "residue_name": res.name, "residue_chain": res.chain,
                        "residue_seqid": res.seqid, "residue_ring": res_ring,
                        "ring_pattern": ring.label + (" (fused)" if ring.fused else ""),
                        "ring_type": ring.label, "fused": ring.fused,
                        "ring_system_smiles": ring.system_smiles or "",
                        "ring_size": len(ring.atoms),
                        "distance": round(geom.distance, 3), "angle": round(geom.angle, 1),
                        "offset": round(geom.offset, 3), "vertical": round(geom.vertical, 3),
                        "stack_type": stacking.classify(geom, cfg.criteria),
                        "polar_contact": polar_d <= cfg.polar_cutoff,
                        "metal_coordination": metal_d <= cfg.metal_cutoff,
                        "residue_polar_min_distance": (
                            round(polar_d, 3) if polar_d != float("inf") else None
                        ),
                        "release_date": inst.release_date, "resolution": inst.resolution,
                    }
                )
    kinds = [s["stack_type"] for s in stacks]
    ligand_row = {
        "heavy_atoms_modelled": mol.GetNumAtoms(),
        "aromatic_rings": len(rings),
        "residues_in_contact": in_contact,
        "polar_contact_residues": polar_residues,
        "metal_contact_residues": metal_residues,
        "face_to_face": kinds.count(FACE_TO_FACE),
        "edge_to_face": kinds.count(EDGE_TO_FACE),
        "status": "ok",
        "smiles": molecules.to_smiles(mol),
    }
    return ligand_row, stacks, contacts


def analyse_entry(client: RcsbClient, cfg: Config, record: dict) -> EntryResult:
    """Fetch coordinates for one entry's candidate ligands and analyse them."""
    result = EntryResult()
    instances = discovery.ligand_instances(
        record, cfg.residue, cfg.min_aromatic_bonds, cfg.loi_only, cfg.exclude_comps
    )
    if not instances:
        return result
    entry_id = instances[0].entry_id
    residues: list[structures.ResidueInstance] = []
    residue_status = "ok"
    try:
        structure = structures.read_cif_text(client.residue_atoms_cif(entry_id, cfg.residue))
        residues = structures.residue_instances(structure, cfg.residue, cfg.atom_scope)
    except RcsbError as exc:
        log.debug("%s residue download failed: %s", entry_id, exc)
        residue_status = "residue-download-failed"
    for inst in instances:
        row = {
            "entry_id": inst.entry_id, "asym_id": inst.asym_id,
            "auth_asym_id": inst.auth_asym_id, "auth_seq_id": inst.auth_seq_id,
            "comp_id": inst.comp_id, "name": inst.name, "is_loi": inst.is_loi,
            "release_date": inst.release_date, "resolution": inst.resolution,
            "method": inst.method, "heavy_atoms_ccd": inst.heavy_atoms,
        }
        if residue_status != "ok":
            result.ligands.append({**row, "status": residue_status})
            continue
        try:
            mol = molecules.mol_from_sdf(client.ligand_sdf(entry_id, inst.asym_id))
        except RcsbError as exc:
            log.debug("%s/%s ligand download failed: %s", entry_id, inst.asym_id, exc)
            result.ligands.append({**row, "status": "ligand-download-failed"})
            continue
        if mol is None:
            result.ligands.append({**row, "status": "unparseable"})
            continue
        lig_row, stacks, contacts = analyse_instance(inst, mol, residues, cfg)
        result.ligands.append({**row, **lig_row})
        result.stacks.extend(stacks)
        result.contacts.extend(contacts)
    return result


def face_to_face_hits(stacks: pd.DataFrame, ligands: pd.DataFrame) -> pd.DataFrame:
    """One row per (entry, ligand instance, residue) with a face-to-face stack (closest ring pair)."""
    f2f = stacks[stacks["stack_type"] == FACE_TO_FACE]
    if f2f.empty:
        return pd.DataFrame(columns=HIT_COLUMNS)
    keys = ["entry_id", "asym_id", "residue_name", "residue_chain", "residue_seqid"]
    patterns = f2f.groupby(keys)["ring_pattern"].agg(lambda s: "; ".join(sorted(set(s))))
    best = f2f.sort_values("distance").drop_duplicates(keys).set_index(keys)
    best["ring_patterns"] = patterns
    best = best.reset_index()
    names = ligands.drop_duplicates(["entry_id", "asym_id"]).set_index(["entry_id", "asym_id"])["name"]
    best["name"] = [names.get((e, a), "") for e, a in zip(best["entry_id"], best["asym_id"])]
    best["residue"] = best["residue_name"] + best["residue_chain"] + best["residue_seqid"].astype(str)
    best["pure_stacking"] = ~(best["polar_contact"] | best["metal_coordination"])
    best["url"] = "https://www.rcsb.org/structure/" + best["entry_id"]
    return best.sort_values(["entry_id", "asym_id", "residue"])[HIT_COLUMNS]


def select_entries(client: RcsbClient, cfg: Config) -> tuple[list[str], dict]:
    """Entries to scan plus search statistics (for ``run_info.json``)."""
    if cfg.entry_ids:
        entries = sorted({e.upper() for e in cfg.entry_ids})
        return entries, {"mode": "explicit", "entries_scanned": len(entries)}
    query = discovery.build_search_query(
        cfg.residue, cfg.min_aromatic_bonds, cfg.loi_only, cfg.max_resolution,
        cfg.released_after, cfg.released_before,
    )
    matching, n_entities = discovery.search_entries(client, query)
    entries = discovery.sample_entries(matching, cfg.max_entries, cfg.seed)
    info = {
        "mode": "search",
        "search_query": query,
        "matching_entries": len(matching),
        "matching_ligand_entities": n_entities,
        "entries_scanned": len(entries),
    }
    return entries, info


def _frame(rows: list[dict], columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=columns)


def run(cfg: Config, client: RcsbClient | None = None) -> dict[str, pd.DataFrame]:
    """Run the full scan and write tables, plots and examples. Returns the tables."""
    cfg.residue = cfg.residue.upper()
    if cfg.residue not in structures.SUPPORTED_RESIDUES:
        raise RuntimeError(
            f"Residue {cfg.residue!r} has no aromatic ring; choose from "
            f"{', '.join(structures.SUPPORTED_RESIDUES)}."
        )
    client = client or RcsbClient(cache_dir=cfg.cache_dir)
    entries, info = select_entries(client, cfg)
    if not entries:
        raise RuntimeError("No PDB entries matched the search criteria.")
    if info["mode"] == "search":
        tqdm.write(
            f"{info['matching_entries']:,} entries ({info['matching_ligand_entities']:,} ligand "
            f"entities) match; scanning {len(entries):,}."
        )

    batches = [entries[i : i + 50] for i in range(0, len(entries), 50)]
    records: dict[str, dict] = {}
    with ThreadPoolExecutor(max_workers=max(1, min(cfg.workers, 4))) as pool:
        for part in tqdm(
            pool.map(lambda b: discovery.fetch_entry_metadata(client, b), batches),
            total=len(batches), desc="Entry metadata", unit="batch",
        ):
            records.update(part)

    ligands: list[dict] = []
    stacks: list[dict] = []
    contacts: list[dict] = []
    failed = 0
    with ThreadPoolExecutor(max_workers=max(1, cfg.workers)) as pool:
        futures = {
            pool.submit(analyse_entry, client, cfg, records[e]): e for e in entries if e in records
        }
        for fut in tqdm(
            as_completed(futures), total=len(futures),
            desc=f"Ligand/{cfg.residue} stacking", unit="entry",
        ):
            try:
                res = fut.result()
            except Exception as exc:  # resilience: one bad entry must not stop the scan
                failed += 1
                tqdm.write(f"[warn] {futures[fut]}: {exc}")
                continue
            ligands.extend(res.ligands)
            stacks.extend(res.stacks)
            contacts.extend(res.contacts)

    ligands_df = _frame(ligands, LIGAND_COLUMNS).sort_values(["entry_id", "asym_id"])
    stacks_df = _frame(stacks, STACK_COLUMNS).sort_values(
        ["entry_id", "asym_id", "residue_chain", "residue_seqid", "distance"]
    )
    contacts_df = _frame(contacts, CONTACT_COLUMNS).sort_values(
        ["entry_id", "asym_id", "residue_chain", "residue_seqid", "pharmacophore"]
    )
    hits_df = face_to_face_hits(stacks_df, ligands_df)

    out = cfg.out_dir
    out.mkdir(parents=True, exist_ok=True)
    ligands_df.to_csv(out / "ligands.csv", index=False)
    stacks_df.to_csv(out / "stacking.csv", index=False)
    contacts_df.to_csv(out / "contacts.csv", index=False)
    hits_df.to_csv(out / "face_to_face_hits.csv", index=False)

    write_examples(
        client, stacks_df, out / "examples", per_group=cfg.examples_per_group,
        max_groups=cfg.max_groups, radius=cfg.pocket_radius,
    )
    plots.make_all(ligands_df, stacks_df, contacts_df, out / "plots", cfg.residue, cfg.criteria,
                   ring_systems=cfg.max_groups)

    info.update(
        {
            "tool_version": __version__,
            "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "entries_without_metadata": len(entries) - len(records),
            "entries_failed": failed,
            "ligand_instances": len(ligands_df),
            "ligand_instances_ok": int((ligands_df["status"] == "ok").sum()),
            "face_to_face_hits": len(hits_df),
            "face_to_face_entries": int(hits_df["entry_id"].nunique()),
            "config": {
                k: (str(v) if isinstance(v, Path) else sorted(v) if isinstance(v, frozenset) else v)
                for k, v in asdict(cfg).items()
            },
        }
    )
    (out / "run_info.json").write_text(json.dumps(info, indent=2, default=str))
    tqdm.write(
        f"Done: {len(ligands_df):,} ligand instances in {len(records):,} entries; "
        f"{len(hits_df):,} face-to-face {cfg.residue} stacks in "
        f"{info['face_to_face_entries']:,} entries -> {out}"
    )
    return {"ligands": ligands_df, "stacks": stacks_df, "contacts": contacts_df, "hits": hits_df}
