"""Example structures for the ring systems most often found stacking face-to-face.

Ring systems are ranked by distinct (entry, ligand) pairs, exactly as in
``plots/face_to_face_ring_systems.png``, and the PDB IDs printed under each structure there are the
files written here. Per ring system, up to ``per_group`` stacks from distinct entries are chosen:
pure stacks (no polar/metal contact to the same residue) first, then by centroid distance.

Each example is the RCSB ModelServer ``residueSurroundings`` pocket (ligand + everything within
``radius`` A, waters included) in the deposited frame, written as PDB with ``REMARK 999`` lines
describing the stack. Entries that PDB format cannot represent (CCD codes longer than 3 characters,
chain IDs longer than 2) are written as mmCIF instead.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import gemmi
import pandas as pd

from .rcsb import RcsbClient, RcsbError
from .stacking import FACE_TO_FACE

log = logging.getLogger(__name__)

_SAFE = re.compile(r"[^A-Za-z0-9]+")

PICK_COLUMNS = [
    "rank", "ring_system_smiles", "n_pairs", "entry_id", "comp_id", "asym_id", "residue",
    "distance", "angle", "offset", "pure_stacking",
]


def slug(label: str) -> str:
    return _SAFE.sub("_", label).strip("_") or "group"


def group_dir_name(rank: int, smiles: str) -> str:
    return f"{rank:02d}_{slug(smiles)[:60]}"


def ring_system_examples(
    stacks: pd.DataFrame, per_group: int = 3, max_groups: int = 30
) -> pd.DataFrame:
    """Rank ring systems in face-to-face stacks and pick example stacks for each (see module doc)."""
    f2f = stacks[
        (stacks["stack_type"] == FACE_TO_FACE) & (stacks["ring_system_smiles"].fillna("") != "")
    ].copy()
    if f2f.empty or per_group <= 0:
        return pd.DataFrame(columns=PICK_COLUMNS)
    f2f["pure_stacking"] = ~(
        f2f["polar_contact"].astype(bool) | f2f["metal_coordination"].astype(bool)
    )
    f2f["residue"] = (
        f2f["residue_name"] + f2f["residue_chain"].astype(str) + f2f["residue_seqid"].astype(str)
    )
    counts = (
        f2f.drop_duplicates(["ring_system_smiles", "entry_id", "comp_id"])
        .groupby("ring_system_smiles").size().rename("n_pairs").reset_index()
        .sort_values(["n_pairs", "ring_system_smiles"], ascending=[False, True])
        .head(max_groups)
    )
    ordered = f2f.sort_values(
        ["pure_stacking", "distance", "entry_id", "asym_id"], ascending=[False, True, True, True]
    )
    picks = []
    for rank, (smiles, n_pairs) in enumerate(counts.itertuples(index=False), start=1):
        chosen = ordered[ordered["ring_system_smiles"] == smiles].drop_duplicates("entry_id")
        chosen = chosen.head(per_group).assign(rank=rank, n_pairs=n_pairs)
        picks.append(chosen)
    return pd.concat(picks)[PICK_COLUMNS].reset_index(drop=True)


def _remarks(pick, radius: float) -> list[str]:
    lines = [
        "face-to-face pi-stack (rcsb-residue-contacts)",
        f"ligand {pick.comp_id} label_asym_id {pick.asym_id}; ring system rank {pick.rank}",
        f"ring system {pick.ring_system_smiles}",
        f"stacks on {pick.residue}: centroid {pick.distance:.2f} A, angle {pick.angle:.1f} deg,"
        f" offset {pick.offset:.2f} A",
        f"polar/metal contact to this residue: {'no' if pick.pure_stacking else 'yes'}",
        f"pocket = residues within {radius:g} A of the ligand (RCSB ModelServer)",
    ]
    return [f"REMARK 999 {line}"[:80] for line in lines]


def to_pdb(cif_text: str, remarks: list[str]) -> str | None:
    """mmCIF pocket -> PDB text with ``remarks`` on top; ``None`` if PDB can't represent it."""
    try:
        doc = gemmi.cif.read_string(cif_text)
        structure = gemmi.make_structure_from_block(doc.sole_block())
    except Exception:
        return None
    model = structure[0] if len(structure) else []
    if any(len(ch.name) > 2 for ch in model) or any(
        len(res.name) > 3 for ch in model for res in ch
    ):
        return None
    pdb = structure.make_pdb_string(gemmi.PdbWriteOptions(minimal=True))
    return "\n".join(remarks) + "\n" + pdb


def write_examples(
    client: RcsbClient,
    stacks: pd.DataFrame,
    out_dir: Path,
    per_group: int = 3,
    max_groups: int = 30,
    radius: float = 6.0,
) -> pd.DataFrame:
    """Download and write example pockets; returns (and writes) ``out_dir/index.csv``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    picks = ring_system_examples(stacks, per_group, max_groups)
    files: list[str | None] = []
    for pick in picks.itertuples(index=False):
        try:
            cif = client.pocket_cif(pick.entry_id, pick.asym_id, radius)
        except RcsbError as exc:
            log.warning("pocket download failed for %s/%s: %s", pick.entry_id, pick.asym_id, exc)
            files.append(None)
            continue
        pdb = to_pdb(cif, _remarks(pick, radius))
        text, ext = (pdb, "pdb") if pdb is not None else (cif, "cif")
        group = out_dir / group_dir_name(pick.rank, pick.ring_system_smiles)
        group.mkdir(parents=True, exist_ok=True)
        fname = f"{pick.entry_id}_{pick.comp_id}_{pick.asym_id}__{pick.residue}.{ext}"
        (group / fname).write_text(text)
        files.append(f"{group.name}/{fname}")
    index = picks.assign(file=files).dropna(subset=["file"]) if len(picks) else picks.assign(file=[])
    index.to_csv(out_dir / "index.csv", index=False)
    return index
