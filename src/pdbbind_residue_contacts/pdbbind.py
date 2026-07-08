"""Discovery of a locally-downloaded PDBBind dataset.

PDBBind is distributed as one directory per complex, e.g. ``data/refined-set/1a4k/`` containing
``1a4k_protein.pdb``, ``1a4k_ligand.sdf``, ``1a4k_ligand.mol2`` and ``1a4k_pocket.pdb``. The
dataset is large (~11 GB) and must be downloaded by the user (see the project README); it is
never bundled with this package.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

PDBBIND_URL = "https://www.pdbbind-plus.org.cn/download?login=true"

# Standard PDBBind subset directory names.
SUBSET_DIRS = {
    "refined": ["refined-set"],
    "general": ["v2020-other-PL", "general-set", "general-set-except-refined"],
    "all": ["refined-set", "v2020-other-PL", "general-set", "general-set-except-refined"],
}

# Directory names that are not complexes.
_NON_ENTRY = {"index", "readme", "readme.txt"}


@dataclass(frozen=True)
class Entry:
    """One PDBBind complex on disk."""

    pdb_id: str
    directory: Path

    @property
    def protein_pdb(self) -> Path:
        return self.directory / f"{self.pdb_id}_protein.pdb"

    @property
    def ligand_sdf(self) -> Path:
        return self.directory / f"{self.pdb_id}_ligand.sdf"

    @property
    def ligand_mol2(self) -> Path:
        return self.directory / f"{self.pdb_id}_ligand.mol2"

    @property
    def pocket_pdb(self) -> Path:
        return self.directory / f"{self.pdb_id}_pocket.pdb"

    def is_valid(self) -> bool:
        return self.protein_pdb.exists() and (
            self.ligand_sdf.exists() or self.ligand_mol2.exists()
        )


def _iter_subset_dirs(data_dir: Path, subset: str) -> list[Path]:
    names = SUBSET_DIRS.get(subset)
    if names is None:
        raise ValueError(f"Unknown subset {subset!r}; choose from {sorted(SUBSET_DIRS)}")
    dirs = [data_dir / name for name in names if (data_dir / name).is_dir()]
    return dirs


def discover_entries(data_dir: Path, subset: str = "all") -> list[Entry]:
    """Return all valid PDBBind entries under ``data_dir`` for the requested subset."""
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        raise FileNotFoundError(
            f"Data directory {data_dir!r} not found. Create a free account and download "
            f"PDBBind from {PDBBIND_URL}, then extract it into that directory."
        )
    subset_dirs = _iter_subset_dirs(data_dir, subset)
    if not subset_dirs:
        raise FileNotFoundError(
            f"No PDBBind subset directories for {subset!r} under {data_dir!r}. "
            f"Expected one of {SUBSET_DIRS[subset]}. Download PDBBind from {PDBBIND_URL}."
        )
    entries: list[Entry] = []
    for sdir in subset_dirs:
        for child in sorted(sdir.iterdir()):
            if not child.is_dir() or child.name.lower() in _NON_ENTRY:
                continue
            entry = Entry(pdb_id=child.name, directory=child)
            if entry.is_valid():
                entries.append(entry)
    return entries


def select_entries(
    data_dir: Path,
    subset: str = "all",
    max_entries: int | None = None,
    seed: int = 0,
) -> list[Entry]:
    """Discover entries and, if capped, take a reproducible random sample of ``max_entries``."""
    entries = discover_entries(data_dir, subset)
    if max_entries is not None and 0 <= max_entries < len(entries):
        rng = random.Random(seed)
        entries = rng.sample(entries, max_entries)
        entries.sort(key=lambda e: e.pdb_id)
    return entries
