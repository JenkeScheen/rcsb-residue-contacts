"""Write exemplary PDB files illustrating each pharmacophore interaction.

For every pharmacophore group we keep a handful of small example structures, each containing the
whole binding pocket (all residues from PDBBind's ``_pocket.pdb``, as ATOM records) plus the
interacting molecule (HETATM records) in their shared coordinate frame -- so the queried residue
is shown in its structural context, not in isolation. When the pocket file is unavailable we fall
back to just the queried residue instance. Examples are buffered in memory and only written at the
end, once the top ``max_groups`` groups by total count are known.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from .contacts import ResidueInstance

_SAFE = re.compile(r"[^A-Za-z0-9]+")


def _slug(label: str) -> str:
    return _SAFE.sub("_", label).strip("_").lower() or "group"


def _atom_line(
    serial: int,
    record: str,
    atom_name: str,
    res_name: str,
    chain: str,
    res_seq: str,
    element: str,
    x: float,
    y: float,
    z: float,
) -> str:
    name = atom_name if len(atom_name) >= 4 else f" {atom_name:<3}"
    name = name[:4]
    res_name = (res_name or "LIG")[:3].rjust(3)
    chain = (chain or "A")[:1]
    try:
        seq = int(re.sub(r"[^0-9-]", "", res_seq) or 0)
    except ValueError:
        seq = 0
    element = element.strip()[:2].rjust(2)
    return (
        f"{record:<6}{serial:>5} {name:<4} {res_name} {chain}{seq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element}"
    )


def render_example_pdb(
    residue: ResidueInstance,
    molecule_resname: str,
    molecule_atoms: list[tuple[str, str, float, float, float]],
    pocket_atoms: list[tuple[str, str, str, str, str, str, float, float, float]] | None = None,
) -> str:
    """Render an example structure as PDB text: binding pocket + interacting molecule.

    ``pocket_atoms`` is the full ``_pocket.pdb`` content as
    ``(record, name, resname, chain, seqid, element, x, y, z)`` tuples; when omitted/empty we fall
    back to just the queried ``residue`` instance. ``molecule_atoms`` is a list of
    ``(atom_name, element, x, y, z)`` tuples, written as HETATM in the shared coordinate frame.
    """
    lines: list[str] = []
    serial = 1
    if pocket_atoms:
        for record, name, resname, chain, seqid, element, x, y, z in pocket_atoms:
            lines.append(
                _atom_line(serial, record, name, resname, chain, seqid, element, x, y, z)
            )
            serial += 1
    else:
        for rec in residue.atoms:
            lines.append(
                _atom_line(
                    serial, "ATOM", rec.name, residue.name, residue.chain,
                    residue.seqid, rec.element, rec.x, rec.y, rec.z,
                )
            )
            serial += 1
    lines.append("TER")
    for name, element, x, y, z in molecule_atoms:
        lines.append(
            _atom_line(
                serial, "HETATM", name, molecule_resname, "L", "1", element, x, y, z
            )
        )
        serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


@dataclass
class _Candidate:
    filename: str
    text: str


@dataclass
class ExampleCollector:
    """Buffers example PDBs per pharmacophore group and writes the top groups at the end."""

    out_dir: Path
    per_group: int = 3
    max_groups: int = 50
    _buffers: dict[str, list[_Candidate]] = field(default_factory=dict)

    def needs(self, label: str) -> bool:
        """True if this group still has room for another example."""
        return len(self._buffers.get(label, [])) < self.per_group

    def offer(self, label: str, entry_id: str, molecule_id: str, pdb_text: str) -> None:
        bucket = self._buffers.setdefault(label, [])
        if len(bucket) >= self.per_group:
            return
        fname = f"{entry_id}_{_slug(molecule_id)}.pdb"
        if any(c.filename == fname for c in bucket):
            return
        bucket.append(_Candidate(filename=fname, text=pdb_text))

    def finalize(self, counts: dict[str, int]) -> Path:
        """Write example folders for the ``max_groups`` largest groups by ``counts``."""
        top = sorted(counts, key=lambda k: counts[k], reverse=True)[: self.max_groups]
        self.out_dir.mkdir(parents=True, exist_ok=True)
        for label in top:
            candidates = self._buffers.get(label)
            if not candidates:
                continue
            group_dir = self.out_dir / _slug(label)
            group_dir.mkdir(parents=True, exist_ok=True)
            for cand in candidates:
                (group_dir / cand.filename).write_text(cand.text)
        return self.out_dir
