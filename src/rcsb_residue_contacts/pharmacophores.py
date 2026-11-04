"""Pharmacophore (functional-group) definitions and assignment.

Each pharmacophore is a SMARTS pattern. For a given ligand component we compute, once, the
set of atom indices covered by every pattern. A pharmacophore is considered to *make* a
contact when at least one of its matched atoms is among the atoms found within the cutoff of
the target residue.
"""

from __future__ import annotations

from rdkit import Chem

from .molecules import HALIDES, METALS

# Ordered from most specific / informative to most generic. Order matters only for the
# element-based fallback below (specific patterns are always all evaluated).
PHARMACOPHORE_SMARTS: list[tuple[str, str]] = [
    ("carboxylate/acid", "[CX3](=O)[OX1H0-,OX2H1]"),
    ("ester", "[CX3](=O)[OX2H0][#6]"),
    ("amide", "[NX3][CX3](=[OX1])"),
    ("sulfonamide", "[SX4](=O)(=O)[NX3]"),
    ("sulfonyl/sulfonate", "[SX4](=O)(=O)"),
    ("phosphate", "[PX4](=O)([O,N])[O,N]"),
    ("nitro", "[NX3](=O)=O"),
    ("nitrile", "[NX1]#[CX2]"),
    ("guanidine", "[NX3][CX3](=[NX2,NX3+])[NX3]"),
    ("imidazole", "c1cnc[nH,n]1"),
    ("pyridine-type aromatic N", "[n]"),
    ("fused aromatic", "c1ccc2ccccc2c1"),
    ("benzene/phenyl", "c1ccccc1"),
    ("aromatic ring (other)", "[a]"),
    ("halogen", "[F,Cl,Br,I]"),
    ("hydroxyl", "[OX2H]"),
    ("ether", "[OX2]([#6])[#6]"),
    ("carbonyl (ketone/aldehyde)", "[CX3]=[OX1]"),
    ("primary amine", "[NX3;H2;!$(NC=O)]"),
    ("secondary/tertiary amine", "[NX3;H1,H0;!$(N=*);!$(NC=O);!$(N~[!#6])]"),
    ("thiol/thioether", "[#16X2]"),
    ("phosphorus", "[#15]"),
]

# Fallback labels keyed by element symbol, used when a contacting atom matched no pattern.
_ELEMENT_FALLBACK = {
    "C": "aliphatic carbon",
    "N": "other nitrogen",
    "O": "other oxygen",
    "S": "other sulfur",
    "P": "phosphorus",
}


def _compile() -> list[tuple[str, Chem.Mol]]:
    compiled = []
    for label, smarts in PHARMACOPHORE_SMARTS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:  # pragma: no cover - static definitions
            raise ValueError(f"Invalid SMARTS for {label!r}: {smarts}")
        compiled.append((label, patt))
    return compiled


_COMPILED = _compile()

# A generic label is dropped from an atom when a more specific label already covers that atom,
# so e.g. a benzene carbon is reported as "benzene/phenyl" rather than also "aromatic ring (other)".
_SUPPRESS_IF_PRESENT: dict[str, set[str]] = {
    "aromatic ring (other)": {
        "benzene/phenyl", "fused aromatic", "pyridine-type aromatic N", "imidazole",
    },
    "carbonyl (ketone/aldehyde)": {"carboxylate/acid", "ester", "amide"},
    "ether": {"ester"},
    "secondary/tertiary amine": {"guanidine", "amide"},
    "primary amine": {"guanidine", "amide"},
}


def atom_label_map(mol: Chem.Mol) -> dict[int, list[str]]:
    """Map each atom index to the list of pharmacophore labels covering it (generics suppressed)."""
    per_atom: dict[int, list[str]] = {}
    for label, patt in _COMPILED:
        for match in mol.GetSubstructMatches(patt, uniquify=True):
            for idx in match:
                per_atom.setdefault(idx, [])
                if label not in per_atom[idx]:
                    per_atom[idx].append(label)
    for idx, labels in per_atom.items():
        present = set(labels)
        per_atom[idx] = [
            l for l in labels if not (_SUPPRESS_IF_PRESENT.get(l, set()) & present)
        ]
    return per_atom


def _fallback_label(atom: Chem.Atom) -> str | None:
    sym = atom.GetSymbol()
    usym = sym.upper()
    if atom.GetAtomicNum() <= 1:
        return None  # hydrogens/dummies are not pharmacophores
    if usym in METALS:
        return "metal ion"
    if usym in HALIDES and atom.GetDegree() == 0:
        return "halide ion"
    return _ELEMENT_FALLBACK.get(sym, f"other ({sym})")


def assign_with_map(
    mol: Chem.Mol, label_map: dict[int, list[str]], contact_indices: set[int]
) -> set[str]:
    """Like :func:`assign` but reusing a precomputed ``atom_label_map`` (faster in loops)."""
    labels: set[str] = set()
    for idx in contact_indices:
        if label_map.get(idx):
            labels.update(label_map[idx])
        elif idx < mol.GetNumAtoms():
            fb = _fallback_label(mol.GetAtomWithIdx(idx))
            if fb is not None:
                labels.add(fb)
    return labels


def assign(mol: Chem.Mol, contact_indices: set[int]) -> set[str]:
    """Return the set of pharmacophore labels made by the given contacting atom indices.

    Every contacting atom is guaranteed to contribute at least one label: if it matched no
    SMARTS pattern it is labelled by an element-based fallback. Requires a correct bond graph,
    which ModelServer SDF (CCD bond orders) provides.
    """
    return assign_with_map(mol, atom_label_map(mol), contact_indices)
