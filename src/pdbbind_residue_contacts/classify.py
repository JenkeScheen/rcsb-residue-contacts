"""Molecule-class helpers.

In the PDBBind workflow the drug-like / solvent split is determined by *source*: the curated
``_ligand.sdf`` is drug-like, while non-water het groups in the protein PDB are solvent-type.
This module holds the small shared vocabulary (class labels, water/ion recognition, and a
reference set of common crystallisation additives) used across the pipeline.
"""

from __future__ import annotations

from .molecules import HALIDES, METALS

DRUG_LIKE = "drug-like"
SOLVENT = "solvent"

# Residue names that are water (never analysed).
WATER_NAMES = {"HOH", "DOD", "WAT", "H2O", "TIP", "SOL"}

# Common crystallisation additives / buffers / cryoprotectants. Informational: any non-water
# het group is treated as solvent regardless, but this set is handy for labelling/reporting.
COMMON_ADDITIVES = {
    "EDO", "GOL", "PEG", "PG4", "PGE", "1PE", "P6G", "MPD", "SO4", "PO4", "IOD",
    "ACT", "ACY", "DMS", "TRS", "EPE", "MES", "FMT", "CIT", "FLC", "BME", "IPA",
    "NO3", "MOH", "BU3", "TLA", "GAL", "NAG", "BMA", "MAN", "SCN", "AZI",
}


def is_water(name: str) -> bool:
    return name.strip().upper() in WATER_NAMES


def is_metal(element: str) -> bool:
    return element.strip().upper() in METALS


def is_halide(element: str) -> bool:
    return element.strip().upper() in HALIDES


def is_common_additive(name: str) -> bool:
    return name.strip().upper() in COMMON_ADDITIVES
