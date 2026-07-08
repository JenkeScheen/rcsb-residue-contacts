"""Classify the *type* of interaction (H-bond, hydrophobic, pi-stacking, ...) with ProLIF.

The rest of this package answers "which pharmacophore of the molecule touches the residue?" by
distance-based SMARTS matching. That says *what group* is nearby but not *how* it interacts. To
name the interaction type we use `ProLIF <https://prolif.readthedocs.io/>`_, which encodes a
protein-ligand interaction fingerprint (Hydrophobic, HBDonor, HBAcceptor, PiStacking, PiCation /
CationPi, Anionic / Cationic, VdWContact, ...).

ProLIF is an **optional** dependency: import failures degrade gracefully (``HAS_PROLIF`` is
``False``) so the core offline pharmacophore scan still runs. Interaction typing needs a correct
bond graph and hydrogens, so it is applied only to the curated drug-like ligand (from the SDF),
never to the hydrogen-less solvent het groups whose bonds are perceived from geometry alone.
"""

from __future__ import annotations

import logging
from pathlib import Path

from rdkit import Chem

log = logging.getLogger(__name__)

try:  # ProLIF (and its MDAnalysis dependency) are optional.
    import prolif as plf

    HAS_PROLIF = True
except Exception:  # pragma: no cover - exercised only when prolif is absent
    plf = None
    HAS_PROLIF = False


def _fingerprint():
    """A ProLIF fingerprint with the default interaction set (created lazily, reused per run)."""
    return plf.Fingerprint()


def make_fingerprint():
    """Return a ProLIF ``Fingerprint`` or ``None`` when ProLIF is unavailable."""
    if not HAS_PROLIF:
        return None
    return _fingerprint()


def ligand_plf(sdf_path: Path, mol2_path: Path):
    """Load the drug-like ligand *with explicit hydrogens* as a ProLIF molecule.

    ProLIF's H-bond and charge rules need hydrogens and a real bond graph, so -- unlike the
    heavy-atom-only path used for distance contacts -- we keep/reperceive hydrogens here. Returns
    ``None`` if nothing usable can be read.
    """
    mol = None
    if sdf_path is not None and Path(sdf_path).exists():
        try:
            supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
            mol = next((m for m in supplier if m is not None), None)
        except Exception as exc:  # pragma: no cover - malformed file
            log.debug("ProLIF ligand SDF read failed for %s: %s", sdf_path, exc)
    if mol is None and mol2_path is not None and Path(mol2_path).exists():
        try:
            mol = Chem.MolFromMol2File(str(mol2_path), removeHs=False, sanitize=True)
        except Exception as exc:  # pragma: no cover - malformed file
            log.debug("ProLIF ligand mol2 read failed for %s: %s", mol2_path, exc)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    try:
        mol = Chem.AddHs(mol, addCoords=True)  # no-op where Hs are already explicit
        return plf.Molecule.from_rdkit(mol)
    except Exception as exc:  # pragma: no cover - unsanitisable ligand
        log.debug("ProLIF ligand build failed for %s: %s", sdf_path, exc)
        return None


def protein_plf(pdb_path: Path):
    """Load a protein PDB as a ProLIF molecule (hydrogens added), or ``None`` on failure.

    PDBBind protein PDBs carry no hydrogens; ProLIF needs them for donor/acceptor perception, so we
    add them at idealised positions after RDKit sanitises the standard-residue bond graph.
    """
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=True, sanitize=True)
    if mol is None:
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=True, sanitize=False)
    if mol is None:
        return None
    try:
        mol = Chem.AddHs(mol, addCoords=True)
        return plf.Molecule.from_rdkit(mol)
    except Exception as exc:  # pragma: no cover - unsanitisable protein
        log.debug("ProLIF protein build failed for %s: %s", pdb_path, exc)
        return None


def target_residues(prot, residue: str) -> list:
    """ProLIF residue ids in ``prot`` whose residue name matches ``residue`` (e.g. all TRPs)."""
    residue = residue.upper()
    return [rid for rid in prot.residues if str(rid.name).upper() == residue]


def residue_interactions(fp, ligand, prot, residue: str) -> list[dict]:
    """Return one row per detected ligand<->target-residue interaction.

    Each row carries the residue instance (chain + seqid) and the ProLIF interaction name. An empty
    list means no interaction with any instance of ``residue`` (or ProLIF failed for this pair).
    """
    residues = target_residues(prot, residue)
    if not residues:
        return []
    try:
        ifp = fp.generate(ligand, prot, residues=residues, metadata=True)
    except Exception as exc:  # pragma: no cover - resilience against odd geometries
        log.debug("ProLIF fingerprint failed: %s", exc)
        return []
    rows: list[dict] = []
    for (_lig_res, prot_res), interactions in ifp.items():
        chain = prot_res.chain
        seqid = prot_res.number
        for name, occurrences in interactions.items():
            # ``occurrences`` is a tuple of metadata dicts -- one per distinct atom grouping.
            for _ in occurrences:
                rows.append(
                    {
                        "residue_chain": "" if chain is None else str(chain),
                        "residue_seqid": str(seqid),
                        "interaction": name,
                    }
                )
    return rows
