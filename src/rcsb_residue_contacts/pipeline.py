"""End-to-end pipeline: scan PDBBind for pharmacophores contacting a residue type."""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from . import contacts, interactions, molecules, pharmacophores, plots
from .classify import DRUG_LIKE, SOLVENT
from .examples import ExampleCollector
from .pdbbind import Entry, select_entries

log = logging.getLogger(__name__)


@dataclass
class Config:
    data_dir: Path
    residue: str = "TRP"
    subset: str = "refined"
    cutoff: float = 4.0
    atom_scope: str = "sidechain"
    max_structures: int | None = 300
    out_dir: Path = Path("results")
    examples_per_group: int = 3
    max_groups: int = 50
    seed: int = 0
    classify_interactions: bool = True


def _mol_atom_tuples(mol: Chem.Mol) -> list[tuple[str, str, float, float, float]]:
    conf = mol.GetConformer()
    out = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        name = atom.GetProp("atomName") if atom.HasProp("atomName") else f"{atom.GetSymbol()}{idx}"
        pos = conf.GetAtomPosition(idx)
        out.append((name, atom.GetSymbol(), pos.x, pos.y, pos.z))
    return out


def _record_molecule(
    cfg: Config,
    mol: Chem.Mol,
    mol_coords,
    mol_class: str,
    molecule_id: str,
    ccd: str,
    residues: list[contacts.ResidueInstance],
    entry_id: str,
    rows: list[dict],
    collector: ExampleCollector,
    connectivity: bool,
    pocket_atoms: list | None = None,
) -> None:
    """Detect contacts to every residue instance and record pharmacophore hits + examples.

    ``connectivity`` selects the labeller: SMARTS (proper bond graph, e.g. SDF ligand) vs the
    connectivity-only graph rules used for hydrogen-less PDB het groups.
    """
    label_map = None if connectivity else pharmacophores.atom_label_map(mol)
    mol_tuples = None
    for res in residues:
        idxs, min_dist = contacts.contacting_atoms(res.coords, mol_coords, cfg.cutoff)
        if not idxs:
            continue
        if connectivity:
            labels = pharmacophores.assign_connectivity(mol, idxs)
        else:
            labels = pharmacophores.assign_with_map(mol, label_map, idxs)
        for label in labels:
            rows.append(
                {
                    "entry_id": entry_id,
                    "ligand_class": mol_class,
                    "molecule_id": molecule_id,
                    "ccd": ccd,
                    "residue_chain": res.chain,
                    "residue_seqid": res.seqid,
                    "pharmacophore": label,
                    "n_contact_atoms": len(idxs),
                    "min_distance": round(min_dist, 3),
                }
            )
            group_key = f"{mol_class}/{label}"
            if collector.needs(group_key):
                if mol_tuples is None:
                    mol_tuples = _mol_atom_tuples(mol)
                pdb_text = _render(res, ccd, mol_tuples, pocket_atoms)
                collector.offer(group_key, entry_id, molecule_id, pdb_text)


def _render(res, ccd, mol_tuples, pocket_atoms=None):
    from .examples import render_example_pdb

    return render_example_pdb(res, ccd, mol_tuples, pocket_atoms)


def _record_interactions(
    cfg: Config, entry: Entry, fp, int_rows: list[dict]
) -> None:
    """Type the drug-like ligand<->residue interactions with ProLIF and append rows.

    Restricted to the curated ligand: ProLIF needs a correct bond graph and hydrogens, which the
    geometry-perceived solvent het groups lack.

    Uses PDBBind's ``_pocket.pdb`` (binding-site residues only) rather than the full protein: any
    residue close enough to interact is in the pocket, so results are identical, but building the
    ProLIF molecule from a whole protein costs seconds/entry vs milliseconds for the pocket. Falls
    back to the full protein only if the pocket file is missing.
    """
    protein_source = entry.pocket_pdb if entry.pocket_pdb.exists() else entry.protein_pdb
    prot = interactions.protein_plf(protein_source)
    if prot is None:
        return
    ligand = interactions.ligand_plf(entry.ligand_sdf, entry.ligand_mol2)
    if ligand is None:
        return
    for row in interactions.residue_interactions(fp, ligand, prot, cfg.residue):
        int_rows.append(
            {
                "entry_id": entry.pdb_id,
                "ligand_class": DRUG_LIKE,
                "molecule_id": entry.pdb_id,
                **row,
            }
        )


def _process_entry(
    cfg: Config,
    entry: Entry,
    rows: list[dict],
    collector: ExampleCollector,
    fp=None,
    int_rows: list[dict] | None = None,
) -> None:
    structure = contacts.read_protein(entry.protein_pdb)
    if structure is None:
        return
    residues = contacts.residue_instances(structure, cfg.residue, cfg.atom_scope)
    if not residues:
        return

    # Full binding pocket for example structures (empty if the file is missing -> per-residue fallback).
    pocket_atoms = (
        contacts.read_pocket_atoms(entry.pocket_pdb) if entry.pocket_pdb.exists() else []
    )

    # Drug-like: the curated PDBBind ligand.
    ligand = molecules.load_ligand(entry.ligand_sdf, entry.ligand_mol2)
    if ligand is not None:
        _record_molecule(
            cfg, ligand, molecules.coords_of(ligand), DRUG_LIKE,
            molecule_id=entry.pdb_id, ccd="LIG", residues=residues,
            entry_id=entry.pdb_id, rows=rows, collector=collector,
            connectivity=False, pocket_atoms=pocket_atoms,
        )
        if fp is not None and int_rows is not None:
            _record_interactions(cfg, entry, fp, int_rows)

    # Solvent: non-water het groups from the protein PDB.
    for het in contacts.solvent_het_groups(structure):
        mol = molecules.build_mol_from_atoms(het.atoms)
        if mol is None or mol.GetNumConformers() == 0:
            continue
        molecule_id = f"{het.ccd}_{het.chain}{het.seqid}"
        _record_molecule(
            cfg, mol, molecules.coords_of(mol), SOLVENT,
            molecule_id=molecule_id, ccd=het.ccd, residues=residues,
            entry_id=entry.pdb_id, rows=rows, collector=collector,
            connectivity=True, pocket_atoms=pocket_atoms,
        )


def run(cfg: Config) -> pd.DataFrame:
    """Run the full scan and write results, plots and example PDBs. Returns the results table."""
    entries = select_entries(cfg.data_dir, cfg.subset, cfg.max_structures, cfg.seed)
    if not entries:
        raise RuntimeError(
            f"No PDBBind entries found under {cfg.data_dir} for subset {cfg.subset!r}."
        )

    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    collector = ExampleCollector(
        out_dir=cfg.out_dir / "examples",
        per_group=cfg.examples_per_group,
        max_groups=cfg.max_groups,
    )

    fp = None
    int_rows: list[dict] = []
    if cfg.classify_interactions:
        if interactions.HAS_PROLIF:
            fp = interactions.make_fingerprint()
        else:
            tqdm.write(
                "[warn] ProLIF not installed; skipping interaction-type classification. "
                "Install with `pip install prolif` (or `pip install '.[interactions]'`)."
            )

    rows: list[dict] = []
    for entry in tqdm(entries, desc=f"Scanning {cfg.residue} contacts", unit="entry"):
        try:
            _process_entry(cfg, entry, rows, collector, fp=fp, int_rows=int_rows)
        except Exception as exc:  # pragma: no cover - resilience
            tqdm.write(f"[warn] {entry.pdb_id}: {exc}")

    df = pd.DataFrame(rows)
    results_csv = cfg.out_dir / "results.csv"
    df.to_csv(results_csv, index=False)

    int_df = pd.DataFrame(
        int_rows,
        columns=["entry_id", "ligand_class", "molecule_id", "residue_chain",
                 "residue_seqid", "interaction"],
    )
    int_df.to_csv(cfg.out_dir / "interactions.csv", index=False)

    # Group counts for examples (per class/pharmacophore) and plots.
    group_counts: Counter = Counter()
    if not df.empty:
        for cls, label in zip(df["ligand_class"], df["pharmacophore"]):
            group_counts[f"{cls}/{label}"] += 1
    collector.finalize(dict(group_counts))

    plots.make_all(df, cfg.out_dir / "plots", residue=cfg.residue, interactions_df=int_df)

    tqdm.write(
        f"Done: {len(entries)} entries, {len(df)} contact records, "
        f"{len(int_df)} interactions -> {cfg.out_dir}"
    )
    return df
