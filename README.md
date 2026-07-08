# pdbbind-residue-contacts

Mine a **local [PDBBind](https://www.pdbbind-plus.org.cn/) dataset** to discover which
**pharmacophores** of non-protein molecules make contact with a given **residue type**, then
count and plot them — split between **drug-like** ligands and **solvent-type** molecules
(ions, buffers, cryoprotectants). Waters are excluded. The default test case is tryptophan.

For every entry the tool:

1. reads the protein PDB and finds each instance of the queried residue (sidechain atoms by default);
2. treats the curated PDBBind `*_ligand.sdf` as the **drug-like** molecule, and every non-water
   het group in the protein PDB (typed as non-polymer by [gemmi](https://gemmi.readthedocs.io/))
   as a **solvent-type** molecule;
3. finds molecule atoms within a distance cutoff of the residue;
4. maps those atoms onto pharmacophore substructures (benzene, nitrile, carboxylate, sulfonamide,
   metal ion, …) via [RDKit](https://www.rdkit.org/) SMARTS matching;
5. classifies the *type* of each drug-like ligand↔residue interaction (hydrophobic, H-bond,
   π-stacking, π-cation, salt bridge, van der Waals, …) with
   [ProLIF](https://prolif.readthedocs.io/);
6. accumulates counts and renders plots, and writes a few exemplary PDB files per pharmacophore group.

Everything runs **fully offline** — no web queries.

## You must download PDBBind yourself

The dataset is large (~11 GB) and licensed, so it is **not** included in this repository and is
`.gitignore`d. Create a free account and download it from:

> **https://www.pdbbind-plus.org.cn/download?login=true**

Extract the subset(s) under `./data` (this is the `--data-dir`, default `./data`) so the layout
looks like this — the per-entry directory and file names are exactly as PDBBind ships them:

```
data/
├── refined-set/                     # scanned by --subset refined (and all)
│   ├── index/                       # PDBBind metadata; ignored by the tool
│   ├── readme/                      # ignored by the tool
│   ├── 10gs/                        # one directory per complex, named by PDB id
│   │   ├── 10gs_protein.pdb         # REQUIRED — protein, source of the queried residue + solvent het groups
│   │   ├── 10gs_ligand.sdf          # REQUIRED — curated drug-like ligand (mol2 used only if the sdf is missing)
│   │   ├── 10gs_ligand.mol2         # optional — fallback for the ligand
│   │   └── 10gs_pocket.pdb          # recommended — binding-site residues; used for interaction typing + example structures
│   ├── 1a4k/
│   │   └── ...
│   └── ...
└── v2020-other-PL/                  # scanned by --subset general (and all); the "general set"
    ├── index/
    ├── readme/
    └── <pdbid>/ ...                 # same four files per entry as above
```

Details:

- An entry is used only if `<pdbid>_protein.pdb` **and** one of `<pdbid>_ligand.sdf` /
  `<pdbid>_ligand.mol2` are present; other files are optional.
- `<pdbid>_pocket.pdb` is not required, but without it interaction typing falls back to parsing the
  whole protein (much slower) and example structures fall back to the queried residue alone instead
  of the full pocket. Keep it.
- `--subset` selects which top-level directories are scanned: `refined` → `refined-set`;
  `general` → `v2020-other-PL` (also `general-set` / `general-set-except-refined` if present);
  `all` → every one of those found. `index`/`readme` directories are skipped automatically.
- `refined-set` alone is enough to try the tool; add `v2020-other-PL` for the general set.

## Installation

Using conda/mamba (recommended):

```bash
conda env create -f environment.yml   # or: mamba env create -f environment.yml
conda activate pdbbind-residue-contacts
```

The environment installs this package in editable mode (including ProLIF for interaction typing).
To install into an existing environment instead:

```bash
pip install -e ".[interactions]"   # or plain `pip install -e .` to skip ProLIF
```

Interaction-type classification needs [ProLIF](https://prolif.readthedocs.io/) (an optional
`interactions` extra). Without it, the tool still runs the full pharmacophore scan; the
interaction-type step is skipped, so `interactions.csv` is written empty and
`interaction_types.png` is a placeholder. Pass `--no-interactions` to skip it deliberately even
when ProLIF is installed.

## Usage

```bash
# Default: tryptophan, sidechain contacts, 4.0 Å, a 300-entry sample of the refined set
pdbbind-residue-contacts

# Explicit options
pdbbind-residue-contacts \
    --residue TRP \
    --data-dir data \
    --subset refined \
    --cutoff 4.0 \
    --atom-scope sidechain \
    --max-structures 300 \
    --out-dir results
```

Key options (`pdbbind-residue-contacts --help` for all):

| Option | Default | Meaning |
| --- | --- | --- |
| `--residue` | `TRP` | Three-letter residue code to query |
| `--cutoff` | `4.0` | Contact distance in Å (heavy-atom to heavy-atom) |
| `--atom-scope` | `sidechain` | `sidechain` or `all` residue atoms |
| `--subset` | `refined` | `refined`, `general`, or `all` |
| `--max-structures` | `300` | Cap on entries scanned (`-1` for all); reproducible random sample |
| `--examples-per-group` | `3` | Exemplary PDBs per pharmacophore group |
| `--max-groups` | `50` | Only write examples for the N largest groups |
| `--no-interactions` | (off) | Skip ProLIF interaction-type classification |

## Output

Written under `--out-dir` (default `results/`):

```
results/
├── results.csv                      # one row per (entry, molecule, pharmacophore) contact
├── interactions.csv                 # one row per drug-like ligand↔residue interaction (ProLIF)
├── plots/
│   ├── druglike_pharmacophores.png  # pharmacophore counts for drug-like ligands
│   ├── solvent_pharmacophores.png   # pharmacophore counts for solvent molecules
│   ├── combined_fraction.png        # normalized drug-like vs solvent composition
│   ├── ligand_class_summary.png     # drug-like vs solvent totals
│   └── interaction_types.png        # counts of interaction types (H-bond, hydrophobic, π-stacking, …)
└── examples/
    └── <class>_<pharmacophore>/     # top-50 groups, up to 3 PDBs each
        └── <pdbid>_<molecule>.pdb   # full binding pocket (_pocket.pdb) + interacting molecule,
                                     # shared coordinate frame
```

`results/` and `examples/` are gitignored.

## Development

```bash
pip install -e ".[dev]"
pytest
```

The unit tests are offline and do not require the PDBBind download.

## License

MIT — see [LICENSE](LICENSE). PDBBind itself is subject to its own license/terms; obtain it
directly from the PDBBind-CN site.
