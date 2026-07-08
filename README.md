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
5. accumulates counts and renders plots, and writes a few exemplary PDB files per pharmacophore group.

Everything runs **fully offline** — no web queries.

## You must download PDBBind yourself

The dataset is large (~11 GB) and licensed, so it is **not** included in this repository and is
`.gitignore`d. Create a free account and download it from:

> **https://www.pdbbind-plus.org.cn/download?login=true**

Extract the subset(s) under `./data` so the layout looks like:

```
data/
├── refined-set/
│   └── <pdbid>/<pdbid>_protein.pdb, <pdbid>_ligand.sdf, <pdbid>_ligand.mol2, <pdbid>_pocket.pdb
└── v2020-other-PL/
    └── <pdbid>/ ...
```

`refined-set` alone is enough to try the tool; add `v2020-other-PL` for the general set.

## Installation

Using conda/mamba (recommended):

```bash
conda env create -f environment.yml   # or: mamba env create -f environment.yml
conda activate pdbbind-residue-contacts
```

The environment installs this package in editable mode. To install into an existing environment
instead:

```bash
pip install -e .
```

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

## Output

Written under `--out-dir` (default `results/`):

```
results/
├── results.csv                      # one row per (entry, molecule, pharmacophore) contact
├── plots/
│   ├── druglike_pharmacophores.png  # pharmacophore counts for drug-like ligands
│   ├── solvent_pharmacophores.png   # pharmacophore counts for solvent molecules
│   ├── combined_fraction.png        # normalized drug-like vs solvent composition
│   └── ligand_class_summary.png     # drug-like vs solvent totals
└── examples/
    └── <class>_<pharmacophore>/     # top-50 groups, up to 3 PDBs each
        └── <pdbid>_<molecule>.pdb   # residue + interacting molecule, shared coordinate frame
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
