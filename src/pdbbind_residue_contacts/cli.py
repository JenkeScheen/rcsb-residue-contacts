"""Command-line interface for the PDBBind pharmacophore-contact scanner."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from tqdm import tqdm

from . import __version__
from .pipeline import Config, run

PDBBIND_URL = "https://www.pdbbind-plus.org.cn/download?login=true"

_DATA_HELP = (
    "PDBBind is not bundled with this tool. Create a free account and download it from\n"
    f"    {PDBBIND_URL}\n"
    "then extract the subsets so that '<data-dir>/refined-set/<id>/<id>_protein.pdb' and/or\n"
    "'<data-dir>/v2020-other-PL/<id>/...' exist (default <data-dir> is ./data)."
)


class _TqdmLoggingHandler(logging.Handler):
    """Route log records through tqdm.write so progress bars are never garbled."""

    def emit(self, record: logging.LogRecord) -> None:  # pragma: no cover - IO
        try:
            tqdm.write(self.format(record))
        except Exception:
            self.handleError(record)


def _configure_logging(verbose: bool) -> None:
    handler = _TqdmLoggingHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
    root = logging.getLogger()
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if verbose else logging.WARNING)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pdbbind-residue-contacts",
        description=(
            "Scan a local PDBBind dataset for the pharmacophores of non-protein molecules "
            "contacting a given residue type, split by drug-like vs solvent molecules."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--residue", default="TRP", help="Three-letter residue code to query.")
    p.add_argument(
        "--data-dir", type=Path, default=Path("data"),
        help="Directory containing the PDBBind subsets (refined-set / v2020-other-PL).",
    )
    p.add_argument(
        "--subset", choices=["refined", "general", "all"], default="refined",
        help="Which PDBBind subset(s) to scan.",
    )
    p.add_argument(
        "--cutoff", type=float, default=4.0,
        help="Contact distance cutoff in angstroms (heavy-atom to heavy-atom).",
    )
    p.add_argument(
        "--atom-scope", choices=["sidechain", "all"], default="sidechain",
        help="Residue atoms considered: sidechain only, or the whole residue.",
    )
    p.add_argument(
        "--max-structures", type=int, default=300,
        help="Cap on number of entries scanned (a reproducible random sample). Use -1 for all.",
    )
    p.add_argument("--out-dir", type=Path, default=Path("results"), help="Output directory.")
    p.add_argument(
        "--examples-per-group", type=int, default=3,
        help="Exemplary PDB files written per pharmacophore group.",
    )
    p.add_argument(
        "--max-groups", type=int, default=50,
        help="Only write example folders for the N largest pharmacophore groups.",
    )
    p.add_argument("--seed", type=int, default=0, help="Random seed for entry sampling.")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _configure_logging(args.verbose)

    if not args.data_dir.is_dir():
        parser.error(f"Data directory '{args.data_dir}' does not exist.\n\n{_DATA_HELP}")

    cfg = Config(
        data_dir=args.data_dir,
        residue=args.residue.upper(),
        subset=args.subset,
        cutoff=args.cutoff,
        atom_scope=args.atom_scope,
        max_structures=None if args.max_structures < 0 else args.max_structures,
        out_dir=args.out_dir,
        examples_per_group=args.examples_per_group,
        max_groups=args.max_groups,
        seed=args.seed,
    )

    try:
        run(cfg)
    except FileNotFoundError as exc:
        print(f"error: {exc}\n\n{_DATA_HELP}", file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
