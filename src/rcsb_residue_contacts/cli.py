"""Command-line interface: ligand <-> aromatic-residue stacking across the PDB via RCSB APIs."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from tqdm import tqdm

from . import __version__, discovery
from .pipeline import Config, run
from .rcsb import RcsbClient, RcsbError
from .stacking import StackCriteria
from .structures import SUPPORTED_RESIDUES


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
        prog="rcsb-residue-contacts",
        description=(
            "Find PDB entries whose ligands pi-stack (face-to-face / edge-to-face) on an aromatic "
            "residue type, using only the RCSB PDB Search, Data and ModelServer APIs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sel = p.add_argument_group("entry selection")
    sel.add_argument("--residue", default="HIS", type=str.upper, choices=SUPPORTED_RESIDUES,
                     help="Aromatic residue type to analyse.")
    sel.add_argument("--max-entries", type=int, default=300,
                     help="Reproducible random sample of matching entries to scan (-1 for all).")
    sel.add_argument("--seed", type=int, default=0, help="Random seed for entry sampling.")
    sel.add_argument("--entries", nargs="+", metavar="PDB_ID",
                     help="Scan these entries instead of searching the archive.")
    sel.add_argument("--released-after", metavar="YYYY-MM-DD",
                     help="Only entries first released on/after this date.")
    sel.add_argument("--released-before", metavar="YYYY-MM-DD",
                     help="Only entries first released before this date.")
    sel.add_argument("--max-resolution", type=float, default=3.0,
                     help="Resolution cutoff in Å (0 disables; NMR entries have no resolution).")
    sel.add_argument("--all-ligands", dest="loi_only", action="store_false",
                     help="Include ligands not flagged by RCSB as the ligand of interest "
                          "(cofactors, additives).")
    sel.add_argument("--min-aromatic-bonds", type=int, default=5,
                     help="Minimum aromatic bonds in the ligand's chemical component.")
    sel.add_argument("--exclude-comps", default="",
                     help="Comma-separated CCD codes to skip, e.g. HEM,FAD,NAD.")
    sel.add_argument("--count-only", action="store_true",
                     help="Only report how many entries/ligands match the search, then exit.")

    geo = p.add_argument_group("stacking geometry (PLIP defaults + clash floor)")
    geo.add_argument("--stack-max-distance", type=float, default=5.5,
                     help="Max ring centroid distance (Å).")
    geo.add_argument("--stack-min-distance", type=float, default=3.0,
                     help="Min ring centroid distance (Å); closer rings overlap (clash).")
    geo.add_argument("--stack-max-angle", type=float, default=30.0,
                     help="Max inter-plane angle for face-to-face (deg).")
    geo.add_argument("--stack-min-t-angle", type=float, default=60.0,
                     help="Min inter-plane angle for edge-to-face (deg).")
    geo.add_argument("--stack-max-offset", type=float, default=2.0,
                     help="Max lateral centroid offset (Å).")
    geo.add_argument("--contact-cutoff", type=float, default=4.0,
                     help="Heavy-atom distance for pharmacophore contacts (Å).")
    geo.add_argument("--atom-scope", choices=["sidechain", "all"], default="sidechain",
                     help="Residue atoms used for pharmacophore contacts.")

    out = p.add_argument_group("output")
    out.add_argument("--out-dir", type=Path, default=Path("results"), help="Output directory.")
    out.add_argument("--cache-dir", type=Path, default=Path(".rcsb_cache"),
                     help="Download cache (reruns are offline-fast).")
    out.add_argument("--no-cache", action="store_true", help="Disable the download cache.")
    out.add_argument("--workers", type=int, default=8, help="Concurrent RCSB requests.")
    out.add_argument("--examples-per-group", type=int, default=3,
                     help="Example pocket structures per ring system.")
    out.add_argument("--max-groups", type=int, default=30,
                     help="Example folders (and ring-systems grid cells) for the N most frequent ring systems.")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def config_from_args(args: argparse.Namespace) -> Config:
    return Config(
        residue=args.residue,
        max_entries=None if args.max_entries < 0 else args.max_entries,
        seed=args.seed,
        entry_ids=args.entries,
        released_after=args.released_after,
        released_before=args.released_before,
        max_resolution=args.max_resolution or None,
        loi_only=args.loi_only,
        min_aromatic_bonds=args.min_aromatic_bonds,
        exclude_comps=frozenset(c.strip().upper() for c in args.exclude_comps.split(",") if c.strip()),
        criteria=StackCriteria(
            max_distance=args.stack_max_distance,
            min_distance=args.stack_min_distance,
            max_parallel_angle=args.stack_max_angle,
            min_perpendicular_angle=args.stack_min_t_angle,
            max_offset=args.stack_max_offset,
        ),
        contact_cutoff=args.contact_cutoff,
        atom_scope=args.atom_scope,
        out_dir=args.out_dir,
        cache_dir=None if args.no_cache else args.cache_dir,
        workers=args.workers,
        examples_per_group=args.examples_per_group,
        max_groups=args.max_groups,
    )


def _count_only(cfg: Config) -> None:
    query = discovery.build_search_query(
        cfg.residue, cfg.min_aromatic_bonds, cfg.loi_only, cfg.max_resolution,
        cfg.released_after, cfg.released_before,
    )
    entries, n_entities = discovery.search_entries(RcsbClient(cache_dir=None), query)
    print(f"{len(entries):,} entries / {n_entities:,} ligand entities match "
          f"(aromatic ligand within 5 Å of {cfg.residue}).")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    _configure_logging(args.verbose)
    cfg = config_from_args(args)
    try:
        if args.count_only:
            _count_only(cfg)
        else:
            run(cfg)
    except RcsbError as exc:
        print(f"error: RCSB request failed: {exc}", file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
