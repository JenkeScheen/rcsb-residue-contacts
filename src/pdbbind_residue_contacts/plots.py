"""Plotting of pharmacophore contact statistics, split by molecule class."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .classify import DRUG_LIKE, SOLVENT  # noqa: E402

log = logging.getLogger(__name__)

_DRUG_COLOR = "#2b6cb0"
_SOLVENT_COLOR = "#c05621"


def _counts(df: pd.DataFrame, cls: str) -> pd.Series:
    subset = df[df["ligand_class"] == cls]
    if subset.empty:
        return pd.Series(dtype=int)
    return subset["pharmacophore"].value_counts()


def _bar(counts: pd.Series, title: str, color: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, max(3, 0.35 * len(counts) + 1)))
    if counts.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
    else:
        counts = counts.sort_values()
        ax.barh(counts.index, counts.values, color=color)
        ax.set_xlabel("contact count")
        for i, v in enumerate(counts.values):
            ax.text(v, i, f" {v}", va="center", fontsize=8)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _combined_fraction(df: pd.DataFrame, residue: str, path: Path) -> None:
    drug = _counts(df, DRUG_LIKE)
    solvent = _counts(df, SOLVENT)
    labels = sorted(set(drug.index) | set(solvent.index))
    if not labels:
        _bar(pd.Series(dtype=int), "no data", _DRUG_COLOR, path)
        return
    drug_frac = np.array([drug.get(l, 0) for l in labels], dtype=float)
    solv_frac = np.array([solvent.get(l, 0) for l in labels], dtype=float)
    if drug_frac.sum():
        drug_frac = drug_frac / drug_frac.sum()
    if solv_frac.sum():
        solv_frac = solv_frac / solv_frac.sum()

    y = np.arange(len(labels))
    h = 0.4
    fig, ax = plt.subplots(figsize=(9, max(3, 0.4 * len(labels) + 1)))
    ax.barh(y - h / 2, drug_frac, height=h, color=_DRUG_COLOR, label="drug-like")
    ax.barh(y + h / 2, solv_frac, height=h, color=_SOLVENT_COLOR, label="solvent")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("fraction of contacts within class")
    ax.set_title(f"Pharmacophore composition contacting {residue}: drug-like vs solvent")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _class_summary(df: pd.DataFrame, residue: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    if df.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
    else:
        counts = df["ligand_class"].value_counts()
        colors = [
            _DRUG_COLOR if c == DRUG_LIKE else _SOLVENT_COLOR for c in counts.index
        ]
        ax.bar(counts.index, counts.values, color=colors)
        for i, v in enumerate(counts.values):
            ax.text(i, v, str(v), ha="center", va="bottom")
        ax.set_ylabel("contact records")
    ax.set_title(f"Contacts with {residue} by molecule class")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def make_all(df: pd.DataFrame, out_dir: Path, residue: str = "TRP") -> Path:
    """Render all summary plots into ``out_dir`` and return the directory."""
    out_dir.mkdir(parents=True, exist_ok=True)
    _bar(
        _counts(df, DRUG_LIKE),
        f"Drug-like pharmacophores contacting {residue}",
        _DRUG_COLOR,
        out_dir / "druglike_pharmacophores.png",
    )
    _bar(
        _counts(df, SOLVENT),
        f"Solvent pharmacophores contacting {residue}",
        _SOLVENT_COLOR,
        out_dir / "solvent_pharmacophores.png",
    )
    _combined_fraction(df, residue, out_dir / "combined_fraction.png")
    _class_summary(df, residue, out_dir / "ligand_class_summary.png")
    return out_dir
