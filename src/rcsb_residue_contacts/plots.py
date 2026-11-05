"""Plots for the ligand <-> residue pi-stacking scan.

Counting unit: a **distinct (entry, ligand component) pair**, so a ligand bound to every chain of
a homo-oligomer counts once. All plots are static PNGs; the CSVs next to them are the table view.

Palette (validated for colour-vision deficiency): slot-1 blue for the series in focus, orange for
the one contrasting series, muted grey for de-emphasised context.
"""

from __future__ import annotations

import io
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from .examples import ring_system_examples  # noqa: E402
from .stacking import EDGE_TO_FACE, FACE_TO_FACE, NO_STACK, StackCriteria  # noqa: E402

log = logging.getLogger(__name__)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"
ORANGE = "#eb6834"

_STYLE = {
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK_2,
    "axes.titlecolor": INK,
    "axes.titlesize": 12,
    "axes.titlelocation": "left",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "grid.linestyle": "-",
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "text.color": INK,
    "legend.frameon": False,
    "legend.labelcolor": INK_2,
    "font.family": "sans-serif",
    "font.size": 9.5,
}

PAIR = ["entry_id", "comp_id"]


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _placeholder(path: Path, title: str, msg: str = "no data") -> None:
    with plt.rc_context(_STYLE):
        fig, ax = plt.subplots(figsize=(6, 2.5))
        ax.text(0.5, 0.5, msg, ha="center", va="center", color=MUTED, transform=ax.transAxes)
        ax.axis("off")
        ax.set_title(title)
        _save(fig, path)


def _xgrid(ax) -> None:
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)


def _bar_end_labels(ax, values, xmax) -> None:
    for i, v in enumerate(values):
        ax.text(v + xmax * 0.01, i, f"{int(v)}", va="center", fontsize=8, color=INK_2)


# --------------------------------------------------------------------------- derivations
def _f2f(stacks: pd.DataFrame) -> pd.DataFrame:
    f2f = stacks[stacks["stack_type"] == FACE_TO_FACE].copy()
    f2f["pure"] = ~(f2f["polar_contact"].astype(bool) | f2f["metal_coordination"].astype(bool))
    return f2f


def interaction_modes(ligands: pd.DataFrame, stacks: pd.DataFrame) -> pd.Series:
    """Classify each (entry, ligand) pair by its strongest interaction mode with the residue."""
    ok = ligands[ligands["status"] == "ok"]
    if ok.empty:
        return pd.Series(dtype=int)
    pairs = ok.groupby(PAIR).agg(
        polar=("polar_contact_residues", lambda s: (s.fillna(0) > 0).any()),
        metal=("metal_contact_residues", lambda s: (s.fillna(0) > 0).any()),
    )
    f2f = _f2f(stacks)
    pure = set(map(tuple, f2f.loc[f2f["pure"], PAIR].values))
    any_f2f = set(map(tuple, f2f[PAIR].values))
    e2f = set(map(tuple, stacks.loc[stacks["stack_type"] == EDGE_TO_FACE, PAIR].values))
    labels = []
    for key, row in pairs.iterrows():
        if key in pure:
            labels.append("face-to-face, no polar/metal contact")
        elif key in any_f2f:
            labels.append("face-to-face + polar/metal contact")
        elif key in e2f:
            labels.append("edge-to-face only")
        elif row["polar"] or row["metal"]:
            labels.append("polar/metal contact only")
        else:
            labels.append("other contacts only")
    order = [
        "face-to-face, no polar/metal contact", "face-to-face + polar/metal contact",
        "edge-to-face only", "polar/metal contact only", "other contacts only",
    ]
    return pd.Series(labels).value_counts().reindex(order, fill_value=0)


# --------------------------------------------------------------------------------- plots
def ring_patterns_plot(stacks: pd.DataFrame, residue: str, path: Path, top: int = 20) -> None:
    """Ranked ligand ring patterns in face-to-face stacks, split by polar co-contact."""
    title = f"Ligand ring patterns stacking face-to-face with {residue}"
    f2f = _f2f(stacks)
    if f2f.empty:
        _placeholder(path, title, "no face-to-face stacks found")
        return
    per = f2f.groupby(["ring_pattern", *PAIR])["pure"].any().reset_index()
    counts = per.groupby("ring_pattern")["pure"].agg(pure="sum", total="count")
    counts["polar"] = counts["total"] - counts["pure"]
    counts = counts.sort_values("total", ascending=False).head(top).iloc[::-1]
    with plt.rc_context(_STYLE):
        fig, ax = plt.subplots(figsize=(8.5, 0.32 * len(counts) + 1.6))
        y = np.arange(len(counts))
        kw = dict(height=0.62, edgecolor=SURFACE, linewidth=1.5)
        ax.barh(y, counts["pure"], color=BLUE, label=f"no polar/metal contact to {residue}", **kw)
        ax.barh(y, counts["polar"], left=counts["pure"], color=ORANGE,
                label=f"also polar/metal contact to the same {residue}", **kw)
        ax.set_yticks(y, counts.index)
        xmax = counts["total"].max()
        _bar_end_labels(ax, counts["total"].values, xmax)
        ax.set_xlim(0, xmax * 1.1)
        ax.set_xlabel("distinct entry–ligand pairs")
        ax.set_title(title)
        ax.legend(loc="lower right")
        _xgrid(ax)
        _save(fig, path)


def ring_systems_plot(
    stacks: pd.DataFrame, residue: str, path: Path, top: int = 30, ids_per_system: int = 3
) -> None:
    """Structure grid of the ring systems most often found stacking face-to-face.

    Each cell shows the rank, the number of distinct entry-ligand pairs and up to
    ``ids_per_system`` example PDB IDs -- the same picks written to ``examples/<rank>_...``.
    """
    title = f"Ring systems stacking face-to-face with {residue} (substituents removed)"
    picks = ring_system_examples(stacks, per_group=ids_per_system, max_groups=top)
    if picks.empty:
        _placeholder(path, title, "no face-to-face stacks found")
        return
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    cells = []
    for (rank, smi, n), group in picks.groupby(["rank", "ring_system_smiles", "n_pairs"], sort=True):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        drawer = rdMolDraw2D.MolDraw2DCairo(300, 200)
        opts = drawer.drawOptions()
        opts.setBackgroundColour((0.988, 0.988, 0.984, 1.0))  # SURFACE
        opts.fixedBondLength = 28  # same scale in every cell
        opts.padding = 0.08
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        img = plt.imread(io.BytesIO(drawer.GetDrawingText()))
        cells.append((rank, n, list(group["entry_id"]), img))

    ncols = min(5, len(cells))
    nrows = -(-len(cells) // ncols)
    cell_w, cell_h, header = 2.3, 2.15, 0.75
    fig_w, fig_h = max(ncols * cell_w, 8.0), nrows * cell_h + header
    with plt.rc_context(_STYLE):
        fig = plt.figure(figsize=(fig_w, fig_h))
        x0 = (fig_w - ncols * cell_w) / 2
        for k, (rank, n, ids, img) in enumerate(cells):
            r, c = divmod(k, ncols)
            left = (x0 + c * cell_w + 0.1) / fig_w
            bottom = (fig_h - header - (r + 1) * cell_h + 0.55) / fig_h
            ax = fig.add_axes([left, bottom, (cell_w - 0.2) / fig_w, (cell_h - 0.6) / fig_h])
            ax.imshow(img)
            ax.axis("off")
            ax.text(0.5, -0.04, f"#{rank}   n = {n}", transform=ax.transAxes, ha="center",
                    va="top", fontsize=10, color=INK)
            ax.text(0.5, -0.2, "  ".join(ids), transform=ax.transAxes, ha="center",
                    va="top", fontsize=9, color=INK_2)
        fig.text(0.1 / fig_w, 1 - 0.1 / fig_h, title, ha="left", va="top", fontsize=12,
                 color=INK)
        fig.text(0.1 / fig_w, 1 - 0.4 / fig_h,
                 "#rank = examples/ folder · n = distinct entry–ligand pairs · PDB IDs = example "
                 "stacks (pure stacking first, then closest)",
                 ha="left", va="top", fontsize=9, color=INK_2)
        fig.savefig(path, dpi=150)
        plt.close(fig)


def geometry_plot(
    stacks: pd.DataFrame, residue: str, criteria: StackCriteria, path: Path
) -> None:
    """Centroid distance vs inter-plane angle and vs offset for every ring pair in the window."""
    title = f"Ligand aromatic rings within {criteria.window:g} Å of {residue} rings"
    if stacks.empty:
        _placeholder(path, title)
        return
    styles = [
        (NO_STACK, MUTED, "not stacked", 0.35, 7),
        (EDGE_TO_FACE, ORANGE, "edge-to-face (T-shaped)", 0.8, 12),
        (FACE_TO_FACE, BLUE, "face-to-face", 0.8, 12),
    ]
    with plt.rc_context(_STYLE):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
        lo, width = criteria.min_distance, criteria.max_distance - criteria.min_distance
        ax1.add_patch(Rectangle((lo, 0), width, criteria.max_parallel_angle,
                                color=BLUE, alpha=0.07, linewidth=0))
        ax1.add_patch(Rectangle((lo, criteria.min_perpendicular_angle), width,
                                90 - criteria.min_perpendicular_angle,
                                color=ORANGE, alpha=0.07, linewidth=0))
        ax2.add_patch(Rectangle((lo, 0), width, criteria.max_offset,
                                color=MUTED, alpha=0.08, linewidth=0))
        handles = []
        for kind, color, label, alpha, size in styles:
            sub = stacks[stacks["stack_type"] == kind]
            kw = dict(s=size, color=color, alpha=alpha, edgecolors=SURFACE, linewidths=0.3)
            ax1.scatter(sub["distance"], sub["angle"], **kw)
            ax2.scatter(sub["distance"], sub["offset"], **kw)
            handles.append(Patch(color=color, label=f"{label} ({len(sub):,})"))
        for ax in (ax1, ax2):
            for x in (criteria.min_distance, criteria.max_distance):
                ax.axvline(x, color=MUTED, linewidth=0.8, linestyle="--")
            ax.set_xlim(left=min(3.0, float(stacks["distance"].min()) - 0.2),
                        right=criteria.window + 0.1)
            ax.set_xlabel("ring centroid distance (Å)")
            ax.grid(True)
            ax.set_axisbelow(True)
        ax1.set_ylim(0, 90)
        ax1.set_yticks([0, 15, 30, 45, 60, 75, 90])
        ax1.set_ylabel("angle between ring planes (°)")
        ax1.set_title("Angle: face-to-face ≤ "
                      f"{criteria.max_parallel_angle:g}°, edge-to-face ≥ "
                      f"{criteria.min_perpendicular_angle:g}°", fontsize=10)
        ax2.axhline(criteria.max_offset, color=MUTED, linewidth=0.8, linestyle="--")
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel("lateral offset (Å)")
        ax2.set_title(f"Either stack needs offset ≤ {criteria.max_offset:g} Å and "
                      f"{criteria.min_distance:g}–{criteria.max_distance:g} Å distance",
                      fontsize=10)
        fig.suptitle(title, x=0.01, ha="left", fontsize=12, color=INK)
        fig.legend(handles=handles, loc="upper right", ncol=3, bbox_to_anchor=(0.99, 0.99))
        _save(fig, path)


def modes_plot(ligands: pd.DataFrame, stacks: pd.DataFrame, residue: str, path: Path) -> None:
    """How ligands that neighbour the residue actually interact with it."""
    title = f"How aromatic ligands next to {residue} interact with it"
    modes = interaction_modes(ligands, stacks)
    if modes.empty or modes.sum() == 0:
        _placeholder(path, title)
        return
    modes = modes.iloc[::-1]
    colors = [BLUE if m.startswith("face-to-face, no") else MUTED for m in modes.index]
    with plt.rc_context(_STYLE):
        fig, ax = plt.subplots(figsize=(8.5, 3.2))
        y = np.arange(len(modes))
        ax.barh(y, modes.values, color=colors, height=0.62, edgecolor=SURFACE, linewidth=1.5)
        ax.set_yticks(y, modes.index)
        xmax = max(modes.max(), 1)
        total = modes.sum()
        for i, v in enumerate(modes.values):
            ax.text(v + xmax * 0.01, i, f"{v} ({100 * v / total:.0f}%)", va="center",
                    fontsize=8, color=INK_2)
        ax.set_xlim(0, xmax * 1.18)
        ax.set_xlabel("distinct entry–ligand pairs (each counted once, first matching mode)")
        ax.set_title(title)
        _xgrid(ax)
        _save(fig, path)


def by_year_plot(stacks: pd.DataFrame, residue: str, path: Path, split_year: int = 2020) -> None:
    """Entries with a face-to-face stack, by initial release year."""
    title = f"Entries with a face-to-face {residue} stack, by release year (scanned entries)"
    f2f = stacks[stacks["stack_type"] == FACE_TO_FACE]
    years = pd.to_numeric(
        f2f.drop_duplicates("entry_id")["release_date"].astype(str).str[:4], errors="coerce"
    ).dropna().astype(int)
    if years.empty:
        _placeholder(path, title, "no face-to-face stacks found")
        return
    counts = years.value_counts().sort_index()
    full = pd.Series(0, index=range(counts.index.min(), counts.index.max() + 1))
    full.update(counts)
    colors = [BLUE if yr >= split_year else MUTED for yr in full.index]
    late = int(full[full.index >= split_year].sum())
    with plt.rc_context(_STYLE):
        fig, ax = plt.subplots(figsize=(9, 3.6))
        ax.bar(full.index, full.values, color=colors, width=0.8, edgecolor=SURFACE, linewidth=1)
        ax.set_ylabel("entries")
        ax.set_title(title)
        ax.grid(axis="y")
        ax.set_axisbelow(True)
        ax.legend(
            handles=[
                Patch(color=MUTED, label=f"released before {split_year}"),
                Patch(color=BLUE, label=f"released {split_year} or later "
                                        f"({late} of {int(full.sum())}, {100 * late / full.sum():.0f}%)"),
            ],
            loc="upper left",
        )
        _save(fig, path)


def pharmacophores_plot(
    contacts: pd.DataFrame, stacks: pd.DataFrame, residue: str, path: Path, top: int = 15
) -> None:
    """Which ligand groups sit within the contact cutoff of the residue: all vs stackers."""
    title = f"Ligand groups in contact with {residue} side chains"
    if contacts.empty:
        _placeholder(path, title)
        return
    stackers = set(map(tuple, stacks.loc[stacks["stack_type"] == FACE_TO_FACE, PAIR].values))
    pairs = contacts.drop_duplicates(["pharmacophore", *PAIR])
    n_all = contacts[PAIR].drop_duplicates().shape[0]
    frac_all = pairs["pharmacophore"].value_counts() / n_all
    is_stacker = [tuple(k) in stackers for k in pairs[PAIR].values]
    sub = pairs[is_stacker]
    n_stack = sub[PAIR].drop_duplicates().shape[0]
    frac_stack = (sub["pharmacophore"].value_counts() / n_stack) if n_stack else frac_all * 0
    labels = frac_all.sort_values(ascending=False).head(top).index[::-1]
    a = frac_all.reindex(labels, fill_value=0).values
    s = frac_stack.reindex(labels, fill_value=0).values
    with plt.rc_context(_STYLE):
        fig, ax = plt.subplots(figsize=(8.5, 0.36 * len(labels) + 1.6))
        y = np.arange(len(labels))
        h = 0.38
        ax.barh(y + h / 2, a, height=h, color=MUTED, edgecolor=SURFACE, linewidth=1,
                label=f"all ligands in contact (n = {n_all})")
        ax.barh(y - h / 2, s, height=h, color=BLUE, edgecolor=SURFACE, linewidth=1,
                label=f"face-to-face stackers (n = {n_stack})")
        ax.set_yticks(y, labels)
        ax.set_xlim(0, 1)
        ax.set_xlabel("fraction of entry–ligand pairs with this group in contact")
        ax.set_title(title)
        ax.legend(loc="lower right")
        _xgrid(ax)
        _save(fig, path)


def make_all(
    ligands: pd.DataFrame,
    stacks: pd.DataFrame,
    contacts: pd.DataFrame,
    out_dir: Path,
    residue: str = "HIS",
    criteria: StackCriteria = StackCriteria(),
    ring_systems: int = 30,
) -> Path:
    """Render every plot into ``out_dir`` and return it.

    ``ring_systems`` should equal the number of example groups written, so that grid cell ``#N``
    and folder ``examples/NN_*`` always describe the same ring system.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    ring_patterns_plot(stacks, residue, out_dir / "face_to_face_ring_patterns.png")
    ring_systems_plot(stacks, residue, out_dir / "face_to_face_ring_systems.png", top=ring_systems)
    geometry_plot(stacks, residue, criteria, out_dir / "stacking_geometry.png")
    modes_plot(ligands, stacks, residue, out_dir / "interaction_modes.png")
    by_year_plot(stacks, residue, out_dir / "face_to_face_by_year.png")
    pharmacophores_plot(contacts, stacks, residue, out_dir / "contact_pharmacophores.png")
    return out_dir
