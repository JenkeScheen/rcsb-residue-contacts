"""Interaction-type classification (ProLIF) and its plot.

The plotting path is offline and prolif-free (it takes a plain DataFrame). The ProLIF module is
only import-guarded here; the heavy fingerprinting is exercised by the pipeline, not unit-tested,
because it needs the licensed PDBBind download.
"""

import pandas as pd

from pdbbind_residue_contacts import interactions, plots


def test_interactions_module_import_guard():
    # Whether or not prolif is installed, the module imports and exposes the flag + factory.
    assert isinstance(interactions.HAS_PROLIF, bool)
    if not interactions.HAS_PROLIF:
        assert interactions.make_fingerprint() is None


def test_interaction_types_plot_with_data(tmp_path):
    df = pd.DataFrame(
        {
            "interaction": [
                "Hydrophobic", "Hydrophobic", "VdWContact", "PiStacking", "HBDonor",
            ]
        }
    )
    out = tmp_path / "interaction_types.png"
    plots._interaction_types(df, "TRP", out)
    assert out.exists() and out.stat().st_size > 0


def test_interaction_types_plot_handles_empty_and_missing(tmp_path):
    # Empty DataFrame (prolif ran, found nothing) and None (prolif not installed) both render.
    for label, value in [("empty", pd.DataFrame(columns=["interaction"])), ("none", None)]:
        out = tmp_path / f"interaction_types_{label}.png"
        plots._interaction_types(value, "TRP", out)
        assert out.exists() and out.stat().st_size > 0


def test_make_all_includes_interaction_plot(tmp_path):
    df = pd.DataFrame(
        {
            "ligand_class": ["drug-like", "solvent"],
            "pharmacophore": ["benzene/phenyl", "sulfate"],
        }
    )
    int_df = pd.DataFrame({"interaction": ["Hydrophobic", "PiStacking"]})
    out_dir = plots.make_all(df, tmp_path, residue="TRP", interactions_df=int_df)
    assert (out_dir / "interaction_types.png").exists()
