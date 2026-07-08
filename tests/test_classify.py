"""Classification helpers."""

from pdbbind_residue_contacts import classify


def test_water():
    assert classify.is_water("HOH")
    assert classify.is_water("hoh")
    assert classify.is_water("DOD")
    assert not classify.is_water("EDO")


def test_metal():
    assert classify.is_metal("ZN")
    assert classify.is_metal("ca")
    assert not classify.is_metal("C")


def test_halide():
    assert classify.is_halide("CL")
    assert not classify.is_halide("N")


def test_common_additive():
    assert classify.is_common_additive("EDO")
    assert classify.is_common_additive("GOL")
    assert not classify.is_common_additive("STI")  # a real drug ligand


def test_class_labels_distinct():
    assert classify.DRUG_LIKE != classify.SOLVENT
