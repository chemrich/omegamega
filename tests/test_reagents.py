"""Tests for reagent_quantities and write_reagent_summary."""

import pytest
import pandas as pd
from pricing import (
    reagent_quantities,
    write_reagent_summary,
    _ENZYME_CONC_U_PER_UL,
    _T4_LIGASE_CONC_U_PER_UL,
    _GG_ENZYME_U_PER_RXN,
    _GG_T4_LIGASE_U_PER_RXN,
    _KAPA_POLY_U_PER_RXN,
    _KAPA_BUFFER_UL_PER_RXN,
    _SPRI_RATIO,
    _PCR_RXN_VOL_UL,
)


def test_schema():
    df = reagent_quantities(3, "BsaI")
    assert list(df.columns) == ["reagent", "quantity_per_pool", "unit", "n_pools", "total_quantity", "notes"]
    assert len(df) == 5


def test_total_quantity_equals_per_pool_times_n_pools():
    for n in [1, 5, 37]:
        df = reagent_quantities(n, "BsaI")
        for _, row in df.iterrows():
            assert row["total_quantity"] == pytest.approx(row["quantity_per_pool"] * n)


def test_n_pools_column_is_uniform():
    df = reagent_quantities(7, "BsaI")
    assert (df["n_pools"] == 7).all()


@pytest.mark.parametrize("enzyme", ["BsaI", "BsmBI", "BbsI"])
def test_enzyme_row_units_and_quantity(enzyme):
    df = reagent_quantities(4, enzyme)
    row = df[df["reagent"] == enzyme].iloc[0]
    assert row["unit"] == "U"
    assert row["quantity_per_pool"] == pytest.approx(_GG_ENZYME_U_PER_RXN)
    assert row["total_quantity"] == pytest.approx(_GG_ENZYME_U_PER_RXN * 4)


@pytest.mark.parametrize("enzyme", ["BsaI", "BsmBI", "BbsI"])
def test_enzyme_notes_contain_correct_volume(enzyme):
    df = reagent_quantities(1, enzyme)
    row = df[df["reagent"] == enzyme].iloc[0]
    expected_ul = _GG_ENZYME_U_PER_RXN / _ENZYME_CONC_U_PER_UL[enzyme]
    assert str(round(expected_ul, 10))[:3] in row["notes"] or f"{expected_ul:.2g}" in row["notes"]


def test_t4_ligase_row():
    df = reagent_quantities(10)
    row = df[df["reagent"] == "T4 DNA Ligase"].iloc[0]
    assert row["unit"] == "U"
    assert row["quantity_per_pool"] == pytest.approx(_GG_T4_LIGASE_U_PER_RXN)
    assert row["total_quantity"] == pytest.approx(_GG_T4_LIGASE_U_PER_RXN * 10)
    assert f"{_T4_LIGASE_CONC_U_PER_UL:.0f}" in row["notes"]


def test_kapa_rows():
    df = reagent_quantities(6)
    poly = df[df["reagent"] == "KAPA HiFi HotStart Polymerase"].iloc[0]
    buf  = df[df["reagent"] == "KAPA HiFi Buffer (5×)"].iloc[0]
    assert poly["unit"] == "U"
    assert poly["quantity_per_pool"] == pytest.approx(_KAPA_POLY_U_PER_RXN)
    assert buf["unit"] == "µL"
    assert buf["quantity_per_pool"] == pytest.approx(_KAPA_BUFFER_UL_PER_RXN)


def test_spri_row():
    df = reagent_quantities(8)
    row = df[df["reagent"] == "SPRI beads"].iloc[0]
    assert row["unit"] == "µL"
    assert row["quantity_per_pool"] == pytest.approx(_SPRI_RATIO * _PCR_RXN_VOL_UL)
    assert row["total_quantity"] == pytest.approx(_SPRI_RATIO * _PCR_RXN_VOL_UL * 8)


def test_invalid_n_pools():
    with pytest.raises(ValueError, match="n_pools"):
        reagent_quantities(0)


def test_invalid_enzyme():
    with pytest.raises(ValueError, match="Unknown enzyme"):
        reagent_quantities(1, "NotAnEnzyme")


def test_write_reagent_summary(tmp_path):
    df = reagent_quantities(5, "BsmBI")
    out = write_reagent_summary(tmp_path, df)
    assert out.exists()
    loaded = pd.read_csv(out)
    assert list(loaded.columns) == list(df.columns)
    assert len(loaded) == len(df)
    assert (loaded["n_pools"] == 5).all()
