"""Unit tests for code/pricing.py.

The two parser fixtures mimic the structure of real Twist OLIGO_POOLS_REGULAR
quote responses captured in scripts/_probe_responses/quote_prod_*.json on
2026-05-09 (Tier 1 / 10-oligo and Tier 5 / 3512-oligo). Only the fields the
parser reads are reproduced; PII (name, address, phone) is omitted.
"""
from __future__ import annotations

import pandas as pd
import pytest

from pricing import (
    DEFAULT_IDT_PRIMER_TABLE,
    DEFAULT_TWIST_TABLE,
    IDTPrimerPricing,
    TwistOligoPoolPricing,
    annotate_pool_stats_with_primer_cost,
    cost_summary,
    parse_oligo_pool_quote,
    wet_lab_steps,
)


@pytest.fixture(scope="module")
def pricing() -> TwistOligoPoolPricing:
    return TwistOligoPoolPricing.from_csv(DEFAULT_TWIST_TABLE)


@pytest.mark.parametrize(
    "n_oligos,max_len,expected_tier,expected_price",
    [
        (10, 300, 1, 1030.0),       # Tier 1, 251–300 nt — matches real prod quote
        (3512, 300, 5, 6181.0),     # Tier 5, 251–300 nt — matches real prod quote
        (10, 350, 1, 1288.0),       # Tier 1, 301–350 nt — matches real prod quote
        (3512, 350, 5, 7727.0),     # Tier 5, 301–350 nt — matches real prod quote
        (2, 120, 1, 400.0),         # exact tier_min, exact length max
        (100, 120, 1, 400.0),       # exact tier_max boundary
        (101, 120, 2, 800.0),       # one over tier 1 -> tier 2
        (50, 250, 1, 689.0),        # mid-tier, mid-length
        (500, 350, 2, 2575.0),      # tier 2, max length bin
        (696000, 350, 28, 143028.0),  # absolute max
    ],
)
def test_offline_quote(pricing, n_oligos, max_len, expected_tier, expected_price):
    q = pricing.quote(n_oligos, max_len)
    assert q.tier == expected_tier
    assert q.pool_price_usd == expected_price


def test_below_min_tier_raises(pricing):
    # Twist's smallest pool is 2 oligos; n_oligos=1 is not orderable.
    with pytest.raises(ValueError, match="not in any Twist tier"):
        pricing.quote(1, 300)


def test_length_bin_picks_smallest_fitting(pricing):
    q = pricing.quote(50, 121)
    assert q.length_bin == "len_121_150"
    q = pricing.quote(50, 200)
    assert q.length_bin == "len_151_200"


def test_length_bin_lower_edge(pricing):
    # Lengths at the bottom of each bin should still pick that bin, not
    # accidentally bleed into the bin below.
    assert pricing.quote(50, 20).length_bin == "len_20_120"   # absolute floor
    assert pricing.quote(50, 121).length_bin == "len_121_150"  # one over floor
    assert pricing.quote(50, 151).length_bin == "len_151_200"
    assert pricing.quote(50, 201).length_bin == "len_201_250"
    assert pricing.quote(50, 251).length_bin == "len_251_300"
    assert pricing.quote(50, 301).length_bin == "len_301_350"


def test_above_max_length_raises(pricing):
    with pytest.raises(ValueError, match="exceeds Twist's max"):
        pricing.quote(50, 351)


def test_above_max_oligos_raises(pricing):
    with pytest.raises(ValueError, match="not in any Twist tier"):
        pricing.quote(696001, 300)


def test_zero_or_negative_raises(pricing):
    with pytest.raises(ValueError):
        pricing.quote(0, 300)
    with pytest.raises(ValueError):
        pricing.quote(10, 0)


# ----- parser tests against minimal mock responses --------------------------

def _mock_quote(unit_price: float, n_oligos: int, length_label: str,
                tier_label: str, subtotal: float, tax: float, total: float,
                business_days: int = 4) -> dict:
    return {
        "tat": {"business_days": business_days},
        "quote": {
            "subtotal": subtotal,
            "tax_total": tax,
            "price": total,
            "quote_lines": [
                {"product_code": "SHIP",
                 "description": "Shipping Costs",
                 "list_unit_price": "35.00", "quantity": "1.00"},
                {"product_code": "Handling",
                 "description": "Handling Costs",
                 "list_unit_price": "25.00", "quantity": "1.00"},
                {"product_code": "104060",
                 "description": f"Oligo Pool {tier_label} ({n_oligos} Oligos) {length_label}",
                 "list_unit_price": str(unit_price), "quantity": "1.00"},
            ],
        },
    }


def test_parse_tier1_shape():
    parsed = parse_oligo_pool_quote(
        _mock_quote(1030.0, 10, "251-300nt", "Tier 1 (2-100 Oligos)",
                    subtotal=1030.0, tax=100.43, total=1190.43)
    )
    assert parsed["pool_subtotal_usd"] == 1030.0
    assert parsed["shipping_usd"] == 35.0
    assert parsed["handling_usd"] == 25.0
    assert parsed["subtotal_usd"] == 1030.0
    assert parsed["total_price_usd"] == 1190.43
    assert parsed["business_days"] == 4
    assert len(parsed["pool_lines"]) == 1
    assert parsed["pool_lines"][0]["unit_price_usd"] == 1030.0


def test_parse_handles_missing_optional_fields():
    minimal = {"quote": {"quote_lines": []}}
    parsed = parse_oligo_pool_quote(minimal)
    assert parsed["pool_subtotal_usd"] == 0.0
    assert parsed["shipping_usd"] == 0.0
    assert parsed["business_days"] is None
    assert parsed["pool_lines"] == []


def test_parse_aggregates_multiple_pool_line_items():
    # If a single quote response carries more than one 104060 line (e.g.
    # multi-pool order), pool_subtotal_usd must be the sum across them and
    # pool_lines must hold all entries.
    response = {
        "tat": {"business_days": 5},
        "quote": {
            "subtotal": 2318.0,
            "tax_total": 0.0,
            "price": 2378.0,
            "quote_lines": [
                {"product_code": "SHIP", "description": "Shipping",
                 "list_unit_price": "35.00", "quantity": "1.00"},
                {"product_code": "Handling", "description": "Handling",
                 "list_unit_price": "25.00", "quantity": "1.00"},
                {"product_code": "104060",
                 "description": "Oligo Pool Tier 1 (10 Oligos) 251-300nt",
                 "list_unit_price": "1030.00", "quantity": "1.00"},
                {"product_code": "104060",
                 "description": "Oligo Pool Tier 1 (10 Oligos) 301-350nt",
                 "list_unit_price": "1288.00", "quantity": "1.00"},
            ],
        },
    }
    parsed = parse_oligo_pool_quote(response)
    assert len(parsed["pool_lines"]) == 2
    assert parsed["pool_subtotal_usd"] == pytest.approx(1030.0 + 1288.0)
    assert parsed["shipping_usd"] == 35.0
    assert parsed["handling_usd"] == 25.0


# ----- cost_summary integration --------------------------------------------
#
# These tests use synthetic Twist + IDT rate cards built into tmp_path
# rather than the bundled CSVs at data/pricing/. The bundled rate cards
# are spot-validated separately in test_offline_quote / test_idt_quote_primer
# (the anchor tests); the cost_summary path is testing aggregation logic,
# not rate-card values, so it shouldn't break when the bundled cards are
# refreshed.

@pytest.fixture
def synthetic_twist_csv(tmp_path):
    path = tmp_path / "twist.csv"
    pd.DataFrame([
        {"tier": 1, "tier_min": 2, "tier_max": 100,
         "len_20_120": 100.0, "len_121_150": 200.0, "len_151_200": 300.0,
         "len_201_250": 400.0, "len_251_300": 500.0, "len_301_350": 600.0},
        {"tier": 2, "tier_min": 101, "tier_max": 1000,
         "len_20_120": 1000.0, "len_121_150": 2000.0, "len_151_200": 3000.0,
         "len_201_250": 4000.0, "len_251_300": 5000.0, "len_301_350": 6000.0},
    ]).to_csv(path, index=False)
    return path


@pytest.fixture
def synthetic_idt_csv(tmp_path):
    path = tmp_path / "idt.csv"
    # $1/bp at 25nmole/STD/plate keeps the arithmetic obvious.
    pd.DataFrame([{
        "scale": "25nmole", "purification": "STD", "format": "plate",
        "length_min": 1, "length_max": 60,
        "price_per_base_usd": 1.00, "plate_setup_usd": 0.0,
        "source": "synthetic-fixture",
    }]).to_csv(path, index=False)
    return path


def test_cost_summary_offline_only(synthetic_twist_csv):
    df = pd.DataFrame({
        "name": [f"oligo_{i}" for i in range(50)],
        "sequence": ["A" * 280] * 50,
    })
    summary = cost_summary(df, table_path=synthetic_twist_csv, twist_quote=False)
    assert summary["n_oligos"] == 50
    assert summary["max_oligo_len_nt"] == 280
    assert summary["offline_tier"] == 1
    assert summary["offline_length_bin"] == "len_251_300"
    assert summary["offline_pool_price_usd"] == 500.0  # synthetic table value
    assert "live_pool_subtotal_usd" not in summary
    assert "primers_total_usd" not in summary  # no pool_stats_df


# ----- IDT primer pricing ---------------------------------------------------

@pytest.fixture(scope="module")
def idt_pricing() -> IDTPrimerPricing:
    return IDTPrimerPricing.from_csv(DEFAULT_IDT_PRIMER_TABLE)


@pytest.mark.parametrize("length,expected", [
    (20, 4.80),    # subramanian primer length — anchored to real IDT quote
    (1, 0.24),     # min length, $0.24/bp
    (60, 14.40),   # max length supported by the bundled row (60 bp)
    (10, 2.40),
])
def test_idt_quote_primer(idt_pricing, length, expected):
    q = idt_pricing.quote_primer(length)
    assert q.price_usd == pytest.approx(expected)
    assert q.scale == "25nmole"
    assert q.purification == "STD"
    assert q.format == "plate"
    assert q.per_base_usd == 0.24
    assert q.plate_setup_usd == 0.0


def test_idt_quote_primer_unknown_config_raises(idt_pricing):
    with pytest.raises(ValueError, match="No IDT pricing row"):
        idt_pricing.quote_primer(20, scale="100nmole")


def test_idt_quote_primer_zero_length_raises(idt_pricing):
    with pytest.raises(ValueError):
        idt_pricing.quote_primer(0)


def test_idt_quote_pool_pair_anchor(idt_pricing):
    # Subramanian primers are all 20 nt → $4.80 each → $9.60 per pair
    pair = idt_pricing.quote_pool_pair(20, 20)
    assert pair["fwd_primer_cost_usd"] == 4.80
    assert pair["rev_primer_cost_usd"] == 4.80
    assert pair["pair_cost_usd"] == 9.60


def test_idt_quote_pool_pair_mismatched_lengths(idt_pricing):
    # Asymmetric pair: 18 nt fwd + 22 nt rev should price each side
    # independently and sum, not e.g. quote both at the average length.
    pair = idt_pricing.quote_pool_pair(18, 22)
    assert pair["fwd_primer_cost_usd"] == pytest.approx(0.24 * 18)
    assert pair["rev_primer_cost_usd"] == pytest.approx(0.24 * 22)
    assert pair["pair_cost_usd"] == pytest.approx(0.24 * 40)


def test_annotate_pool_stats_with_primer_cost():
    df = pd.DataFrame({
        "pool": [0, 1],
        "pfwd_sequence": ["A" * 20, "G" * 20],
        "prev_sequence": ["T" * 20, "C" * 20],
    })
    out = annotate_pool_stats_with_primer_cost(df)
    assert (out["fwd_primer_cost_usd"] == 4.80).all()
    assert (out["rev_primer_cost_usd"] == 4.80).all()
    assert (out["primer_pair_cost_usd"] == 9.60).all()
    # original columns preserved
    assert "pool" in out.columns
    assert "pfwd_sequence" in out.columns


def test_cost_summary_with_pool_stats(synthetic_twist_csv, synthetic_idt_csv):
    oligo_df = pd.DataFrame({
        "name": [f"o{i}" for i in range(60)],
        "sequence": ["A" * 280] * 60,
    })
    pool_stats_df = pd.DataFrame({
        "pool": [0, 1, 2],
        "pfwd_sequence": ["A" * 20] * 3,
        "prev_sequence": ["T" * 20] * 3,
    })
    summary = cost_summary(
        oligo_df,
        table_path=synthetic_twist_csv,
        idt_table_path=synthetic_idt_csv,
        pool_stats_df=pool_stats_df,
    )
    assert summary["n_pools"] == 3
    assert summary["n_primer_pairs"] == 3
    # synthetic IDT card: $1/bp × 20 nt × 2 primers = $40/pair
    assert summary["primers_per_pool_avg_usd"] == pytest.approx(40.00)
    assert summary["primers_total_usd"] == pytest.approx(120.00)
    assert summary["idt_scale"] == "25nmole"
    # wet-lab counts piggy-back on pool_stats_df availability
    assert summary["wetlab_pcrs"] == 3
    assert summary["wetlab_pcr_cleanups"] == 3
    assert summary["wetlab_quants"] == 3
    assert summary["wetlab_assembly_reactions"] == 3
    assert summary["wetlab_final_cleanups"] == 1
    assert summary["wetlab_transformations"] == 1
    assert summary["wetlab_total_steps"] == 14  # 4*3 + 2


# ----- wet-lab step counts --------------------------------------------------

@pytest.mark.parametrize("n_pools,total", [
    (1, 6),     # 4 + 2
    (7, 30),    # 7-pool gfp library
    (37, 150),  # 37-pool fpbase library (avgFP corpus, 350 nt default)
])
def test_wet_lab_steps_total(n_pools, total):
    s = wet_lab_steps(n_pools)
    assert s["wetlab_total_steps"] == total
    assert s["wetlab_pcrs"] == n_pools
    assert s["wetlab_pcr_cleanups"] == n_pools
    assert s["wetlab_quants"] == n_pools
    assert s["wetlab_assembly_reactions"] == n_pools
    assert s["wetlab_final_cleanups"] == 1
    assert s["wetlab_transformations"] == 1


def test_wet_lab_steps_zero_raises():
    with pytest.raises(ValueError):
        wet_lab_steps(0)
    with pytest.raises(ValueError):
        wet_lab_steps(-1)
