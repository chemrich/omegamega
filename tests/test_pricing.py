"""Unit tests for code/pricing.py.

The two parser fixtures mimic the structure of real Twist OLIGO_POOLS_REGULAR
quote responses captured in scripts/_probe_responses/quote_prod_*.json on
2026-05-09 (Tier 1 / 10-oligo and Tier 5 / 3512-oligo). Only the fields the
parser reads are reproduced; PII (name, address, phone) is omitted.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))

from pricing import (  # noqa: E402
    DEFAULT_IDT_PRIMER_TABLE,
    DEFAULT_TWIST_TABLE,
    IDTPrimerPricing,
    TwistOligoPoolPricing,
    annotate_pool_stats_with_primer_cost,
    cost_summary,
    parse_oligo_pool_quote,
)


@pytest.fixture(scope="module")
def pricing() -> TwistOligoPoolPricing:
    return TwistOligoPoolPricing.from_csv(DEFAULT_TWIST_TABLE)


@pytest.mark.parametrize(
    "n_oligos,max_len,expected_tier,expected_price",
    [
        (10, 300, 1, 1030.0),       # Tier 1, 251–300 nt — matches real prod quote
        (3512, 300, 5, 6181.0),     # Tier 5, 251–300 nt — matches real prod quote
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


# ----- cost_summary integration --------------------------------------------

def test_cost_summary_offline_only(tmp_path):
    df = pd.DataFrame({
        "name": [f"oligo_{i}" for i in range(50)],
        "sequence": ["A" * 280] * 50,
    })
    summary = cost_summary(df, twist_quote=False)
    assert summary["n_oligos"] == 50
    assert summary["max_oligo_len_nt"] == 280
    assert summary["offline_tier"] == 1
    assert summary["offline_length_bin"] == "len_251_300"
    assert summary["offline_pool_price_usd"] == 1030.0
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


def test_cost_summary_with_pool_stats():
    oligo_df = pd.DataFrame({
        "name": [f"o{i}" for i in range(60)],
        "sequence": ["A" * 280] * 60,
    })
    pool_stats_df = pd.DataFrame({
        "pool": [0, 1, 2],
        "pfwd_sequence": ["A" * 20] * 3,
        "prev_sequence": ["T" * 20] * 3,
    })
    summary = cost_summary(oligo_df, pool_stats_df=pool_stats_df)
    assert summary["n_pools"] == 3
    assert summary["n_primer_pairs"] == 3
    assert summary["primers_per_pool_avg_usd"] == pytest.approx(9.60)
    assert summary["primers_total_usd"] == pytest.approx(28.80)
    assert summary["idt_scale"] == "25nmole"
