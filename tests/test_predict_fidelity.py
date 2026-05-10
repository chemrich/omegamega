"""Unit tests for code/predict_fidelity.py.

The module is hot-loop code that was rewritten from a pandas/Bio.Seq path to
a numpy-cache path; CLAUDE.md flags it as performance-critical with a
non-obvious cache invalidation rule (the cache is keyed on `id(data)`).
These tests pin numerical correctness against a hand-computed 4-site
ligation matrix and guard the cross-instance cache invariant.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from predict_fidelity import (
    _TABLE_CACHE,
    correct_ligations,
    geneset_fidelity,
    predict_average,
    predict_fidelity,
    predict_minimum,
    predict_minimum_site,
    site_probability,
    total_ligations,
)


SITES = ["AAAA", "TTTT", "GGGG", "CCCC"]
# Hand-built ligation matrix.
# Rows/cols index SITES in order. Diagonal = self-ligation noise; the
# AAAA/TTTT and GGGG/CCCC pairs are reverse-complement partners with
# strong on-target counts (100 and 200 respectively).
_RAW = np.array(
    [
        [10, 100, 1, 2],     # AAAA
        [100, 5, 3, 4],      # TTTT (RC of AAAA)
        [1, 3, 20, 200],     # GGGG
        [2, 4, 200, 30],     # CCCC (RC of GGGG)
    ],
    dtype=np.int64,
)


def _matrix() -> pd.DataFrame:
    return pd.DataFrame(_RAW.copy(), index=SITES, columns=SITES)


# ----- correctness against hand-computed ground truth ----------------------

def test_correct_ligations_pulls_off_diagonal_pair():
    m = _matrix()
    assert correct_ligations("AAAA", "TTTT", m) == 100
    assert correct_ligations("TTTT", "AAAA", m) == 100
    assert correct_ligations("GGGG", "CCCC", m) == 200


def test_total_ligations_sums_site_against_pool_and_rcs():
    m = _matrix()
    # site AAAA (row 0) against ["AAAA", "GGGG"] + their RCs ["TTTT", "CCCC"]
    # = M[0, [0,2,1,3]].sum() = 10 + 1 + 100 + 2 = 113
    assert total_ligations("AAAA", ["AAAA", "GGGG"], m) == 113


def test_site_probability_matches_hand_computation():
    m = _matrix()
    # AAAA vs ["AAAA","GGGG"]: correct = M[0,1]+M[1,0] = 200,
    # total = (10+1+100+2)+(100+3+5+4) = 113+112 = 225 → 200/225
    p = site_probability("AAAA", ["AAAA", "GGGG"], m)
    assert p == pytest.approx(200 / 225)


def test_predict_fidelity_is_product_of_per_site_probabilities():
    m = _matrix()
    expected = (200 / 225) * (400 / 460)
    assert predict_fidelity(["AAAA", "GGGG"], m) == pytest.approx(expected)


def test_predict_minimum_site_picks_lowest_per_site_probability():
    m = _matrix()
    # 200/225 ≈ 0.8889; 400/460 ≈ 0.8696 → min is GGGG
    assert predict_minimum_site(["AAAA", "GGGG"], m) == pytest.approx(400 / 460)


def test_geneset_fidelity_returns_per_gene_products():
    m = _matrix()
    fids = geneset_fidelity([["AAAA"], ["GGGG"]], m)
    assert fids == pytest.approx([200 / 225, 400 / 460])


def test_predict_minimum_and_average_aggregate_across_genes():
    m = _matrix()
    fids = [200 / 225, 400 / 460]
    assert predict_minimum([["AAAA"], ["GGGG"]], m) == pytest.approx(min(fids))
    assert predict_average([["AAAA"], ["GGGG"]], m) == pytest.approx(sum(fids) / 2)


def test_accepts_numpy_array_of_sites():
    # SAPool.optimize passes site selections as numpy arrays of dtype object;
    # _to_idx must round-trip through .tolist().
    m = _matrix()
    sites = np.array(["AAAA", "GGGG"], dtype=object)
    assert predict_fidelity(sites, m) == pytest.approx((200 / 225) * (400 / 460))


# ----- cache behavior ------------------------------------------------------

def test_distinct_dataframe_instances_with_identical_contents_yield_same_result():
    """The cache is keyed on id(data). Two distinct DataFrames with identical
    contents should still produce identical fidelity (separate cache slots,
    same numeric output)."""
    df1 = _matrix()
    df2 = _matrix()
    assert id(df1) != id(df2)

    f1 = predict_fidelity(["AAAA", "GGGG"], df1)
    f2 = predict_fidelity(["AAAA", "GGGG"], df2)
    assert f1 == f2

    # Both DataFrames should now be cached separately.
    assert id(df1) in _TABLE_CACHE
    assert id(df2) in _TABLE_CACHE


def test_cache_reuses_view_across_calls():
    """Repeated calls with the same DataFrame instance must hit the cache,
    not rebuild the numpy view. We assert this by snapshotting the cached
    dict's identity."""
    df = _matrix()
    predict_fidelity(["AAAA", "GGGG"], df)
    cached_first = _TABLE_CACHE[id(df)]
    predict_fidelity(["TTTT", "CCCC"], df)
    cached_second = _TABLE_CACHE[id(df)]
    assert cached_first is cached_second
