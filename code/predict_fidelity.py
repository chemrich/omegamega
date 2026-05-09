"""Fidelity prediction for Golden Gate site sets.

The original implementation used pandas .loc with string indices and
Bio.Seq.reverse_complement inside the SA inner loop, which dominated runtime.
This module precomputes a numpy-indexed view of the ligation matrix once per
input DataFrame and serves all queries from that cached view. Public
signatures are unchanged.
"""

from __future__ import annotations

from itertools import chain
from typing import Iterable

import numpy as np
import pandas as pd

_RC_TABLE = str.maketrans("ACGT", "TGCA")


def _rc(site: str) -> str:
    return site.translate(_RC_TABLE)[::-1]


_TABLE_CACHE: dict[int, dict] = {}


def _get_table(data: pd.DataFrame) -> dict:
    """Return a cached numpy-indexed view of `data`, keyed by DataFrame id.

    The view exposes:
        idx_of: dict[str, int] — site string -> canonical index
        matrix: np.ndarray (N,N) of ligation counts in canonical ordering
        rc_idx: np.ndarray (N,) where rc_idx[i] = canonical idx of reverse complement of site i
    """
    key = id(data)
    cached = _TABLE_CACHE.get(key)
    if cached is not None:
        return cached

    cols = list(data.columns)
    idx_of = {s: i for i, s in enumerate(cols)}
    matrix = data.reindex(index=cols, columns=cols).to_numpy(dtype=np.int64, copy=True)
    rc_idx = np.array([idx_of[_rc(s)] for s in cols], dtype=np.int64)

    cached = {"idx_of": idx_of, "matrix": matrix, "rc_idx": rc_idx}
    _TABLE_CACHE[key] = cached
    return cached


def _to_idx(sites: Iterable[str], idx_of: dict[str, int]) -> np.ndarray:
    if isinstance(sites, np.ndarray):
        sites_iter = sites.tolist()
    else:
        sites_iter = list(sites)
    return np.fromiter((idx_of[s] for s in sites_iter), dtype=np.int64, count=len(sites_iter))


def correct_ligations(site: str, wc_site: str, data_matrix: pd.DataFrame) -> int:
    """Number of correct ligations observed for a GG site."""
    table = _get_table(data_matrix)
    return int(table["matrix"][table["idx_of"][site], table["idx_of"][wc_site]])


def total_ligations(site: str, sites: Iterable[str], data_matrix: pd.DataFrame) -> int:
    """Total ligation events involving `site` against `sites` and their RCs."""
    table = _get_table(data_matrix)
    s_idx = table["idx_of"][site]
    sites_idx = _to_idx(sites, table["idx_of"])
    targets = np.concatenate([sites_idx, table["rc_idx"][sites_idx]])
    return int(table["matrix"][s_idx, targets].sum())


def site_probability(site: str, sites: Iterable[str], data_matrix: pd.DataFrame) -> float:
    """Probability of site orthogonality vs. proposed site set."""
    table = _get_table(data_matrix)
    matrix = table["matrix"]
    rc_idx = table["rc_idx"]
    s_idx = table["idx_of"][site]
    wc_idx = rc_idx[s_idx]

    sites_idx = _to_idx(sites, table["idx_of"])
    targets = np.concatenate([sites_idx, rc_idx[sites_idx]])

    correct = matrix[s_idx, wc_idx] + matrix[wc_idx, s_idx]
    total = matrix[s_idx, targets].sum() + matrix[wc_idx, targets].sum()
    return float(correct / total)


def _vectorized_probabilities(sites_idx: np.ndarray, table: dict) -> np.ndarray:
    """Return per-site orthogonality probabilities for an integer-indexed site array."""
    matrix = table["matrix"]
    rc_idx = table["rc_idx"]
    wc_idx = rc_idx[sites_idx]
    targets = np.concatenate([sites_idx, wc_idx])

    correct = matrix[sites_idx, wc_idx] + matrix[wc_idx, sites_idx]
    total = matrix[sites_idx][:, targets].sum(axis=1) + matrix[wc_idx][:, targets].sum(axis=1)
    return correct / total


def predict_fidelity(sites: Iterable[str], data: pd.DataFrame) -> float:
    """Predicted fidelity for a set of GG sites (product of per-site probabilities)."""
    table = _get_table(data)
    sites_idx = _to_idx(sites, table["idx_of"])
    probs = _vectorized_probabilities(sites_idx, table)
    return float(np.prod(probs))


def predict_minimum_site(sites: Iterable[str], data: pd.DataFrame) -> float:
    """Least orthogonal site in the proposed set."""
    table = _get_table(data)
    sites_idx = _to_idx(sites, table["idx_of"])
    probs = _vectorized_probabilities(sites_idx, table)
    return float(probs.min())


def geneset_fidelity(gene_sites: Iterable[Iterable[str]], data: pd.DataFrame) -> list[float]:
    """Per-gene fidelities given a set of genes (typically a pool)."""
    gene_sites_list = [list(g) for g in gene_sites]
    table = _get_table(data)
    pool_sites = list(chain.from_iterable(gene_sites_list))
    pool_idx = _to_idx(pool_sites, table["idx_of"])
    pool_probs = _vectorized_probabilities(pool_idx, table)

    out: list[float] = []
    cursor = 0
    for gene in gene_sites_list:
        n = len(gene)
        out.append(float(np.prod(pool_probs[cursor : cursor + n])))
        cursor += n
    return out


def predict_minimum(gene_sites: Iterable[Iterable[str]], data: pd.DataFrame) -> float:
    """Lowest fidelity for any gene in a pool."""
    return min(geneset_fidelity(gene_sites, data))


def predict_average(gene_sites: Iterable[Iterable[str]], data: pd.DataFrame) -> float:
    """Average fidelity across genes in a pool."""
    fids = geneset_fidelity(gene_sites, data)
    return sum(fids) / len(fids)
