"""Unit tests for code/helpers.py.

Helpers are pure functions but several have subtle invariants:
- dna_contains_seq with reverse_complement=True must catch RC matches
  (CLAUDE.md: listing 'ATA' must also forbid 'TAT').
- count_sequence_element with overlapped=True must count all overlapping
  hits, not just non-overlapping ones.
- random_dna(0) must return '' rather than raising or producing a stray base.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from helpers import (
    _rc,
    count_sequence_element,
    dna_contains_seq,
    dynamic_chunker,
    flatten,
    index_array,
    overlap_split,
    random_dna,
    shuffler,
    unique_orthogonal,
)


# ----- _rc -----------------------------------------------------------------

@pytest.mark.parametrize("seq,expected", [
    ("ACGT", "ACGT"),       # palindrome
    ("AAAA", "TTTT"),
    ("ATCG", "CGAT"),
    ("", ""),
    ("acgt", "acgt"),       # case preserved (lowercase palindrome)
    ("AaCcGgTt", "aAcCgGtT"),
])
def test_rc(seq, expected):
    assert _rc(seq) == expected


def test_rc_round_trip():
    seq = "ATCGATCGGGCCAATT"
    assert _rc(_rc(seq)) == seq


# ----- dna_contains_seq ----------------------------------------------------

def test_dna_contains_seq_direct_match():
    assert dna_contains_seq("AAAGGGCCCTTT", "GGG") is True


def test_dna_contains_seq_rc_match_when_enabled():
    # 'ATA' RC is 'TAT'. A sequence containing only TAT should match when
    # reverse_complement=True, and miss when False. CLAUDE.md calls this
    # invariant out explicitly: listing 'ATA' as illegal forbids 'TAT' too.
    assert dna_contains_seq("CCCTATCCC", "ATA", reverse_complement=True) is True
    assert dna_contains_seq("CCCTATCCC", "ATA", reverse_complement=False) is False


def test_dna_contains_seq_case_insensitive():
    assert dna_contains_seq("aaagggccc", "GGG") is True


def test_dna_contains_seq_no_match():
    assert dna_contains_seq("AAAA", "GGG", reverse_complement=True) is False


def test_dna_contains_seq_multiple_elements():
    # any match returns True
    assert dna_contains_seq("AAAACCCC", "TTT", "GGG", reverse_complement=True) is True


# ----- count_sequence_element ---------------------------------------------

def test_count_sequence_element_counts_overlapping_matches():
    # 'AA' in 'AAAA' overlaps three times: positions 0,1,2.
    # The function also folds in the RC ('TT'), which does not appear here.
    assert count_sequence_element("AAAA", "AA") == 3


def test_count_sequence_element_includes_reverse_complement():
    # 'GGG' RC is 'CCC'. Both appear once each → 2.
    assert count_sequence_element("AAAGGGCCCTTT", "GGG") == 2


def test_count_sequence_element_zero_when_absent():
    assert count_sequence_element("AAAA", "GGG") == 0


# ----- unique_orthogonal --------------------------------------------------

def test_unique_orthogonal_passes_for_distinct_non_palindromic_sites():
    sites = pd.Series(["AAAA", "GGGG", "ATCG"])  # distinct + RCs distinct
    assert unique_orthogonal(sites) is True


def test_unique_orthogonal_fails_when_set_contains_a_site_and_its_rc():
    # 'AAAA' and 'TTTT' are RC partners → not orthogonal
    sites = pd.Series(["AAAA", "TTTT", "GGGG"])
    assert unique_orthogonal(sites) is False


def test_unique_orthogonal_fails_for_palindromic_site():
    # 'ACGT' is its own RC; including it folds the set on itself.
    sites = pd.Series(["AAAA", "ACGT"])
    assert unique_orthogonal(sites) is False


# ----- random_dna ---------------------------------------------------------

def test_random_dna_zero_returns_empty_string():
    assert random_dna(0) == ""


def test_random_dna_negative_returns_empty_string():
    assert random_dna(-5) == ""


@pytest.mark.parametrize("size", [1, 10, 100])
def test_random_dna_length_and_alphabet(size):
    seq = random_dna(size)
    assert len(seq) == size
    assert set(seq) <= set("ATCG")


# ----- dynamic_chunker, flatten, index_array, overlap_split, shuffler -----

def test_dynamic_chunker_partitions_by_chunk_sizes():
    out = list(dynamic_chunker(range(10), [2, 3, 5]))
    assert out == [[0, 1], [2, 3, 4], [5, 6, 7, 8, 9]]


def test_dynamic_chunker_raises_when_iterable_too_short():
    # PEP 479: StopIteration inside a generator surfaces as RuntimeError.
    with pytest.raises(RuntimeError):
        list(dynamic_chunker(range(3), [2, 3]))


def test_flatten_one_level():
    assert flatten([[1, 2], [3], [4, 5]]) == [1, 2, 3, 4, 5]


def test_index_array_marks_indices():
    arr = index_array([1, 3], length=5)
    assert isinstance(arr, np.ndarray)
    assert arr.tolist() == [0, 1, 0, 1, 0]


def test_overlap_split_shape():
    # overlap=1 returns sliding pairs (each wrapped in a single-element list)
    out = overlap_split([1, 2, 3, 4], overlap=1)
    assert out == [[[1, 2]], [[2, 3]], [[3, 4]]]


def test_shuffler_does_not_mutate_input():
    original = [1, 2, 3, 4, 5]
    snapshot = original.copy()
    shuffled = shuffler(original)
    assert original == snapshot
    assert sorted(shuffled) == snapshot
