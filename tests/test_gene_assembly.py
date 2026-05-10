"""Round-trip test for Gene fragmentation + reassembly.

The full pipeline is end-to-end-tested only by the test_install.yml smoke run,
which is slow and not part of `uv run pytest`. This test pins the central
correctness property (a Gene's oligos reassemble back to the original
sequence and translate to the same protein) on a tiny synthetic gene with
a manually pinned GG site, so we don't depend on the SA optimizer.
"""
from __future__ import annotations

import pandas as pd
import pytest

from data_classes import EnzymeTypes, Primer, define_enzyme
from library_classes import Gene


# 21 bp = 7 codons: M K L G T K *
SEQ = "ATGAAACTGGGTACCAAATAA"
# At pos=5 the in-sequence 4-mer is "ACTG" — not in the disallowed-palindrome
# list and not the RC of either backbone site.
GGSITE = "ACTG"
GGPOS = 5
UPSTREAM_BB = "GAGA"
DOWNSTREAM_BB = "GTGT"


def _build_gene(oligo_len: int = 100) -> Gene:
    enzyme = define_enzyme(EnzymeTypes.BsaI)
    # Primers can be any non-enzyme-site sequence; their content doesn't
    # affect fragmentation or assembly correctness.
    fwd = Primer("fwd", "ATCG" * 5, True)
    rev = Primer("rev", "TGCA" * 5, False)
    gene = Gene(
        name="test_gene",
        sequence=SEQ,
        upstream_bbsite=UPSTREAM_BB,
        downstream_bbsite=DOWNSTREAM_BB,
        other_used_sites=[],
        enzyme=enzyme,
        forward_primer=fwd,
        reverse_primer=rev,
        oligo_len=oligo_len,
        illegal_dna_sequences=tuple(),
    )
    # Pin the GG site assignment that the SA optimizer would otherwise pick.
    gene.assigned_sites = pd.DataFrame({"ggsite": [GGSITE], "pos": [GGPOS]})
    return gene


def test_gene_assembles_returns_true_for_valid_assignment():
    gene = _build_gene()
    assert gene.gene_assembles() is True


def test_oligos_reassemble_to_original_sequence_with_backbones():
    gene = _build_gene()
    enzyme = gene.enzyme
    enz_flank = len(enzyme.seq) + enzyme.padding  # 7 nt at each end
    site_size = enzyme.site_size  # 4

    # Skip primers and padding to keep the assembly assertion exact; both
    # are bookkeeping around the assembly substring, not part of the
    # round-trip property under test.
    oligos = gene.get_oligos(add_primers=False, pad_oligo=False)
    assert len(oligos) == 2

    # Strip the BsaI flanking site + 1nt pad from each oligo.
    bare = [o[enz_flank:-enz_flank] for o in oligos]

    # Reassemble: consecutive fragments overlap by `site_size` (the GG site).
    assembled = bare[0] + bare[1][site_size:]

    # Strip the upstream/downstream backbone GG sites.
    assert assembled.startswith(UPSTREAM_BB)
    assert assembled.endswith(DOWNSTREAM_BB)
    assert assembled[site_size:-site_size] == SEQ


def test_gene_assembles_raises_when_assigned_site_does_not_match_sequence():
    """If assigned_sites contains a 4-mer that isn't actually at that
    position in the gene, the boundary check between consecutive fragments
    fails and gene_assembles raises."""
    gene = _build_gene()
    # 'TTTT' is not at position 5 in SEQ, so the GG site overlap won't match.
    gene.assigned_sites = pd.DataFrame({"ggsite": ["TTTT"], "pos": [GGPOS]})

    with pytest.raises(ValueError, match="Golden Gate sites do not match"):
        gene.gene_assembles()
