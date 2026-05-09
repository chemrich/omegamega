# Subramanian primer set notes

OMEGA uses 20-mer orthogonal primers from Subramanian, Russ & Ranganathan (2018) for sub-pool amplification ("subra_*" primers). This note records what's in this repo vs. the upstream paper.

## Source

> Subramanian, S. K., Russ, W. P. & Ranganathan, R. *A set of experimentally validated, mutually orthogonal primers for combinatorially specifying genetic components.* Synth. Biol. 3(1), ysx008 (2018). [DOI: 10.1093/synbio/ysx008](https://doi.org/10.1093/synbio/ysx008)

The paper: 185 candidate primers were designed and tested. After interaction-matrix analysis at a 0.95 dissimilarity threshold, the authors identified a validated **orthogonal subset** to use in combinatorial assemblies. Supplementary Table S1 (xlsx) gives the keep/drop call per primer.

The xlsx flags **165 primers as `Yes` (keep)** and **20 as `-` (drop)**. The paper text says "166 mutually orthogonal" — there is a one-off mismatch between the manuscript and the supplemental table; the table is the authoritative list.

## Reference materials in this repo

- `docs/Subramanian A set of experimentally validated, mutually orthogonal primers for combinatorially specifying genetic components.pdf` — the paper
- `docs/subramanian_supp.pdf` — supplementary material (template oligo design + full tree figure)
- `docs/subramanian_supp_st_1.xlsx` — Supplementary Table S1, the keep/drop list

## Files derived from these references

- `data/subramanian_185_primers.csv` — all 185 candidates (`name,sequence,is_orthogonal`)
- `data/subramanian_orthogonal.csv` — the 165-primer validated set in OMEGA's paired `fwd_name,fwd_sequence,rev_name,rev_sequence` format (82 fwd + 83 rev). Drop-in for `--primers`.

## What's in `data/test_primers.csv` (upstream)

The bundled OMEGA primer file ships **184 of the 185** Subramanian candidates: 92 in the `fwd_*` columns (`subra_0`–`subra_91`) and 92 in the `rev_*` columns (`subra_93`–`subra_184`). It omits `subra_92` (`ATGCTAGCTGGAACTATCGG`).

Cross-referencing against Supplementary Table S1:

- **All 20 non-orthogonal primers (`-` in S1) are present** in `test_primers.csv`. 10 are in the fwd column, 10 in the rev column.
- **The single primer `test_primers.csv` omits — `subra_92` — is one of the validated orthogonal `Yes` primers**, not one of the dropped 20.

So the upstream file simultaneously (a) ships every primer Subramanian flagged as cross-reactive, and (b) drops one of the validated orthogonal ones. Net composition: 184 = 164 orthogonal + 20 non-orthogonal.

### The 20 non-orthogonal primers shipped in `test_primers.csv`

| Subramanian id | name in repo | column | sequence |
|---|---|---|---|
| 7 | subra_6 | fwd | AACTCCATCGGACTAATGCG |
| 8 | subra_7 | fwd | AACTCGATAACAGGGAACCC |
| 9 | subra_8 | fwd | AACTGACTTGTAAAACGCGC |
| 10 | subra_9 | fwd | AACTTGTAGATAGACGCCGG |
| 19 | subra_18 | fwd | AATAGAGCACGGAAGCCAAC |
| 66 | subra_65 | fwd | AGTACGGACACCAAGATTGC |
| 72 | subra_71 | fwd | AGTCCATTGAACACAGGAGC |
| 76 | subra_75 | fwd | AGTGCAAAGCCAGACAGTAG |
| 82 | subra_81 | fwd | ATAAGCGAGAGTTTCCCTCC |
| 87 | subra_86 | fwd | ATATCTCGAAGGAACTGCGG |
| 114 | subra_113 | rev | CTAGGAAAAGGCACCTGACC |
| 118 | subra_117 | rev | CTCACGCAGATAGTACGGTG |
| 138 | subra_137 | rev | TAACAGATTGCCAGCTCGAC |
| 143 | subra_142 | rev | TAATGATCTCACGCCTGACG |
| 155 | subra_154 | rev | TATCCAGAACAGGCATTGGC |
| 161 | subra_160 | rev | TCAGGGATGACCGAAAAGTC |
| 168 | subra_167 | rev | TCTCCTTGGTGAAAAGAGCG |
| 176 | subra_175 | rev | TGGACCAATAAGCCACAGTG |
| 179 | subra_178 | rev | TTAAGGGGCAGTACGAATCC |
| 182 | subra_181 | rev | TTCATAAGGGAATCCACGCC |

## Recommendation

For new library designs, use `data/subramanian_orthogonal.csv` (165 validated primers). The existing `data/test_primers.csv` is preserved as-is to avoid silently changing behavior for runs that pin to it.

## Verification

The diff was produced by parsing `docs/subramanian_supp_st_1.xlsx` directly and matching sequences against `data/test_primers.csv` with a CSV parser (after stripping CRLF line endings). Every dropped primer matched by exact sequence; every match also lined up with the expected `subra_(id-1)` name.
