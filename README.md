# omegamega

> Fork of [romerolab/omega](https://github.com/romerolab/omega) — adds cost accounting, validated primers, faster simulated annealing, and uv-based packaging. GPL-3.0-or-later.

[![DOI](https://zenodo.org/badge/964209533.svg)](https://doi.org/10.5281/zenodo.17637682)

OMEGA designs oligopools that assemble into custom gene libraries via Golden Gate, presented in [Freschlin et al. (bioRxiv 2025)](https://www.biorxiv.org/content/10.1101/2025.03.22.644747v1). This fork keeps the design algorithm and assembly protocol intact and layers on the things you need to actually order and execute a library.

## What this fork adds

- **Cost accounting.** Per-run `cost_summary.csv` with Twist oligo-pool list price (offline tier table), IDT primer-pair cost (offline rate card), and an optional live `OLIGO_POOLS_REGULAR` quote against the Twist API. See [docs/PRICING_NOTES.md](docs/PRICING_NOTES.md).
- **Wet-lab step counts.** Same summary reports the number of PCRs, cleanups, Golden Gate assemblies, and transformations to instantiate the library (3·N + 2 steps for N subpools).
- **Validated primer set.** [`data/subramanian_orthogonal.csv`](data/subramanian_orthogonal.csv) ships the 165 primers Subramanian et al. flagged as orthogonal in their Supplementary Table S1. Upstream's `data/test_primers.csv` ships all 20 primers they flagged as cross-reactive and is missing one validated primer (`subra_92`); we keep both for parity but recommend the validated set. See [docs/PRIMER_NOTES.md](docs/PRIMER_NOTES.md).
- **Faster simulated annealing.** `predict_fidelity` was the bottleneck; the inner loop now uses a numpy-indexed view of the ligation matrix cached per DataFrame.
- **FPbase benchmark corpus.** [`code/build_fpbase_corpus.py`](code/build_fpbase_corpus.py) pulls FPbase, codon-optimizes for E. coli (avoiding BsaI/BsmBI/BbsI sites), and bins by length for size-controlled benchmarking.
- **uv-based packaging.** `uv sync` instead of conda. Locked in [`uv.lock`](uv.lock).

## Install

```
uv sync
```

This creates `.venv/` from `pyproject.toml` + `uv.lock`. All commands below use `uv run` to execute inside that environment. If you prefer Colab, the `omega_google_colab.ipynb` notebook still works.

Verify with the bundled smoke test (~seconds):

```
uv run python ./code/omega.py genes --config configs/test_install.yml
```

## Designing a library

Two CLI subcommands do real work, plus one for re-pricing:

| command | purpose |
|---|---|
| `genes` | Full library design from a FASTA of codon-optimized genes. |
| `junctions` | Optimize a standalone GG junction set without library sequences (used for Fig. 2 of the paper). |
| `costs` | Re-price an existing output dir. |

The flow: copy [`configs/template.yml`](configs/template.yml), fill in your inputs, and run `genes`. Any config key can be overridden on the command line as `--<key> <value>`.

```
uv run python ./code/omega.py genes --config configs/template.yml \
    --input_seqs path/to/library.fasta \
    --njunctions 50 \
    --upstream_bbsite AATG --downstream_bbsite TTAG \
    --primers ./data/subramanian_orthogonal.csv \
    --nopt_steps 1000 --nopt_runs 5 --njobs 8
```

For 50 junctions, 5–10 independent optimizations of 1,000 steps per subpool is usually sufficient. Bumping `nopt_runs` beats bumping `nopt_steps` for marginal fidelity gains. Longer constructs need more fragments per gene, which means fewer genes per pool: `genes_per_pool ≈ (njunctions − 2) / (nfrags − 1)`.

A small full-library run (7 subpools, 50 junctions, 5 opt seeds):

```
uv run python ./code/omega.py genes --config configs/genes_test.yml --njobs 8
```

## Cost + wet-lab summary

`pricing_enabled` is on by default; the `genes` flow finishes with something like:

```
Twist offline list price: $6,181.00 (tier 5, len_251_300)
IDT primer pairs: $528.00 (55 pools at $9.60/pool)
Wet-lab steps: 167 total (55 PCRs + 55 PCR cleanups + 55 GG assemblies + 1 final cleanup + 1 transformation)
Cost summary saved to output/<run>/cost_summary.csv
```

Re-price an existing run without re-optimizing:

```
uv run python ./code/omega.py costs --output_dir output/<run>
```

Add a live Twist quote (requires `TWIST_JWT_TOKEN`, `TWIST_END_USER_TOKEN`, `TWIST_USER_EMAIL` env vars and a usable shipping address on your account):

```
uv run python ./code/omega.py costs --output_dir output/<run> --twist_quote true
```

The live quote was spot-validated against two production quotes on 2026-05-09 — Tier 1 (10 oligos) and Tier 5 (3,512 oligos) — both subtotals matched the offline table exactly.

For pricing schema, table sources, and the OLIGO_POOLS_REGULAR API specifics, see [docs/PRICING_NOTES.md](docs/PRICING_NOTES.md).

## Output files

Every `genes` run writes to the directory passed as `--output_dir`:

| file | contents |
|---|---|
| `oligo_order.csv` | Just `name,sequence` for the designed oligos. Submit directly to your oligo vendor. |
| `optimization_results.csv` | Per-gene record: name, submitted DNA, fragment oligos, fwd + rev primers, fidelity. |
| `pool_stats.csv` | Per-subpool record: fidelities, gene/site counts, optimization seed, enzyme, primer pair, primer cost. |
| `cost_summary.csv` | One-row library-level summary: oligo counts, Twist tier price, IDT primer total, wet-lab step counts, optional live-quote columns. |
| `experiment_details.txt` | Plain-text dump of the ligation-data experimental conditions used for fidelity scoring. |

Generating an IDT bulk-quote upload from your primer file:

```
uv run python scripts/generate_idt_quote_request.py \
    --input data/subramanian_orthogonal.csv \
    --output-prefix data/pricing/idt_quote_request_subramanian_orthogonal
```

Writes one `<prefix>_plate<N>.xls` per 96-well plate in IDT's plate-upload format.

## Fidelity scoring

OMEGA reports three fidelity metrics per pool:

- **`fidelity`** — the Pryor et al. ligation-fidelity calculation across all GG sites in the pool. Assumes a single sequential assembly; doesn't reflect OMEGA's combinatorial conditions. Used as the optimization objective.
- **`min_gene_fidelity`** — the lowest per-gene fidelity in the pool, accounting for the complex assembly background and the relevant gene length. Most relevant metric for OMEGA-style libraries.
- **`min_site_fidelity`** — the least orthogonal site selected.

Default ligation data is `T4_18h_37C` (Potapov et al.); the bundled BsaI cycling data and other Pryor et al. datasets are also available via `--ligation_data`. For the OMEGA paper conditions, **stick with BsaI + T4_18h_37C** unless you have a reason to deviate — it generated the highest-fidelity sites in the original characterization.

## Assembly protocol

Adapted from [Twist's oligopool amplification guidelines](https://www.twistbioscience.com/resources/protocol/twist-oligo-pool-amplification-guidelines) (FRM-001034 REV 8). The wet-lab counts in `cost_summary.csv` enumerate steps 2, 4, 5, 10, and 11.

1. Resuspend the lyophilized oligopool in 10 mM Tris pH 8.0 to ≥20 ng/µL (total yield in ng is printed on the shipping tube).
2. Set up one PCR reaction per subpool (Table 1).
3. Amplify with the protocol in Table 2. **Notes:** (a) annealing temperature must be optimized for your primer sequences — 61 °C is a tested starting point for Subramanian primers; (b) cycle count is length-dependent per Twist: 6–10 cycles for 20–100 nt oligos, 10–12 for 100–150 nt, 12–14 for 151–350 nt. Minimizing cycles avoids overamplification bias. (c) Expected yield assuming 1.8× per cycle: 10 ng × 1.8¹² ≈ 12 µg per 25 µL reaction at 12 cycles (theoretical; primer concentration typically caps practical yield at ~1–2 µg). Either way, far in excess of the ng-scale insert mass needed for cleanup (step 4) and assembly (Table 3).
4. Column cleanup of each PCR reaction individually.
5. Set up Golden Gate assembly per subpool (Table 3).
6. Digest 2 hr at 37 °C.
7. Add 1,000 U T4 Ligase ([NEB M0202T](https://www.neb.com/en-us/products/m0202-t4-dna-ligase)) — use the higher-concentration stock to keep volumes constant.
8. Ligate 18 hr at 37 °C.
9. Heat-inactivate T4 ligase 15 min at 65 °C.
10. Combine all subpool assemblies into a single tube; final column cleanup of the pooled library.
11. Transform.

If a downstream step needs PCR on the assembled library, run another digest afterward — empty vectors are smaller than complete assemblies and will be enriched in the PCR product.

**Table 1 — PCR setup** (per 25 µL reaction; KAPA HiFi HotStart PCR Kit, Roche)

| Component | Final concentration | Volume |
|--|--|--|
| 5× KAPA HiFi Fidelity Buffer | 1× | 5 µL |
| 10 mM dNTP mix | 0.3 mM each | 0.75 µL |
| Forward primer (10 µM) | 0.3 µM | 0.75 µL |
| Reverse primer (10 µM) | 0.3 µM | 0.75 µL |
| Oligopool (20 ng/µL) | 0.4 ng/µL (10 ng total) | 0.5 µL |
| KAPA HiFi HotStart DNA Polymerase (1 U/µL) | 0.5 U/reaction | 0.5 µL |
| Nuclease-free water | — | to 25 µL |

**Table 2 — PCR protocol**

| Step | Temperature | Time |
|--|--|--|
| Initial denaturation | 95 °C | 3 min |
| Denaturation | 98 °C | 20 sec |
| Annealing | optimum (~61 °C for Subramanian primers) | 15 sec |
| Extension | 72 °C | 15 sec |
| Repeat steps 2–4 | — | 6–14 cycles (length-dependent; see step 3) |
| Final extension | 72 °C | 1 min |

**Table 3 — Golden Gate assembly** (per 20 µL reaction; T4 ligase added after a 2-hr digest)

| Component | Amount |
|--|--|
| PCR product | 18:1 insert:vector molar ratio |
| Destination vector | 75 ng |
| BsaI (15 U/µL) | 15 U |
| T4 Ligase buffer (10×) | 2 µL |
| Nuclease-free water | to 20 µL |

## Options reference

Set in a YAML config or override on the command line.

**Required**

- `input_seqs` — FASTA of codon-optimized library sequences. Example in `data/fastas/`.
- `primers` — paired-primer CSV (`fwd_name,fwd_sequence,rev_name,rev_sequence`). Use `data/subramanian_orthogonal.csv` for new designs.
- `upstream_bbsite`, `downstream_bbsite` — vector ligation sites (e.g. `AATG`, `TTAG`).
- `njunctions` — number of GG sites per subpool, *including* the 2 backbone sites.

**Optional**

- `output_dir` (default `output`) — where to write outputs; created if missing.
- `enzyme` (default `BsaI`) — one of `BsaI`, `BsmBI`, `BbsI`.
- `ligation_data` (default `T4_18h_37C`) — `T4_{01h,18h}_{25C,37C}` (Potapov et al.) or `{BsaI,BbsI,BsmBI,Esp3I}_cycling` (Pryor et al.).
- `nopt_steps` (default 1000), `nopt_runs` (default 5) — opt budget per subpool. `opt_seeds` (mutually exclusive with `nopt_runs`) sets explicit seeds for reproducibility.
- `njobs` (default 1) — parallel optimization runs via joblib. Increase for faster pool optimization.
- `oligo_len` (default 300) — max oligo length.
- `add_primers` (default `true`), `pad_oligos` (default `true`) — controls whether primers are appended and whether random DNA pads oligos to uniform length.
- `illegal_dna_sequences` — sequences excluded from random padding (with reverse-complement awareness, so `'ATA'` also forbids `'TAT'`).
- `other_used_sites` — extra GG sites in your assembly that aren't backbone-vector sites.
- `pricing_enabled` (default `true`) — write `cost_summary.csv` and per-pool primer costs in `pool_stats.csv`. Set `false` to skip pricing entirely.
- `twist_quote` (default `false`) — also file a live Twist API quote and merge the parsed numbers into `cost_summary.csv`. Requires `TWIST_*` env vars (see [docs/PRICING_NOTES.md](docs/PRICING_NOTES.md)).

## Citing

If you use this code, please cite the OMEGA paper:

> Freschlin, C. R., Yang, K. K., Romero, P. A. *Scalable and cost-efficient custom gene library assembly from oligopools.* bioRxiv (2025). [doi:10.1101/2025.03.22.644747](https://doi.org/10.1101/2025.03.22.644747)

This fork has its own Zenodo DOI (above) for fork-specific changes (cost accounting, validated primers, refactored SA, uv packaging).

## References

- Pryor, J. M. et al. [Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized assembly design.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0238592) *PLoS One* 15, e0238592 (2020).
- Subramanian, S. K., Russ, W. P. & Ranganathan, R. [A set of experimentally validated, mutually orthogonal primers for combinatorially specifying genetic components.](https://academic.oup.com/synbio/article/3/1/ysx008/4817474) *Synth. Biol.* 3, ysx008 (2018).
- Potapov, V. et al. [Comprehensive profiling of four base overhang ligation fidelity by T4 DNA ligase and application to DNA assembly.](https://pubs.acs.org/doi/10.1021/acssynbio.8b00333) *ACS Synth. Biol.* 7, 2665–2674 (2018).
