# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository

OMEGAMEGA is a fork of [romerolab/omega](https://github.com/romerolab/omega). It designs oligopools for scalable Golden Gate gene library assembly (Freschlin et al. bioRxiv 2025). Licensed GPL-3.0-or-later.

## Environment & commands

This fork uses **uv** (not conda). All Python is run through `uv run` inside `.venv/`.

```
uv sync                                                           # install/refresh deps from uv.lock
uv run python ./code/omega.py genes --config configs/test_install.yml   # fast smoke test
uv run python ./code/omega.py genes --config configs/genes_test.yml     # small full run
uv run pytest                                                     # there are no committed tests yet
uv run ruff check code/                                           # lint
```

The CLI is built with `jsonargparse.CLI` over two functions: `genes` (full library design) and `junctions` (design a set of GG sites without library sequences, used for Fig. 2 of the paper). Any config key can be overridden by `--<key> <value>` on the command line.

`opt_seeds` and `nopt_runs` are mutually exclusive — set exactly one. The default config sets `opt_seeds: null` and uses `nopt_runs`.

## Architecture

Entry point `code/omega.py` orchestrates two flows. The `genes` flow is the main one:

1. **Parse inputs** — `data_classes.define_enzyme` resolves the Type IIS enzyme to an `Enzyme` dataclass; `define_ligation_data` loads a Potapov/Pryor ligation-frequency CSV from `data/ligation_data/` (controlled by `constants_v.LIGATION_DATA`); `PrimerIterator` reads the primer CSV and filters out any primer that contains an enzyme recognition site.
2. **Build a `Library`** (`code/library_classes.py`) from the FASTA records. The library estimates `nfrags` (how many oligos each gene needs to be split into to fit in `oligo_len` after primers + enzyme padding), partitions genes into pools sized so `(nfrags-1) * genes_per_pool + bbsites + other_used_sites <= njunctions`, and assigns one primer pair per pool.
3. **Optimize each pool × seed in parallel** via `joblib.Parallel(n_jobs=njobs)`. Each job constructs a `SAPool` (default, simulated annealing) or `Pool` (greedy), randomizes a starting set of orthogonal GG sites for every gene, then runs `nopt_steps` of optimization. Fidelity is scored by `predict_fidelity.predict_fidelity` against the loaded ligation matrix. Best run per pool is kept.
4. **Package outputs** — `package_library`, `package_oligos`, and a per-pool stats DataFrame are written to `output_dir` as `optimization_results.csv`, `oligo_order.csv`, `pool_stats.csv`, and `experiment_details.txt`.

Key per-module responsibilities:

- `library_classes.py` — `Library` (top-level coordinator), `Pool` (greedy), `SAPool` (simulated annealing; the production path), `Gene` (per-construct fragmentation, padding, primer attachment, and assembly verification). `SAPool.optimize` is the hot loop.
- `predict_fidelity.py` — fidelity scoring. **Performance critical.** The original pandas/Bio.Seq path was the bottleneck; this module precomputes a numpy-indexed view of the ligation DataFrame keyed by `id(data)` in `_TABLE_CACHE` and serves all queries from that cached view. Public signatures are unchanged. If you change the ligation DataFrame in place, the cache will go stale — pass a new DataFrame instance instead.
- `helpers.py` — reverse-complement (`_rc` via `str.translate`), orthogonality check, regex-based illegal-sequence detection (`dna_contains_seq` with `lru_cache`d compiled patterns).
- `junctions.py` — standalone optimizer for raw GG junction sets (no library sequences); the second CLI subcommand.
- `data_classes.py` — `Enzyme`, `LigationData`, `PrimerIterator`, and the `EnzymeTypes` / `LigationDataOpt` enums consumed by jsonargparse.
- `constants_v.py` — registry mapping `LigationDataOpt` values to ligation-data CSV paths and experimental conditions. Add new ligation datasets here.

## Data layout

- `data/ligation_data/pmid_30335370/` — Potapov et al. T4 datasets (default `T4_18h_37C`).
- `data/ligation_data/pmid_32877448/` — Pryor et al. cycling datasets per enzyme.
- `data/fastas/` — example/benchmark FASTAs. `fpbase_*.fasta` are produced by `code/build_fpbase_corpus.py` (pulls FPbase, codon-optimizes for E. coli with `dnachisel` while avoiding BsaI/BsmBI/BbsI sites) and binned by `code/segment_fpbase.py`.
- `data/test_primers.csv` — upstream primer file. **Has known issues**: ships all 20 primers Subramanian et al. flagged as cross-reactive and is missing one validated orthogonal primer (`subra_92`). Preserved as-is for parity with upstream. Prefer `data/subramanian_orthogonal.csv` (165 validated primers, paired fwd/rev) for new designs. See `docs/PRIMER_NOTES.md`.

## Conventions

- Default ligation data is `T4_18h_37C` and default optimizer is `simulated_annealing`. The README strongly recommends BsaI + T4_18h_37C; only deviate with reason.
- `padding` adds random DNA between primers and enzyme sites to equalize oligo lengths. It rejects any draw that would create an enzyme site or a user-listed `illegal_dna_sequences` element across the padding/oligo/primer boundary, and raises after `max_attempts=100000` if the constraints are unsatisfiable. The bundled `fpbase_avgfp_bench.yml` deliberately omits `'ATA'` from `illegal_dna_sequences` because at ~5% per-position match probability it makes some FPbase genes unsatisfiable; `test_install.yml` and `genes_test.yml` keep `'ATA'` for parity with upstream.
- `illegal_dna_sequences` is checked with reverse-complement awareness, so listing `'ATA'` also forbids `'TAT'`.
