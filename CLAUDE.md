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
uv run python ./code/omega.py costs --output_dir output/<run>           # re-price an existing run
uv run pytest                                                     # tests live under tests/
uv run ruff check code/                                           # lint
```

The CLI is built with `jsonargparse.CLI` over three functions: `genes` (full library design), `junctions` (design a set of GG sites without library sequences, used for Fig. 2 of the paper), and `costs` (re-price an existing output dir from `oligo_order.csv` + `pool_stats.csv`). Any config key can be overridden by `--<key> <value>` on the command line.

`opt_seeds` and `nopt_runs` are mutually exclusive — set exactly one. The default config sets `opt_seeds: null` and uses `nopt_runs`.

## Architecture

Entry point `code/omega.py` orchestrates three flows (`genes`, `junctions`, `costs`). The `genes` flow is the main one:

1. **Parse inputs** — `data_classes.define_enzyme` resolves the Type IIS enzyme to an `Enzyme` dataclass; `define_ligation_data` loads a Potapov/Pryor ligation-frequency CSV from `data/ligation_data/` (controlled by `constants_v.LIGATION_DATA`); `PrimerIterator` reads the primer CSV and filters out any primer that contains an enzyme recognition site.
2. **Build a `Library`** (`code/library_classes.py`) from the FASTA records. The library estimates `nfrags` (how many oligos each gene needs to be split into to fit in `oligo_len` after primers + enzyme padding), partitions genes into pools sized so `(nfrags-1) * genes_per_pool + bbsites + other_used_sites <= njunctions`, and assigns one primer pair per pool.
3. **Optimize each pool × seed in parallel** via `joblib.Parallel(n_jobs=njobs)`. Each job constructs a `SAPool` (default, simulated annealing) or `Pool` (greedy), randomizes a starting set of orthogonal GG sites for every gene, then runs `nopt_steps` of optimization. Fidelity is scored by `predict_fidelity.predict_fidelity` against the loaded ligation matrix. Best run per pool is kept.
4. **Package outputs** — `package_library`, `package_oligos`, and a per-pool stats DataFrame are written to `output_dir` as `optimization_results.csv`, `oligo_order.csv`, `pool_stats.csv`, and `experiment_details.txt`.
5. **Cost + wet-lab annotation** (default-on, opt-out via `--pricing_enabled false`) — `code/pricing.py` adds `fwd_primer_cost_usd`/`rev_primer_cost_usd`/`primer_pair_cost_usd` columns to `pool_stats.csv` from the IDT rate card and writes a single-row `cost_summary.csv` with the Twist pool subtotal, IDT primer total, and wet-lab step counts (`wetlab_pcrs`, `wetlab_pcr_cleanups`, `wetlab_quants`, `wetlab_assembly_reactions`, `wetlab_final_cleanups`, `wetlab_transformations`, `wetlab_total_steps`) derived from the README's assembly protocol. `--twist_quote true` also files a live `OLIGO_POOLS_REGULAR` quote against the Twist API.

Key per-module responsibilities:

- `library_classes.py` — `Library` (top-level coordinator), `Pool` (greedy), `SAPool` (simulated annealing; the production path), `Gene` (per-construct fragmentation, padding, primer attachment, and assembly verification). `SAPool.optimize` is the hot loop.
- `predict_fidelity.py` — fidelity scoring. **Performance critical.** The original pandas/Bio.Seq path was the bottleneck; this module precomputes a numpy-indexed view of the ligation DataFrame keyed by `id(data)` in `_TABLE_CACHE` and serves all queries from that cached view. Public signatures are unchanged. If you change the ligation DataFrame in place, the cache will go stale — pass a new DataFrame instance instead.
- `helpers.py` — reverse-complement (`_rc` via `str.translate`), orthogonality check, regex-based illegal-sequence detection (`dna_contains_seq` with `lru_cache`d compiled patterns).
- `junctions.py` — standalone optimizer for raw GG junction sets (no library sequences); the second CLI subcommand.
- `data_classes.py` — `Enzyme`, `LigationData`, `PrimerIterator`, and the `EnzymeTypes` / `LigationDataOpt` enums consumed by jsonargparse.
- `constants_v.py` — registry mapping `LigationDataOpt` values to ligation-data CSV paths and experimental conditions. Add new ligation datasets here.
- `pricing.py` — `TwistOligoPoolPricing` (offline tier-table lookup), `IDTPrimerPricing` (offline per-base lookup), `parse_oligo_pool_quote` (parser anchored to real Twist responses), `live_oligo_pool_quote` (full submit→score→quote flow), and `cost_summary` / `annotate_pool_stats_with_primer_cost` helpers used by the `genes` and `costs` flows. Loaders cache by `Path` in module-level `_TABLE_CACHE`.
- `vendors/twist.py` — vendored Twist TAPI client (see `code/vendors/NOTICE.md`). Used only by the optional `--twist_quote true` live-quote path.

## Data layout

- `data/ligation_data/pmid_30335370/` — Potapov et al. T4 datasets (default `T4_18h_37C`).
- `data/ligation_data/pmid_32877448/` — Pryor et al. cycling datasets per enzyme.
- `data/fastas/` — example/benchmark FASTAs. `fpbase_*.fasta` are produced by `code/build_fpbase_corpus.py` (pulls FPbase, codon-optimizes for E. coli with `dnachisel` while avoiding BsaI/BsmBI/BbsI sites) and binned by `code/segment_fpbase.py`.
- `data/test_primers.csv` — upstream primer file. **Has known issues**: ships all 20 primers Subramanian et al. flagged as cross-reactive and is missing one validated orthogonal primer (`subra_92`). Preserved as-is for parity with upstream. Prefer `data/subramanian_orthogonal.csv` (165 validated primers, paired fwd/rev) for new designs. See `docs/PRIMER_NOTES.md`.
- `data/pricing/` — rate cards consumed by `code/pricing.py`. `twist_oligo_pools.csv` (28 size tiers × 6 length bins, flat per-pool prices) and `idt_primers.csv` (per-base + per-plate-setup keyed on scale/purification/format/length) each ship with a sibling `*.meta.json` recording source, retrieval date, and anchor numbers. The `idt_quote_request_*.xls` workbooks are user-facing IDT bulk-upload artifacts produced by `scripts/generate_idt_quote_request.py`. `example_plate-file-upload.xls` is IDT's reference template. See `docs/PRICING_NOTES.md`.

## Conventions

- Default ligation data is `T4_18h_37C` and default optimizer is `simulated_annealing`. The README strongly recommends BsaI + T4_18h_37C; only deviate with reason.
- `padding` adds random DNA between primers and enzyme sites to equalize oligo lengths. It rejects any draw that would create an enzyme site or a user-listed `illegal_dna_sequences` element across the padding/oligo/primer boundary, and raises after `max_attempts=100000` if the constraints are unsatisfiable. The bundled `fpbase_avgfp_bench.yml` deliberately omits `'ATA'` from `illegal_dna_sequences` because at ~5% per-position match probability it makes some FPbase genes unsatisfiable; `test_install.yml` and `genes_test.yml` keep `'ATA'` for parity with upstream.
- `illegal_dna_sequences` is checked with reverse-complement awareness, so listing `'ATA'` also forbids `'TAT'`.
- `code/vendors/` is **vendored** from chemrich/construct_compiler — see `code/vendors/NOTICE.md` for the source commit. Don't edit in place expecting upstream to know; resync by re-copying.
- `code/pricing.py` writes `cost_summary.csv` to `output_dir` after `genes` runs. Offline lookup against `data/pricing/twist_oligo_pools.csv` is the always-on path (cached by `Path` in `pricing._TABLE_CACHE` like `predict_fidelity._TABLE_CACHE`); `--twist_quote true` adds a live OLIGO_POOLS_REGULAR quote via `vendors.twist.TwistVendor`. See `docs/PRICING_NOTES.md` for env-var requirements and the spot-validated price-table-vs-live-quote check (10/3,512 oligos × 300 nt → Tier 1/Tier 5 prices match exactly).
- For OLIGO_POOLS_REGULAR quotes, the delivery format goes on the **container** (`type: TUBE`, `fill_method: Vertical`), not in `order_settings`. This differs from gene-fragment quotes (which use `order_settings: [{name: "Delivery Format", product_code: "SER_PKG_TUBE"}]`).
- **Filing ad-hoc live Twist quotes** (e.g., to validate a new rate-card column): four credentials are required, all PII-scoped — never hardcode any of them in tracked files. `TWIST_JWT_TOKEN` lives in `~/.zshrc`. `TWIST_END_USER_TOKEN` is exported in the user's interactive shell but not in any rc file (inherited into Bash-tool subshells but not grep-able from rc files). `TWIST_USER_EMAIL` and `TWIST_USER_PHONE` should also be exported in `~/.zshrc`; if either is missing, ask the user and pass it inline (`TWIST_USER_EMAIL=… TWIST_USER_PHONE=… uv run python -c "..."` or `phone=…` kwarg to `live_oligo_pool_quote`). The public entry point is `pricing.live_oligo_pool_quote(oligos, name=..., phone=..., delivery_type='TUBE', fill_method='Vertical', allow_pending_address=...)`, returning a `TwistLiveQuote` (subtotal, shipping, total, business_days, per-pool line items). Each oligo is `(name, sequence)`; pad to the target length with random non-enzyme-site DNA to land in the intended length bin. Quote latency is dominated by Twist's scoring step (~2–5 min); the `score_timeout` (300 s) and `quote_timeout` (1800 s) defaults are usually fine. Validation pattern: anchor each length bin in `data/pricing/twist_oligo_pools.csv` against Tier 1 (10 oligos) + Tier 5 (~3,512 oligos) live quotes, then record the anchor numbers in the sibling `*.meta.json`. Quote responses contain shipping addresses and `owner_contact` — write only into `scripts/_probe_responses/` (gitignored), never into tracked files.
- IDT primer pricing is table-only (no API path — IDT's SciTools Plus doesn't expose pricing). `data/pricing/idt_primers.csv` is keyed on `(scale, purification, format, length range)`; the bundled row anchors to a real IDT plate quote (25nmole/STD/plate = $0.24/bp, $0 setup, validated by 96+69 primer plates totalling $792). `pool_stats.csv` gains `fwd_primer_cost_usd`/`rev_primer_cost_usd`/`primer_pair_cost_usd` columns; `cost_summary.csv` gains `primers_total_usd`, `n_pools`, `primers_per_pool_avg_usd`, and `idt_*` config columns.
