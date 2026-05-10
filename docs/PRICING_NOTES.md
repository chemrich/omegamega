# Pricing notes

OMEGA optionally estimates the cost of ordering a designed library. Two paths
are available; both are wired into the `genes` flow and the standalone
`costs` subcommand.

## Twist oligo pools — offline tier table

`data/pricing/twist_oligo_pools.csv` encodes Twist's 28-tier × 6-length-bin
flat-rate table for oligo pools. Pricing is keyed by `(n_oligos, max_oligo_len_nt)`:

* Tier rows: `tier_min`–`tier_max` oligo counts (e.g., Tier 5 = 2,001–6,000).
* Length bins: 20–120, 121–150, 151–200, 201–250, 251–300, 301–350 nt.
* Cell value: total pool list price in USD.

This table was supplied by the user from their Twist eCommerce account portal
on 2026-05-09 (see `data/pricing/twist_oligo_pools.meta.json`). Account-specific
quotes may differ; treat the table as a published list-price estimate.

## Twist oligo pools — live API quote

Pass `--twist_quote true` to `genes` or `costs` to file a real
`OLIGO_POOLS_REGULAR` quote against your Twist account and surface the parsed
numbers (subtotal, shipping, handling, tax, total, business-day TAT).

Required environment variables (read by `code/vendors/twist.py` — vendored
from chemrich/construct_compiler):

```
TWIST_JWT_TOKEN          # Authorization header (JWT prefix added)
TWIST_END_USER_TOKEN     # X-End-User-Token header
TWIST_USER_EMAIL         # scopes /v1/users/{email}/ paths
TWIST_USER_PHONE         # required by create_quote if your profile lacks one
```

Account requirements:

* A shipping address (Twist accepts `PENDING_REVIEW` for quoting; OMEGA's
  helper has `allow_pending_address` but the `--twist_quote` CLI flag does
  not yet expose it — call `pricing.live_oligo_pool_quote(...)` directly to
  override).
* A name and phone on either the profile or via env var.
* Twist whitelists the requesting IP — contact `b2b-support@twistbioscience.com`
  to register if a fresh request returns 403/500.

Offline-table prices have been spot-validated against the live Twist API
at the corners that bracket OMEGA's design space (smallest pool / largest
pool × shortest practical oligo / longest practical oligo):

| pool | tier × length bin | offline price | live `subtotal` | match | validated |
|---|---|---|---|---|---|
| 10 oligos × 300 nt    | Tier 1 × len_251_300 | $1,030.00 | $1,030.00 | ✅ | 2026-05-09 |
| 3,512 oligos × 300 nt | Tier 5 × len_251_300 | $6,181.00 | $6,181.00 | ✅ | 2026-05-09 |
| 10 oligos × 350 nt    | Tier 1 × len_301_350 | $1,288.00 | $1,288.00 | ✅ | 2026-05-10 |
| 3,512 oligos × 350 nt | Tier 5 × len_301_350 | $7,727.00 | $7,727.00 | ✅ | 2026-05-10 |

Live `total_price_usd` adds shipping ($35), handling ($25), and tax (which
varies by destination state).

### Choosing oligo length: the cost trade-off isn't flat

The per-cell flat rate is higher at longer lengths, but the same library
needs fewer oligos and fewer subpools when oligos are longer (each oligo
covers more gene content; `nfrags` decreases). Net library cost depends on
which side of a tier boundary you land on, and on whether your IDT/labor
savings outweigh the Twist increase. Concrete worked example from the
bundled `configs/fpbase_avgfp_bench.yml` (878 FPbase genes, 642–804 bp):

| metric | `oligo_len: 300` | `oligo_len: 350` | Δ |
|---|---|---|---|
| subpools | 55 | 37 | −33% |
| total oligos | 3,512 | 2,634 | −25% |
| Twist tier × bin | Tier 5 × len_251_300 | Tier 5 × len_301_350 | same tier |
| Twist cost | $6,181.00 | $7,727.00 | +$1,546 |
| IDT primer pairs | $528.00 | $355.20 | −$172.80 |
| **Reagent total** | **$6,709.00** | **$8,082.20** | **+$1,373** |
| wet-lab steps (4·N + 2) | 222 | 150 | −72 steps |

For this library, both lengths land in Tier 5 (the 2001–6000-oligo
bracket), so going longer just pays the higher per-cell rate without
escaping the tier. A library with a smaller oligo count near a tier
boundary could move *down* a tier at 350 nt and pay *less* Twist
overall — the right call is empirical, not formulaic.

## Output

When pricing runs, OMEGA writes `cost_summary.csv` next to `oligo_order.csv`.
Offline-only runs have these columns:

```
n_oligos,max_oligo_len_nt,
offline_tier,offline_tier_min,offline_tier_max,
offline_length_bin,offline_pool_price_usd,offline_source
```

If `pool_stats.csv` is available (always written by `genes`; loaded
opportunistically by `costs`), the IDT primer columns and wet-lab step
counts are added too:

```
n_pools,n_primer_pairs,
primers_per_pool_avg_usd,primers_total_usd,
idt_scale,idt_purification,idt_format,
wetlab_pcrs,wetlab_pcr_cleanups,wetlab_quants,wetlab_assembly_reactions,
wetlab_final_cleanups,wetlab_transformations,wetlab_total_steps
```

`--twist_quote true` adds:

```
live_construct_id,live_quote_id,
live_pool_subtotal_usd,live_shipping_usd,live_handling_usd,
live_tax_total_usd,live_total_price_usd,live_business_days,
live_pool_lines
```

## Repricing past runs

```
uv run python ./code/omega.py costs --output_dir output/<existing_run>
uv run python ./code/omega.py costs --output_dir output/<existing_run> --twist_quote true
```

## Updating the Twist table

Re-export the price grid from your Twist eCommerce account, replace
`data/pricing/twist_oligo_pools.csv`, bump `retrieved` in the sibling
`.meta.json`. The cache in `pricing._TABLE_CACHE` is keyed by `Path`, so a
new path or a fresh process picks up changes automatically.

## IDT primer pricing (per-pool primer pairs)

`data/pricing/idt_primers.csv` encodes a per-base + per-plate-setup table
keyed on `(scale, purification, format, length range)`. The bundled row
was anchored to a real IDT cart snapshot for the 165 Subramanian
orthogonal primers (`GetCartSnapShot 36128990`, retrieved 2026-05-10):

| config | per-base | plate setup | example: 20 nt | example: 1 pair |
|---|---|---|---|---|
| 25nmole / STD / plate | $0.24 | $0.00 | $4.80 | $9.60 |

Plate 1 = 96 primers × $4.80 = $460.80; plate 2 = 69 × $4.80 = $331.20;
subtotal $792.00 (the rate per primer doesn't change between plates).

IDT's API does **not** expose primer pricing — SciTools Plus only covers
codon optimization, oligo analysis, and gBlock screening — so this is a
table-only path (no `--idt_quote` live mode). Add new rows to the CSV for
other scales (`100nmole`, `250nmole`), purifications (`HPLC`, `PAGE`),
formats (`tube`), or length ranges as you collect quotes.

When the bundled rate card is present, OMEGA does two things automatically:

* `pool_stats.csv` gains `fwd_primer_cost_usd`, `rev_primer_cost_usd`,
  `primer_pair_cost_usd` columns (one row per subpool).
* `cost_summary.csv` gains `n_pools`, `n_primer_pairs`,
  `primers_per_pool_avg_usd`, `primers_total_usd`, `idt_scale`,
  `idt_purification`, `idt_format`.

For libraries that reuse a primer set you've already ordered, treat
`primers_total_usd` as a sunk cost rather than recurring.

## Generating an IDT bulk-quote upload

`scripts/generate_idt_quote_request.py` flattens an OMEGA paired-primer
CSV (`fwd_name,fwd_sequence,rev_name,rev_sequence`) into IDT's plate-upload
.xls format. Each output workbook holds up to 96 primers in row-major
order (A1, A2, ..., A12, B1, ...). Use `--wells-per-plate 384` for
384-well plates.

```
uv run python scripts/generate_idt_quote_request.py \
    --input data/subramanian_orthogonal.csv \
    --output-prefix data/pricing/idt_quote_request_subramanian_orthogonal
```

At the IDT bulk-upload page, set Scale = 25 nm and Purification = Standard
Desalt — those aren't template fields. The IDT-supplied template lives at
`data/pricing/example_plate-file-upload.xls` for reference.

## Wet-lab step counts

Per the README's assembly protocol, instantiating an N-pool library takes:

| step | per-library count |
|---|---|
| PCR amplifications  | N (one per subpool) |
| PCR cleanups        | N |
| Golden Gate assembly reactions | N |
| Final pool-and-cleanup | 1 |
| Transformations | 1 |
| **Total**          | **3·N + 2** |

These are emitted as `wetlab_*` columns in `cost_summary.csv` and printed
to stdout. Optional/downstream steps (plating, colony picking, sequencing,
re-PCR + re-digest after transformation per the README's note) aren't
counted because they depend on what you do with the library. The
"transformation = 1" assumption follows the README's combined-pool
protocol; if you instead transform each subpool separately to characterize
per-pool fidelity before pooling, multiply by N.

Worked examples:

| library | N pools | total wet-lab steps |
|---|---|---|
| `test_install` (60 oligos) | 1  | 5   |
| `genes_test` (7 subpools)  | 7  | 23  |
| `fpbase_avgfp_bench` (3,512 oligos) | 55 | 167 |
