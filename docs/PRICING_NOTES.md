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

The live quote was spot-validated on 2026-05-09 against two production quotes:

| pool | tier (offline) | offline price | live `subtotal` | match |
|---|---|---|---|---|
| 10 oligos × 300 nt   | Tier 1 | $1,030.00 | $1,030.00 | ✅ |
| 3,512 oligos × 300 nt | Tier 5 | $6,181.00 | $6,181.00 | ✅ |

Live `total_price_usd` adds shipping ($35), handling ($25), and tax (which
varies by destination state).

## Output

When pricing runs, OMEGA writes `cost_summary.csv` next to `oligo_order.csv`.
Offline-only runs have these columns:

```
n_oligos,max_oligo_len_nt,
offline_tier,offline_tier_min,offline_tier_max,
offline_length_bin,offline_pool_price_usd,offline_source
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

## Updating the table

Re-export the price grid from your Twist eCommerce account, replace
`data/pricing/twist_oligo_pools.csv`, bump `retrieved` in the sibling
`.meta.json`. The cache in `pricing._TABLE_CACHE` is keyed by `Path`, so a
new path or a fresh process picks up changes automatically.
