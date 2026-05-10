"""Pricing helpers for OMEGA library designs.

Two paths:

* :class:`TwistOligoPoolPricing` — offline tier lookup against a CSV-encoded
  Twist price table (``data/pricing/twist_oligo_pools.csv``). Always available.
* :func:`live_oligo_pool_quote` — submits the designed oligos to Twist's TAPI
  as an ``OLIGO_POOL`` construct, files an ``OLIGO_POOLS_REGULAR`` quote, and
  parses the response. Requires ``TWIST_JWT_TOKEN``, ``TWIST_END_USER_TOKEN``,
  ``TWIST_USER_EMAIL`` env vars and a usable shipping address on the account.

Spot-validated against four prod quotes covering the corners of OMEGA's
design space (smallest pool / largest pool × shortest practical oligo /
longest practical oligo), all matching the offline table exactly:
Tier 1 × 300 nt → $1,030 and Tier 5 × 300 nt → $6,181 on 2026-05-09;
Tier 1 × 350 nt → $1,288 and Tier 5 × 350 nt → $7,727 on 2026-05-10.
"""
from __future__ import annotations

import os
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd

DEFAULT_TWIST_TABLE = (
    Path(__file__).resolve().parent.parent / "data" / "pricing" / "twist_oligo_pools.csv"
)
DEFAULT_IDT_PRIMER_TABLE = (
    Path(__file__).resolve().parent.parent / "data" / "pricing" / "idt_primers.csv"
)


_LENGTH_BINS: list[tuple[str, int, int]] = [
    ("len_20_120", 20, 120),
    ("len_121_150", 121, 150),
    ("len_151_200", 151, 200),
    ("len_201_250", 201, 250),
    ("len_251_300", 251, 300),
    ("len_301_350", 301, 350),
]

_TABLE_CACHE: dict[Path, pd.DataFrame] = {}


@dataclass(frozen=True)
class TwistOfflineQuote:
    """Result of an offline tier lookup."""
    tier: int
    tier_min: int
    tier_max: int
    length_bin: str
    length_min: int
    length_max: int
    n_oligos: int
    max_oligo_len_nt: int
    pool_price_usd: float
    source: str = "list_table"


@dataclass(frozen=True)
class TwistLiveQuote:
    """Parsed result of a live OLIGO_POOLS_REGULAR quote."""
    quote_id: str
    construct_id: str
    pool_subtotal_usd: float
    shipping_usd: float
    handling_usd: float
    tax_total_usd: float
    total_price_usd: float
    business_days: Optional[int]
    pool_lines: list[dict] = field(default_factory=list)


class TwistOligoPoolPricing:
    """Offline tier-based pricing for Twist oligo pools."""

    def __init__(self, table: pd.DataFrame, source_path: Optional[Path] = None):
        self.table = table
        self.source_path = source_path

    @classmethod
    def from_csv(cls, path: str | Path) -> "TwistOligoPoolPricing":
        path = Path(path)
        cached = _TABLE_CACHE.get(path)
        if cached is None:
            cached = pd.read_csv(path)
            _TABLE_CACHE[path] = cached
        return cls(cached, source_path=path)

    def quote(self, n_oligos: int, max_oligo_len_nt: int) -> TwistOfflineQuote:
        if n_oligos < 1:
            raise ValueError(f"n_oligos must be >= 1, got {n_oligos}")
        if max_oligo_len_nt < 1:
            raise ValueError(f"max_oligo_len_nt must be >= 1, got {max_oligo_len_nt}")

        max_supported_len = _LENGTH_BINS[-1][2]
        if max_oligo_len_nt > max_supported_len:
            raise ValueError(
                f"max_oligo_len_nt {max_oligo_len_nt} exceeds Twist's max "
                f"({max_supported_len} nt). Consider --twist-quote for a custom quote."
            )

        bin_key, bin_lo, bin_hi = next(
            b for b in _LENGTH_BINS if max_oligo_len_nt <= b[2]
        )

        rows = self.table[
            (self.table.tier_min <= n_oligos) & (self.table.tier_max >= n_oligos)
        ]
        if rows.empty:
            min_supported = int(self.table.tier_min.min())
            max_supported = int(self.table.tier_max.max())
            raise ValueError(
                f"n_oligos {n_oligos} not in any Twist tier "
                f"(supported range {min_supported}–{max_supported} oligos per pool)."
            )
        row = rows.iloc[0]
        return TwistOfflineQuote(
            tier=int(row.tier),
            tier_min=int(row.tier_min),
            tier_max=int(row.tier_max),
            length_bin=bin_key,
            length_min=bin_lo,
            length_max=bin_hi,
            n_oligos=int(n_oligos),
            max_oligo_len_nt=int(max_oligo_len_nt),
            pool_price_usd=float(row[bin_key]),
        )

    def get_tier_efficiency(self, n_oligos: int) -> dict:
        """Calculate how close we are to the tier boundaries."""
        rows = self.table[
            (self.table.tier_min <= n_oligos) & (self.table.tier_max >= n_oligos)
        ]
        if rows.empty:
            return {}
        row = rows.iloc[0]
        
        tier_max = int(row.tier_max)
        tier_min = int(row.tier_min)
        
        # Check if we are just above a tier boundary
        prev_tier_max = 0
        if int(row.tier) > 1:
            prev_tier = self.table[self.table.tier == int(row.tier) - 1]
            if not prev_tier.empty:
                prev_tier_max = int(prev_tier.iloc[0].tier_max)
        
        return {
            'current_n': n_oligos,
            'tier_min': tier_min,
            'tier_max': tier_max,
            'prev_tier_max': prev_tier_max,
            'over_prev_tier': n_oligos - prev_tier_max if prev_tier_max > 0 else 0,
            'percent_of_tier': (n_oligos - tier_min) / (tier_max - tier_min) if tier_max > tier_min else 1.0
        }


def parse_oligo_pool_quote(response: dict) -> dict:
    """Pull pool/shipping/handling/total out of an OLIGO_POOLS_REGULAR quote.

    Returns a dict with the keyed numeric fields plus ``pool_lines`` (one
    entry per Twist quote line whose description starts with ``Oligo Pool``).
    Anchored to the response shape captured in
    ``scripts/_probe_responses/quote_prod_*.json`` on 2026-05-09.
    """
    inner = response.get("quote") or {}
    lines = inner.get("quote_lines") or []

    def _line_total(line: dict) -> float:
        unit = float(line.get("list_unit_price") or 0.0)
        qty = float(line.get("quantity") or 1.0)
        return unit * qty

    pool_lines_raw = [
        l for l in lines if (l.get("description") or "").startswith("Oligo Pool")
    ]
    shipping = sum(_line_total(l) for l in lines if l.get("product_code") == "SHIP")
    handling = sum(
        _line_total(l) for l in lines if l.get("product_code") == "Handling"
    )
    pool_subtotal = sum(_line_total(l) for l in pool_lines_raw)

    return {
        "pool_subtotal_usd": pool_subtotal,
        "shipping_usd": float(shipping),
        "handling_usd": float(handling),
        "subtotal_usd": float(inner.get("subtotal") or 0.0),
        "tax_total_usd": float(inner.get("tax_total") or 0.0),
        "total_price_usd": float(inner.get("price") or 0.0),
        "business_days": (response.get("tat") or {}).get("business_days"),
        "pool_lines": [
            {
                "product_code": l.get("product_code"),
                "description": l.get("description"),
                "unit_price_usd": float(l.get("list_unit_price") or 0.0),
            }
            for l in pool_lines_raw
        ],
    }


def live_oligo_pool_quote(
    oligos: Iterable[tuple[str, str]],
    *,
    name: Optional[str] = None,
    vendor=None,
    user_email: Optional[str] = None,
    sandbox: bool = False,
    phone: Optional[str] = None,
    recipient_address_id: Optional[str] = None,
    delivery_type: str = "TUBE",
    fill_method: str = "Vertical",
    allow_pending_address: bool = False,
    score_timeout: float = 300.0,
    quote_timeout: float = 1800.0,
) -> TwistLiveQuote:
    """File a real OLIGO_POOLS_REGULAR quote and return parsed pricing.

    See module docstring for env-var requirements. ``vendor`` may be passed
    in to reuse an authenticated TwistVendor; otherwise one is constructed
    from env vars + ``user_email`` (defaults to the value baked into the
    end-user token).
    """
    if vendor is None:
        from vendors.twist import TwistVendor  # local import to keep optional
        vendor = TwistVendor(user_email=user_email or os.environ.get("TWIST_USER_EMAIL"),
                             sandbox=sandbox)
    if not vendor.authenticated:
        raise RuntimeError(
            "Twist credentials missing. Set TWIST_JWT_TOKEN, "
            "TWIST_END_USER_TOKEN, TWIST_USER_EMAIL."
        )

    oligo_list = list(oligos)
    if not oligo_list:
        raise ValueError("No oligos to quote")
    sequences = [seq.upper() for _, seq in oligo_list]

    construct_name = (name or f"omega-pool-{int(time.time())}")[:32]
    submit_url = vendor._user_url("constructs")
    submit_resp = vendor._request("POST", submit_url, json={
        "sequences": sequences,
        "name": construct_name,
        "type": "OLIGO_POOL",
        "adapters_on": False,
    })
    construct_id = submit_resp["id"] if isinstance(submit_resp, dict) else None
    if not construct_id:
        raise RuntimeError(f"Construct creation returned no id: {submit_resp!r}")

    deadline = time.time() + score_timeout
    while time.time() < deadline:
        described = vendor.describe_constructs([construct_id])
        if described and described[0].get("scored"):
            scored = described[0]
            if scored.get("score") != "BUILDABLE":
                raise RuntimeError(
                    f"Construct {construct_id} scored {scored.get('score')!r}; "
                    "cannot quote"
                )
            break
        time.sleep(3.0)
    else:
        raise TimeoutError(f"Construct {construct_id} not scored in {score_timeout:.0f}s")

    user = vendor.get_user()
    first_name = user.get("first_name") or ""
    last_name = user.get("last_name") or ""
    if phone is None:
        phone = user.get("phone_number") or os.environ.get("TWIST_USER_PHONE")
    if not phone:
        raise RuntimeError(
            "No phone number — pass `phone=` or set TWIST_USER_PHONE."
        )

    if recipient_address_id is None:
        addrs = vendor.list_addresses()
        accepted = {"VERIFIED", "APPROVED"}
        if allow_pending_address:
            accepted.add("PENDING_REVIEW")
        ship = [
            a for a in addrs
            if a.get("address_type") == "Shipping"
            and a.get("verification_status") in accepted
        ]
        if not ship:
            raise RuntimeError(
                "No usable shipping address on the Twist account "
                "(pass allow_pending_address=True to accept PENDING_REVIEW)."
            )
        recipient_address_id = next(
            (a["id"] for a in ship if a.get("is_default")), ship[0]["id"]
        )

    quote_resp = vendor.create_quote(
        external_id=f"omega-pool-{int(time.time())}",
        containers=[{
            "constructs": [{"id": construct_id, "index": 1}],
            "type": delivery_type,
            "fill_method": fill_method,
        }],
        order_sub_product_type="OLIGO_POOLS_REGULAR",
        recipient_address_id=recipient_address_id,
        first_name=first_name,
        last_name=last_name,
        phone=phone,
    )
    quote_id = quote_resp.get("id")
    if not quote_id:
        raise RuntimeError(f"create_quote returned no id: {quote_resp!r}")

    deadline = time.time() + quote_timeout
    while time.time() < deadline:
        q = vendor.get_quote(quote_id) or {}
        status = (q.get("status_info") or {}).get("status", "")
        if status == "SUCCESS":
            parsed = parse_oligo_pool_quote(q)
            return TwistLiveQuote(
                quote_id=quote_id,
                construct_id=construct_id,
                pool_subtotal_usd=parsed["pool_subtotal_usd"],
                shipping_usd=parsed["shipping_usd"],
                handling_usd=parsed["handling_usd"],
                tax_total_usd=parsed["tax_total_usd"],
                total_price_usd=parsed["total_price_usd"],
                business_days=parsed["business_days"],
                pool_lines=parsed["pool_lines"],
            )
        if status == "ERROR":
            err = (q.get("status_info") or {}).get("error", "<no detail>")
            raise RuntimeError(f"Quote {quote_id} ended in ERROR: {err}")
        time.sleep(5.0)
    raise TimeoutError(f"Quote {quote_id} did not reach SUCCESS in {quote_timeout:.0f}s")


@dataclass(frozen=True)
class IDTPrimerQuote:
    """Result of an offline IDT primer-pricing lookup."""
    length_nt: int
    scale: str
    purification: str
    format: str
    price_usd: float
    per_base_usd: float
    plate_setup_usd: float
    source: str


class IDTPrimerPricing:
    """Per-primer pricing for IDT, encoded as a tiered per-base + setup table.

    The bundled table at ``data/pricing/idt_primers.csv`` encodes 25 nmole DNA
    Plate Oligo / Standard Desalting at $0.24/bp with no per-plate setup
    fee, anchored to a real IDT cart snapshot for the Subramanian primer set.
    Other scales/purifications can be added as additional rows.
    """

    def __init__(self, table: pd.DataFrame, source_path: Optional[Path] = None):
        self.table = table
        self.source_path = source_path

    @classmethod
    def from_csv(cls, path: str | Path) -> "IDTPrimerPricing":
        path = Path(path)
        cached = _TABLE_CACHE.get(path)
        if cached is None:
            cached = pd.read_csv(path)
            _TABLE_CACHE[path] = cached
        return cls(cached, source_path=path)

    def quote_primer(
        self,
        length_nt: int,
        *,
        scale: str = "25nmole",
        purification: str = "STD",
        format: str = "plate",
    ) -> IDTPrimerQuote:
        if length_nt < 1:
            raise ValueError(f"length_nt must be >= 1, got {length_nt}")
        rows = self.table[
            (self.table.scale == scale)
            & (self.table.purification == purification)
            & (self.table.format == format)
            & (self.table.length_min <= length_nt)
            & (self.table.length_max >= length_nt)
        ]
        if rows.empty:
            raise ValueError(
                f"No IDT pricing row for length={length_nt} scale={scale!r} "
                f"purification={purification!r} format={format!r}. "
                f"Add a row to {self.source_path}."
            )
        row = rows.iloc[0]
        per_base = float(row.price_per_base_usd)
        return IDTPrimerQuote(
            length_nt=int(length_nt),
            scale=scale,
            purification=purification,
            format=format,
            price_usd=per_base * length_nt,
            per_base_usd=per_base,
            plate_setup_usd=float(row.plate_setup_usd),
            source=str(row.source),
        )

    def quote_pool_pair(
        self,
        fwd_len: int,
        rev_len: int,
        **kwargs,
    ) -> dict:
        """Cost of one (forward, reverse) primer pair for a single subpool."""
        f = self.quote_primer(fwd_len, **kwargs)
        r = self.quote_primer(rev_len, **kwargs)
        return {
            "fwd_primer_cost_usd": f.price_usd,
            "rev_primer_cost_usd": r.price_usd,
            "pair_cost_usd": f.price_usd + r.price_usd,
            "per_base_usd": f.per_base_usd,
            "scale": f.scale,
            "purification": f.purification,
            "format": f.format,
            "source": f.source,
        }


# NEB stock concentrations (U/µL) used to convert units → pipette volumes.
_ENZYME_CONC_U_PER_UL: dict[str, float] = {
    "BsaI":  20.0,  # BsaI-HFv2, R3733
    "BsmBI": 10.0,  # BsmBI-v2, R0580
    "BbsI":  10.0,  # BbsI-HF, R3539
}
_ENZYME_NEB_CAT: dict[str, str] = {
    "BsaI":  "NEB R3733",
    "BsmBI": "NEB R0580",
    "BbsI":  "NEB R3539",
}
_T4_LIGASE_CONC_U_PER_UL = 400.0   # NEB M0202T high-concentration stock
_KAPA_POLY_CONC_U_PER_UL  = 1.0    # KAPA HiFi HotStart Polymerase
_PCR_RXN_VOL_UL            = 25.0   # Table 1 reaction volume
_GG_ENZYME_U_PER_RXN       = 15.0   # Table 3
_GG_T4_LIGASE_U_PER_RXN    = 1000.0 # Table 3 (added post-digest)
_KAPA_POLY_U_PER_RXN       = 0.5    # Table 1
_KAPA_BUFFER_UL_PER_RXN    = 5.0    # Table 1 (5× buffer)
_SPRI_RATIO                = 1.0    # Table 4 default for 350 nt amplicons


def reagent_quantities(n_pools: int, enzyme_name: str = "BsaI") -> pd.DataFrame:
    """Reagent quantities to physically instantiate a designed library.

    Covers the four consumables with meaningful volume per run: KAPA HiFi
    (PCR), the Type IIS enzyme (GG digest), T4 DNA Ligase (GG ligation), and
    SPRI beads (PCR cleanup). Quantities follow Tables 1–4 of the README.

    Args:
        n_pools: Number of subpools (PCR reactions = GG reactions = SPRI
            cleanups).
        enzyme_name: Type IIS enzyme used; one of BsaI, BsmBI, BbsI.

    Returns:
        DataFrame with columns reagent, quantity_per_pool, unit, n_pools,
        total_quantity, notes.
    """
    if n_pools < 1:
        raise ValueError(f"n_pools must be >= 1, got {n_pools}")
    if enzyme_name not in _ENZYME_CONC_U_PER_UL:
        raise ValueError(
            f"Unknown enzyme {enzyme_name!r}. "
            f"Supported: {sorted(_ENZYME_CONC_U_PER_UL)}"
        )

    enz_conc = _ENZYME_CONC_U_PER_UL[enzyme_name]
    enz_cat  = _ENZYME_NEB_CAT[enzyme_name]
    enz_ul   = _GG_ENZYME_U_PER_RXN / enz_conc

    t4_ul    = _GG_T4_LIGASE_U_PER_RXN / _T4_LIGASE_CONC_U_PER_UL
    kapa_ul  = _KAPA_POLY_U_PER_RXN / _KAPA_POLY_CONC_U_PER_UL
    spri_ul  = _SPRI_RATIO * _PCR_RXN_VOL_UL

    rows = [
        {
            "reagent":           "KAPA HiFi HotStart Polymerase",
            "quantity_per_pool": _KAPA_POLY_U_PER_RXN,
            "unit":              "U",
            "n_pools":           n_pools,
            "total_quantity":    _KAPA_POLY_U_PER_RXN * n_pools,
            "notes":             f"{kapa_ul:.2g} µL/pool at {_KAPA_POLY_CONC_U_PER_UL:.0f} U/µL (KAPA HiFi HotStart PCR Kit)",
        },
        {
            "reagent":           "KAPA HiFi Buffer (5×)",
            "quantity_per_pool": _KAPA_BUFFER_UL_PER_RXN,
            "unit":              "µL",
            "n_pools":           n_pools,
            "total_quantity":    _KAPA_BUFFER_UL_PER_RXN * n_pools,
            "notes":             "PCR step (Table 1); bundled with KAPA HiFi HotStart PCR Kit",
        },
        {
            "reagent":           enzyme_name,
            "quantity_per_pool": _GG_ENZYME_U_PER_RXN,
            "unit":              "U",
            "n_pools":           n_pools,
            "total_quantity":    _GG_ENZYME_U_PER_RXN * n_pools,
            "notes":             f"{enz_ul:.2g} µL/pool at {enz_conc:.0f} U/µL ({enz_cat}); GG assembly digest",
        },
        {
            "reagent":           "T4 DNA Ligase",
            "quantity_per_pool": _GG_T4_LIGASE_U_PER_RXN,
            "unit":              "U",
            "n_pools":           n_pools,
            "total_quantity":    _GG_T4_LIGASE_U_PER_RXN * n_pools,
            "notes":             f"{t4_ul:.2g} µL/pool at {_T4_LIGASE_CONC_U_PER_UL:.0f} U/µL (NEB M0202T); GG assembly ligation",
        },
        {
            "reagent":           "SPRI beads",
            "quantity_per_pool": spri_ul,
            "unit":              "µL",
            "n_pools":           n_pools,
            "total_quantity":    spri_ul * n_pools,
            "notes":             f"{_SPRI_RATIO:.1f}× of {_PCR_RXN_VOL_UL:.0f} µL PCR; AMPure XP / Mag-Bind / equivalent",
        },
    ]
    return pd.DataFrame(rows)


def write_reagent_summary(output_dir: str | Path, df: pd.DataFrame) -> Path:
    """Write the DataFrame from :func:`reagent_quantities` to ``reagent_summary.csv``."""
    out = Path(output_dir) / "reagent_summary.csv"
    df.to_csv(out, index=False)
    return out


def wet_lab_steps(n_pools: int) -> dict:
    """Count of wet-lab manipulations to physically instantiate a designed
    library, following the README's assembly protocol.

    Per subpool: 1 PCR amplification, 1 SPRI cleanup, 1 PicoGreen quant
    (with normalization), 1 Golden Gate assembly. Then 1 final
    pool-and-cleanup, and 1 transformation of the combined library.

    Optional/downstream steps (plating, colony picking, sequencing,
    re-PCR + re-digest after transformation) aren't counted because they
    depend on what you do with the library.
    """
    if n_pools < 1:
        raise ValueError(f"n_pools must be >= 1, got {n_pools}")
    return {
        "wetlab_pcrs": n_pools,
        "wetlab_pcr_cleanups": n_pools,
        "wetlab_quants": n_pools,
        "wetlab_assembly_reactions": n_pools,
        "wetlab_final_cleanups": 1,
        "wetlab_transformations": 1,
        "wetlab_total_steps": 4 * n_pools + 2,
    }


def cost_summary(
    oligo_order_df: pd.DataFrame,
    *,
    table_path: Path = DEFAULT_TWIST_TABLE,
    twist_quote: bool = False,
    twist_kwargs: Optional[dict] = None,
    pool_stats_df: Optional[pd.DataFrame] = None,
    idt_table_path: Path = DEFAULT_IDT_PRIMER_TABLE,
    idt_scale: str = "25nmole",
    idt_purification: str = "STD",
    idt_format: str = "plate",
) -> dict:
    """Build a single-row cost summary for an OMEGA oligo_order DataFrame.

    Always includes Twist offline tier-table pricing. With ``twist_quote=True``
    also submits a live OLIGO_POOLS_REGULAR quote. With ``pool_stats_df``
    provided, also computes IDT primer costs from the pool's primer pairs.
    """
    n_oligos = int(len(oligo_order_df))
    max_len = int(oligo_order_df.sequence.str.len().max())

    pricing = TwistOligoPoolPricing.from_csv(table_path)
    offline = pricing.quote(n_oligos, max_len)
    summary: dict[str, Any] = {
        "n_oligos": n_oligos,
        "max_oligo_len_nt": max_len,
        "offline_tier": offline.tier,
        "offline_tier_min": offline.tier_min,
        "offline_tier_max": offline.tier_max,
        "offline_length_bin": offline.length_bin,
        "offline_pool_price_usd": offline.pool_price_usd,
        "offline_source": offline.source,
    }

    if twist_quote:
        oligos = list(zip(oligo_order_df.name, oligo_order_df.sequence))
        live = live_oligo_pool_quote(oligos, **(twist_kwargs or {}))
        summary.update({
            "live_construct_id": live.construct_id,
            "live_quote_id": live.quote_id,
            "live_pool_subtotal_usd": live.pool_subtotal_usd,
            "live_shipping_usd": live.shipping_usd,
            "live_handling_usd": live.handling_usd,
            "live_tax_total_usd": live.tax_total_usd,
            "live_total_price_usd": live.total_price_usd,
            "live_business_days": live.business_days,
            "live_pool_lines": "; ".join(
                f"{p['product_code']}={p['unit_price_usd']:.2f} ({p['description']})"
                for p in live.pool_lines
            ),
        })

    if pool_stats_df is not None and len(pool_stats_df) > 0:
        n_pools = int(len(pool_stats_df))
        idt = IDTPrimerPricing.from_csv(idt_table_path)
        per_pool: list[float] = []
        for _, row in pool_stats_df.iterrows():
            pair = idt.quote_pool_pair(
                len(row["pfwd_sequence"]),
                len(row["prev_sequence"]),
                scale=idt_scale,
                purification=idt_purification,
                format=idt_format,
            )
            per_pool.append(pair["pair_cost_usd"])
        summary.update({
            "n_pools": n_pools,
            "n_primer_pairs": n_pools,
            "primers_per_pool_avg_usd": sum(per_pool) / len(per_pool),
            "primers_total_usd": float(sum(per_pool)),
            "idt_scale": idt_scale,
            "idt_purification": idt_purification,
            "idt_format": idt_format,
        })
        summary.update(wet_lab_steps(n_pools))

    return summary


def write_cost_summary(output_dir: str | Path, summary: dict) -> Path:
    """Write the dict from :func:`cost_summary` to ``cost_summary.csv``."""
    out = Path(output_dir) / "cost_summary.csv"
    pd.DataFrame([summary]).to_csv(out, index=False)
    return out


def annotate_pool_stats_with_primer_cost(
    pool_stats_df: pd.DataFrame,
    *,
    idt_table_path: Path = DEFAULT_IDT_PRIMER_TABLE,
    scale: str = "25nmole",
    purification: str = "STD",
    format: str = "plate",
) -> pd.DataFrame:
    """Return ``pool_stats_df`` with IDT per-pool primer cost columns added."""
    idt = IDTPrimerPricing.from_csv(idt_table_path)
    out = pool_stats_df.copy()
    fwd, rev, pair = [], [], []
    for _, row in out.iterrows():
        q = idt.quote_pool_pair(
            len(row["pfwd_sequence"]),
            len(row["prev_sequence"]),
            scale=scale, purification=purification, format=format,
        )
        fwd.append(q["fwd_primer_cost_usd"])
        rev.append(q["rev_primer_cost_usd"])
        pair.append(q["pair_cost_usd"])
    out["fwd_primer_cost_usd"] = fwd
    out["rev_primer_cost_usd"] = rev
    out["primer_pair_cost_usd"] = pair
    return out
