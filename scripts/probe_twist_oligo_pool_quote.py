"""Probe the Twist TAPI for an OLIGO_POOLS_REGULAR quote.

Stages:
  1. auth     — call get_user() and confirm name/email/phone come back
  2. address  — call list_addresses() and pick a verified shipping addr
  3. submit   — submit oligos as an OLIGO_POOL construct and poll scoring
  4. quote    — file create_quote with order_sub_product_type=OLIGO_POOLS_REGULAR
                and poll until SUCCESS, dumping the full response to disk

The script is intended to be run interactively with --stage. Default sandbox.

Auth via env:
  TWIST_JWT_TOKEN, TWIST_END_USER_TOKEN, TWIST_USER_EMAIL (must match the
  end-user-token payload; required, no default).

Usage:
  uv run python scripts/probe_twist_oligo_pool_quote.py --stage auth
  uv run python scripts/probe_twist_oligo_pool_quote.py --stage address
  uv run python scripts/probe_twist_oligo_pool_quote.py --stage submit \
      --oligos output/fpbase_avgfp_bench/oligo_order.csv --limit 5
  uv run python scripts/probe_twist_oligo_pool_quote.py --stage quote \
      --construct-id <uuid>
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))

from vendors.twist import TwistVendor  # noqa: E402

RESPONSES_DIR = ROOT / "scripts" / "_probe_responses"


def _now() -> str:
    return time.strftime("%Y%m%dT%H%M%S")


def _vendor(sandbox: bool) -> TwistVendor:
    email = os.environ.get("TWIST_USER_EMAIL")
    if not email:
        raise SystemExit("TWIST_USER_EMAIL must be set (Twist account email)")
    return TwistVendor(user_email=email, sandbox=sandbox)


def stage_auth(args) -> int:
    v = _vendor(args.sandbox)
    if not v.authenticated:
        print("Missing TWIST_JWT_TOKEN or TWIST_END_USER_TOKEN", file=sys.stderr)
        return 1
    user = v.get_user()
    out = {
        "id": user.get("id"),
        "first_name": user.get("first_name"),
        "last_name": user.get("last_name"),
        "email": user.get("email"),
        "phone_number": user.get("phone_number"),
        "default_shipping_address": user.get("default_shipping_address"),
        "default_billing_address": user.get("default_billing_address"),
    }
    print(json.dumps(out, indent=2))
    return 0


def stage_address(args) -> int:
    v = _vendor(args.sandbox)
    addrs = v.list_addresses()
    print(f"Returned {len(addrs)} addresses")
    summary = []
    for a in addrs:
        summary.append({
            "id": a.get("id"),
            "address_type": a.get("address_type"),
            "verification_status": a.get("verification_status"),
            "is_default": a.get("is_default"),
            "city": a.get("city"),
            "state": a.get("state"),
        })
    print(json.dumps(summary, indent=2))
    return 0


def _read_oligos(path: Path, limit: int | None) -> list[tuple[str, str]]:
    """Read (name, sequence) tuples from oligo_order.csv. CSV has columns
    `,name,sequence` (an unnamed index column from pandas)."""
    import csv
    rows: list[tuple[str, str]] = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((row["name"], row["sequence"]))
            if limit and len(rows) >= limit:
                break
    return rows


def stage_submit(args) -> int:
    v = _vendor(args.sandbox)
    oligos = _read_oligos(Path(args.oligos), args.limit)
    print(f"Submitting {len(oligos)} oligos as OLIGO_POOL construct...")
    if not oligos:
        print("No oligos to submit", file=sys.stderr)
        return 1

    sequences = [seq.upper() for _, seq in oligos]
    name = args.name or f"omega-probe-{_now()}"
    body: dict = {
        "sequences": sequences,
        "name": name[:32],
        "type": "OLIGO_POOL",
        "adapters_on": False,
    }
    url = v._user_url("constructs")
    resp = v._request("POST", url, json=body)
    out_path = RESPONSES_DIR / f"submit_{'sandbox' if args.sandbox else 'prod'}_{_now()}.json"
    out_path.write_text(json.dumps(resp, indent=2) + "\n")
    print(f"Saved POST /constructs/ response to {out_path}")
    print(json.dumps(resp, indent=2)[:400])
    return 0


def stage_describe(args) -> int:
    v = _vendor(args.sandbox)
    result = v.describe_constructs([args.construct_id])
    out_path = RESPONSES_DIR / f"describe_{'sandbox' if args.sandbox else 'prod'}_{_now()}.json"
    out_path.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Saved describe response to {out_path}")
    print(json.dumps(result, indent=2)[:1200])
    return 0


def stage_quote(args) -> int:
    v = _vendor(args.sandbox)
    user = v.get_user()
    first_name = user.get("first_name") or ""
    last_name = user.get("last_name") or ""
    phone = user.get("phone_number") or os.environ.get("TWIST_USER_PHONE", "")
    if not phone:
        print("No phone on profile or TWIST_USER_PHONE env var", file=sys.stderr)
        return 1
    addrs = v.list_addresses()
    accepted = {"VERIFIED", "APPROVED"}
    if args.allow_pending:
        accepted.add("PENDING_REVIEW")
    candidates = [
        a for a in addrs
        if a.get("address_type") == "Shipping"
        and a.get("verification_status") in accepted
    ]
    if not candidates:
        print("No usable shipping address", file=sys.stderr)
        return 1
    addr = next((a for a in candidates if a.get("is_default")), candidates[0])
    address_id = addr["id"]
    print(f"Using shipping address {address_id} "
          f"({addr.get('city')}, {addr.get('state')}) [{addr.get('verification_status')}]")

    # Oligo-pool containers carry the delivery format on the container itself
    # (per OpenAPI example), not in `order_settings`.
    containers = [{
        "constructs": [{"id": args.construct_id, "index": 1}],
        "type": "TUBE",
        "fill_method": "Vertical",
    }]
    external_id = f"omega-pool-probe-{int(time.time())}"
    print(f"Filing quote external_id={external_id} for construct {args.construct_id}...")
    quote = v.create_quote(
        external_id=external_id,
        containers=containers,
        order_sub_product_type="OLIGO_POOLS_REGULAR",
        recipient_address_id=address_id,
        first_name=first_name,
        last_name=last_name,
        phone=phone,
    )
    quote_id = quote.get("id")
    if not quote_id:
        print(f"create_quote returned no id: {quote}", file=sys.stderr)
        return 1
    print(f"Quote {quote_id}; polling until SUCCESS...")

    deadline = time.time() + 1800.0
    last_status = ""
    while time.time() < deadline:
        q = v.get_quote(quote_id) or {}
        status = (q.get("status_info") or {}).get("status", "")
        if status != last_status:
            print(f"  status: {status}", flush=True)
            last_status = status
        if status == "SUCCESS":
            out_path = RESPONSES_DIR / f"quote_{'sandbox' if args.sandbox else 'prod'}_{_now()}.json"
            out_path.write_text(json.dumps(q, indent=2) + "\n")
            print(f"Saved quote response to {out_path}")
            return 0
        if status == "ERROR":
            err = (q.get("status_info") or {}).get("error", "<no detail>")
            print(f"Quote ended in ERROR: {err}", file=sys.stderr)
            return 1
        time.sleep(5.0)
    print("Timeout waiting for quote", file=sys.stderr)
    return 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True,
                    choices=["auth", "address", "submit", "describe", "quote"])
    sb = ap.add_mutually_exclusive_group()
    sb.add_argument("--sandbox", action="store_true", default=True,
                    help="Use Twist sandbox (default)")
    sb.add_argument("--prod", action="store_false", dest="sandbox",
                    help="Use Twist production")
    ap.add_argument("--oligos", type=Path,
                    default=ROOT / "output" / "fpbase_avgfp_bench" / "oligo_order.csv")
    ap.add_argument("--limit", type=int, default=None,
                    help="Submit only the first N oligos (probe-friendly)")
    ap.add_argument("--name", default=None)
    ap.add_argument("--construct-id", default=None,
                    help="UUID returned by --stage submit")
    ap.add_argument("--allow-pending", action="store_true",
                    help="Accept PENDING_REVIEW shipping addresses")
    ap.add_argument("--delivery-format", default="SER_PKG_TUBE",
                    help="Twist delivery-format product_code (default SER_PKG_TUBE)")
    args = ap.parse_args()

    if args.stage == "auth":
        return stage_auth(args)
    if args.stage == "address":
        return stage_address(args)
    if args.stage == "submit":
        return stage_submit(args)
    if args.stage == "describe":
        if not args.construct_id:
            print("--construct-id required", file=sys.stderr)
            return 2
        return stage_describe(args)
    if args.stage == "quote":
        if not args.construct_id:
            print("--construct-id required", file=sys.stderr)
            return 2
        return stage_quote(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
