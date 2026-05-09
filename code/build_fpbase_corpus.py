"""Pull all proteins from FPbase, codon-optimize for E. coli, write a FASTA.

Output is a deterministic FASTA suitable as `--input_seqs` for OMEGA. Each
sequence is reverse-translated with dnachisel under constraints that avoid
the Type IIS restriction sites OMEGA may use (BsaI, BsmBI, BbsI), keeps
GC% in a sensible range, and matches E. coli codon usage.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import requests
from dnachisel import (
    AvoidHairpins,
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceTranslation,
    reverse_translate,
)
from joblib import Parallel, delayed
from tqdm import tqdm

FPBASE_URL = "https://www.fpbase.org/api/proteins/?format=json"
FORBIDDEN_SITES = ["GGTCTC", "CGTCTC", "GAAGAC"]  # BsaI, BsmBI, BbsI
AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def fetch_fpbase(cache: Path) -> list[dict]:
    if cache.exists():
        return json.loads(cache.read_text())
    cache.parent.mkdir(parents=True, exist_ok=True)
    print(f"Fetching {FPBASE_URL}", file=sys.stderr)
    r = requests.get(FPBASE_URL, headers={"User-Agent": "omegamega-corpus-builder"}, timeout=60)
    r.raise_for_status()
    cache.write_text(r.text)
    return r.json()


def filter_proteins(records: list[dict], min_len: int, max_len: int) -> list[dict]:
    keep = []
    for p in records:
        seq = (p.get("seq") or "").strip().upper()
        if not seq or not AA_RE.fullmatch(seq):
            continue
        if not (min_len <= len(seq) <= max_len):
            continue
        keep.append({"slug": p["slug"], "name": p["name"], "seq": seq})
    return keep


def codon_optimize(record: dict, organism: str) -> tuple[str, str, str | None]:
    """Return (slug, dna, error). dna is empty if optimization failed."""
    aa = record["seq"]
    seed_dna = reverse_translate(aa)
    constraints = [EnforceTranslation()]
    for site in FORBIDDEN_SITES:
        constraints.append(AvoidPattern(site))
    constraints.append(EnforceGCContent(mini=0.3, maxi=0.7, window=80))
    constraints.append(AvoidHairpins(stem_size=20, hairpin_window=200))
    problem = DnaOptimizationProblem(
        sequence=seed_dna,
        constraints=constraints,
        objectives=[CodonOptimize(species=organism)],
        logger=None,
    )
    try:
        problem.resolve_constraints()
        problem.optimize()
    except Exception as e:
        return record["slug"], "", f"{type(e).__name__}: {e}"
    return record["slug"], problem.sequence, None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--cache",
        type=Path,
        default=Path("data/fpbase_raw.json"),
        help="Path to cache the raw FPbase JSON dump",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("data/fastas/fpbase_codon_optimized.fasta"),
    )
    ap.add_argument("--organism", default="e_coli")
    ap.add_argument("--min-len", type=int, default=150)
    ap.add_argument("--max-len", type=int, default=400)
    ap.add_argument("--limit", type=int, default=0, help="0 = no limit (whole FPbase)")
    ap.add_argument("--njobs", type=int, default=-1)
    args = ap.parse_args()

    raw = fetch_fpbase(args.cache)
    records = filter_proteins(raw, args.min_len, args.max_len)
    if args.limit:
        records = records[: args.limit]
    print(f"FPbase total={len(raw)} usable={len(records)}", file=sys.stderr)

    t0 = time.time()
    results = Parallel(n_jobs=args.njobs, backend="loky")(
        delayed(codon_optimize)(r, args.organism)
        for r in tqdm(records, desc="codon-optimize")
    )
    elapsed = time.time() - t0

    by_slug = {slug: (dna, err) for slug, dna, err in results}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    skipped: list[tuple[str, str]] = []
    with args.out.open("w") as f:
        for r in records:
            dna, err = by_slug[r["slug"]]
            if err or not dna:
                skipped.append((r["slug"], err or "empty"))
                continue
            f.write(f">{r['slug']} {r['name']}\n{dna}\n")
            written += 1

    print(
        f"wrote {written} sequences to {args.out} "
        f"(skipped {len(skipped)}) in {elapsed:.1f}s",
        file=sys.stderr,
    )
    if skipped:
        log = args.out.with_suffix(".skipped.tsv")
        log.write_text("\n".join(f"{s}\t{e}" for s, e in skipped) + "\n")
        print(f"  skip log: {log}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
