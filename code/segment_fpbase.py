"""Segment the full FPbase corpus into size bins (FASTA per bin).

Each bin is written to data/fastas/fpbase_<label>.fasta. Bins are defined
relative to a reference protein (default avGFP, 238 aa) so the corpus can
be benchmarked on a size-controlled subset.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from Bio import SeqIO


# (label, lower_aa, upper_aa) — inclusive bounds. AvGFP = 238 aa.
DEFAULT_BINS: list[tuple[str, int, int]] = [
    ("avgfp_size", 190, 286),  # 238 ± 20%
    ("small", 150, 189),
    ("large", 287, 400),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input",
        type=Path,
        default=Path("data/fastas/fpbase_codon_optimized.fasta"),
    )
    ap.add_argument("--out-dir", type=Path, default=Path("data/fastas"))
    ap.add_argument(
        "--prefix",
        default="fpbase_",
        help="Output filenames are <prefix><label>.fasta",
    )
    args = ap.parse_args()

    records = list(SeqIO.parse(args.input, "fasta"))
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"input: {args.input} ({len(records)} records)", file=sys.stderr)
    for label, lo_aa, hi_aa in DEFAULT_BINS:
        lo_nt, hi_nt = lo_aa * 3, hi_aa * 3  # codons w/o stop
        kept = [r for r in records if lo_nt <= len(r.seq) <= hi_nt]
        out = args.out_dir / f"{args.prefix}{label}.fasta"
        with out.open("w") as f:
            for r in kept:
                f.write(f">{r.id} {r.description.split(' ', 1)[1] if ' ' in r.description else ''}\n{str(r.seq)}\n")
        print(
            f"  bin {label!r}: {lo_aa}-{hi_aa} aa ({lo_nt}-{hi_nt} bp) -> "
            f"{len(kept)} records -> {out}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
