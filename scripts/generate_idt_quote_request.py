"""Flatten an OMEGA paired-primer CSV into IDT's bulk-upload format.

OMEGA primer CSVs use the paired layout
``fwd_name,fwd_sequence,rev_name,rev_sequence`` (one row per primer pair).
IDT's bulk-order tool expects one primer per row with a Name + Sequence
plus optional Scale and Purification columns.

Usage:
    uv run python scripts/generate_idt_quote_request.py \\
        --input data/subramanian_orthogonal.csv \\
        --output data/pricing/idt_quote_request_subramanian_orthogonal.csv \\
        --scale '25 nm' --purification STD
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def flatten_primer_csv(input_path: Path) -> list[tuple[str, str]]:
    """Yield (name, sequence) for every non-empty primer in a paired CSV.

    Skips rows where a fwd or rev cell is empty (the orthogonal set has one
    row with an empty fwd, since it ships 82 fwd + 83 rev).
    """
    primers: list[tuple[str, str]] = []
    seen: set[str] = set()
    with input_path.open() as f:
        reader = csv.DictReader(f)
        required = {"fwd_name", "fwd_sequence", "rev_name", "rev_sequence"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(
                f"{input_path} missing required columns; got {reader.fieldnames}"
            )
        for row in reader:
            for name_col, seq_col in [("fwd_name", "fwd_sequence"),
                                      ("rev_name", "rev_sequence")]:
                name, seq = row[name_col].strip(), row[seq_col].strip()
                if not name or not seq:
                    continue
                if name in seen:
                    continue
                seen.add(name)
                primers.append((name, seq.upper()))
    return primers


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, required=True,
                    help="Paired-primer CSV (fwd_name,fwd_sequence,rev_name,rev_sequence)")
    ap.add_argument("--output", type=Path, required=True,
                    help="Output CSV in IDT bulk-upload format")
    ap.add_argument("--scale", default="25 nm",
                    help="IDT synthesis scale label (default '25 nm')")
    ap.add_argument("--purification", default="STD",
                    help="IDT purification grade (default STD = standard desalt)")
    args = ap.parse_args()

    primers = flatten_primer_csv(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Sequence Name", "Sequence", "Scale", "Purification"])
        for name, seq in primers:
            writer.writerow([name, seq, args.scale, args.purification])

    lengths = [len(s) for _, s in primers]
    print(f"Wrote {len(primers)} primers to {args.output}")
    print(f"  scale={args.scale}  purification={args.purification}")
    print(f"  length range: {min(lengths)}-{max(lengths)} nt "
          f"(median {sorted(lengths)[len(lengths)//2]} nt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
