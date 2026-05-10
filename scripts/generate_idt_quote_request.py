"""Flatten an OMEGA paired-primer CSV into IDT's plate-upload .xls format.

OMEGA primer CSVs use the paired layout
``fwd_name,fwd_sequence,rev_name,rev_sequence`` (one row per primer pair).
IDT's bulk-quote tool ingests an Excel workbook with three columns
``Well Position``, ``Name``, ``Sequence`` (one row per oligo, one workbook
per plate). Template at ``data/pricing/example_plate-file-upload.xls``.

Usage:
    uv run python scripts/generate_idt_quote_request.py \\
        --input data/subramanian_orthogonal.csv \\
        --output-prefix data/pricing/idt_quote_request_subramanian_orthogonal

Writes one ``<prefix>_plate<N>.xls`` per 96-well plate. Scale (25 nm)
and Purification (STD) are not template fields — set them when
uploading via IDT's web tool.
"""
from __future__ import annotations

import argparse
import csv
import string
from pathlib import Path

import xlwt


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


def well_positions(rows: int = 8, cols: int = 12):
    """Yield 'A1', 'A2', ..., 'A12', 'B1', ... — row-major plate layout."""
    for r in string.ascii_uppercase[:rows]:
        for c in range(1, cols + 1):
            yield f"{r}{c}"


def write_plate(out_path: Path, plate_primers: list[tuple[str, str]],
                wells_per_plate: int) -> None:
    rows = 8 if wells_per_plate == 96 else 16
    cols = wells_per_plate // rows
    wells = list(well_positions(rows=rows, cols=cols))

    book = xlwt.Workbook(encoding="utf-8")
    sheet = book.add_sheet("Sheet1")
    sheet.write(0, 0, "Well Position")
    sheet.write(0, 1, "Name")
    sheet.write(0, 2, "Sequence")
    for i, (name, seq) in enumerate(plate_primers, start=1):
        sheet.write(i, 0, wells[i - 1])
        sheet.write(i, 1, name)
        sheet.write(i, 2, seq)
    book.save(str(out_path))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, required=True,
                    help="Paired-primer CSV (fwd_name,fwd_sequence,rev_name,rev_sequence)")
    ap.add_argument("--output-prefix", type=Path, required=True,
                    help="Output path stem; '_plate<N>.xls' is appended per plate")
    ap.add_argument("--wells-per-plate", type=int, default=96,
                    choices=[96, 384],
                    help="Plate format (default 96-well)")
    args = ap.parse_args()

    primers = flatten_primer_csv(args.input)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    n_plates = (len(primers) + args.wells_per_plate - 1) // args.wells_per_plate
    for plate_idx in range(n_plates):
        chunk = primers[plate_idx * args.wells_per_plate:
                        (plate_idx + 1) * args.wells_per_plate]
        out_path = Path(f"{args.output_prefix}_plate{plate_idx + 1}.xls")
        write_plate(out_path, chunk, args.wells_per_plate)
        print(f"  Plate {plate_idx + 1}: {len(chunk)} primers -> {out_path}")

    lengths = [len(s) for _, s in primers]
    print(f"\nTotal: {len(primers)} primers across {n_plates} plate(s) "
          f"({args.wells_per_plate}-well)")
    print(f"  length range: {min(lengths)}-{max(lengths)} nt "
          f"(median {sorted(lengths)[len(lengths)//2]} nt)")
    print("  IDT upload notes: set Scale=25 nm, Purification=Standard Desalt "
          "in the web tool")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
