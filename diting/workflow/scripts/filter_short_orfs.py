#!/usr/bin/env python3
"""Remove ORF DNA sequences shorter than a requested minimum length."""

from __future__ import annotations

import argparse
from pathlib import Path

from fasta_io import read_fasta, write_record


def filter_short_orfs(
    input_fasta: str | Path,
    output_fasta: str | Path,
    minimum_length: int = 100,
) -> None:
    if minimum_length < 1:
        raise ValueError("minimum_length must be at least 1")

    seen_ids: set[str] = set()
    retained = 0
    with open(output_fasta, "w", encoding="utf-8") as output:
        for record in read_fasta(input_fasta):
            if record.identifier in seen_ids:
                raise ValueError(f"Duplicate FASTA ID: {record.identifier}")
            seen_ids.add(record.identifier)
            if len(record.sequence) >= minimum_length:
                write_record(output, record)
                retained += 1

    if not seen_ids:
        raise ValueError("Input FASTA contains no sequences")
    if retained == 0:
        raise ValueError(
            f"No ORFs remain after filtering sequences shorter than {minimum_length} nt"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--minimum-length", type=int, default=100)
    args = parser.parse_args()
    filter_short_orfs(args.input, args.output, args.minimum_length)


if __name__ == "__main__":
    main()
