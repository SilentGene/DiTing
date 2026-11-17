#!/usr/bin/env python3
"""Normalize CoverM dense TPM output to DiTing's stable abundance interfaces."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

from fasta_io import read_fasta


def normalize_coverm_tpm(
    coverm_table: str | Path,
    catalogue_fasta: str | Path,
    samples: list[str],
    output_matrix: str | Path,
    abundance_dir: str | Path,
) -> None:
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("Samples must be a non-empty list of unique names")

    catalogue_ids = [record.identifier for record in read_fasta(catalogue_fasta)]
    if len(catalogue_ids) != len(set(catalogue_ids)):
        raise ValueError("Catalogue contains duplicate gene IDs")

    values_by_gene: dict[str, list[str]] = {}
    with open(coverm_table, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if len(header) != len(samples) + 1:
            raise ValueError(
                f"CoverM output has {len(header) - 1} sample columns; expected {len(samples)}"
            )
        for line_number, raw_line in enumerate(handle, start=2):
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != len(samples) + 1:
                raise ValueError(f"Malformed CoverM row at line {line_number}")
            gene_id, raw_values = fields[0], fields[1:]
            if gene_id in values_by_gene:
                raise ValueError(f"Duplicate gene in CoverM output: {gene_id}")
            normalized_values = []
            for value in raw_values:
                try:
                    numeric_value = float(value)
                except ValueError as exc:
                    raise ValueError(
                        f"Non-numeric TPM value {value!r} at line {line_number}"
                    ) from exc
                normalized_values.append(value if math.isfinite(numeric_value) else "0")
            values_by_gene[gene_id] = normalized_values

    unknown = sorted(set(values_by_gene) - set(catalogue_ids))
    if unknown:
        raise ValueError(f"CoverM output contains unknown catalogue gene: {unknown[0]}")

    rows = [
        (gene_id, values_by_gene.get(gene_id, ["0"] * len(samples)))
        for gene_id in catalogue_ids
    ]
    with open(output_matrix, "w", encoding="utf-8") as output:
        output.write("gene_id\t" + "\t".join(samples) + "\n")
        for gene_id, values in rows:
            output.write(gene_id + "\t" + "\t".join(values) + "\n")

    os.makedirs(abundance_dir, exist_ok=True)
    for sample_index, sample in enumerate(samples):
        output_path = Path(abundance_dir) / f"{sample}.abundance"
        with open(output_path, "w", encoding="utf-8") as output:
            output.write("#ID\tgene_abundance\n")
            for gene_id, values in rows:
                output.write(f"{gene_id}\t{values[sample_index]}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--catalogue", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--output-matrix", required=True)
    parser.add_argument("--abundance-dir", required=True)
    args = parser.parse_args()
    normalize_coverm_tpm(
        args.input,
        args.catalogue,
        args.samples,
        args.output_matrix,
        args.abundance_dir,
    )


if __name__ == "__main__":
    main()
