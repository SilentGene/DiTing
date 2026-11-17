#!/usr/bin/env python3
"""Merge per-sample Prodigal FFN/FAA files with synchronized unique IDs."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path
from typing import Iterable
from urllib.parse import quote

from fasta_io import FastaRecord, read_fasta, write_record


def merge_orfs(
    sample_files: Iterable[tuple[str, str | Path, str | Path]],
    output_ffn: str | Path,
    output_faa: str | Path,
    id_map: str | Path,
) -> None:
    used_ids: set[str] = set()
    base_occurrences: dict[str, int] = {}
    record_count = 0

    with (
        open(output_ffn, "w", encoding="utf-8") as ffn_out,
        open(output_faa, "w", encoding="utf-8") as faa_out,
        open(id_map, "w", encoding="utf-8") as map_out,
    ):
        map_out.write("sample\toriginal_gene_id\tcatalogue_gene_id\n")

        for sample, ffn_path, faa_path in sorted(sample_files, key=lambda item: item[0]):
            safe_sample = quote(sample, safe="._-")
            ffn_records = read_fasta(ffn_path)
            faa_records = read_fasta(faa_path)
            for position, pair in enumerate(
                itertools.zip_longest(ffn_records, faa_records), start=1
            ):
                ffn_record, faa_record = pair
                if ffn_record is None or faa_record is None:
                    raise ValueError(
                        f"FFN/FAA record count mismatch for sample {sample!r} at record {position}"
                    )
                if ffn_record.identifier != faa_record.identifier:
                    raise ValueError(
                        f"FFN/FAA ID mismatch for sample {sample!r} at record {position}: "
                        f"{ffn_record.identifier!r} != {faa_record.identifier!r}"
                    )

                original_id = ffn_record.identifier
                base_id = f"{safe_sample}__{original_id}"
                occurrence = base_occurrences.get(base_id, 0) + 1
                base_occurrences[base_id] = occurrence
                new_id = base_id if occurrence == 1 else f"{base_id}__dup{occurrence}"
                while new_id in used_ids:
                    occurrence += 1
                    base_occurrences[base_id] = occurrence
                    new_id = f"{base_id}__dup{occurrence}"
                used_ids.add(new_id)

                write_record(
                    ffn_out,
                    FastaRecord(new_id + ffn_record.description, ffn_record.sequence),
                )
                write_record(
                    faa_out,
                    FastaRecord(new_id + faa_record.description, faa_record.sequence),
                )
                map_out.write(f"{sample}\t{original_id}\t{new_id}\n")
                record_count += 1

    if record_count == 0:
        raise ValueError("No ORF records were found in the supplied FFN/FAA files")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orf-dir", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--output-ffn", required=True)
    parser.add_argument("--output-faa", required=True)
    parser.add_argument("--id-map", required=True)
    args = parser.parse_args()

    orf_dir = Path(args.orf_dir)
    sample_files = [
        (sample, orf_dir / f"{sample}.ffn", orf_dir / f"{sample}.faa")
        for sample in args.samples
    ]
    merge_orfs(sample_files, args.output_ffn, args.output_faa, args.id_map)


if __name__ == "__main__":
    main()
