"""Small, dependency-free FASTA helpers used by the workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def identifier(self) -> str:
        return self.header.split(None, 1)[0]

    @property
    def description(self) -> str:
        parts = self.header.split(None, 1)
        return "" if len(parts) == 1 else " " + parts[1]


def read_fasta(path: str | Path) -> Iterator[FastaRecord]:
    header: str | None = None
    sequence_parts: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence_parts)
                    if not sequence:
                        raise ValueError(f"Empty FASTA sequence for {header!r} in {path}")
                    yield FastaRecord(header, sequence)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before first FASTA header at {path}:{line_number}")
                sequence_parts.append(line)

    if header is not None:
        sequence = "".join(sequence_parts)
        if not sequence:
            raise ValueError(f"Empty FASTA sequence for {header!r} in {path}")
        yield FastaRecord(header, sequence)


def write_record(handle: TextIO, record: FastaRecord, width: int = 60) -> None:
    handle.write(f">{record.header}\n")
    for start in range(0, len(record.sequence), width):
        handle.write(record.sequence[start:start + width] + "\n")


def write_fasta(path: str | Path, records: Iterable[FastaRecord]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for record in records:
            write_record(handle, record)
