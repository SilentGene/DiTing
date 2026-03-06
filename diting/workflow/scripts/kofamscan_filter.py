#!/usr/bin/env python3
"""
Improved filter for KOFamScan detail-tsv output.

Outputs a tab-delimited table with:
gene_ID, KO_ID, KO_function, above_thrshld

Rules:
- Do NOT discard hits solely because score < threshold.
- Optionally filter hits by a global maximum e-value (--evalue / -E).
- For each gene, keep only the best KO hit by:
    1) smaller e-value wins
    2) if e-value ties, larger score wins
    3) if still ties, keep the first encountered
- above_thrshld is:
    - "yes" if score >= thrshld
    - "no"  if score <  thrshld
    - "NA"  if thrshld is empty/unavailable

Sample input:
#	gene name	KO	thrshld	score	E-value	KO definition
#	---------	------	-------	------	---------	-------------
*	Inre_C_10_12_Ctg7576_5	K00399	775.53	965.1	1.40E-292	methyl-coenzyme M reductase alpha subunit [EC:2.8.4.1]
	Inre_C_10_12_Ctg7576_5	K00401	504.3	11.5	0.0012	methyl-coenzyme M reductase beta subunit [EC:2.8.4.1]
*	Inre_C_10_12_Ctg8269_2	K00399	775.53	1009.8	0	methyl-coenzyme M reductase alpha subunit [EC:2.8.4.1]
*	Inre_C_10_12_Ctg15638_1	K00399	775.53	969.5	6.20E-294	methyl-coenzyme M reductase alpha subunit [EC:2.8.4.1]
	Inre_C_10_12_Ctg442556_1	K00399	775.53	154.3	5.80E-47	methyl-coenzyme M reductase alpha subunit [EC:2.8.4.1]
	Inre_C_10_12_Ctg442556_1	K21050		13.2	0.00064	HCMV unique short US2 glycoprotein
	Inre_C_10_12_Ctg442556_1	K12813	1207.5	9.9	0.0021	pre-mRNA-splicing factor ATP-dependent RNA helicase DHX16 [EC:5.6.2.6]
	Inre_C_10_12_Ctg547855_1	K00399	775.53	154.3	5.60E-47	methyl-coenzyme M reductase alpha subunit [EC:2.8.4.1]

"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Dict, Optional


@dataclass
class BestHit:
    ko_id: str
    ko_function: str
    evalue: float
    score: float
    thrshld: Optional[float]  # None if missing
    order: int  # line order for stable tie-breaking


def parse_evalue(s: str) -> float:
    """
    Parse e-values like '6.7e-11' or '0.0035' to float.
    If parsing fails, return +inf so it will never be selected as best.
    """
    try:
        return float(s)
    except Exception:
        return math.inf


def parse_float_optional(s: str) -> Optional[float]:
    s = s.strip()
    if s == "":
        return None
    try:
        return float(s)
    except Exception:
        return None


def above_threshold_label(score: float, thrshld: Optional[float]) -> str:
    if thrshld is None:
        return "NA"
    return "yes" if score >= thrshld else "no"


def is_better_hit(new: BestHit, old: BestHit) -> bool:
    # 1) smaller evalue
    if new.evalue < old.evalue:
        return True
    if new.evalue > old.evalue:
        return False
    # 2) higher score
    if new.score > old.score:
        return True
    if new.score < old.score:
        return False
    # 3) keep first => only replace if new.order < old.order (won't happen in streaming)
    return False


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Improved filter for KOFamScan 'detail-tsv' output. "
            "Keeps one best KO per gene; optionally filters by e-value; "
            "adds above_thrshld column."
        )
    )

    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input-file",
        help="Input KOFamScan detail-tsv table",
        metavar="<FILE>",
        required=True
    )
    parser.add_argument(
        "-o", "--output-file",
        help='Output table filename (default: "output.tsv")',
        metavar="<FILE>",
        default="output.tsv"
    )
    parser.add_argument(
        "-s", "--sample",
        help="Sample name to prepend to output rows",
        metavar="<STR>",
        required=True
    )
    parser.add_argument(
        "-E", "--evalue",
        help="Maximum allowed E-value to consider a hit (e.g. 1e-5). "
             "Hits with evalue > this will be ignored for choosing the best KO.",
        metavar="<FLOAT>",
        type=float,
        default=None
    )

    args = parser.parse_args()

    # Track all genes encountered (even those without hits) to ensure output rows exist.
    seen_genes = set()

    best_by_gene: Dict[str, BestHit] = {}

    order = 0
    with open(args.input_file, "r", encoding="utf-8") as annots:
        for raw_line in annots:
            if raw_line.startswith("#"):
                continue

            clean = raw_line.lstrip("*").lstrip().rstrip("\n")
            line = clean.split("\t")
            if len(line) == 0:
                continue

            gene_id = line[0].strip()
            if gene_id == "":
                continue

            seen_genes.add(gene_id)

            # Some lines may contain only the gene_id (unannotated)
            if len(line) == 1:
                continue

            # Expected columns:
            # 0 gene, 1 KO, 2 thrshld, 3 score, 4 evalue, 5 definition
            # Guard against malformed lines:
            if len(line) < 6:
                continue

            ko_id = line[1].strip() or "NA"
            thrshld = parse_float_optional(line[2])
            score = parse_float_optional(line[3])
            evalue = parse_evalue(line[4])
            ko_def = line[5].strip().strip('"') if line[5] is not None else "NA"

            if score is None:
                # can't rank by score reliably
                continue

            # Global e-value filtering (optional)
            if args.evalue is not None and evalue > args.evalue:
                continue

            candidate = BestHit(
                ko_id=ko_id,
                ko_function=ko_def,
                evalue=evalue,
                score=float(score),
                thrshld=thrshld,
                order=order
            )
            order += 1

            if gene_id not in best_by_gene:
                best_by_gene[gene_id] = candidate
            else:
                if is_better_hit(candidate, best_by_gene[gene_id]):
                    best_by_gene[gene_id] = candidate

    with open(args.output_file, "w", encoding="utf-8") as out:
        for gene_id in sorted(seen_genes):
            if gene_id not in best_by_gene:
                continue

            hit = best_by_gene[gene_id]
            if hit.ko_id == "NA":
                continue
                
            out.write(f"{args.sample}\t{gene_id}\t{hit.ko_id}\n")


if __name__ == "__main__":
    main()