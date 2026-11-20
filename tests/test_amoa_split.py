"""Tests for handling KEGG's split of amoA (K28504) from pmoA (K10944)."""

from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path


WORKFLOW = Path(__file__).parents[1] / "diting" / "workflow"
sys.path.insert(0, str(WORKFLOW))

from scripts.func import (  # noqa: E402
    detect_amoa_ko,
    hierarchical_ko_abundance_among_samples,
)
from scripts.pathway import Pathway  # noqa: E402

KO_LIST_HEADER = "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tdefinition\n"
CYCLE_TABLE = WORKFLOW / "table" / "KO_affilated_to_biogeochemical_cycle.tab"
AMMONIA_STEP = "Nitrification, ammonia -> hydroxylamine (amoABC)"
METHANE_STEP = "Methane oxidation, methane -> methanol (mmoBCDXYZ, amoABC)"


class AmoaSplitTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def write_ko_list(self, kos: list[str], name: str = "ko_list") -> Path:
        path = self.root / name
        body = KO_LIST_HEADER + "".join(f"{ko}\t100\tfull\tall\t1.0\tdesc\n" for ko in kos)
        if name.endswith(".gz"):
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(body)
        else:
            path.write_text(body, encoding="utf-8")
        return path

    def test_detect_legacy_release(self):
        ko_list = self.write_ko_list(["K10944", "K10945", "K10946"])
        self.assertEqual(detect_amoa_ko(ko_list), "K10944")

    def test_detect_split_release(self):
        ko_list = self.write_ko_list(["K10944", "K10945", "K10946", "K28504"])
        self.assertEqual(detect_amoa_ko(ko_list), "K28504")

    def test_detect_split_release_gzipped(self):
        ko_list = self.write_ko_list(["K10944", "K28504"], name="ko_list.gz")
        self.assertEqual(detect_amoa_ko(ko_list), "K28504")

    def test_detect_ignores_k28504_outside_first_column(self):
        path = self.root / "ko_list"
        path.write_text(KO_LIST_HEADER + "K10944\t1\tfull\tall\t1.0\tsee K28504\n", encoding="utf-8")
        self.assertEqual(detect_amoa_ko(path), "K10944")

    def solve(self, amoa_ko: str) -> dict:
        abundance = {"s_K10944": 30.0, "s_K10945": 60.0, "s_K10946": 90.0, "s_K28504": 9.0}
        return Pathway(abundance, "s", amoa_ko).solve_pathway()

    def test_legacy_formula_uses_k10944_for_ammonia_oxidation(self):
        result = self.solve("K10944")
        self.assertAlmostEqual(result[AMMONIA_STEP], (30 + 60 + 90) / 3)

    def test_split_formula_uses_k28504_for_ammonia_oxidation(self):
        result = self.solve("K28504")
        self.assertAlmostEqual(result[AMMONIA_STEP], (9 + 60 + 90) / 3)
        self.assertAlmostEqual(result["Nitrification"], ((9 + 60 + 90) / 3) / 3)

    def test_methane_oxidation_always_uses_k10944(self):
        expected = (30 + 60 + 90) / 3
        self.assertAlmostEqual(self.solve("K10944")[METHANE_STEP], expected)
        self.assertAlmostEqual(self.solve("K28504")[METHANE_STEP], expected)

    def test_default_amoa_ko_is_legacy(self):
        abundance = {"s_K10944": 30.0, "s_K10945": 60.0, "s_K10946": 90.0, "s_K28504": 9.0}
        self.assertAlmostEqual(Pathway(abundance, "s").solve_pathway()[AMMONIA_STEP], 60.0)

    def hierarchy_rows(self, amoa_ko: str) -> list[list[str]]:
        ko_table = self.root / "ko_abundance.tab"
        ko_table.write_text(
            "k_number\tS1\nK10944\t1.5\nK10945\t2.5\nK10946\t3.5\nK28504\t7.5\n",
            encoding="utf-8",
        )
        output = self.root / f"hierarchy_{amoa_ko}.tab"
        hierarchical_ko_abundance_among_samples(ko_table, CYCLE_TABLE, output, amoa_ko)
        return [line.split("\t") for line in output.read_text(encoding="utf-8").splitlines()]

    def ammonia_step_rows(self, rows: list[list[str]]) -> list[list[str]]:
        start = next(i for i, row in enumerate(rows) if row[1].strip('"').startswith("Nitrification, ammonia"))
        section = [rows[start]]
        for row in rows[start + 1:]:
            if row[1]:
                break
            section.append(row)
        return section

    def test_hierarchy_keeps_k10944_for_legacy_release(self):
        section = self.ammonia_step_rows(self.hierarchy_rows("K10944"))
        self.assertEqual([row[2] for row in section], ["K10944", "K10945", "K10946"])
        self.assertEqual(section[0][-1], "1.5")

    def test_hierarchy_swaps_in_k28504_for_split_release(self):
        rows = self.hierarchy_rows("K28504")
        section = self.ammonia_step_rows(rows)
        self.assertEqual([row[2] for row in section], ["K28504", "K10945", "K10946"])
        self.assertTrue(section[0][1].strip('"').startswith("Nitrification, ammonia"))
        self.assertEqual(section[0][-1], "7.5")
        methane_kos = [row[2] for row in rows if row[2] == "K10944"]
        self.assertEqual(len(methane_kos), 1, "pmoA should remain under methane oxidation")


if __name__ == "__main__":
    unittest.main()
