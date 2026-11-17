from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


WORKFLOW = Path(__file__).parents[1] / "diting" / "workflow"
SCRIPTS = WORKFLOW / "scripts"
sys.path.insert(0, str(WORKFLOW))
sys.path.insert(0, str(SCRIPTS))

from fasta_io import read_fasta  # noqa: E402
from filter_short_orfs import filter_short_orfs  # noqa: E402
from merge_orfs import merge_orfs  # noqa: E402
from normalize_coverm_tpm import normalize_coverm_tpm  # noqa: E402
from scripts.func import merge_abun_ko  # noqa: E402


class CatalogueScriptTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def write(self, name: str, contents: str) -> Path:
        path = self.root / name
        path.write_text(contents, encoding="utf-8")
        return path

    def test_merge_prefixes_ids_and_handles_duplicates(self):
        s1_ffn = self.write("s1.ffn", ">gene_1 # dna\nAAA\n>gene_1\nCCC\n")
        s1_faa = self.write("s1.faa", ">gene_1 # protein\nK\n>gene_1\nP\n")
        s2_ffn = self.write("s2.ffn", ">gene_1\nGGG\n")
        s2_faa = self.write("s2.faa", ">gene_1\nG\n")
        out_ffn, out_faa, id_map = (
            self.root / "all.ffn",
            self.root / "all.faa",
            self.root / "ids.tsv",
        )

        merge_orfs(
            [("s2", s2_ffn, s2_faa), ("s1", s1_ffn, s1_faa)],
            out_ffn,
            out_faa,
            id_map,
        )

        self.assertEqual(
            [record.identifier for record in read_fasta(out_ffn)],
            ["s1__gene_1", "s1__gene_1__dup2", "s2__gene_1"],
        )
        self.assertIn(">s1__gene_1 # dna", out_ffn.read_text(encoding="utf-8"))
        self.assertEqual(
            [record.identifier for record in read_fasta(out_faa)],
            ["s1__gene_1", "s1__gene_1__dup2", "s2__gene_1"],
        )

    def test_merge_rejects_ffn_faa_mismatch(self):
        ffn = self.write("bad.ffn", ">gene_1\nAAA\n")
        faa = self.write("bad.faa", ">gene_2\nK\n")
        with self.assertRaisesRegex(ValueError, "ID mismatch"):
            merge_orfs(
                [("sample", ffn, faa)],
                self.root / "out.ffn",
                self.root / "out.faa",
                self.root / "map.tsv",
            )

    def test_merge_encodes_whitespace_in_sample_name(self):
        ffn = self.write("spaced.ffn", ">gene_1\nAAA\n")
        faa = self.write("spaced.faa", ">gene_1\nK\n")
        output = self.root / "out.ffn"
        merge_orfs(
            [("sample one", ffn, faa)],
            output,
            self.root / "out.faa",
            self.root / "map.tsv",
        )
        self.assertEqual(next(read_fasta(output)).identifier, "sample%20one__gene_1")

    def test_filter_keeps_100_nt_and_longer(self):
        source = self.write(
            "all.ffn",
            f">short\n{'A' * 99}\n>exact\n{'C' * 100}\n>long\n{'G' * 101}\n",
        )
        output = self.root / "filtered.ffn"
        filter_short_orfs(source, output, 100)
        self.assertEqual(
            [record.identifier for record in read_fasta(output)], ["exact", "long"]
        )

    def test_filter_rejects_duplicate_ids(self):
        source = self.write("duplicate.ffn", ">gene\nAAA\n>gene\nCCC\n")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            filter_short_orfs(source, self.root / "filtered.ffn", 1)

    def test_normalize_coverm_adds_missing_zero_rows(self):
        catalogue = self.write("nr.ffn", ">gene1\nAAA\n>gene2\nCCC\n")
        coverm = self.write("coverm.tsv", "Contig\treads TPM\ngene1\t1000000\n")
        matrix = self.root / "matrix.tsv"
        abundance_dir = self.root / "Abundance"
        normalize_coverm_tpm(coverm, catalogue, ["sampleA"], matrix, abundance_dir)
        self.assertEqual(
            matrix.read_text(encoding="utf-8"),
            "gene_id\tsampleA\ngene1\t1000000\ngene2\t0\n",
        )
        self.assertIn(
            "gene2\t0", (abundance_dir / "sampleA.abundance").read_text(encoding="utf-8")
        )

    def test_normalize_coverm_converts_non_finite_values_to_zero(self):
        catalogue = self.write("nr.ffn", ">gene1\nAAA\n")
        coverm = self.write("coverm.tsv", "Contig\treads TPM\ngene1\tNaN\n")
        matrix = self.root / "matrix.tsv"
        normalize_coverm_tpm(
            coverm, catalogue, ["sampleA"], matrix, self.root / "Abundance"
        )
        self.assertIn("gene1\t0", matrix.read_text(encoding="utf-8"))

    def test_merge_abundance_and_catalogue_annotations(self):
        matrix = self.write(
            "matrix.tsv",
            "gene_id\tsampleA\tsampleB\ngene1\t750000\t0\ngene2\t250000\t1000000\n",
        )
        annotations = self.write(
            "ko_merged.txt", "#gene_id\tk_number\ngene1\tK00001\ngene2\tK00002\n"
        )
        output = self.root / "ko_abun.txt"
        merge_abun_ko(matrix, annotations, output)
        self.assertEqual(
            output.read_text(encoding="utf-8"),
            "#sample\tk_number\trelative_abundance\tgene_id\n"
            "sampleA\tK00001\t750000\tgene1\n"
            "sampleB\tK00001\t0\tgene1\n"
            "sampleA\tK00002\t250000\tgene2\n"
            "sampleB\tK00002\t1000000\tgene2\n",
        )


if __name__ == "__main__":
    unittest.main()
